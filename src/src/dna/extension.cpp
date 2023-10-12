/****
DIAMOND protein aligner
Copyright (C) 2022 Dimitrios Koutsogiannis

Code developed by Dimitrios Koutsogiannis <dimitrios.koutsogiannis@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "../lib/ksw2/ksw2.h"
#include "seed_set_dna.h"
#include "../util/util.h"
#include "bindings/cpp/WFAligner.hpp"
#include "extension.h"
#include "../util/sequence/translate.h"
#include "timer.h"
#include "../search/finger_print.h"




const EMap<DNAExtensionAlgo> EnumTraits<DNAExtensionAlgo>::to_string = { {DNAExtensionAlgo::KSW,"ksw"}, { DNAExtensionAlgo::WFA,"wfa"} };
const SEMap<DNAExtensionAlgo> EnumTraits<DNAExtensionAlgo>::from_string = { {"ksw",DNAExtensionAlgo::KSW}, {"wfa",DNAExtensionAlgo::WFA} };


namespace Dna {


    struct ExtendedSeed {
        ExtendedSeed(int i_min, int i_max, int j_min, int j_max) :
                i_min_(i_min),
                i_max_(i_max),
                j_min_(j_min),
                j_max_(j_max) {}

        int i_min_extended() const { return i_min_; };

        int i_max_extended() const { return i_max_; };

        int j_min_extended() const { return j_min_; };

        int j_max_extended() const { return j_max_; };

    private:
        int i_min_, i_max_, j_min_, j_max_;
    };

    struct Cigar {
    public:

        Cigar &operator+(const Cigar &other) {
            this->cigar_data.insert(this->cigar_data.end(), other.cigar_data.begin(), other.cigar_data.end());
            this->score_ += other.score_;
            return *this;
        }

        int32_t score() const { return score_; }

        int max_query() const { return max_query_; }

        int max_target() const { return max_target_; }

        std::vector<std::pair<int, char>> cigar_data;

    protected:
        Cigar() = default;

        int32_t score_{0};
        int max_query_{0}, max_target_{0};
    };

    struct KswCigar : public Cigar {
        KswCigar(const Sequence &tseq, const Sequence &qseq, const Search::Config &cfg, int flag) {


            int a = cfg.score_builder->reward(), b = cfg.score_builder->penalty() < 0 ? cfg.score_builder->penalty()
                                                                                      : -cfg.score_builder->penalty(); // a>0 and b<0
            int8_t mat[NUCLEOTIDE_COUNT * NUCLEOTIDE_COUNT] = {static_cast<int8_t>(a), static_cast<int8_t>(b),
                                                               static_cast<int8_t>(b),
                                                               static_cast<int8_t>(b), 0, static_cast<int8_t>(b),
                                                               static_cast<int8_t>(a),
                                                               static_cast<int8_t>(b), static_cast<int8_t>(b), 0,
                                                               static_cast<int8_t>(b),
                                                               static_cast<int8_t>(b), static_cast<int8_t>(a),
                                                               static_cast<int8_t>(b), 0,
                                                               static_cast<int8_t>(b), static_cast<int8_t>(b),
                                                               static_cast<int8_t>(b),
                                                               static_cast<int8_t>(a), 0, 0, 0, 0, 0, 0};
            ksw_extz_t ez;
            memset(&ez, 0, sizeof(ksw_extz_t));


            ksw_extz2_sse(nullptr, qseq.length(), reinterpret_cast<const uint8_t *>(qseq.data()), tseq.length(),
                          reinterpret_cast<const uint8_t *>(tseq.data()), NUCLEOTIDE_COUNT, mat,
                          (int8_t) cfg.score_builder->gap_open(), (int8_t) cfg.score_builder->gap_extend(), KSW2_BAND,
                          config.zdrop, KSW2_END_BONUS, flag, &ez);


            this->score_ = ez.max;
            this->max_query_ = ez.max_q;
            this->max_target_ = ez.max_t;
            cigar_data.reserve(ez.n_cigar);
            for (int i = 0; i < ez.n_cigar; ++i)
                cigar_data.emplace_back(ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);


            free(ez.cigar);

        }
    };


    struct WfaCigar : public Cigar {
        WfaCigar(const Sequence &tseq, const Sequence &qseq, const Search::Config &cfg, bool left = false) {

            std::string tseq2 = tseq.to_string();
            std::string qseq2 = qseq.to_string();

            int tl = tseq2.length(), ql = qseq2.length();

            thread_local std::unique_ptr<wfa::WFAlignerGapAffine> aligner;
            if (!aligner)
                aligner.reset(new wfa::WFAlignerGapAffine(0, -config.mismatch_penalty, config.gap_open, config.gap_extend,
                                                    wfa::WFAligner::Alignment));
            aligner->setHeuristicNone();
            aligner->setHeuristicZDrop(config.zdrop, WFA_CUTOFF_STEPS);
            aligner->setHeuristicWFadaptive(10,50,1);

            aligner->alignExtension(tseq2.c_str(), tl, qseq2.c_str(), ql);

            auto cigar = aligner->getCIGAR(true);


            max_query_ = -1;
            max_target_ = -1;
            int steps = 0;
            for (char c: cigar) {
                if (isdigit(c)) {
                    steps = steps * 10 + (c - '0');
                    continue;
                }
                cigar_data.emplace_back(steps, c);
                switch (c) {
                    case '=':
                        this->score_ += cfg.score_builder->reward() * steps;
                        max_query_ += steps;
                        max_target_ += steps;
                        break;
                    case 'X':
                        this->score_ += cfg.score_builder->penalty() * steps;
                        max_query_ += steps;
                        max_target_ += steps;
                        break;
                    case 'I':
                        this->score_ -= cfg.score_builder->gap_open() + steps * cfg.score_builder->gap_extend();
                        max_query_ += steps;
                        break;
                    case 'D':
                        this->score_ -= cfg.score_builder->gap_open() + steps * cfg.score_builder->gap_extend();
                        max_target_ += steps;
                        break;
                    default:
                        throw std::runtime_error(std::string("WFA Cigar: Invalid Cigar Symbol ") + c);

                }
                steps = 0;
            }
            if (left)
                std::reverse(cigar_data.begin(), cigar_data.end());
        }
    };


    bool intersection(SeedMatch &hit, const std::vector<ExtendedSeed> &extended, int window) {
        return std::any_of(extended.begin(), extended.end(), [hit, window](ExtendedSeed s) {
            return hit.i() + 1 >= s.i_min_extended() && hit.i() + window <= s.i_max_extended() &&
                   hit.j() + 1 >= s.j_min_extended() && hit.j() + window <= s.j_max_extended();
        });
    }

 void calculate_ungapped_scores(SeedMatch &hit, const Sequence &target, const Sequence &query) {

#ifdef __SSE2__
        int i = 0;
        int score=0;
        alignas(16) __m128i q, t;
        while(hit.i() - i*16 >= 0 && hit.j() - i*16 >=0){
            q = letter_mask(_mm_loadu_si128((__m128i const*)(query.data()+hit.i() - (i+1)*16)));
            t = letter_mask(_mm_loadu_si128((__m128i const*)(target.data()+hit.j() - (i+1)*16)));
            auto equal = popcount64((uint64_t)_mm_movemask_epi8(_mm_cmpeq_epi8(q, t)));
            if(equal != 16) break;
            score+= 16;
            i+=16;
        }
        while (query[hit.i() - i] == target[hit.j() - i]) {
            score++;
            i++;
        }
        hit.l_ungapped(score);

        i = 0;

        while(hit.i() + i*16 < query.length() && hit.j() + i*16 < target.length()){
            q = letter_mask(_mm_loadu_si128((__m128i const*)(query.data()+hit.i()+1 + i*16)));
            t = letter_mask(_mm_loadu_si128((__m128i const*)(target.data()+hit.j()+1 + i*16)));
            auto equal = popcount64((uint64_t)_mm_movemask_epi8(_mm_cmpeq_epi8(q, t)));
            if(equal!=16) break;
            score+= 16;
            i+=16;
        }
        ++i;
        while (query[hit.i() + i] == target[hit.j() + i]) {
            score++;
            i++;
        }
        hit.score(score);


#else
        int score = 0;
        int i = 0;

        while (query[hit.i() - i] == target[hit.j() - i]) {
            score++;
            i++;
        }
        hit.l_ungapped(score);

        i = 1;
        while (query[hit.i() + i] == target[hit.j() + i]) {
            score++;
            i++;
        }
        hit.score(score);

#endif
    }


    void cigar_to_hsp(const Cigar &cigar, const Sequence &target, const Sequence &query, const SeedMatch &hit, Hsp &out,
                      bool reverse) {

        int pattern_pos = hit.i() - cigar.max_query() - 1;
        int text_pos = hit.j() - cigar.max_target() - 1;
        out.query_range.begin_ = pattern_pos;
        out.subject_range.begin_ = text_pos;

        for (auto operation: cigar.cigar_data) {
            switch (operation.second) {
                case 'M':
                case '=':
                case 'X':
                    for (int j = 0; j < operation.first; ++j) {
                        out.push_match(target[text_pos], query[pattern_pos], true);
                        pattern_pos++;
                        text_pos++;
                    }
                    break;
                case 'D':
                    out.push_gap(op_deletion, operation.first, target.data() + operation.first + text_pos);
                    text_pos += operation.first;
                    break;
                case 'I':
                    out.push_gap(op_insertion, operation.first, query.data() + operation.first + pattern_pos);
                    pattern_pos += operation.first;
                    break;
                default:
                    break;
            }
        }

        out.query_range.end_ = pattern_pos;
        out.subject_range.end_ = text_pos;
        out.transcript.push_terminator();
        out.target_seq = target;
        out.query_source_range = out.query_range;
        out.subject_source_range = reverse ? Interval(out.subject_range.end_, out.subject_range.begin_) : Interval(
                out.subject_range.begin_, out.subject_range.end_);
        out.frame = reverse - 1;
    }

    void cigar_to_hsp(const Sequence &target, const Sequence &query, const SeedMatch &hit, Hsp &out, bool reverse) {
        int pattern_pos = hit.i() - hit.l_ungapped();
        int text_pos = hit.j() - hit.l_ungapped();
        out.query_range.begin_ = pattern_pos;
        out.subject_range.begin_ = text_pos;

        for (int i = 0; i <= hit.ungapped_score(); ++i) {
            out.push_match(target[text_pos], query[pattern_pos], true);
            pattern_pos++;
            text_pos++;
        }

        out.query_range.end_ = pattern_pos;
        out.subject_range.end_ = text_pos;
        out.transcript.push_terminator();
        out.target_seq = target;
        out.query_source_range = out.query_range;
        out.subject_source_range = reverse ? Interval(out.subject_range.end_, out.subject_range.begin_) : Interval(
                out.subject_range.begin_, out.subject_range.end_);
        out.frame = reverse - 1;

    }

    Extension::Match
    target_extension(const Search::Config &cfg, const BlockId id, const Sequence &query, SeedMatch *begin,
                     SeedMatch *end, bool reverse) {
        std::vector<ExtendedSeed> extended{};

        Extension::Match m = Extension::Match(id, cfg.target->seqs()[id], ::Stats::TargetMatrix(), 0, 0);
        for (auto hit = begin; hit != end; ++hit) {
            if ((intersection(*hit, extended, cfg.minimizer_window)))
                continue;
            auto q = Byte_finger_print_48(query.data()+hit->i());
            const Sequence &target = cfg.target->seqs()[id];
            //Hamming Filter
            if(q.match({target.data()+hit->j()}) < 40) continue;
                if (hit->ungapped_score() == query.length() - 1) {
                Hsp out = Hsp();

                out.score = hit->ungapped_score() * cfg.score_builder->reward();
                out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
                out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());
                if (out.evalue >= config.max_evalue)
                    continue;

                cigar_to_hsp(cfg.target->seqs()[id], query, *hit, out, reverse);
                m.hsp.push_back(out);
                extended.emplace_back(hit->i() - hit->l_ungapped(),
                                      out.query_range.end_,
                                      hit->j() - hit->l_ungapped(),
                                      out.subject_range.end_);

            } else {
                Sequence query_right = query.subseq(hit->i(), query.length());
                Sequence target_right = target.subseq(hit->j(), std::min(target.length(),
                                                                         hit->j() + (int) query_right.length() * 2));

                std::vector<Letter> query_left = query.subseq(0, hit->i()).reverse();
                std::vector<Letter> target_left = target.subseq(std::max(0, (int) (hit->j() - (query_left.size() * 2))),
                                                                hit->j()).reverse();

                ExtensionTimer timer;
                auto start_extend = std::chrono::high_resolution_clock::now();
                Cigar extension = config.dna_extension == DNAExtensionAlgo::WFA ?
                                  static_cast<Cigar>(WfaCigar(target_left, query_left, cfg, 1) +
                                                     WfaCigar(target_right, query_right, cfg))
                                                                                :
                                  KswCigar(target_left, query_left, cfg, KSW_FLAG_L) +
                                  KswCigar(target_right, query_right, cfg, KSW_FLAG_R);
                auto end_extend = std::chrono::high_resolution_clock::now();
                timer.update(1,end_extend-start_extend);


                Hsp out = Hsp();
                out.score = extension.score();
                out.bit_score = cfg.score_builder->blast_bit_Score(out.score);
                out.evalue = cfg.score_builder->blast_eValue(out.score, query.length());
                if (out.evalue >= config.max_evalue)
                    continue;

                cigar_to_hsp(extension, cfg.target->seqs()[id], query, *hit, out, reverse);

                m.hsp.push_back(out);
                extended.emplace_back(hit->i() - extension.max_query(),
                                      out.query_range.end_,
                                      hit->j() - extension.max_target(),
                                      out.subject_range.end_);

            }
        }
        return m;
    }

    std::vector<Extension::Match> query_extension(const Search::Config &cfg, const Sequence &query, bool reverse) {
        std::vector<Extension::Match> matches;

        auto hits = seed_lookup(query, cfg.target->seqs(), cfg.dna_ref_index.get(), cfg.minimizer_window);
        std::vector<int> tmp;
        std::for_each(hits.begin(), hits.end(), [&](SeedMatch &hit) {
            calculate_ungapped_scores(hit, cfg.target->seqs()[hit.id()], query);
        });

        std::sort(hits.begin(), hits.end(), std::greater<>());
        auto it0 = merge_keys(hits.begin(), hits.end(), [](const SeedMatch &hit1) { return hit1.id(); });
        while(it0.good()) {
            for(auto a = it0.begin(); a < it0.end(); ++a) {
                std::cout << a->id()<< "  ||  " << a->j() << "\n";
            }
            ++it0;
            std::cout << "------------------------------------------\n";
        }

        std::sort(hits.begin(), hits.end(), std::greater<>());

        auto it = merge_keys(hits.begin(), hits.end(), [](const SeedMatch &hit1) { return hit1.id(); });
        while (it.good()){
                Extension::Match m = target_extension(cfg, it.key(), query,it.begin().operator->(),it.end().operator->(),reverse);
                if (!m.hsp.empty()) {
                    matches.push_back(m);
                }
            ++it;
        }

        return matches;
    }
    std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,const Sequence &query) {

        std::vector<Extension::Match> matches = query_extension(cfg,query, false);
        std::vector<Extension::Match> reverse = query_extension(cfg, Translator::reverse(query), true);

        matches.insert( matches.end(), reverse.begin(), reverse.end());

        return {matches, Extension::Stats()};
    }
}






