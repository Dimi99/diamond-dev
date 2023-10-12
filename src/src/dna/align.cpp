#include "align.h"
#include <cstring>
#include "../lib/ksw2/ksw2.h"
#include "seed_set_dna.h"
#include "../util/util.h"
#include "bindings/cpp/WFAligner.hpp"
#include "extension.h"
#include "../util/sequence/translate.h"
#include <queue>
#include <mutex>
#include <condition_variable>
#include <iostream>
#include "../align/target.h"


std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,const Sequence &query) {
    task_timer timer(flag_any(flags, DP::Flags::PARALLEL) ? config.target_parallel_verbosity : UINT_MAX);
//Loading Seed Hits
    timer.go("Loading seed hits");
    Extension::SeedHitList l = load_hits(begin, end, cfg.target->seqs());
    auto hits = seed_lookup(query, cfg.target->seqs(), cfg.dna_ref_index.get(), cfg.minimizer_window);
    stat.inc(Statistics::TARGET_HITS0, l.target_block_ids.size());
    stat.inc(Statistics::TIME_LOAD_HIT_TARGETS, timer.microseconds());
    timer.finish();

    const int64_t target_count = (int64_t)l.target_block_ids.size();
    if (target_count == 0 && !config.swipe_all) {
        if (add_self_aln(cfg))
            return { {Match::self_match(query_id, cfg.query->seqs()[query_id])}, Stats() };
        return { {}, Stats() };
    }
    const int64_t chunk_size = ranking_chunk_size(target_count, cfg.target->seqs().letters(), cfg.max_target_seqs);

    if (chunk_size < target_count || config.global_ranking_targets > 0) {
        timer.go("Sorting targets by score");
        std::sort(l.target_scores.begin(), l.target_scores.end());
        stat.inc(Statistics::TIME_SORT_TARGETS_BY_SCORE, timer.microseconds());
        timer.finish();
        if (config.global_ranking_targets > 0)
            return make_pair(GlobalRanking::ranking_list(query_id, l.target_scores.begin(), l.target_scores.end(), l.target_block_ids.begin(), l.seed_hits, cfg), Stats());
    }

    pair<vector<Match>, Stats> r = extend(query_id, cfg, stat, flags, l);
    return r;
}