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

#include "../data/block/block.h"
#include "../align/extend.h"
#include "../lib/ksw2/ksw2.h"



#pragma once


namespace Dna{
const int KSW2_END_BONUS = 5;
const int KSW2_BAND = 64;
const int WFA_CUTOFF_STEPS = 10;
const int KSW_FLAG_R = KSW_EZ_EXTZ_ONLY | KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
const int KSW_FLAG_L = KSW_FLAG_R | KSW_EZ_REV_CIGAR;

std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg, const Sequence &query);//, ExtensionTimer& extend_time);
}
