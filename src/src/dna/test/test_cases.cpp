/****
DIAMOND protein aligner
Copyright (C) 2019-2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include "test.h"

using std::vector;

namespace DnaTest {

    const vector<TestCase> test_cases = {
            { "blastn (default)", "blastn -p1" },
            { "blastn (multithreaded)", "blastn -p4" },
            { "blastn (blocked)", "blastn -c1 -b0.00002 -p4" },
            { "blastn (target seqs)", "blastn -k3 -c1 -p4" },
            { "blastn (evalue)", "blastn -e10000 -c1 -p4" },
            { "blastn (pairwise format)", "blastn -c1 -f0 -p4" },
            { "blastn (XML format)", "blastn -c1 -f xml -p4" },
            {"blastn (Benchmark Configuration)", "blastn --band 70 -z 70 -p 50"},
            { "blastn (PAF format)", "blastn -c1 -f paf -p1" }
    };

    const vector<uint64_t> ref_hashes = {
            0x36bf16afef49c7ad,
            0x36bf16afef49c7ad,
            0x36bf16afef49c7ad,
            0x7ed13391c638dc2e,
            0x61ac7ee1bb73d36d,
            0xd62b1c97fb27608f,
            0x2dd4b2985c1bebd2,
            0x7ed13391c638dc2e,
            0x98b810039a4e6e1b,
            0x9a20976998759371,
            0xa67de9d0530d5968,
            0xa67de9d0530d5968,
            0x3d593e440ca8eb97,
            0x487a213a131d4958,
            0x201a627d0d128fd5,
            0xe787dcb23cc5b120,
            0x5aa4baf48a888be9,
            0x21f14583e88a13ac,
            0xe5eb09929f0cc182,
            0x713deb9a5ae4b9e,
    };

}