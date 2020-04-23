/*
 * Copyright (C) 2012 Kengo Sato
 *
 * This file is part of DAFS.
 *
 * DAFS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DAFS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DAFS.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __INC_TYPEDEFS_H__
#define __INC_TYPEDEFS_H__

#include <vector>

typedef unsigned int uint;

typedef std::vector<float> VF;
typedef std::vector<VF> VVF;

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;

typedef std::vector<uint> VU;
typedef std::vector<VU> VVU;

typedef std::vector<char> VC;
typedef std::vector<VC> VVC;
  
typedef std::vector<std::pair<uint,float> > SV; // sparse vectors
typedef std::vector<SV> MP;
typedef std::vector<SV> BP;

typedef std::vector< std::pair<uint, std::vector<bool> > > ALN; // alignments

#endif

// Local Variables:
// mode: C++
// End:
