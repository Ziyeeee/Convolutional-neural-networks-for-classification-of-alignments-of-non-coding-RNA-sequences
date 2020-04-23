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

#ifndef __INC_ALIFOLD_H__
#define __INC_ALIFOLD_H__

#include <vector>
#include <string>
#include "fa.h"
#include "typedefs.h"

class Alifold
{
public:
  Alifold(float th) : th_(th) { }
  void fold(const ALN& aln, const std::vector<Fasta>& fa, BP& bp) const;
  void fold(const ALN& aln, const std::vector<Fasta>& fa,
            const std::string& str, BP& bp) const;

  float energy_of_struct(const ALN& aln, const std::vector<Fasta>& fa,
                         const std::string& str, float& cv) const;
  
private:
  static char** alloc_aln(const ALN& aln, const std::vector<Fasta>& fa);
  static void free_aln(char** seqs);
  
private:
  float th_;
};

#endif  //  __INC_ALIFOLD_H__

// Local Variables:
// mode: C++
// End:
