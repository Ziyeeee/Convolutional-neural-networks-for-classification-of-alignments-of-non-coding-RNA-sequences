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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "alifold.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
};
};

#define FOREACH(itr, i, v) for (itr i=(v).begin(); i!=(v).end(); ++i)

#ifdef HAVE_VIENNA18
typedef Vienna::plist pair_info;
#else
typedef Vienna::pair_info pair_info;
#endif

void
Alifold::
fold(const ALN& aln, const std::vector<Fasta>& fa, BP& bp) const
{
  const uint L=aln.front().second.size();
  //const uint N=aln.size();

  char** seqs = alloc_aln(aln, fa);
    
  // scaling parameters to avoid overflow
  std::string res(L+1, ' ');
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, &res[0]);
#else
  double min_en = Vienna::alifold(seqs, &res[0]);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

  pair_info* pi;
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, NULL, &pi);
#else
  Vienna::alipf_fold(seqs, NULL, &pi);
#endif

  bp.resize(L);
  for (uint k=0; pi[k].i!=0; ++k)
    if (pi[k].p>th_)
      bp[pi[k].i-1].push_back(std::make_pair(pi[k].j-1, pi[k].p));

  free(pi);
  Vienna::free_alipf_arrays();
  free_aln(seqs);
}

void
Alifold::
fold(const ALN& aln, const std::vector<Fasta>& fa, const std::string& str, BP& bp) const
{
  const uint L=aln.front().second.size();
  std::string p(str);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  char** seqs = alloc_aln(aln, fa);

  // scaling parameters to avoid overflow
  std::string res(p);
#ifdef HAVE_VIENNA20
  double min_en = Vienna::alifold((const char**)seqs, &res[0]);
#else
  double min_en = Vienna::alifold(seqs, &res[0]);
#endif
  double kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
  Vienna::pf_scale = exp(-(1.07*min_en)/kT/L);
  Vienna::free_alifold_arrays();

  pair_info* pi;
#ifdef HAVE_VIENNA20
  Vienna::alipf_fold((const char**)seqs, &p[0], &pi);
#else
  Vienna::alipf_fold(seqs, &p[0], &pi);
#endif

  bp.resize(L);
  for (uint k=0; pi[k].i!=0; ++k)
    if (pi[k].p>th_)
      bp[pi[k].i-1].push_back(std::make_pair(pi[k].j-1, pi[k].p));

  free(pi);
  Vienna::free_alipf_arrays();
  free_aln(seqs);
  Vienna::fold_constrained = bk;
}

float
Alifold::
energy_of_struct(const ALN& aln, const std::vector<Fasta>& fa,
                 const std::string& str, float& cv) const
{
  const uint L=aln.front().second.size();
  const uint N=aln.size();

  std::string res(L+1, ' ');
  char** seqs = alloc_aln(aln, fa);
    
#ifdef HAVE_VIENNA20
  float min_en = Vienna::energy_of_alistruct((const char **)seqs, str.c_str(), N, &cv);
#else
  std::vector<float> cv_temp(2);
  float min_en = Vienna::energy_of_alistruct(seqs, str.c_str(), N, &cv_temp[0]);
  cv = cv_temp[1];
#endif
  free_aln(seqs);
  return min_en;
}

//static
char**
Alifold::
alloc_aln(const ALN& aln, const std::vector<Fasta>& fa)
{
  const uint L=aln.front().second.size();
  const uint N=aln.size();

  char** seqs = new char*[N+1];
  seqs[N] = NULL;
  char** s = seqs;
  FOREACH (ALN::const_iterator, it, aln)
  {
    *s = new char[L+1];
    for (uint i=0, j=0; i!=L; ++i)
      (*s)[i] = it->second[i] ? fa[it->first].seq()[j++] : '-';
    (*s)[L] = 0;
    ++s;
  }
  return seqs;
}

//static
void
Alifold::free_aln(char** seqs)
{
  for (char** s=seqs; *s!=NULL; ++s)
    delete[] *s;
  delete[] seqs;
}
