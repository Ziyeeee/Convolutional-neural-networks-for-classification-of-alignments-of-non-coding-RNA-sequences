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

#include "fold.h"
#include <cassert>
#include <cstring>
#include <cmath>
#include <sys/errno.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

namespace Vienna {
extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/alifold.h>
#include <ViennaRNA/aln_util.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/PS_dot.h>
  extern void read_parameter_file(const char fname[]);
};
};

#ifndef FLT_OR_DBL
typedef Vienna::FLT_OR_DBL FLT_OR_DBL;
#endif

extern "C" {
#include "boltzmann_param.h"
};

// some constants
const uint Fold::Decoder::n_support_brackets=4+26;
const char* Fold::Decoder::left_brackets ="([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets=")]}>abcdefghijklmnopqrstuvwxyz";

void
Fold::Model::
calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp)
{
  const uint N=fa.size();
  bp.resize(N);
  for (uint i=0; i!=N; ++i)
    this->calculate(fa[i].seq(), bp[i]);
}

RNAfold::
RNAfold(bool bl, const char* param, float th)
  : Fold::Model(th)
{
  if (bl) copy_boltzmann_parameters();
  if (param) Vienna::read_parameter_file(param);
}

void
RNAfold::
calculate(const std::string& seq, BP& bp)
{
  uint L=seq.size();
  bp.resize(L);
  Vienna::pf_scale = -1;

  // scaling parameters to avoid overflow
  if (1 /*L>1600*/)
  {
    std::string res(L+1, ' ');
    float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &res[0]);
    float sfact = 1.07;
    float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
    Vienna::free_arrays();
  }
#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), NULL);

#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      const float& p = pr[iindx[i+1]-(j+1)];
      if (p>threshold())
        bp[i].push_back(std::make_pair(j, p));
    }
  Vienna::free_pf_arrays();
}

void
RNAfold::
calculate(const std::string& seq, const std::string& str, BP& bp)
{
  assert(seq.size()==str.size());
  std::string p(str);
  std::replace(p.begin(), p.end(), '.', 'x');
  std::replace(p.begin(), p.end(), '?', '.');

  int bk = Vienna::fold_constrained;
  Vienna::fold_constrained = 1;

  uint L=seq.size();
  bp.resize(L);
  Vienna::pf_scale = -1;

  // scaling parameters to avoid overflow
  if (1 /*L>1600*/)
  {
    std::string res(p);
    float min_en = Vienna::fold(const_cast<char*>(seq.c_str()), &res[0]);
    float sfact = 1.07;
    float kT = (Vienna::temperature+273.15)*1.98717/1000.; /* in Kcal */
    Vienna::pf_scale = exp(-(sfact*min_en)/kT/seq.size());
    Vienna::free_arrays();
  }

#ifndef HAVE_VIENNA20
  Vienna::init_pf_fold(L);
#endif
  Vienna::pf_fold(const_cast<char*>(seq.c_str()), &p[0]);

#ifdef HAVE_VIENNA20
  FLT_OR_DBL* pr = Vienna::export_bppm();
  int* iindx = Vienna::get_iindx(seq.size());
#else
  FLT_OR_DBL* pr = Vienna::pr;
  int* iindx = Vienna::iindx;
#endif
  for (uint i=0; i!=L-1; ++i)
    for (uint j=i+1; j!=L; ++j)
    {
      const float& p = pr[iindx[i+1]-(j+1)];
      if (p>threshold())
        bp[i].push_back(std::make_pair(j, p));
    }
  Vienna::free_pf_arrays();
  Vienna::fold_constrained = bk;
}

CONTRAfold::
CONTRAfold(float th)
  : Fold::Model(th), CONTRAFOLD::CONTRAfold<float>()
{
}

void
CONTRAfold::
calculate(const std::string& seq, BP& bp)
{
  bp.resize(seq.size());
  std::vector<float> posterior;
  ComputePosterior(seq, posterior);
  for (uint i=0, k=0; i!=seq.size()+1; ++i)
  {
    for (uint j=i; j!=seq.size()+1; ++j, ++k)
    {
      if (i!=0 && posterior[k]>threshold())
        bp[i-1].push_back(std::make_pair(j-1, posterior[k]));
    }
  }
}

void
CONTRAfold::
calculate(const std::string& seq, const std::string& str, BP& bp)
{
  bp.resize(seq.size());
  std::vector<float> posterior;
  SetConstraint(str);
  ComputePosterior(seq, posterior);
  for (uint i=0, k=0; i!=seq.size()+1; ++i)
  {
    for (uint j=i; j!=seq.size()+1; ++j, ++k)
    {
      if (i!=0 && posterior[k]>threshold())
        bp[i-1].push_back(std::make_pair(j-1, posterior[k]));
    }
  }
}

AUXFold::
AUXFold(const std::string& file, float th)
  : Fold::Model(th),
    file_(file)
{
}

void
AUXFold::
calculate(const std::string& seq, BP& bp)
{
  throw "not supported";
}

void
AUXFold::
calculate(const std::string& seq, const std::string& str, BP& bp)
{
  throw "not supported";
}

static
void
load_bp(std::istream& is, std::vector<BP>& bp)
{
  std::string s, t;
  uint x, i, j;
  float p;
  while (std::getline(is, s))
  {
    std::istringstream ss(s);
    if (s[0]=='>')
    {
      ss >> t >> x;
      assert(x-1<bp.size());
    }
    else
    {
      ss >> i;
      if (i-1>=bp[x-1].size()) bp[x-1].resize(i);
      while (ss >> t)
      {
        if (sscanf(t.c_str(), "%u:%f", &j, &p)==2)
        {
          assert(i<j);
          bp[x-1][i-1].push_back(std::make_pair(j-1, p));
        }
      }
    }
  }
}

void
AUXFold::
calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp)
{
  const uint N=fa.size();
  bp.resize(N);
  std::ifstream is(file_.c_str());
  if (is.is_open())
    load_bp(is, bp);
  else
    throw strerror(errno);
#ifndef NDEBUG
  for (uint i=0; i!=fa.size(); ++i)
  {
    assert(bp[i].size()==fa[i].size());
  }
#endif
}
  
