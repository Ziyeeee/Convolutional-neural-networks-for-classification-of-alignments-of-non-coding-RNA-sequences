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

#include "ipknot.h"
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include "ip.h"

static
void
decompose_plevel(const VU& ss, VU& plevel);

static
std::string
make_brackets(const VU& ss, const VU& plevel);

IPknot::
IPknot(const VF& th, int n_th /*=1*/)
  : th_(th),
    alpha_(th_.size(), 1.0),
    levelwise_(true),
    stacking_constraints_(true),
    n_th_(n_th)
{
}

float
IPknot::
decode(float w, const VVF& p, const VVF& q, VU& ss)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, w, p, q);
  make_constraints(ip);
  return solve(ip, ss);
}

float
IPknot::
decode(const VVF& p, VU& ss, std::string& str)
{
  IP ip(IP::MAX, n_th_);
  make_objective(ip, p);
  make_constraints(ip);
  float s=solve(ip, ss);
  str=::make_brackets(ss, plevel_);
  return s;
}

void
IPknot::
make_brackets(const VU& ss, std::string& str) const
{
  VU plevel;
  decompose_plevel(ss, plevel);
  str=::make_brackets(ss, plevel);
}

void
IPknot::
make_objective(IP& ip, float w, const VVF& p, const VVF& q) 
{
  const uint L=p.size();
  const uint P=th_.size();
  
  v_.clear(); v_.resize(P, VVI(L, VI(L, -1)));
  w_.clear(); w_.resize(P, VVI(L));
    
  // make objective variables with their weights
  for (uint j=1; j!=L; ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      for (uint lv=0; lv!=P; ++lv)
      {
        float s=w*(p[i][j]-th_[lv])-q[i][j];
        if (s>0.0)
        {
          v_[lv][i][j] = ip.make_variable(s*alpha_[lv]);
          w_[lv][i].push_back(j);
        }
      }
    }
  }
  ip.update();
}

void
IPknot::
make_objective(IP& ip, const VVF& p) 
{
  const uint L=p.size();
  const uint P=th_.size();

  v_.clear(); v_.resize(P, VVI(L, VI(L, -1)));
  w_.clear(); w_.resize(P, VVI(L));
    
  // make objective variables with their weights
  for (uint j=1; j!=L; ++j)
  {
    for (uint i=j-1; i!=-1u; --i)
    {
      for (uint lv=0; lv!=P; ++lv)
      {
        float s=p[i][j]-th_[lv];
        if (s>0.0)
        {
          v_[lv][i][j] = ip.make_variable(s*alpha_[lv]);
          w_[lv][i].push_back(j);
        }
      }
    }
  }
  ip.update();
}

void
IPknot::
make_constraints(IP& ip)
{
  const uint L=v_[0].size();
  const uint P=th_.size();

  // constraint 1: each s_i is paired with at most one base
  for (uint i=0; i!=L; ++i)
  {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint lv=0; lv!=P; ++lv)
    {
      for (uint j=0; j<i; ++j)
        if (v_[lv][j][i]>=0)
          ip.add_constraint(row, v_[lv][j][i], 1);
      for (uint j=i+1; j<L; ++j)
        if (v_[lv][i][j]>=0)
          ip.add_constraint(row, v_[lv][i][j], 1);
    }
  }

  if (levelwise_)
  {
    // constraint 2: disallow pseudoknots in x[lv]
    for (uint lv=0; lv!=P; ++lv)
      for (uint i=0; i<w_[lv].size(); ++i)
        for (uint p=0; p<w_[lv][i].size(); ++p)
        {
          uint j=w_[lv][i][p];
          for (uint k=i+1; k<j; ++k)
            for (uint q=0; q<w_[lv][k].size(); ++q)
            {
              uint l=w_[lv][k][q];
              if (j<l)
              {
                int row = ip.make_constraint(IP::UP, 0, 1);
                ip.add_constraint(row, v_[lv][i][j], 1);
                ip.add_constraint(row, v_[lv][k][l], 1);
              }
            }
        }

    // constraint 3: any x[t]_kl must be pseudoknotted with x[u]_ij for t>u
    for (uint lv=1; lv!=P; ++lv)
      for (uint k=0; k<w_[lv].size(); ++k)
        for (uint q=0; q<w_[lv][k].size(); ++q)
        {
          uint l=w_[lv][k][q];
          for (uint plv=0; plv!=lv; ++plv)
          {
            int row = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(row, v_[lv][k][l], -1);
            for (uint i=0; i<k; ++i)
              for (uint p=0; p<w_[plv][i].size(); ++p)
              {
                uint j=w_[plv][i][p];
                if (k<j && j<l)
                  ip.add_constraint(row, v_[plv][i][j], 1);
              }
            for (uint i=k+1; i<l; ++i)
              for (uint p=0; p<w_[plv][i].size(); ++p)
              {
                uint j=w_[plv][i][p];
                if (l<j)
                  ip.add_constraint(row, v_[plv][i][j], 1);
              }
          }
        }
  }

  if (stacking_constraints_)
  {
    for (uint lv=0; lv!=P; ++lv)
    {
      // upstream
      for (uint i=0; i<L; ++i)
      {
        int row = ip.make_constraint(IP::LO, 0, 0);
        for (uint j=0; j<i; ++j)
          if (v_[lv][j][i]>=0)
            ip.add_constraint(row, v_[lv][j][i], -1);
        if (i>0)
          for (uint j=0; j<i-1; ++j)
            if (v_[lv][j][i-1]>=0)
              ip.add_constraint(row, v_[lv][j][i-1], 1);
        if (i+1<L)
          for (uint j=0; j<i+1; ++j)
            if (v_[lv][j][i+1]>=0)
              ip.add_constraint(row, v_[lv][j][i+1], 1);
      }

      // downstream
      for (uint i=0; i<L; ++i)
      {
        int row = ip.make_constraint(IP::LO, 0, 0);
        for (uint j=i+1; j<L; ++j)
          if (v_[lv][i][j]>=0)
            ip.add_constraint(row, v_[lv][i][j], -1);
        if (i>0)
          for (uint j=i; j<L; ++j)
            if (v_[lv][i-1][j]>=0)
              ip.add_constraint(row, v_[lv][i-1][j], 1);
        if (i+1<L)
          for (uint j=i+2; j<L; ++j)
            if (v_[lv][i+1][j]>=0)
              ip.add_constraint(row, v_[lv][i+1][j], 1);
      }
    }
  }
}

float
IPknot::
solve(IP& ip, VU& ss)
{
  const uint L=v_[0].size();
  const uint P=th_.size();

  // execute optimization
  float s=ip.solve();

  // build the result
  ss.resize(L);
  std::fill(ss.begin(), ss.end(), -1u);
  plevel_.resize(L);
  std::fill(plevel_.begin(), plevel_.end(), -1u);
  for (uint lv=0; lv!=P; ++lv)
  {
    for (uint i=0; i<L; ++i)
      for (uint j=i+1; j<L; ++j)
        if (v_[lv][i][j]>=0 && ip.get_value(v_[lv][i][j])>0.5)
        {
          ss[i]=j; //ss[j]=i;
          plevel_[i]=plevel_[j]=lv;
        }
  }

  if (!levelwise_) decompose_plevel(ss, plevel_);

  return s;
}

struct cmp_by_degree : public std::less<int>
{
  cmp_by_degree(const VVU& g) : g_(g) {}
  bool operator()(int x, int y) const { return g_[y].size()<g_[x].size(); }
  const VVU& g_;
};

struct cmp_by_count : public std::less<int>
{
  cmp_by_count(const VU& count) : count_(count) { }
  bool operator()(int x, int y) const { return count_[y]<count_[x]; }
  const VU& count_;
};

static void
decompose_plevel(const VU& ss, VU& plevel)
{
  // resolve the symbol of brackets by the graph coloring problem
  uint L=ss.size();
    
  // make an adjacent graph, in which pseudoknotted base-pairs are connected.
  VVU g(L);
  for (uint i=0; i!=L; ++i)
  {
    if (ss[i]==-1u || ss[i]<=i) continue;
    uint j=ss[i];
    for (uint k=i+1; k!=L; ++k)
    {
      if (ss[k]==-1u || ss[k]<=k) continue;
      uint l=ss[k];
      if (k<j && j<l)
      {
        g[i].push_back(k);
        g[k].push_back(i);
      }
    }
  }
  // vertices are indexed by the position of the left base
  VU v;
  for (uint i=0; i!=ss.size(); ++i)
    if (ss[i]!=-1u && i<ss[i]) 
      v.push_back(i);
  // sort vertices by degree
  std::sort(v.begin(), v.end(), cmp_by_degree(g));

  // determine colors
  VU c(L, -1u);
  uint max_color=0;
  for (uint i=0; i!=v.size(); ++i)
  {
    // find the smallest color that is unused
    VU used;
    for (uint j=0; j!=g[v[i]].size(); ++j)
      if (c[g[v[i]][j]]!=-1u) used.push_back(c[g[v[i]][j]]);
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    uint j=0;
    for (j=0; j!=used.size(); ++j)
      if (used[j]!=j) break;
    c[v[i]]=j;
    max_color=std::max(max_color, j);
  }

  // renumber colors in decentant order by the number of base-pairs for each color
  VU count(max_color+1, 0);
  for (uint i=0; i!=c.size(); ++i)
    if (c[i]!=-1u) count[c[i]]++;
  VU idx(count.size());
  for (uint i=0; i!=idx.size(); ++i) idx[i]=i;
  sort(idx.begin(), idx.end(), cmp_by_count(count));
  VU rev(idx.size());
  for (uint i=0; i!=rev.size(); ++i) rev[idx[i]]=i;
  plevel.resize(L);
  for (uint i=0; i!=c.size(); ++i)
    plevel[i]= c[i]!=-1u ? rev[c[i]] : -1u;
}

static
std::string
make_brackets(const VU& ss, const VU& plevel) 
{
  std::string r(ss.size(), '.');
  for (uint i=0; i!=ss.size(); ++i)
  {
    if (ss[i]!=-1u && i<ss[i])
    {
      uint j=ss[i];
      assert(plevel[i]!=-1u);
      if (plevel[i]<Fold::Decoder::n_support_brackets)
      {
        r[i]=Fold::Decoder::left_brackets[plevel[i]];
        r[j]=Fold::Decoder::right_brackets[plevel[i]];
      }
    }
  }
  return r;
}
