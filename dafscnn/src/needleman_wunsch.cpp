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

#include "needleman_wunsch.h"
#include <cassert>
#include <algorithm>

float
NeedlemanWunsch::
decode(const VVF& p, const VVF& q, VU& al) const
{
  const uint L1 = p.size();
  const uint L2 = p[0].size();

  VVF dp(L1+1, VF(L2+1, 0.0));
  VVC tr(L1+1, VC(L2+1, ' '));
  for (uint i=1; i!=L1+1; ++i) tr[i][0] = 'X';
  for (uint k=1; k!=L2+1; ++k) tr[0][k] = 'Y';

  for (uint i=1; i!=L1+1; ++i)
  {
    for (uint k=1; k!=L2+1; ++k)
    {
      float v = dp[i-1][k-1]+p[i-1][k-1]-th_+q[i-1][k-1];
      char t = 'M';
      if (v<dp[i-1][k])
      {
        v = dp[i-1][k];
        t = 'X';
      }
      if (v<dp[i][k-1])
      {
        v = dp[i][k-1];
        t = 'Y';
      }
      dp[i][k] = v;
      tr[i][k] = t;
    }
  }


  // traceback
  std::string rpath;
  int i = L1, k = L2;
  while (i>0 || k>0)
  {
    rpath.push_back(tr[i][k]);
    switch (tr[i][k])
    {
      case 'M': --i; --k; break;
      case 'X': --i; break;
      case 'Y': --k; break;
      default: assert(!"unreachable"); break;
    }
  }
  std::string vpath(rpath.size(), ' ');
  std::reverse_copy(rpath.begin(), rpath.end(), vpath.begin());

  // decode the resultant alignment
  al.resize(L1, -1u);
  for (uint i=0, k=0, p=0; p!=vpath.size(); ++p)
  {
    switch (vpath[p])
    {
      case 'M':
        assert(i<L1); assert(k<L2);
        al[i++]=k++;
        break;
      case 'X':
        assert(i<L1); assert(k<=L2);
        al[i++]=-1u;
        break;
      case 'Y':
        assert(i<=L1); assert(k<L2);
        k++;
        break;
      default: break;
    }
  }

  return dp[L1][L2];
}

float
NeedlemanWunsch::
decode(const VVF& p, VU& al) const
{
  const uint L1 = p.size();
  const uint L2 = p[0].size();

  VVF dp(L1+1, VF(L2+1, 0.0));
  VVC tr(L1+1, VC(L2+1, ' '));
  for (uint i=1; i!=L1+1; ++i) tr[i][0] = 'X';
  for (uint k=1; k!=L2+1; ++k) tr[0][k] = 'Y';

  for (uint i=1; i!=L1+1; ++i)
  {
    for (uint k=1; k!=L2+1; ++k)
    {
      float v = dp[i-1][k-1]+p[i-1][k-1]-th_;
      char t = 'M';
      if (v<dp[i-1][k])
      {
        v = dp[i-1][k];
        t = 'X';
      }
      if (v<dp[i][k-1])
      {
        v = dp[i][k-1];
        t = 'Y';
      }
      dp[i][k] = v;
      tr[i][k] = t;
    }
  }


  // traceback
  std::string rpath;
  int i = L1, k = L2;
  while (i>0 || k>0)
  {
    rpath.push_back(tr[i][k]);
    switch (tr[i][k])
    {
      case 'M': --i; --k; break;
      case 'X': --i; break;
      case 'Y': --k; break;
      default: assert(!"unreachable"); break;
    }
  }
  std::string vpath(rpath.size(), ' ');
  std::reverse_copy(rpath.begin(), rpath.end(), vpath.begin());

  // decode the resultant alignment
  al.resize(L1, -1u);
  for (uint i=0, k=0, p=0; p!=vpath.size(); ++p)
  {
    switch (vpath[p])
    {
      case 'M':
        assert(i<L1); assert(k<L2);
        al[i++]=k++;
        break;
      case 'X':
        assert(i<L1); assert(k<=L2);
        al[i++]=-1u;
        break;
      case 'Y':
        assert(i<=L1); assert(k<L2);
        k++;
        break;
      default: break;
    }
  }

  return dp[L1][L2];
}

void
SparseNeedlemanWunsch::
initialize(const VVF& p)
{
  const uint L1 = p.size();
  const uint L2 = p[0].size();
  env_.resize(L1+1);
  std::fill(env_.begin(), env_.end(), std::make_pair(0, 0));

  for (uint i=1; i!=L1+1; ++i)
  {
    // find the first alignable point
    for (uint k=1; k!=L2+1; ++k)
    {
      if (p[i-1][k-1]-th_>=0.0)
      {
        env_[i-1].first = std::min(env_[i-1].first, k-1);
        env_[i].first = k;
        break;
      }
    }
    // no alignable point
    if (env_[i].first==0)          
    {
      env_[i].first = env_[i-1].first;
      env_[i].second = env_[i-1].second;
      continue;
    }

    // find the last alignable point
    for (uint k=L2; k!=0; --k)
    {
      if (p[i-1][k-1]-th_>=0.0)
      {
        env_[i-1].second = std::max(env_[i-1].second, k-1);
        env_[i].second = k;
        break;
      }
    }
  }
  assert(env_[0].first==0);
  env_[L1].second=L2;

  // force the envelope to be monotonic
  for (uint i=L1, v=L2; i!=0; --i)
    env_[i].first = v = std::min(v, env_[i].first);
  for (uint i=0, v=0; i!=L1+1; ++i)
    env_[i].second = v = std::max(v, env_[i].second);

  // ensure the connectivity
  for (uint i=1; i!=L1+1; ++i)
  {
    if (env_[i-1].second<env_[i].first)
      env_[i].first = env_[i-1].second;
  }
}

float
SparseNeedlemanWunsch::
decode(const VVF& p, const VVF& q, VU& al) const
{
  const uint L1 = p.size();
  const uint L2 = p[0].size();

  VVF dp(L1+1, VF(L2+1, 0.0));
  VVC tr(L1+1, VC(L2+1, ' '));
  for (uint i=1; i!=L1+1; ++i) tr[i][0] = 'X';
  for (uint k=1; k!=L2+1; ++k) tr[0][k] = 'Y';
  for (uint i=1; i!=L1+1; ++i)
  {
    for (uint k=env_[i].first; k<=env_[i].second; ++k)
    {
      if (k==0) continue;
      float v = dp[i-1][k-1]+p[i-1][k-1]-th_+q[i-1][k-1];
      char t = 'M';
      if (v<dp[i-1][k])
      {
        v = dp[i-1][k];
        t = 'X';
      }
      if (v<dp[i][k-1])
      {
        v = dp[i][k-1];
        t = 'Y';
      }
      dp[i][k] = v;
      tr[i][k] = t;
    }
  }

  // traceback
  std::string rpath;
  int i = L1, k = L2;
  while (i>0 || k>0)
  {
    rpath.push_back(tr[i][k]);
    switch (tr[i][k])
    {
      case 'M': --i; --k; break;
      case 'X': --i; break;
      case 'Y': --k; break;
      default: assert(!"unreachable"); break;
    }
  }
  std::string vpath(rpath.size(), ' ');
  std::reverse_copy(rpath.begin(), rpath.end(), vpath.begin());

  // decode the resultant alignment
  al.resize(L1, -1u);
  for (uint i=0, k=0, p=0; p!=vpath.size(); ++p)
  {
    switch (vpath[p])
    {
      case 'M':
        assert(i<L1); assert(k<L2);
        al[i++]=k++;
        break;
      case 'X':
        assert(i<L1); assert(k<=L2);
        al[i++]=-1u;
        break;
      case 'Y':
        assert(i<=L1); assert(k<L2);
        k++;
        break;
      default: break;
    }
  }

  return dp[L1][L2];
}

float
SparseNeedlemanWunsch::
decode(const VVF& p, VU& al) const
{
  const uint L1 = p.size();
  const uint L2 = p[0].size();

  VVF dp(L1+1, VF(L2+1, 0.0));
  VVC tr(L1+1, VC(L2+1, ' '));
  for (uint i=1; i!=L1+1; ++i) tr[i][0] = 'X';
  for (uint k=1; k!=L2+1; ++k) tr[0][k] = 'Y';
  for (uint i=1; i!=L1+1; ++i)
  {
    for (uint k=env_[i].first; k<=env_[i].second; ++k)
    {
      if (k==0) continue;
      float v = dp[i-1][k-1]+p[i-1][k-1]-th_;
      char t = 'M';
      if (v<dp[i-1][k])
      {
        v = dp[i-1][k];
        t = 'X';
      }
      if (v<dp[i][k-1])
      {
        v = dp[i][k-1];
        t = 'Y';
      }
      dp[i][k] = v;
      tr[i][k] = t;
    }
  }

  // traceback
  std::string rpath;
  int i = L1, k = L2;
  while (i>0 || k>0)
  {
    rpath.push_back(tr[i][k]);
    switch (tr[i][k])
    {
      case 'M': --i; --k; break;
      case 'X': --i; break;
      case 'Y': --k; break;
      default: assert(!"unreachable"); break;
    }
  }
  std::string vpath(rpath.size(), ' ');
  std::reverse_copy(rpath.begin(), rpath.end(), vpath.begin());

  // decode the resultant alignment
  al.resize(L1, -1u);
  for (uint i=0, k=0, p=0; p!=vpath.size(); ++p)
  {
    switch (vpath[p])
    {
      case 'M':
        assert(i<L1); assert(k<L2);
        al[i++]=k++;
        break;
      case 'X':
        assert(i<L1); assert(k<=L2);
        al[i++]=-1u;
        break;
      case 'Y':
        assert(i<=L1); assert(k<L2);
        k++;
        break;
      default: break;
    }
  }

  return dp[L1][L2];
}
