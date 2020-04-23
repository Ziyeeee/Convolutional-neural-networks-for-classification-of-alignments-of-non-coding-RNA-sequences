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

#include "config.h"
#include <cassert>
#include <iostream>
#include <queue>
#include "tree.h"

GuideTree::
GuideTree(uint n)
  : n_(n),
    tree_(2*n-1, std::make_pair(0.0, std::make_pair(-1u, -1u)))
{
}

void
GuideTree::
build(const std::vector<float>& dist)
{
  uint n=n_;
  std::vector<std::vector<float> > d(n, std::vector<float>(n, 0.0));
  std::vector<uint> idx(2*n-1, -1u);
  for (uint i=0; i!=n; ++i) idx[i]=i;

  std::priority_queue<node_t> pq;
  for (uint i=0, k=0; i!=n-1; ++i)
  {
    for (uint j=i+1; j!=n; ++j, ++k)
    {
      d[i][j] = d[j][i] = dist[k];
      pq.push(std::make_pair(dist[k], std::make_pair(i,j)));
    }
  }

  while (!pq.empty())
  {
    node_t t=pq.top(); pq.pop();
    if (idx[t.second.first]!=-1u && idx[t.second.second]!=-1u)
    {
      assert(n<tree_.size());
      uint l = idx[t.second.first];
      uint r = idx[t.second.second];
      idx[t.second.first] = idx[t.second.second] = -1u;
      for (uint i=0; i!=n; ++i)
      {
        if (idx[i]!=-1u)
        {
          uint ii=idx[i];
          d[ii][l] = d[l][ii] = (d[ii][l]+d[ii][r])*t.first/2;
          pq.push(std::make_pair(d[ii][l], std::make_pair(i,n)));
        }
      }
      tree_[n] = t;
      idx[n++] = l;
    }
  }
  assert(n==tree_.size());
}

void
GuideTree::
print(std::ostream& os, int i) const
{
  if (tree_[i].second.first!=-1u)
  {
    os << "[ " << tree_[i].first << " ";
    print(os, tree_[i].second.first);
    os << " ";
    print(os, tree_[i].second.second);
    os << " ]";
  }
  else
  {
    assert(tree_[i].second.second==-1u);
    os << i;
  }
}
