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

#ifndef __INC_GUIDE_TREE_H__
#define __INC_GUIDE_TREE_H__

#include <vector>
#include <iostream>
#include "typedefs.h"

class GuideTree
{
private:
  typedef std::pair<float, std::pair<uint,uint> > node_t;

public:
  GuideTree(uint n);
  void build(const std::vector<float>& dist);

  void print(std::ostream& os, int i) const;
  void print(std::ostream& os) const
  {
    print(os, tree_.size()-1);
    os << std::endl;
  }

  template < class F, class T >
  void process(F& func, T& t) const
  {
    proces(func, t, tree_.size()-1);
  }
  
private:
  template < class F, class T >
  void process(F& func, T& t, uint pos) const
  {
    if (tree_[pos].second.first!=-1u)
    {
      assert(tree_[pos].second.second!=-1u);
      T l, r;
      process(func, l, tree_[pos].second.first);
      process(func, r, tree_[pos].second.second);
      func(t, l, r);
    }
    else
    {
      func(t, pos);
    }
  }  

private:
  uint n_;
  std::vector<node_t> tree_;
};

#endif  // __INC_GUIDE_TREE_H__
