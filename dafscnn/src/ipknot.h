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

#ifndef __INC_IPKNOT_H__
#define __INC_IPKNOT_H__

#include "fold.h"
#include "ip.h"

class IPknot : public Fold::Decoder
{
public:
  IPknot(const VF& th, int n_th = 1);
  float decode(float w, const VVF& p, const VVF& q, VU& ss);
  float decode(const VVF& p, VU& ss, std::string& str);
  void make_brackets(const VU& ss, std::string& str) const;

private:
  void make_objective(IP& ip, float w, const VVF& p, const VVF& q);
  void make_objective(IP& ip, const VVF& p);
  void make_constraints(IP& ip);
  float solve(IP& ip, VU& ss);
  
private:
  VF th_;
  VF alpha_;
  bool levelwise_;
  bool stacking_constraints_;
  std::vector<VVI> v_;
  std::vector<VVI> w_;
  uint n_th_;
  VU plevel_;
};

#endif  //  __INC_IPKNOT_H__

// Local Variables:
// mode: C++
// End:
