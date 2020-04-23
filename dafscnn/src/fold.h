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

#ifndef __INC_FOLD_H__
#define __INC_FOLD_H__

#include <string>
#include <vector>

#include "typedefs.h"
#include "contrafold/contrafold.h"
#include "fa.h"

namespace Fold 
{
  // base class for modeling probability distribution of secondary structures
  class Model
  {
  public:
    Model(float th) : th_(th) { }
    virtual ~Model() { }
    virtual void calculate(const std::string& seq, BP& bp) = 0;
    virtual void calculate(const std::string& seq, const std::string& str, BP& bp) = 0;
    virtual void calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp);
    float threshold() const { return th_; }

  private:
    float th_;
  };

  class Decoder
  {
  public:
    static const uint n_support_brackets;
    static const char* left_brackets;
    static const char* right_brackets;

  public:
    Decoder() { }
    virtual ~Decoder() { }
    virtual float decode(float w, const VVF& p, const VVF& q, VU& ss) = 0;
    virtual float decode(const VVF& p, VU& ss, std::string& str) = 0;
    virtual void make_brackets(const VU& ss, std::string& str) const = 0;
  };
}

class RNAfold : public Fold::Model
{
public:
  RNAfold(bool bl, const char* param, float th);
  void calculate(const std::string& seq, BP& bp);
  void calculate(const std::string& seq, const std::string& str, BP& bp);
};

class CONTRAfold : public Fold::Model, CONTRAFOLD::CONTRAfold<float>
{
public:
  CONTRAfold(float th);
  void calculate(const std::string& seq, BP& bp);
  void calculate(const std::string& seq, const std::string& str, BP& bp);
};

class AUXFold : public Fold::Model
{
public:
  AUXFold(const std::string& file, float th);
  ~AUXFold() { }
  void calculate(const std::string& seq, BP& bp);
  void calculate(const std::string& seq, const std::string& str, BP& bp);
  void calculate(const std::vector<Fasta>& fa, std::vector<BP>& bp);

private:
  std::string file_;
};

#endif  //  __INC_FOLD_H__

// Local Variables:
// mode: C++
// End:
