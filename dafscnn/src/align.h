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

#ifndef __INC_ALIGN_H__
#define __INC_ALIGN_H__

#include <string>
#include <vector>

#include "typedefs.h"
#include "fa.h"
#include "probconsRNA/probcons.h"
#include "contralign/contralign.h"
#if 0
#include "partalign/partalign.h"
#endif

namespace Align
{
  // base class for modeling probability distirubtion of alignments
  class Model
  {
  public:
    Model(float th) : th_(th) { }
    virtual ~Model() { }
    virtual void calculate(const std::string& seq1, const std::string& seq2, MP& mp) = 0;
#if 0
    virtual void calculate(const std::string& seq1, const std::string& seq2,
                           const BP& bp1, const BP& bp2, MP& mp)
    {
      return calculate(seq1, seq2, mp);
    }
#endif
    virtual void calculate(const std::vector<Fasta>& fa, std::vector<std::vector<MP> >& mp);
    float threshold() const { return th_; };

  private:
    float th_;
  };

  class Decoder
  {
  public:
    Decoder() { }
    virtual ~Decoder() { }
    virtual void initialize(const VVF& p) { }
    virtual float decode(const VVF& p, const VVF& q, VU& al) const = 0;
    virtual float decode(const VVF& p, VU& al) const = 0;
  };
}

class ProbCons : public Align::Model, PROBCONS::Probcons
{
public:
  ProbCons(float th);
  ~ProbCons() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
};

class CONTRAlign : public Align::Model, CONTRALIGN::CONTRAlign<float>
{
public:
  CONTRAlign(float th);
  ~CONTRAlign() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
};

#if 0
class PartAlign : public Align::Model, PARTALIGN::PartAlign<LogValue<float> >
{
public:
  PartAlign(float th, const std::string& arg);
  ~PartAlign() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
  void calculate(const std::string& seq1, const std::string& seq2,
                 const BP& bp1, const BP& bp2, MP& mp);
};
#endif

class AUXAlign : public Align::Model
{
public:
  AUXAlign(const std::string& file, float th);
  ~AUXAlign() { }
  void calculate(const std::string& seq1, const std::string& seq2, MP& mp);
  void calculate(const std::vector<Fasta>& fa, std::vector<std::vector<MP> >& mp);

private:
  std::string file_;
};

#endif  //  __INC_ALIGN_H__

// Local Variables:
// mode: C++
// End:
