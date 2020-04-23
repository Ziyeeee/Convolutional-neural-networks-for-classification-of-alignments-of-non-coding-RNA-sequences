// $Id:$

#ifndef __INC_OPTIMIZER_H__
#define __INC_OPTIMIZER_H__

#include <vector>
#include "../typedefs.h"

namespace PARTALIGN
{
  class Optimizer
  {
  public:
    Optimizer(bool alpha, bool beta, bool gap, bool ext, bool sm,
              float c=1.0, float eta0=0.5, uint t_max=100, float eps=1e-5);
    void set_parameters(float alpha, float beta, float gap, float ext, const VVF& sm);
    void set_parameters(float alpha, float beta, float gap, float ext, const float sm[10]);
    void get_parameters(float& alpha, float& beta, float& gap, float& ext, VVF& sm);
    void get_parameters(float& alpha, float& beta, float& gap, float& ext, float sm[10]);

    template < class V, class Seq >
    void sgd(const std::vector<std::pair<Seq,Seq> >& seq,
             const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln);
    template < class V, class Seq >
    void lbfgs(const std::vector<std::pair<Seq,Seq> >& seq,
               const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln);

  private:    
    float alpha_;
    float beta_;
    float gap_;
    float ext_;
    VVF sm_;
    float c_;
    float eta0_;
    uint t_max_;
    float eps_;
    bool update_alpha_;
    bool update_beta_;
    bool update_gap_;
    bool update_ext_;
    bool update_sm_;
  };
};

#endif  //  __INC_OPTIMIZER_H__
