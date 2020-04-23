// $Id:$

#ifndef __INC_PART_ALIGN_H__
#define __INC_PART_ALIGN_H__

#include <vector>
#include <string>
#include "../typedefs.h"
#include "log_value.h"

namespace PARTALIGN
{
  typedef std::pair<std::string, BP> BPSeq;

  template < class V >
  class PartAlign
  {
  public:

  private:
    typedef std::vector<V> VV;
    typedef std::vector<VV> VVV;

  public:
    PartAlign() { }

    void load_sequences(const std::string& x, const std::string& y);
    void load_sequences(const std::string& x, const std::string& y,
                        const BP& px, const BP& py);
    void load_sequences(const BPSeq& x, const BPSeq& y)
    {
      load_sequences(x.first, y.first, x.second, y.second);
    }

    void set_parameters(float alpha, float beta, float gap, float ext)
    {
      alpha_=alpha;
      beta_=beta;
      gap_=gap;
      ext_=ext;
    }
    void set_scoring_matrix(const VVF& sm) { sm_=sm; }
    void set_scoring_matrix(const float sm[10]);

    float score(uint i, uint j) const { return alpha_*bp_score(i,j)+nobp_score(i,j); }
    float bp_score(uint i, uint j) const { return px_l_[i]*py_l_[i] + px_l_[i]*py_l_[i]; }
    float nobp_score(uint i, uint j) const { return sm_[x_[i]][y_[j]]*qx_[i]*qy_[j]; }

    V compute_forward();
    V compute_backward();
    void compute_posterior(VVF& posterior);
    void compute_expected_counts(float& ex_gap, float& ex_ext, VVF& posterior);
    float compute_gradients(const std::vector<bool>& bx, const std::vector<bool>& by,
                            float& g_alpha, float& g_beta, float& g_gap, float& g_ext,
                            VVF& g_s);
  private:
    VI x_, y_;
    VVF sm_;
    VF px_l_, px_r_, qx_, py_l_, py_r_, qy_;
    float alpha_, beta_, gap_, ext_;
    VVV Fm_, Fx_, Fy_;
    VVV Bm_, Bx_, By_;
  };
};

#endif  //  __INC_PART_ALIGN_H__
