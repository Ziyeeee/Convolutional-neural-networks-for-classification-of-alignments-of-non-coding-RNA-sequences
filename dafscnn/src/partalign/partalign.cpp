// $Id:$

#include "config.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "partalign.h"

using namespace PARTALIGN;

#define FOREACH(itr, i, v) for (itr i=(v).begin(); i!=(v).end(); ++i)

enum { A=0, C=1, G=2, U=3, N=4 };

template < class V >
void
PartAlign<V>::
set_scoring_matrix(const float sm[10])
{
  sm_.clear();
  sm_.resize(4, VF(4, 0.0));
  for (uint i=0; i<4; ++i)
  {
    sm_[i][i] = *sm++;
    for (uint j=i+1; j<4; ++j)
      sm_[i][j] = sm_[j][i] = *sm++;
  }
}

template < class V >
V
PartAlign<V>::
compute_forward()
{
  const uint L1 = x_.size();
  const uint L2 = y_.size();
  V beta_gap = ExpOf(beta_*gap_);
  V beta_ext = ExpOf(beta_*ext_);

  Fm_.clear(); Fm_.resize(L1+1, VV(L2+1, 0.0));
  Fx_.clear(); Fx_.resize(L1+1, VV(L2+1, 0.0));
  Fy_.clear(); Fy_.resize(L1+1, VV(L2+1, 0.0));

  Fm_[0][0]=1.0;
  for (uint i=1; i!=L1+1; ++i)
  {
    Fx_[i][0] += beta_gap*Fm_[i-1][0];
    Fx_[i][0] += beta_ext*Fx_[i-1][0];
  }
  for (uint j=1; j!=L2+1; ++j)
  {
    Fy_[0][j] += beta_gap*Fm_[0][j-1];
    Fy_[0][j] += beta_ext*Fy_[0][j-1];
  }
  for (uint i=1; i!=L1+1; ++i)
    for (uint j=1; j!=L2+1; ++j)
    {
      V beta_s = ExpOf(beta_*score(i-1,j-1));
      Fm_[i][j] += beta_s*Fm_[i-1][j-1];
      Fm_[i][j] += beta_s*Fx_[i-1][j-1];
      Fm_[i][j] += beta_s*Fy_[i-1][j-1];

      Fx_[i][j] += beta_gap*Fm_[i-1][j];
      Fx_[i][j] += beta_ext*Fx_[i-1][j];

      Fy_[i][j] += beta_gap*Fm_[i][j-1];
      Fy_[i][j] += beta_gap*Fx_[i][j-1];
      Fy_[i][j] += beta_ext*Fy_[i][j-1];
    }

  return Fm_[L1][L2]+Fx_[L1][L2]+Fy_[L1][L2];
}

template < class V >
V
PartAlign<V>::
compute_backward()
{
  const uint L1 = x_.size();
  const uint L2 = y_.size();
  V beta_gap = ExpOf(beta_*gap_);
  V beta_ext = ExpOf(beta_*ext_);

  Bm_.clear(); Bm_.resize(L1+1, VV(L2+1, 0.0));
  Bx_.clear(); Bx_.resize(L1+1, VV(L2+1, 0.0));
  By_.clear(); By_.resize(L1+1, VV(L2+1, 0.0));

  Bm_[L1][L2]=Bx_[L1][L2]=By_[L1][L2]=1.0;
  for (uint i=L1; i!=0; --i)
    for (uint j=L2; j!=0; --j)
    {
      V beta_s = ExpOf(beta_*score(i-1,j-1));
      Bm_[i-1][j-1] += beta_s*Bm_[i][j];
      Bx_[i-1][j-1] += beta_s*Bm_[i][j];
      By_[i-1][j-1] += beta_s*Bm_[i][j];

      Bm_[i-1][j] += beta_gap*Bx_[i][j];
      Bx_[i-1][j] += beta_ext*Bx_[i][j];

      Bm_[i][j-1] += beta_gap*By_[i][j];
      Bx_[i][j-1] += beta_gap*By_[i][j];
      By_[i][j-1] += beta_ext*By_[i][j];
    }
  for (uint j=1; j!=L2+1; ++j)
  {
    Bm_[0][j-1] += beta_gap*By_[0][j];
    By_[0][j-1] += beta_ext*By_[0][j];
  }
  for (uint i=1; i!=L1+1; ++i)
  {
    Bm_[i-1][0] += beta_gap*Bx_[i][0];
    Bx_[i-1][0] += beta_ext*Bx_[i][0];
  }

  return Bm_[0][0];
}

template < class V >
void
PartAlign<V>::
compute_posterior(VVF& posterior)
{
  const uint L1 = x_.size();
  const uint L2 = y_.size();
  V Z = Fm_[L1][L2]+Fx_[L1][L2]+Fy_[L1][L2];

  posterior.clear();
  posterior.resize(L1, VF(L2, 0.0));

  for (uint i=1; i!=L1+1; ++i)
    for (uint j=1; j!=L2+1; ++j)
    {
      posterior[i-1][j-1] = Fm_[i][j]*Bm_[i][j] / Z;
    }
}

template < class V >
void
PartAlign<V>::
compute_expected_counts(float& ex_gap, float& ex_ext, VVF& posterior)
{
  const uint L1 = x_.size();
  const uint L2 = y_.size();
  V beta_gap = ExpOf(beta_*gap_);
  V beta_ext = ExpOf(beta_*ext_);
  V Z = Fm_[L1][L2]+Fx_[L1][L2]+Fy_[L1][L2];
  ex_gap=ex_ext=0.0;
  posterior.clear();
  posterior.resize(L1, VF(L2, 0.0));
  
  for (uint i=1; i!=L1+1; ++i)
  {
    ex_gap += Fm_[i-1][0]*beta_gap*Bx_[i][0] / Z;
    ex_ext += Fx_[i-1][0]*beta_ext*Bx_[i][0] / Z;
  }
  for (uint j=1; j!=L2+1; ++j)
  {
    ex_gap += Fm_[0][j-1]*beta_gap*By_[0][j] / Z;
    ex_ext += Fy_[0][j-1]*beta_ext*By_[0][j] / Z;
  }
  for (uint i=1; i!=L1+1; ++i)
    for (uint j=1; j!=L2+1; ++j)
    {
      posterior[i-1][j-1] = Fm_[i][j]*Bm_[i][j] / Z;

      ex_gap += Fm_[i-1][j]*beta_gap*Bx_[i][j] / Z;
      ex_ext += Fx_[i-1][j]*beta_ext*Bx_[i][j] / Z;

      ex_gap += Fm_[i][j-1]*beta_gap*By_[i][j] / Z;
      ex_gap += Fx_[i][j-1]*beta_gap*By_[i][j] / Z;
      ex_ext += Fy_[i][j-1]*beta_ext*By_[i][j] / Z;
    }
}

static
void
transform(VI& xx, const std::string& x)
{
  xx.resize(x.size());
  for (uint i=0; i!=xx.size(); ++i)
  {
    switch (x[i])
    {
      case 'a': case 'A':
        xx[i]=A;
        break;        
      case 'c': case 'C':
        xx[i]=C;
        break;        
      case 'g': case 'G':
        xx[i]=G;
        break;        
      case 't': case 'T': case 'u': case 'U':
        xx[i]=U;
        break;
      default:
        break;
    }
  }
}

static
void
init_bp(VF& pl, VF& pr, VF& q, uint sz)
{
  pl.clear(); pl.resize(sz, 0.0);
  pr.clear(); pr.resize(sz, 0.0);
  q.clear(); q.resize(sz, 1.0);
}

static
void
calc_bp(VF& pl, VF& pr, VF& q, const BP& px)
{
  for (uint i=0; i!=px.size(); ++i)
    FOREACH(SV::const_iterator, it, px[i])
    {
      pl[i] += it->second;
      q[i] -= it->second;
      pl[it->first] += it->second;
      q[it->first] -= it->second;
    }
}

template < class V >
void
PartAlign<V>::
load_sequences(const std::string& x, const std::string& y)
{
  transform(x_, x);
  init_bp(px_l_, px_r_, qx_, x.size());
  transform(y_, y);
  init_bp(py_l_, py_r_, qy_, y.size());
  alpha_=0.0;
}

template < class V >
void
PartAlign<V>::
load_sequences(const std::string& x, const std::string& y, const BP& px, const BP& py)
{
  transform(x_, x);
  init_bp(px_l_, px_r_, qx_, x.size());
  calc_bp(px_l_, px_r_, qx_, px);

  transform(y_, y);
  init_bp(py_l_, py_r_, qy_, y.size());
  calc_bp(py_l_, py_r_, qy_, py);
}

template < class V >
float
PartAlign<V>::
compute_gradients(const std::vector<bool>& bx, const std::vector<bool>& by,
                  float& g_alpha, float& g_beta, float& g_gap, float& g_ext,
                  VVF& g_sm)
{
  assert(bx.size()==by.size());
  const uint L1 = x_.size();
  const uint L2 = y_.size();
  enum { M, IX, IY };

  uint s = M;
  g_sm.clear(); g_sm.resize(4, VF(4, 0.0));
  g_gap=g_ext=g_alpha=g_beta=0.0;
  for (uint p=0, i=0, j=0; p!=bx.size(); ++p)
  {
    if (bx[p] && by[p])         // match
    {
      const float nobp_score = this->nobp_score(i,j);
      const float bp_score = this->bp_score(i,j);
      g_sm[x_[i]][y_[j]] += beta_*qx_[i]*qy_[j];
      g_alpha += beta_*bp_score;
      g_beta += alpha_*bp_score+nobp_score;
      s = M;
      i++; j++;
    }
    else if (bx[p])             // insert to x
    {
      if (s==IX)
      {
        g_beta += ext_;
        g_ext += beta_;
      }
      else
      {
        g_beta += gap_;
        g_gap += beta_;
      }
      s = IX;
      i++;
    }
    else if (by[p])             // insert to y
    {
      if (s==IY)
      {
        g_beta += ext_;
        g_ext += beta_;
      }
      else
      {
        g_beta += gap_;
        g_gap += beta_;
      }
      s = IY;
      j++;
    }
    else
      assert(!"unreachable");
  }
  float v = beta_*g_beta;
  
  V Z = compute_forward();
  compute_backward();
  float ex_gap, ex_ext;
  VVF posterior;
  compute_expected_counts(ex_gap, ex_ext, posterior);

  for (uint i=0; i!=L1; ++i)
    for (uint j=0; j!=L2; ++j)
    {
      const float ex = posterior[i][j];
      const float nobp_score = this->nobp_score(i,j);
      const float bp_score = this->bp_score(i,j);
      g_sm[x_[i]][y_[j]] -= beta_*qx_[i]*qy_[j]*ex;
      g_alpha -= beta_*bp_score*ex;
      g_beta -= (alpha_*bp_score+nobp_score)*ex;
    }
  g_beta -= ex_gap*gap_ + ex_ext*ext_;
  g_gap -= beta_*ex_gap;
  g_ext -= beta_*ex_ext;

  return v-log(Z);
}

template
class PartAlign<float>;

#include "log_value.h"

template
class PartAlign<LogValue<float> >;

#if 0
int
main()
{
  PartAlign<float> a;
  std::string x("augcaugc"), y("auggaugg");
  a.load_sequences(x, y);
  a.compute_forward();
  a.compute_backward();
  a.compute_posterior();

  return 0;
}

#endif
