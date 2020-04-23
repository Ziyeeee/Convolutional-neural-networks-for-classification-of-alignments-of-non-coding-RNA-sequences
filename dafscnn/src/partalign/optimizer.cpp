// $Id:$

#include "config.h"
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include "partalign.h"
#include "optimizer.h"
#include "lbfgsb.h"

using namespace PARTALIGN;
using std::isinf;

Optimizer::
Optimizer(bool alpha, bool beta, bool gap, bool ext, bool sm,
          float c, float eta0, uint t_max, float eps)
  : c_(c), eta0_(eta0), t_max_(t_max), eps_(eps),
    update_alpha_(alpha), update_beta_(beta),
    update_gap_(gap), update_ext_(ext), update_sm_(sm)
{
}

void
Optimizer::
set_parameters(float alpha, float beta, float gap, float ext, const VVF& sm)
{
  alpha_ = alpha;
  beta_ = beta;
  gap_ = gap;
  ext_ = ext;
  sm_ = sm;
}

void
Optimizer::
set_parameters(float alpha, float beta, float gap, float ext, const float sm[10])
{
  VVF m(4, VF(4, 0.0));
  for (uint i=0; i<4; ++i)
  {
    m[i][i] = *sm++;
    for (uint j=i+1; j<4; ++j)
    {
      m[i][j] = m[j][i] = *sm++;
    }
  }
  set_parameters(alpha, beta, gap, ext, m);
}

void
Optimizer::
get_parameters(float& alpha, float& beta, float& gap, float& ext, VVF& sm)
{
  alpha = alpha_;
  beta = beta_;
  gap = gap_;
  ext = ext_;
  sm = sm_;
}

void
Optimizer::
get_parameters(float& alpha, float& beta, float& gap, float& ext, float sm[10])
{
  VVF m(4, VF(4, 0.0));
  get_parameters(alpha, beta, gap, ext, m);
  for (uint i=0; i<4; ++i)
  {
    *sm++ = m[i][i];
    for (uint j=i+1; j<4; ++j)
      *sm++ = m[i][j];
  }
}

template < class V, class Seq >
void
Optimizer::
sgd(const std::vector<std::pair<Seq,Seq> >& seq,
    const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln)
{
  assert(seq.size()==aln.size());
  float eta=eta0_;               // step width
  PartAlign<V> a;
  VF f_prev(seq.size(), 0.0);
  std::vector<uint> idx(seq.size());
  for (uint s=0; s!=seq.size(); ++s) idx[s]=s;
  uint n=0;
  float C=c_/seq.size();

  std::cout << seq.size() << std::endl;
  for (uint t=0; t!=t_max_; ++t)
  {
    float alpha_prev=alpha_, beta_prev=beta_, gap_prev=gap_, ext_prev=ext_;
    VVF sm_prev=sm_;
    float f_sum=0.0, f_diff=0.0;
    // randomize the order of the optimization
    std::random_shuffle(idx.begin(), idx.end());
    
    for (uint s=0; s!=idx.size(); ++s)
    {
      // set the current parameters
      a.set_parameters(alpha_, beta_, gap_, ext_);
      a.set_scoring_matrix(sm_);

      // calculate gradients
      float g_alpha, g_beta, g_gap, g_ext;
      VVF g_sm(4, VF(4, 0.0));
      a.load_sequences(seq[idx[s]].first, seq[idx[s]].second);
      float f = a.compute_gradients(aln[idx[s]].first, aln[idx[s]].second,
                                    g_alpha, g_beta, g_gap, g_ext, g_sm);

      if (isinf(f)) continue;
#if 0
      std::cout << "\tf:" << f << " eta:" << eta << std::endl;
      std::cout << "\t" << alpha_ << ", " << beta_ << ", "
                << gap_ << ", " << ext_ << std::endl;
      std::cout << "\t";
      for (uint i=0; i<4; ++i)
      {
        std::cout << sm_[i][i] << ", ";
        for (uint j=i+1; j<4; ++j)
          std::cout << sm_[i][j] << ", ";
      }
      std::cout << std::endl;
      std::cout << "\t" << g_alpha << ", " << g_beta << ", "
                << g_gap << ", " << g_ext << std::endl;
      std::cout << "\t";
      for (uint i=0; i<4; ++i)
      {
        std::cout << g_sm[i][i] << ", ";
        for (uint j=i+1; j<4; ++j)
          std::cout << g_sm[i][j] << ", ";
      }
      std::cout << std::endl << std::endl;
#endif
      // update parameters
      if (update_alpha_) alpha_ = std::min(std::max(0.0f, alpha_+eta*(g_alpha-alpha_/C)), alpha_+1.0f);
      if (update_beta_) beta_ = std::min(std::max(1e-3f, beta_+eta*(g_beta-beta_/C)), beta_+1.0f);
      if (update_gap_) gap_ = std::min(std::max(gap_-1.0f, gap_+eta*(g_gap-gap_/C)), gap_+1.0f);
      if (update_ext_) ext_ = std::min(std::max(ext_-1.0f, ext_+eta*(g_ext-ext_/C)), ext_+1.0f);
      if (update_sm_)
        for (uint i=0; i<4; ++i)
        {
          sm_[i][i] = std::min(std::max(sm_[i][i]-1.0f, sm_[i][i]+eta*(g_sm[i][i]-sm_[i][i]/C)), sm_[i][i]+1.0f);
          for (uint j=i+1; j<4; ++j)
            sm_[i][j] = sm_[j][i] = std::min(std::max(sm_[i][j]-1.0f, sm_[i][j]+eta*(g_sm[i][j]+g_sm[j][i]-sm_[i][j]/C)), sm_[i][j]+1.0f); // assume the symmetry
        }

      f_sum += f;
      f_diff += std::abs(f-f_prev[idx[s]]);
      // update the step width
      if (f_prev[idx[s]]<f || t==0)
      {
        n++;
        eta = eta0_ * seq.size() / (seq.size()+n);
      }
      f_prev[idx[s]] = f;
    }

    float rmsd=0.0;
    rmsd += (alpha_-alpha_prev)*(alpha_-alpha_prev);
    rmsd += (beta_-beta_prev)*(beta_-beta_prev);
    rmsd += (gap_-gap_prev)*(gap_-gap_prev);
    rmsd += (ext_-ext_prev)*(ext_-ext_prev);
    for (uint i=0; i<4; ++i)
    {
      rmsd += (sm_[i][i]-sm_prev[i][i])*(sm_[i][i]-sm_prev[i][i]);
      for (uint j=i+1; j<4; ++j)
        rmsd += (sm_[i][j]-sm_prev[i][j])*(sm_[i][j]-sm_prev[i][j]);
    }
    rmsd = sqrt(rmsd);
    std::cout << "t=" << t << " f=" << f_sum << " (" << f_diff << ") "
              << "RMSD=" << rmsd << " "
              << "eta=" << eta
              << std::endl;
#if 0
    std::cout << " " << alpha_ << ", " << beta_ << ", "
              << gap_ << ", " << ext_ << std::endl;
    std::cout << " ";
    for (uint i=0; i<4; ++i)
    {
      std::cout << sm_[i][i] << ", ";
      for (uint j=i+1; j<4; ++j)
        std::cout << sm_[i][j] << ", ";
    }
    std::cout << std::endl;
#endif

    // check the convergence
    //std::cout << "f: " << t << ": " << f_sum << " (" << f_diff << ")" <<std::endl;
    //if (f_diff/std::abs(f_sum)<eps_) break;
    if (rmsd<eps_) break;
  }
}

template < class V, class Seq >
void
Optimizer::
lbfgs(const std::vector<std::pair<Seq,Seq> >& seq,
      const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln)
{
  assert(seq.size()==aln.size());

  uint k;
  // set the initial parameters
  std::vector<double> x;
  x.push_back(alpha_);
  x.push_back(beta_);
  x.push_back(gap_);
  x.push_back(ext_);
  for (uint i=0; i<4; ++i)
  {
    x.push_back(sm_[i][i]);
    for (uint j=i+1; j<4; ++j)    // assume the symmetric matrix
      x.push_back(sm_[i][j]);
  }
  // initialize the optimizer
  std::vector<double> lbd(x.size(), 0.0);
  std::vector<double> ubd(x.size(), 0.0);
  std::vector<long> nbd(x.size(), LBFGSB::UNBOUND);
  k=0;
  lbd[k]=0.0; nbd[k]=LBFGSB::LOWER_BOUND; k++; // alpha
  lbd[k]=1e-3; nbd[k]=LBFGSB::LOWER_BOUND; k++; // beta
  LBFGSB lbfgsb(/*factr, pgtol*/ 1e1, eps_);
  lbfgsb.initialize(x.size(), 5, &lbd[0], &ubd[0], &nbd[0]);

  int iflag=1;
  PartAlign<V> a;
  for (uint t=0; t!=t_max_ && iflag>0; ++t)
  {
    // set the current parameters
    k=0;
    if (update_alpha_) alpha_ = x[k]; k++;
    if (update_beta_) beta_ = x[k]; k++;
    if (update_gap_) gap_ = x[k]; k++;
    if (update_ext_) ext_ = x[k]; k++;
    a.set_parameters(alpha_, beta_, gap_, ext_);
    if (update_sm_)
      for (uint i=0; i<4; ++i)
      {
        sm_[i][i]=x[k++];
        for (uint j=i+1; j<4; ++j)
          sm_[i][j]=sm_[j][i]=x[k++];
      }
    a.set_scoring_matrix(sm_);

    // calculate gradients for all alignments
    double f=0.0;
    std::vector<double> g(x.size(), 0.0);
    for (uint s=0; s!=seq.size(); ++s)
    {
      float g_alpha, g_beta, g_gap, g_ext;
      VVF g_sm(4, VF(4, 0.0));
      a.load_sequences(seq[s].first, seq[s].second);
      double ff = a.compute_gradients(aln[s].first, aln[s].second, g_alpha, g_beta, g_gap, g_ext, g_sm);
      if (isinf(ff)) continue;
      f -= ff;
      k=0;
      g[k++] -= g_alpha;
      g[k++] -= g_beta;
      g[k++] -= g_gap;
      g[k++] -= g_ext;
      for (uint i=0; i<4; ++i)
      {
        g[k++] -= g_sm[i][i];
        for (uint j=i+1; j<4; ++j)
          g[k++] -= g_sm[i][j]+g_sm[j][i];
      }
    }

    // add L2 priors
    for (uint k=0; k!=x.size(); ++k)
    {
      f += x[k]*x[k]/(2*c_);
      g[k] += x[k]/c_;
    }

    std::cout << "t=" << t << " f=" << -f << std::endl;
#if 1
    std::copy(x.begin(), x.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
    std::copy(g.begin(), g.end(), std::ostream_iterator<double>(std::cout, ", "));
    std::cout << std::endl;
#endif

    iflag = lbfgsb.update(&x[0], f, &g[0]);
  }

  // set the current parameters
  k=0;
  if (update_alpha_) alpha_ = x[k]; k++;
  if (update_beta_) beta_ = x[k]; k++;
  if (update_gap_) gap_ = x[k]; k++;
  if (update_ext_) ext_ = x[k]; k++;
  if (update_sm_)
    for (uint i=0; i<4; ++i)
    {
      sm_[i][i]=x[k++];
      for (uint j=i+1; j<4; ++j)
        sm_[i][j]=sm_[j][i]=x[k++];
    }
}

template
void
Optimizer::
sgd<LogValue<float>,BPSeq>(const std::vector<std::pair<BPSeq,BPSeq> >& seq,
                           const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln);

template
void
Optimizer::
sgd<LogValue<float>,std::string>(const std::vector<std::pair<std::string,std::string> >& seq,
                                 const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln);

template
void
Optimizer::
lbfgs<LogValue<float>,BPSeq>(const std::vector<std::pair<BPSeq,BPSeq> >& seq,
                             const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln);

template
void
Optimizer::
lbfgs<LogValue<float>,std::string>(const std::vector<std::pair<std::string,std::string> >& seq,
                             const std::vector<std::pair<std::vector<bool>,std::vector<bool> > >& aln);
