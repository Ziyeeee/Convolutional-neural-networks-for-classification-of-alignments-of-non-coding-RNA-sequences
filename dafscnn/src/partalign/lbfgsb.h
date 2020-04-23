/* $Id: lbfgsb.h 272 2008-08-13 03:21:46Z satoken $ */

#ifndef __INC_LBFGSB_H__
#define __INC_LBFGSB_H__

#include <vector>
#include <algorithm>
#include <cfloat>
#include <cstring>

extern "C" {
  int setulb_(long *n, long *m, double *x, double *l, double *u, long *nbd,
	      double *f, double *g, double *factr, double *pgtol,
	      double *wa, long *iwa, char *task, long *iprint,
	      char *csave, long *lsave, long *isave, double *dsave,
	      short task_len, short csave_len);
};

class LBFGSB
{
public:
  enum { UNBOUND=0,
	 LOWER_BOUND=1,
	 LOWER_UPPER_BOUND=2,
	 UPPER_BOUND=3 };
  
private:
  long n_;
  long m_;
  std::vector<double> l_;
  std::vector<double> u_;
  std::vector<long> nbd_;
  std::vector<double> wa_;
  std::vector<long> iwa_;
  char task_[60];
  char csave_[60];
  long lsave_[4];
  long isave_[44];
  double dsave_[29];
  long iprint_;
  double factr_;
  double pgtol_;
    
  
public:
  LBFGSB(double factr=1e1, double pgtol=1e-5)
    : n_(0), m_(0), l_(), u_(), nbd_(), wa_(), iwa_(),
      iprint_(1), factr_(factr), pgtol_(pgtol)
  {
  }
  
  void initialize(long n, long m)
  {
    std::vector<double> l(n);
    std::vector<double> u(n);
    std::vector<long> nbd(n);
    std::fill(nbd.begin(), nbd.end(), 0);
    initialize(n, m, &l[0], &u[0], &nbd[0]);
  }

  void initialize(long n, long m,
		  const double *l, const double *u, const long *nbd)
  {
    n_ = n;
    m_ = m;
    wa_.resize(2*m*n+4*n+11*m*m+8*m);
    iwa_.resize(3*n);
    l_.resize(n);
    u_.resize(n);
    nbd_.resize(n);
    std::copy(l, l+n, l_.begin());
    std::copy(u, u+n, u_.begin());
    std::copy(nbd, nbd+n, nbd_.begin());
    strcpy(task_, "START");
  }
    
  int update(double *x, double f, const double *g)
  {
    if (strncmp(task_, "START", 5)==0) {
      setulb_(&n_, &m_, x, &l_[0], &u_[0], &nbd_[0], &f, const_cast<double*>(g),
	      &factr_, &pgtol_, &wa_[0], &iwa_[0],
	      task_, &iprint_, csave_, lsave_, isave_, dsave_, 60, 60);
      if (strncmp(task_, "FG", 2)!=0) return -1;
    }
    
    while (1) {
      setulb_(&n_, &m_, x, &l_[0], &u_[0], &nbd_[0], &f, const_cast<double*>(g),
	      &factr_, &pgtol_, &wa_[0], &iwa_[0],
	      task_, &iprint_, csave_, lsave_, isave_, dsave_, 60, 60);
      //std::cout << task_ << std::endl;
      if (strncmp(task_, "FG", 2)==0) return 1;
      else if (strncmp(task_, "NEW_X", 5)==0) continue;
      else if (strncmp(task_, "CONV", 4)==0) return 0;
      else {
	std::cout << task_ << std::endl;
	return -1;
      }
    }
  }
};

#endif /* __INC_LBFGSB_H__ */

// Local Variables:
// mode:C++
// End:
