// $Id:$

#include "config.h"

#include <cstring>
#include <iostream>
#include <vector>

#include "optimizer.h"
#include "log_value.h"
#include "../fa.h"
#include "../fold.h"
#include "partalign.h"

#define FOREACH(itr, i, v) for (itr i=(v).begin(); i!=(v).end(); ++i)
#define CUTOFF 0.01

class PATrain
{
public:
  PATrain()
    : c_(1.0),
      eta0_(0.5),
      t_max_(100),
      eps_(1e-5),
      alpha_(0.0),
      beta_(1.0),
      gap_(-10.0),
      ext_(-5.0),
      use_lbfgs_(false),
      use_alpha_(true),
      use_beta_(true),
      use_gap_(true),
      use_ext_(true),
      use_sm_(true),
      en_s_(NULL),
      use_stdin(false)
  {
    static float sm[10] = {     // RIBOSUM85-60
      2.22, -1.86, -1.46, -1.39,
      /*  */ 1.16, -2.48, -1.05,
      /*         */ 1.03, -1.74,
      /*                */ 1.65
    };
    std::copy(sm, sm+10, sm_);
  }

  ~PATrain()
  {
    if (en_s_) delete en_s_;
  }

  PATrain& parse_options(int& argc, char**& argv)
  {
    int ch;
    const char* s_en = NULL;
    while ((ch=getopt(argc, argv, "hc:a:b:g:e:t:d:lABGESs:f"))!=-1)
    {
      switch (ch)
      {
        case 'c':
          c_ = atof(optarg);
          break;
        case 'a':
          alpha_ = atof(optarg);
          break;
        case 'b':
          beta_ = atof(optarg);
          break;
        case 'g':
          gap_ = atof(optarg);
          break;
        case 'e':
          ext_ = atof(optarg);
          break;
        case 't':
          t_max_ = atoi(optarg);
          break;
        case 'd':
          eta0_ = atof(optarg);
          break;
        case 'l':
          use_lbfgs_ = true;
          break;
        case 'A':
          use_alpha_ = false;
          break;
        case 'B':
          use_beta_ = false;
          break;
        case 'G':
          use_gap_ = false;
          break;
        case 'E':
          use_ext_ = false;
          break;
        case 'S':
          use_sm_ = false;
          break;
        case 's':
          s_en = optarg;
          break;
        case 'f':
          use_stdin = true;
          break;
        case 'h': case '?': default:
          exit(0);
          break;
      }
    }
    argc -= optind;
    argv += optind;

    if (s_en)
    {
      if (strcasecmp(s_en, "Boltzmann")==0)
        en_s_ = new RNAfold(true, NULL, CUTOFF);
      if (strcasecmp(s_en, "Vienna")==0)
        en_s_ = new RNAfold(false, NULL, CUTOFF);
      else if (strcasecmp(s_en, "CONTRAfold")==0)
        en_s_ = new CONTRAfold(CUTOFF);
    }

    return *this;
  }

  bool read_data(const char* file,
                 std::pair<std::string,std::string>& s,
                 std::pair<std::vector<bool>,std::vector<bool> >& a)
  {
    std::vector<Fasta> fa;
    if (Fasta::load(fa, file)<2) return false;
    return read_data(fa[0], fa[1], s, a);
  }

  bool read_data(const Fasta& fa0, const Fasta& fa1,
                 std::pair<std::string,std::string>& s,
                 std::pair<std::vector<bool>,std::vector<bool> >& a)
  {
    assert(fa0.size()==fa1.size());
    const uint L = fa0.size();
    for (uint j=0; j!=L; ++j)
    {
      if (fa0.seq()[j]=='-' && fa1.seq()[j]=='-') continue;
      if (fa0.seq()[j]!='-') s.first.push_back(fa0.seq()[j]);
      if (fa1.seq()[j]!='-') s.second.push_back(fa1.seq()[j]);
      a.first.push_back(fa0.seq()[j]!='-');
      a.second.push_back(fa1.seq()[j]!='-');
    }
    return true;
  }

  int run(uint n, char** files)
  {
    std::cout << "loading ..." << std::flush;
    std::vector< std::pair<std::string,std::string> > seq;
    std::vector<std::pair<std::vector<bool>,std::vector<bool> > > aln;
    if (use_stdin)
    {
      std::string l;
      while (std::cin >> l)
      {
        std::vector<Fasta> fa;
        if (Fasta::load(fa, l.c_str())<2) continue;
        for (std::vector<Fasta>::iterator x=fa.begin(); x!=fa.end(); ++x)
          if (x->name()=="SS_cons") { fa.erase(x); break; }
        for (uint j=0; j!=fa.size()-1; ++j)
          for (uint k=j+1; k!=fa.size(); ++k)
          {
            std::pair<std::string,std::string> s;
            std::pair<std::vector<bool>,std::vector<bool> > a;
            if (!read_data(fa[j], fa[k], s, a)) continue;
            seq.push_back(s);
            aln.push_back(a);
          }
      }
    }
    else
    {
      for (uint i=0; i!=n; ++i) 
      {
        std::vector<Fasta> fa;
        if (Fasta::load(fa, files[i])<2) continue;
        for (std::vector<Fasta>::iterator x=fa.begin(); x!=fa.end(); ++x)
          if (x->name()=="SS_cons") { fa.erase(x); break; }
        for (uint j=0; j!=fa.size()-1; ++j)
          for (uint k=j+1; k!=fa.size(); ++k)
          {
            std::pair<std::string,std::string> s;
            std::pair<std::vector<bool>,std::vector<bool> > a;
            if (!read_data(fa[j], fa[k], s, a)) continue;
            seq.push_back(s);
            aln.push_back(a);
          }
      }
    }
    std::cout << " done. (" << seq.size() << " seqs)" << std::endl;

    PARTALIGN::Optimizer opt(use_alpha_, use_beta_, use_gap_, use_ext_, use_sm_,
                             c_, eta0_, t_max_, eps_);
    opt.set_parameters(alpha_, beta_, gap_, ext_, sm_);

    if (en_s_)
    {
      using PARTALIGN::BPSeq;
      std::vector< std::pair<BPSeq,BPSeq> > bpseq;
      std::cout << "calculating BP ..." << std::flush;
      for (uint i=0; i!=seq.size(); ++i)
      {
        BP bp1, bp2;
        en_s_->calculate(seq[i].first, bp1);
        en_s_->calculate(seq[i].second, bp2);
        bpseq.push_back(std::make_pair(std::make_pair(seq[i].first, bp1),
                                       std::make_pair(seq[i].second, bp2)));
      }
      std::cout << " done." << std::endl;
      if (!use_lbfgs_)
        opt.sgd<LogValue<float>,BPSeq>(bpseq, aln);
      else
        opt.lbfgs<LogValue<float>,BPSeq>(bpseq, aln);
      opt.get_parameters(alpha_, beta_, gap_, ext_, sm_);
    }
    else
    {
      if (!use_lbfgs_)
        opt.sgd<LogValue<float>,std::string>(seq, aln);
      else
        opt.lbfgs<LogValue<float>,std::string>(seq, aln);
      opt.get_parameters(alpha_, beta_, gap_, ext_, sm_);
    }

    std::cout << alpha_ << ", "
              << beta_ << ", "
              << gap_ << ", "
              << ext_ << std::endl;
    std::copy(sm_, sm_+10,
              std::ostream_iterator<float>(std::cout, ", "));
    std::cout << std::endl;

    return 0;
  }

private:
  float c_;
  float eta0_;
  uint t_max_;
  float eps_;
  float alpha_;
  float beta_;
  float gap_;
  float ext_;
  float sm_[10];
  bool use_lbfgs_;
  bool use_alpha_;
  bool use_beta_;
  bool use_gap_;
  bool use_ext_;
  bool use_sm_;
  Fold::Model* en_s_;                  // folding engine
  bool use_stdin;
};

int
main(int argc, char* argv[])
{
  PATrain pa;
  pa.parse_options(argc, argv);
  return pa.run(argc, argv);
}
