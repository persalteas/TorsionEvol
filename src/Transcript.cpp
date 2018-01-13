#include "Transcript.h"

Transcript::Transcript(uint TUindex, DNApos TSS, DNApos TTS, int strand,
                       double rate, int dx)
    : TUindex_(TUindex), TSS_(TSS), TTS_(TTS), s_(strand), r_(rate) {
  size_ =
      1 + ::abs(int(TTS_) - int(TSS_)); // need to allow negative diffference
  size_n_ = size_ / dx;
  start_ = TSS_ / dx;
  end_ = TTS_ / dx;
  expr_count_ = 0;
}

Transcript::Transcript(const Transcript &tr)
    : TUindex_(tr.TUindex_), TSS_(tr.TSS_), TTS_(tr.TTS_), size_(tr.size_),
      start_(tr.start_), end_(tr.end_), size_n_(tr.size_n_), s_(tr.s_),
      r_(tr.r_) {
  // This copy constructor preserves the position end rate but resets the
  // expression count:
  expr_count_ = 0;
}

/* calculate the initiation rate (math formula) */
double Transcript::f_init_rate(double sigma_t, double epsilon, double m,
                               vector<DNApos> &Barr_pos,
                               vector<double> &Barr_sigma) {
  // get the topological domain of this TU
  DNApos start = start_;
  uint domain_index =
      find_if(Barr_pos.begin(), Barr_pos.end(),
              [start](DNApos_n barr) -> bool { return barr > start; }) -
      Barr_pos.begin();
  // The torsion at this TU
  double sigma_ = Barr_sigma[domain_index - 1];
  // The init rate, function of the torsion
  init_rate_ = r_ * exp(m / (1 + exp((sigma_ - sigma_t) / epsilon)));
  return init_rate_;
}

/* normalization of init rate to get a probability */
double Transcript::f_prob_init_rate(double sum_init_rates, int DELTA_T) {
  return (1.0 - exp(-sum_init_rates * double(DELTA_T))) * init_rate_ /
         sum_init_rates;
}

void Transcript::shift(int d, int d_n) {
  TSS_ += d;
  TTS_ += d;
  start_ += d_n;
  end_ += d_n;
}
