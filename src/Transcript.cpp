#include "Transcript.h"

Transcript::Transcript(uint TUindex, DNApos TSS, DNApos TTS, int strand,
                       double rate, int dx)
    : TUindex_(TUindex), TSS_(TSS), TTS_(TTS), s_(strand), r_(rate) {
  size_ = 1 + ::abs(int(TTS_) - int(TSS_));
  size_n_ = size_ / dx;
  start_ = 1 + TSS_ / dx;
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
                               const vector<DNApos> &Barr_pos,
                               const vector<double> &Barr_sigma) {
  // get the topological domain of this TU
  DNApos start = start_;
  uint domain_index;
  if (start < Barr_pos.back())
    domain_index =
        find_if(Barr_pos.begin(), Barr_pos.end(),
                [start](DNApos_n barr) -> bool { return barr > start; }) -
        Barr_pos.begin();
  else
    domain_index = 0;
  // The torsion at this TU
  double sigma_ = Barr_sigma[domain_index];
  // The init rate, function of the torsion
  init_rate_ = r_ * exp(m / (1 + exp((sigma_ - sigma_t) / epsilon)));
  if (init_rate_ != init_rate_) {
    cout << "## erreur, init_rate = nan ##" << endl;
  }
  return init_rate_;
}

/* normalization of init rate to get a probability */
double Transcript::f_prob_init_rate(double sum_init_rates, int DELTA_T) {
  return (1.0 - exp(-sum_init_rates * double(DELTA_T))) * init_rate_ /
         sum_init_rates;
}

void Transcript::shift(int d_n) {
  int d_ = d_n * int(size_ / size_n_);
  TSS_ += d_;
  TTS_ += d_;
  start_ += d_n;
  end_ += d_n;
}

void Transcript::revert(DNApos posA, DNApos posB) {
  start_ = posA + posB - start_;
  end_ = posA + posB - end_;
  uint bc = (posA + posB) * size_ / size_n_;
  TSS_ = bc - TSS_;
  TTS_ = bc - TTS_;
  s_ = -s_;
}