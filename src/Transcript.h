#ifndef TRANSCRIPT_H_
#define TRANSCRIPT_H_

#include "utils.h"

class Transcript {
public:
  Transcript(uint TUindex, DNApos TSS, DNApos TTS, int strand, double rate,
             int dx);
  Transcript(const Transcript &tr); // customized copy constructor
  uint TUindex_;
  DNApos TSS_;
  DNApos TTS_;
  size_t size_;
  DNApos_n start_;
  DNApos_n end_;
  size_t size_n_;
  int s_;
  double r_;
  uint expr_count_;
  double f_init_rate(double sigma_t, double epsilon, double m,
                     const std::vector<DNApos> &Barr_pos,
                     const std::vector<double> &Barr_sigma);
  double f_prob_init_rate(double sum_init_rates, int DELTA_T);
  void shift(int d_n);
  void revert(DNApos posA, DNApos posB);

private:
  double init_rate_;
};

#endif // TRANSCRIPT_H_
