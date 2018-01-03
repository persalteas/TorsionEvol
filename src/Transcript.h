#ifndef TRANSCRIPT_H_
#define TRANSCRIPT_H_

#include "utils.h"

typedef unsigned int DNApos;
typedef unsigned int DNApos_n;

class Transcript 
{
    public:
        Transcript(uint TUindex, DNApos TSS, DNApos TTS, int strand, double rate, int dx);
        uint TUindex_;
        DNApos TSS_;
        DNApos TTS_;
        size_t size_;
        DNApos_n start_;
        DNApos_n end_;
        size_t size_n_;
        int s_;
        double r_;
        double sigma_;
        double init_rate_;
        double prob_init_rate_;
        void f_init_rate(double sigma_t, double epsilon, double m, vector<DNApos>& Barr_pos, vector<double>& Barr_sigma);
        void f_prob_init_rate(double sum_init_rates, int DELTA_T);
        uint domain( vector<DNApos>& Barr_pos);

};

#endif //TRANSCRIPT_H_