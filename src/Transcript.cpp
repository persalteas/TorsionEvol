#include "Transcript.h"

Transcript::Transcript(uint TUindex, DNApos TSS, DNApos TTS, int strand, double rate, int dx)
                      : TUindex_(TUindex), TSS_(TSS), TTS_(TTS), s_(strand), r_(rate)
{
    size_ = 1 + ::abs(int(TTS_) - int(TSS_)); // need to allow negative diffference
    size_n_ = size_/dx;
    start_ = TSS_/dx;
    end_ = TTS_/dx;
}

/* calculate the initiation rate (math formula) */
void Transcript::f_init_rate( double sigma_t, double epsilon, double m, 
                                vector<DNApos>& Barr_pos, vector<double>& Barr_sigma)
{
    // get the topological domain of this TU
    DNApos start = start_;
    uint domain_index = find_if(Barr_pos.begin(), Barr_pos.end(), 
						    [start](DNApos_n barr)->bool{ return barr>start; }) - Barr_pos.begin();
    // The torsion at this TU
    sigma_ = Barr_sigma[domain_index-1];
    // The init rate, function of the torsion
	init_rate_ = r_ * exp(m/(1+exp((sigma_-sigma_t)/epsilon)));
}

/* normalization of init rate to get a probability */
void Transcript::f_prob_init_rate(double sum_init_rates, int DELTA_T)
{
    prob_init_rate_ = (1.0 - exp(-sum_init_rates*double(DELTA_T)))*init_rate_/sum_init_rates;
	// ( 1 - np.exp(-sum_init_rate*DELTA_T)) * (init_rate/sum_init_rate)
}

