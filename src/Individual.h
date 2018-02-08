#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "RNAP.h"
#include "Random.h"
#include "Transcript.h"
#include "utils.h"
// #include <boost/move/unique_ptr.hpp>

// namespace boost_mov = boost::movelib;

class Individual {
private:
  static double _mean; // poisson distribution parameter for indels
  static double _p; // probability for an inversion to occur during replication
  static float _RNAPs_genSC;           // from Meyer simulation
  static Params *_params;              // params from Meyer simulation
  static std::vector<double> _target_envir; // optimal genes proportions to reach
  static Random _rand_generator;
  uint _genome_size{};
  double _fitness{};
  std::vector<Transcript> _genes;
  std::vector<DNApos_n> _barr_fix;
  std::vector<uint> Barr_pos;
  std::vector<int> Barr_type;
  std::vector<double> Barr_sigma;
  std::vector<int> Dom_size;
  std::vector<size_t> Dom_size_1;
  std::vector<RNAP> RNAPs_hooked;
  std::vector<uint> rm_RNAPs_idx;
  std::vector<int> tss_and_unhooked_RNAPs;
  std::vector<int> picked_tr;
  std::vector<double> init_rates;
  std::vector<double> prob_init_rates;
  std::vector<double> prob_unhooked_rates;
  std::vector<double> all_prob;
  std::vector<DNApos> genes_pos;
  void estimate_exression(void);
  DNApos get_rnd_pos_in_domain(int dom);
  int get_rnd_dom_btwn_genes(void);
  void display_state(void);

public:
  Individual(const Individual &indiv); // copy constructor
  Individual(unsigned genome_size, // custom constructor with genes and barriers
             std::vector<Transcript> &tr, std::vector<DNApos> &Barr_fix);
  static double get_mean(void) { return _mean; }
  static double get_p(void) { return _p; }
  static void set_mutation(double mean, double p);
  static void set_simulation_params(Params *params);
  static void set_target_envir(std::vector<double> &env);
  unsigned get_size(void) const { return _genome_size; }
  unsigned get_n_genes(void) const { return _genes.size(); }
  void update_fitness(void);
  unsigned get_n_barriers(void) const { return _barr_fix.size(); }
  double get_fitness(void) const { return _fitness; }
  void mutate(void);
};

#endif
