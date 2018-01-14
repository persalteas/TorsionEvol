#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "Random.h"
#include "Transcript.h"
#include "utils.h"

class Individual {
private:
  static double _mean; // poisson distribution parameter for indels
  static double _p; // probability for an inversion to occur during replication
  static float _RNAPs_genSC;           // from Meyer simulation
  static Params *_params;              // params from Meyer simulation
  static vector<double> _target_envir; // optimal genes proportions to reach
  static Random _rand_generator;
  vector<Transcript> *_genes;
  vector<DNApos> *_barr_fix;
  unsigned _genome_size{};
  float _fitness{};
  void estimate_exression(void);
  DNApos get_rnd_pos_btwn_genes(vector<DNApos> &gene_pos,
                                vector<size_t> &Dom_size);
  DNApos get_rnd_pos_in_domain(uint dom, vector<DNApos> &gene_pos,
                               vector<size_t> &Dom_size);
  uint get_rnd_dom_btwn_genes(const vector<size_t> &Dom_size);
  void display_state(vector<DNApos> &Barr_pos, vector<int> &Barr_type);
  void display_state(void);

public:
  Individual(void);                    // default constructor
  Individual(const Individual &indiv); // copy constructor
  Individual(Individual &&indiv);      // move constructor
  Individual(unsigned genome_size, // custom constructor with genes and barriers
             vector<Transcript> &tr, vector<DNApos> &Barr_fix);
  Individual &operator=(const Individual &indiv); // copy assignment operator
  Individual &operator=(Individual &&indiv);      // move assignment operator
  static double get_mean(void) { return _mean; }
  static double get_p(void) { return _p; }
  static void set_mutation(double mean, double p);
  static void set_simulation_params(Params *params);
  static void set_target_envir(vector<double> &env);
  unsigned get_size(void) const { return _genome_size; }
  unsigned get_n_genes(void) const { return _genes->size(); }
  vector<Transcript> &get_genes(void) const { return *_genes; }
  void update_fitness(void);
  unsigned get_n_barriers(void) const { return _barr_fix->size(); }
  vector<DNApos> &get_barriers(void) const { return *_barr_fix; }
  float get_fitness(void) const { return _fitness; }
  void mutate(void);
  ~Individual(void); // destructor
};

#endif
