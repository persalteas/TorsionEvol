#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

class Individual
{
  private:
    static double _mean; // poisson distribution parameter for indels
    static double _p; // probability for an inversion to occur during replication
    unsigned _size{};
    unsigned _n_genes{};
    unsigned* _genes{}; // even: start; uneven: end
    unsigned _n_barriers{};
    unsigned* _barriers{};
    float _fitness{};
  public:
    Individual(void) {} // default constructor
    Individual(const Individual& indiv); // copy constructor
    Individual(Individual&& indiv); // move constructor
    Individual(unsigned size,
               unsigned n_genes,
               unsigned n_barriers,
               unsigned* genes,
               unsigned* barriers,
               float fitness);
    Individual& operator=(const Individual& indiv); // copy assignment operator
    Individual& operator=(Individual&& indiv); // move assignment operator
    static double get_mean(void) {return _mean;}
    static double get_p(void) {return _p;}
    static void set_mutation(double mean, double p);
    unsigned get_size(void) const {return _size;}
    unsigned get_n_genes(void) const {return _n_genes;}
    unsigned* get_genes(void) const {return _genes;}
    unsigned gene_start(unsigned index) const {return _genes[2*index];}
    unsigned gene_end(unsigned index) const {return _genes[2*index+1];}
    unsigned get_n_barriers(void) const {return _n_barriers;}
    unsigned* get_barriers(void) const {return _barriers;}
    unsigned barrier(unsigned index) const {return _barriers[index];}
    float get_fitness(void) const {return _fitness;}
    void mutate(void);
    void reset(void);
    ~Individual(void); // destructor
};

#endif
