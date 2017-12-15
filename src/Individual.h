#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

class Individual
{
  private:
    static double mean; // poisson distribution parameter for indels
    static double p; // probability for an inversion to occur during replication
    unsigned size{};
    unsigned n_genes{};
    unsigned* genes{}; // even: start; uneven: end
    unsigned n_barriers{};
    unsigned* barriers{};
    float fitness{};
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
    static double mean(void) {return mean;}
    static double p(void) {return p;}
    static void set_mutation(double indel_mean, double inv_p);
    unsigned size(void) {return size;}
    unsigned n_genes(void) {return n_genes;}
    unsigned* genes(void) {return genes;}
    unsigned gene_start(unsigned index) {return genes[2*index];}
    unsigned gene_end(unsigned index) {return genes[2*index+1];}
    unsigned n_barriers(void) {return n_barriers;}
    unsigned* barriers(void) {return barriers;}
    unsigned barrier(unsigned index) {return barriers[index];}
    float fitness(void) {return fitness;}
    void mutate(void);
    void reset(void);
    ~Individual(void); // destructor
};

#endif
