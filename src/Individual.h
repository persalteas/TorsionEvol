#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

class Individual
{
  private:
    unsigned n_genes{};
    unsigned n_barriers{};
    unsigned* gene_start{};
    unsigned* gene_end{};
    unsigned* barrier{};
    float fitness{};
  public:
    Individual(void); // default constructor
    Individual(const Individual& indiv); // copy constructor
    Individual(Individual&& indiv); // move constructor
    Individual(unsigned n_genes,
               unsigned n_barriers,
               unsigned* gene_start,
               unsigned* gene_end,
               unsigned* barrier,
               float fitness);
    Individual& operator=(const Individual& indiv); // copy assignment operator
    Individual& operator=(Individual&& indiv); // move assignment operator
    Individual* mutate();
    unsigned n_genes(void) {return n_genes;}
    unsigned n_barriers(void) {return n_barriers;}
    unsigned* gene_start(unsigned index);
    unsigned* gene_end(unsigned index);
    unsigned* barrier(unsigned index);
    float fitness(void) {return fitness;}
    void reset(void);
    ~Individual(void); // destructor
};

#endif
