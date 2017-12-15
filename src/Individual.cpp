#include "Individual.cpp"
#include <random>

double Individual::mean = 5;
double Individual::p = 0.5;

Individual::Individual(const Individual& indiv)
{
  size = indiv.size();
  n_genes = indiv.n_genes();
  genes = new unsigned[2 * n_genes];
  for (unsigned* i = genes, j = indiv.genes(); i < genes + 2 * n_genes; ++i, ++j)
    *i = *j;
  n_barriers = indiv.n_barriers();
  barriers = new unsigned[n_barriers];
  for (unsigned* i = barriers, j = indiv.barriers(); i < barriers + n_barriers; ++i, ++j)
    *i = *j;
  fitness = indiv.fitness();
}

Individual::Individual(Individual&& indiv)
{
  size = indiv.size();
  n_genes = indiv.n_genes();
  genes = new unsigned[2 * n_genes];
  for (unsigned* i = genes, j = indiv.genes(); i < genes + 2 * n_genes; ++i, ++j)
    *i = *j;
  n_barriers = indiv.n_barriers();
  barriers = new unsigned[n_barriers];
  for (unsigned* i = barriers, j = indiv.barriers(); i < barriers + n_barriers; ++i, ++j)
    *i = *j;
  fitness = indiv.fitness();
	indiv.reset();
}

Individual::Individual(unsigned size,
                       unsigned n_genes,
                       unsigned n_barriers,
                       unsigned* genes,
                       unsigned* barriers,
                       float fitness)
{
  this->size = size;
  this->n_genes = n_genes;
  this->genes = new unsigned[2 * n_genes];
  for (unsigned* i = this->genes, j = genes; i < genes + 2 * n_genes; ++i, ++j)
    *i = *j;
  this->n_barriers = n_barriers;
  barriers = new unsigned[n_barriers];
  for (unsigned* i = this->barriers, j = barriers; i < barriers + n_barriers; ++i, ++j)
    *i = *j;
  this->fitness = fitness;
}

Individual& Individual::operator=(const Individual& indiv)
{
	if (this != &indiv)
	{
	  size = indiv.size();
    n_genes = indiv.n_genes();
    delete[] genes;
    genes = new unsigned[2 * n_genes];
    for (unsigned* i = genes, j = indiv.genes(); i < genes + 2 * n_genes; ++i, ++j)
      *i = *j;
    n_barriers = indiv.n_barriers();
    delete[] barriers;
    barriers = new unsigned[n_barriers];
    for (unsigned* i = barriers, j = indiv.barriers(); i < barriers + n_barriers; ++i, ++j)
      *i = *j;
    fitness = indiv.fitness();
	}
	return *this;
}

Individual& Individual::operator=(Individual&& indiv)
{
	if (this != &indiv)
	{
	  size = indiv.size();
    n_genes = indiv.n_genes();
    delete[] genes;
    genes = new unsigned[2 * n_genes];
    for (unsigned* i = genes, j = indiv.genes(); i < genes + 2 * n_genes; ++i, ++j)
      *i = *j;
    n_barriers = indiv.n_barriers();
    delete[] barriers;
    barriers = new unsigned[n_barriers];
    for (unsigned* i = barriers, j = indiv.barriers(); i < barriers + n_barriers; ++i, ++j)
      *i = *j;
    fitness = indiv.fitness();
		indiv.reset();
	}
	return *this;
}

void Individual::set_mutation(double indel_nb, double inv_p)
{
  mean = indel_nb;
  p = inv_p;
}

void Individual::mutate(void)
{
  
  return void;
}

void Individual::reset(void)
{
  size = 0;
  n_genes = 0;
  n_barriers = 0;
  delete genes;
  genes = nullptr;
  delete barriers
  barriers = nullptr;
  fitness = 0;
  return void;
}

Individual::~Individual(void)
{
  delete[] genes;
  delete[] barriers
}
