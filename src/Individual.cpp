#include "Individual.h"
#include <random>

double Individual::_mean = 5;
double Individual::_p = 0.5;

Individual::Individual(const Individual& indiv)
{
  _size = indiv.get_size();
  _n_genes = indiv.get_n_genes();
  _genes = new unsigned[2 * _n_genes];
  for (unsigned *i = _genes, *j = indiv.get_genes(); i < _genes + 2 * _n_genes; ++i, ++j)
    *i = *j;
  _n_barriers = indiv.get_n_barriers();
  _barriers = new unsigned[_n_barriers];
  for (unsigned *i = _barriers, *j = indiv.get_barriers(); i < _barriers + _n_barriers; ++i, ++j)
    *i = *j;
  _fitness = indiv.get_fitness();
}

Individual::Individual(Individual&& indiv)
{
  _size = indiv.get_size();
  _n_genes = indiv.get_n_genes();
  _genes = new unsigned[2 * _n_genes];
  for (unsigned *i = _genes, *j = indiv.get_genes(); i < _genes + 2 * _n_genes; ++i, ++j)
    *i = *j;
  _n_barriers = indiv.get_n_barriers();
  _barriers = new unsigned[_n_barriers];
  for (unsigned *i = _barriers, *j = indiv.get_barriers(); i < _barriers + _n_barriers; ++i, ++j)
    *i = *j;
  _fitness = indiv.get_fitness();
	indiv.reset();
}

Individual::Individual(unsigned size,
                       unsigned n_genes,
                       unsigned n_barriers,
                       unsigned* genes,
                       unsigned* barriers,
                       float fitness)
{
  _size = size;
  _n_genes = n_genes;
  _genes = new unsigned[2 * _n_genes];
  for (unsigned *i = _genes, *j = genes; i < _genes + 2 * _n_genes; ++i, ++j)
    *i = *j;
  _n_barriers = n_barriers;
  _barriers = new unsigned[_n_barriers];
  for (unsigned *i = _barriers, *j = barriers; i < _barriers + _n_barriers; ++i, ++j)
    *i = *j;
  _fitness = fitness;
}

Individual& Individual::operator=(const Individual& indiv)
{
	if (this != &indiv)
	{
    _size = indiv.get_size();
    _n_genes = indiv.get_n_genes();
    delete[] _genes;
    _genes = new unsigned[2 * _n_genes];
    for (unsigned *i = _genes, *j = indiv.get_genes(); i < _genes + 2 * _n_genes; ++i, ++j)
      *i = *j;
    _n_barriers = indiv.get_n_barriers();
    delete[] _barriers;
    _barriers = new unsigned[_n_barriers];
    for (unsigned *i = _barriers, *j = indiv.get_barriers(); i < _barriers + _n_barriers; ++i, ++j)
      *i = *j;
    _fitness = indiv.get_fitness();
	}
	return *this;
}

Individual& Individual::operator=(Individual&& indiv)
{
	if (this != &indiv)
	{
    _size = indiv.get_size();
    _n_genes = indiv.get_n_genes();
    delete[] _genes;
    _genes = new unsigned[2 * _n_genes];
    for (unsigned *i = _genes, *j = indiv.get_genes(); i < _genes + 2 * _n_genes; ++i, ++j)
      *i = *j;
    _n_barriers = indiv.get_n_barriers();
    delete[] _barriers;
    _barriers = new unsigned[_n_barriers];
    for (unsigned *i = _barriers, *j = indiv.get_barriers(); i < _barriers + _n_barriers; ++i, ++j)
      *i = *j;
    _fitness = indiv.get_fitness();
		indiv.reset();
	}
	return *this;
}

void Individual::set_mutation(double indel_nb, double inv_p)
{
  _mean = indel_nb;
  _p = inv_p;
}

void Individual::mutate(void)
{
  
}

void Individual::reset(void)
{
  _size = 0;
  _n_genes = 0;
  _n_barriers = 0;
  delete _genes;
  _genes = nullptr;
  delete _barriers;
  _barriers = nullptr;
  _fitness = 0;
}

Individual::~Individual(void)
{
  delete[] _genes;
  delete[] _barriers;
}
