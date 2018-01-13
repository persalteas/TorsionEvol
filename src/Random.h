#ifndef RANDOM_H
#define RANDOM_H

#include <random>

class Random {
private:
  std::default_random_engine generator{};
  std::bernoulli_distribution bernoulli;
  std::poisson_distribution<int> poisson;
  std::uniform_int_distribution<int> distribution{};

public:
  Random(void) : Random(0.5, 4){};
  Random(double bernoulli_probability, double poisson_mean) {
    bernoulli = std::bernoulli_distribution(bernoulli_probability);
    poisson = std::poisson_distribution<int>(poisson_mean);
  }
  bool inversion_occurs(void) { return bernoulli(generator); }
  int indels_nb(void) { return poisson(generator); }
};

#endif