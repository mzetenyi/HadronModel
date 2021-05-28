#ifndef RANDOMNUMBERGENERATOR_HPP
#define RANDOMNUMBERGENERATOR_HPP

#include <random>

class RandomNumberGenerator {
public:
  RandomNumberGenerator(unsigned long int seed = 0) :
    generator(seed),
    distribution(0.,1.) {
  }

  ~RandomNumberGenerator() {
  }

  double operator()() {
    return distribution(generator);
  }

  double operator()(double min, double max) {
    return min + (max-min)*distribution(generator);
  }

  void setSeed(unsigned long int seed) {
    generator.seed(seed);
  }

  static RandomNumberGenerator rn;
private:
  std::mt19937_64 generator;
  std::uniform_real_distribution<double> distribution;
};

extern RandomNumberGenerator& rn;

#endif
