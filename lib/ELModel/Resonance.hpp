#ifndef RESONANCE_HPP
#define RESONANCE_HPP

#include <vector>
#include <string>
#include "halfint.hpp"

std::vector<std::string> getResonances();

std::vector<int> getIndexRange(std::string resonance, int i);

std::vector<halfint> getHelicityRange(std::string resonance);

#endif // RESONANCE_HPP