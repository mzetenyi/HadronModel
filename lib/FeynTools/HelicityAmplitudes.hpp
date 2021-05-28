#ifndef HELICITYAMPLITUDES_HPP
#define HELICITYAMPLITUDES_HPP

#include <vector>
#include <string>
#include "halfint.hpp"
#include "utils.hpp"
#include "MultiArray.hpp"

class HelicityAmplitudes {
 public:
  HelicityAmplitudes(std::vector<halfint> spins, std::vector< std::vector<std::string> > channels);
  dcomplex value(std::vector<std::string> names, std::vector<halfint> spins) const;
  double magnitude(std::vector<std::string> names, std::vector<halfint> spins) const;
  double phase(std::vector<std::string> names, std::vector<halfint> spins) const;
  //  double relativeError(std::vector<std::string> names) const;
  //  void setRelativeError(std::string name, double val);
  //  MultiArray<double> get_relative_errors() const;
  //  void set_relative_errors(const MultiArray<double>& relerrs);
  //  dcomplex add(std::string name, std::vector<halfint> spins, dcomplex amplitude);
  dcomplex& operator()(std::vector<std::string> names, std::vector<halfint> spins);
  const dcomplex& operator()(std::vector<std::string> names, std::vector<halfint> spins) const;
  //  MultiArray<dcomplex> densityMatrix(uint i, std::vector<std::string> channels) const;
 private:
  std::vector<NamedIndex> NIs;
  MultiArray<dcomplex> amplitudes;
  //  MultiArray<double> relative_errors;
  static std::vector<NamedIndex> get_named_indices(std::vector< std::vector<std::string> > channels);
  static std::vector<idx> get_idx_list(std::vector<NamedIndex> NIs, std::vector<halfint> spins);
  static std::vector<idx> get_idx_list(std::vector<NamedIndex> NIs);
  std::vector<halfint> get_indices(std::vector<std::string> names, std::vector<halfint> spins) const;
  std::vector<halfint> get_indices(std::vector<std::string> names) const;
};





#endif // HELICITYAMPLITUDES_HPP
