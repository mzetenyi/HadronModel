#include "HelicityAmplitudes.hpp"

using namespace std;

HelicityAmplitudes::HelicityAmplitudes(vector<halfint> spins, vector< vector<string> > channels) :
  NIs(get_named_indices(channels)),
  amplitudes(get_idx_list(NIs,spins)) {
  //  relative_errors(get_idx_list(NIs)) {
  //  cerr << "HelicityAmplitudes::constructor" << endl;
  //  cerr << "   spins: ";
  //  for (auto spin : spins) cerr << spin << " ";
  //  cerr << endl;
  //  cerr << "   channels: ";
  //  for (auto channel : channels) cerr << channel << " ";
  //  cerr << endl;
}

dcomplex HelicityAmplitudes::value(vector<string> names, vector<halfint> spins) const {
  return amplitudes(get_indices(names,spins));
}

double HelicityAmplitudes::magnitude(vector<string> names, vector<halfint> spins) const {
  return abs(amplitudes(get_indices(names,spins)));
}

double HelicityAmplitudes::phase(vector<string> names, vector<halfint> spins) const {
  return arg(amplitudes(get_indices(names,spins)));
}

dcomplex& HelicityAmplitudes::operator()(vector<string> names, vector<halfint> spins) {
  return amplitudes(get_indices(names,spins));
}

const dcomplex& HelicityAmplitudes::operator()(vector<string> names, vector<halfint> spins) const {
  return amplitudes(get_indices(names,spins));
}

vector<NamedIndex> HelicityAmplitudes::get_named_indices(vector< vector<string> > channels) {
  vector<NamedIndex> NIs;
  for (auto namelist : channels) {
    NIs.push_back(NamedIndex(namelist));
  }
  return NIs;
}

vector<idx> HelicityAmplitudes::get_idx_list(vector<NamedIndex> NIs, vector<halfint> spins) {
  vector<idx> ret;
  for (auto NI : NIs) {
    ret.push_back(NI.getIdx());
  }
  for (auto spin : spins) {
    ret.push_back(idx_spin(spin));
  }
  return ret;
}

vector<idx> HelicityAmplitudes::get_idx_list(vector<NamedIndex> NIs) {
  vector<idx> ret;
  for (auto NI : NIs) {
    ret.push_back(NI.getIdx());
  }
  return ret;
}

vector<halfint> HelicityAmplitudes::get_indices(vector<string> names, vector<halfint> spins) const {
  vector<halfint> indices;
  if (names.size() != NIs.size()) {
    cerr << "Incorrect number of chanel names, [ ";
    for (string name : names) cerr << name << " ";
    cerr << "] given for " << NIs.size() << "names";
    exit(0);
  }
  for (uint i(0); i<NIs.size(); i++) {
    indices.push_back(halfint(NIs[i](names[i])));
  }
  indices.insert(indices.end(),spins.begin(),spins.end());
  return indices;
}

vector<halfint> HelicityAmplitudes::get_indices(vector<string> names) const {
  vector<halfint> indices;
  if (names.size() != NIs.size()) {
    cerr << "Incorrect number of chanel names, [ ";
    for (string name : names) cerr << name << " ";
    cerr << "] given for " << NIs.size() << "names";
    exit(0);
  }
  for (uint i(0); i<NIs.size(); i++) {
    indices.push_back(halfint(NIs[i](names[i])));
  }
  return indices;
}

