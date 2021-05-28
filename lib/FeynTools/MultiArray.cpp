#include "MultiArray.hpp"
#include <cstdlib>
#include <stdexcept>
#include <iomanip>

idx::idx() : _min(0), _max(0), _delta(0), _dim(1) {}
idx::idx(int _dim) : _min(0), _max(_dim-1), _delta(1), _dim(_dim) {
  // cerr << "idx::constr(int) - _delta = " << _delta << endl;
  // PR(_min); PR(_max); PR(_dim);
}
idx::idx(halfint _min, halfint _max, halfint delta) : _min(_min), _max(_max), _delta(delta) {
  // cerr << "idx::constr - _delta = " << _delta << "  " << delta << endl;
  // PR(_min); PR(_max);
  assert(twice(_max-_min)%twice(_delta) == 0);
  _dim = twice(_max-_min)/twice(_delta) + 1;
  /*
  cerr << "idx::constr: ";
  PL(_min); PL(_max); PL(_delta); PR(_dim);
  PL(half); PR(3*half);
  //*/
}
halfint idx::min() const { return _min; }
halfint idx::max() const { return _max; }
halfint idx::delta() const { return _delta; }
int idx::dim() const { return _dim; }
int idx::operator()(halfint h) const {
  // cerr << "idx::operator()" << endl;
  // cerr << "_min = " << _min << endl;
  // cerr << "_max = " << _max << endl;
  // cerr << "_delta = " << _delta << endl;
  // cerr << "_dim = " << _dim << endl;
  // cerr << "h = " << h << endl;
  if (!(match(h))) {
    cerr << "idx::operator()(halfint h) const" << endl;
    PL(h); PL(_min); PL(_max); PL(_delta); PR(_dim);
    throw out_of_range("invalid argument in idx::operator()(halfint h)");
  }
  //  assert(match(h));
  return twice(h-_min)/twice(_delta);
}

bool idx::match(halfint h) const {
  return ((_min<=h) and (h<=_max) and (twice(h-_min)%twice(_delta) == 0));
}

bool operator==(const idx& lhs, const idx& rhs) {
  return (lhs._min==rhs._min and lhs._max==rhs._max and lhs._delta==rhs._delta);
}

bool operator!=(const idx& lhs, const idx& rhs) {
  return !(lhs == rhs);
}

ostream& operator<<(ostream& out, const idx& IDX) {
  out << "idx object:" << endl;
  out << "  _min =   " << IDX._min << endl;
  out << "  _max =   " << IDX._max << endl;
  out << "  _delta = " << IDX._delta << endl;
  out << "  _dim =   " << IDX._dim << endl;
  return out;
}

const idx idx::idx_null(0);
const idx idx::idx_0(_0,_0,_1);
const idx idx::idx_lor(_0,_3,_1);
const idx idx::idx_s1(-_1,_1,_1);
const idx idx::idx_s10(-_1,_1,_2);
const idx idx::idx_s1h(-half,half,_1);
const idx idx::idx_s3h(-3*half,3*half,_1);
const idx idx::idx_s5h(-5*half,5*half,_1);

const idx& idx_null(idx::idx_null);
const idx& idx_0(idx::idx_0);
const idx& idx_lor(idx::idx_lor);
const idx& idx_s1(idx::idx_s1);
const idx& idx_s10(idx::idx_s10);
const idx& idx_s1h(idx::idx_s1h);
const idx& idx_s3h(idx::idx_s3h);
const idx& idx_s5h(idx::idx_s5h);



const idx& idx_spin(halfint spin) {
  if (spin==half) {
    return idx_s1h;
  } else if (spin==3*half) {
    return idx_s3h;
  } else if (spin==5*half) {
    return idx_s5h;
  } else if (spin==_1) {
    return idx_s1;
  } 
  cerr << "Unimplemented spin (" << spin << ") in idx_spin(halfint spin)" << endl;
  exit(0);
  // never get here:
  return idx_null;
}


//class NamedIndex {
//public:
NamedIndex::NamedIndex(const vector<string>& names) :
  names(names) {
  for (uint i(0); i<names.size(); i++) {
    name_to_int.insert(pair<string,int>(names[i],i));
    int_to_name.insert(pair<int,string>(i,names[i]));
  }
}


halfint NamedIndex::operator()(const string& name) const {
  auto it = name_to_int.find(name);
  if (it == name_to_int.end()) {
    throw out_of_range("invalid argument in NamedIndex::operator()(const string& name) - key \"" + name + "\" does not exist");
  }
  return halfint(it->second);
}

string NamedIndex::operator()(int i) const {
  auto it = int_to_name.find(i);
  if (it == int_to_name.end()) {
    throw out_of_range("invalid argument in NamedIndex::operator()(int i) - key \"" + to_string(i) + "\" does not exist");
  }
  return it->second;
}

vector<string> NamedIndex::getNames() const {
  return names;
}

idx NamedIndex::getIdx() const {
  return idx(names.size());
}

vector<idx> NamedIndex::getIdxList(const vector<NamedIndex> NIs) {
  vector<idx> ids;
  for (auto NI : NIs) {
    ids.push_back(NI.getIdx());
  }
  return ids; 
}

ostream& operator<<(ostream& out, const NamedIndex& NI) {
  out << "NamedIndex{ ";
  for (auto name : NI.names) out << name << ", ";
  out << "}" << endl;
  out << "name_to_int map:" << endl;
  for (const auto& item : NI.name_to_int) {
    out << setw(15) << item.first << " -> " << setw(5) << item.second << endl;
  }
  out << "int_to_name map:" << endl;
  for (const auto& item : NI.int_to_name) {
    out << setw(5) << item.first << " -> " << setw(15) << item.second << endl;
  }
  return out;
}

