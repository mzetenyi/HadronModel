#include "Config.hpp"

#include <iostream>

using namespace std;

void remove_ws(string& str) {
  str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
}

template <>
int string_to<int>(const string& s) {
  return stoi(s);
}
template <>
long string_to<long>(const string& s) {
  return stol(s);
}
template <>
long long string_to<long long>(const string& s) {
  return stoll(s);
}
template <>
unsigned int string_to<unsigned int>(const string& s) {
  return stoul(s);
}  // stou does not exist!
template <>
unsigned long string_to<unsigned long>(const string& s) {
  return stoul(s);
}
template <>
unsigned long long string_to<unsigned long long>(const string& s) {
  return stoull(s);
}
template <>
float string_to<float>(const string& s) {
  return stof(s);
}
template <>
double string_to<double>(const string& s) {
  return stod(s);
}
template <>
long double string_to<long double>(const string& s) {
  return stold(s);
}

template <>
halfint string_to<halfint>(const string& s) {
  string::size_type sz;
  int i = stoi(s, &sz);
  string rest = s.substr(sz);
  if (rest.size() >= 2 and rest.at(0) == '/' and rest.at(1) == '2') {
    return i * half;
  }
  return 2 * i * half;
}

void Params::add_(const string& key, const string& val) {
  size_t pos;
  if ((pos = key.find('.')) != string::npos) {
    groups[key.substr(0, pos)].add(key.substr(pos + 1), val);
  } else {
    data[key] = val;
  }
}

void Params::add(const string& key, const string& val) {
  if (val.empty()) {
    // check for load(filename):
    if (key.compare(0, 5, "load[") == 0) {
      string filename = key.substr(5);
      filename.erase(filename.size() - 1);
      ifstream in(filename);
      assure(in, filename);
      load(in);
    } else {
      add_(key, val);
    }
  } else {
    if (val.compare(0, 5, "load[") == 0) {
      string filename = val.substr(5);
      filename.erase(filename.size() - 1);
      ifstream in(filename);
      assure(in, filename);
      groups[key].load(in);
    } else {
      add_(key, val);
    }
  }
}

void Params::add(const string& line) {
  size_t eqpos = line.find('=');
  add(line.substr(0, eqpos),
      (eqpos == string::npos) ? "" : line.substr(eqpos + 1));
}

Params& Params::addGroup(const string& key) {
  groups[key] = Params();
  return groups.at(key);
}

void Params::load(uint argc, char** argv) {
  for (uint i(1); i < argc; i++) {
    add(argv[i]);
  }
}

void Params::load(istream& in) {
  string line;
  while (getline(in, line) and (line.find('}') == string::npos)) {
    remove_ws(line);
    if (not line.empty()) {
      size_t pos;
      if ((pos = line.find('{')) != string::npos) {
        Params& newGroup = addGroup(line.substr(0, pos));
        newGroup.load(in);
      } else {
        add(line);
      }
    }
  }
}

const string& Params::get_(const string& key) const {
  return data.at(key);
}

const string& Params::get(const string& key) const {
  size_t pos = key.find('.');
  if (pos == string::npos) {
    return get_(key);
  }
  return getGroup(key.substr(0, pos)).get(key.substr(pos + 1));
}

bool Params::groupExists(const string& key) const {
  return (groups.find(key) != groups.end());
}

bool Params::exists(const string& key) const {
  size_t pos = key.find('.');
  if (pos == string::npos) {
    return (data.find(key) != data.end());
  }
  return groupExists(key.substr(0, pos)) and
         getGroup(key.substr(0, pos)).exists(key.substr(pos + 1));
}

bool Params::hasValue(const string& key) const {
  return (exists(key) and not get(key).empty());
}

const Params& Params::getGroup(const string& key) const {
  return groups.at(key);
}

void Params::list(ostream& out, string head) const {
  for (auto it = data.begin(); it != data.end(); it++) {
    out << head << it->first;
    if (not it->second.empty()) {
      out << " = " << it->second;
    }
    out << endl;
  }
  for (auto it = groups.begin(); it != groups.end(); it++) {
    out << head << it->first << " {" << endl;
    it->second.list(out, head + "  ");
    out << head << "}" << endl;
  }
}

template <>
void Params::add<string>(const string& key, const string& val) {
  data[key] = val;
}

template <>
string Params::get_<string>(const string& key) const {
  return get_(key);
}

template <>
bool Params::get_<bool>(const string& key) const {
  return exists(key);
}

void Config::load(uint argc, char** argv) {
  Params& params = getParams();
  params.load(argc, argv);
}

const string& Config::get(const string& key) {
  return getParams().get(key);
}

bool Config::groupExists(const string& key) {
  return getParams().groupExists(key);
}

bool Config::exists(const string& key) { return getParams().exists(key); }

bool Config::hasValue(const string& key) {
  return getParams().hasValue(key);
}

const Params& Config::getGroup(const string& key) {
  return getParams().getGroup(key);
}

void Config::list(ostream& out, string head) {
  return getParams().list(out, head);
}

Params& Config::getParams() {
  static Params params;
  return params;
}

bool isSet(const string& key) { return Config::exists(key); }

vector<double> getRange(const string& key, double min, double max,
                             double inc) {
  min = getParam<double>(key + "min", min);
  max = getParam<double>(key + "max", max);
  inc = getParam<double>("d" + key, inc);
  vector<double> ret;
  for (double val(min); val<max; val+=inc) {
    ret.push_back(val);
  }
  return ret;
}