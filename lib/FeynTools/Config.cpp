#include "Config.hpp"
#include <iostream>

void remove_ws(std::string& str) {
  str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
}


template <> int string_to<int>(const std::string& s) { return std::stoi(s); }
template <> long string_to<long>(const std::string& s) { return std::stol(s); }
template <> long long string_to<long long>(const std::string& s) { return std::stoll(s); }
template <> unsigned int string_to<unsigned int>(const std::string& s) { return std::stoul(s); } // std::stou does not exist!
template <> unsigned long string_to<unsigned long>(const std::string& s) { return std::stoul(s); }
template <> unsigned long long string_to<unsigned long long>(const std::string& s) { return std::stoull(s); }
template <> float string_to<float>(const std::string& s) { return std::stof(s); }
template <> double string_to<double>(const std::string& s) { return std::stod(s); }
template <> long double string_to<long double>(const std::string& s) { return std::stold(s); }

template <> halfint string_to<halfint>(const std::string& s) {
  std::string::size_type sz;
  int i = std::stoi(s,&sz);
  std::string rest = s.substr(sz);
  if (rest.size()>=2 and rest.at(0)=='/' and rest.at(1)=='2') { 
    return i*half;
  }
  return 2*i*half;
}

void Params::add_(const std::string& key, const std::string& val) {
  size_t pos;
  if ((pos = key.find('.')) != std::string::npos) {
    groups[key.substr(0,pos)].add(key.substr(pos+1),val);
  } else {
    data[key] = val;
  }
}

void Params::add(const std::string& key, const std::string& val) {
  if (val.empty()) {
    //check for load(filename):
    if (key.compare(0,5,"load[") == 0) {
      std::string filename = key.substr(5);
      filename.erase(filename.size()-1);
      std::ifstream in(filename);
      assure(in,filename);
      load(in);
    } else {
      add_(key,val);
    }
  } else {
    if (val.compare(0,5,"load[") == 0) {
      std::string filename = val.substr(5);
      filename.erase(filename.size()-1);
      std::ifstream in(filename);
      assure(in,filename);
      groups[key].load(in);
    } else {
      add_(key,val);
    }
  }
}

void Params::add(const std::string& line) {
  size_t eqpos = line.find('=');
  add(line.substr(0,eqpos), (eqpos==std::string::npos) ? "" : line.substr(eqpos+1));
}

Params& Params::addGroup(const std::string& key) {
  groups[key] = Params();
  return groups.at(key);
}

void Params::load(uint argc, char** argv) {
  for (uint i(1); i<argc; i++) { add(argv[i]); }
}

void Params::load(std::istream& in) {
  std::string line;
  while ( getline(in,line) and (line.find('}')==std::string::npos) ) {
    remove_ws(line);
    if (not line.empty()) {
      size_t pos;
      if ((pos = line.find('{')) != std::string::npos) {
        Params& newGroup = addGroup(line.substr(0,pos));
        newGroup.load(in);
      } else {
        add(line);
      }
    }
  }
}

const std::string& Params::get_(const std::string& key) const {
  return data.at(key);
}

const std::string& Params::get(const std::string& key) const {
  size_t pos = key.find('.');
  if (pos==std::string::npos) {
    return get_(key);
  }
  return getGroup(key.substr(0,pos)).get(key.substr(pos+1));
}

bool Params::groupExists(const std::string& key) const {
  return (groups.find(key) != groups.end());
}

bool Params::exists(const std::string& key) const {
  size_t pos = key.find('.');
  if (pos==std::string::npos) {
    return (data.find(key) != data.end());
  }
  return groupExists(key.substr(0,pos)) and getGroup(key.substr(0,pos)).exists(key.substr(pos+1));
}

bool Params::hasValue(const std::string& key) const {
  return (exists(key) and not get(key).empty());
}

const Params& Params::getGroup(const std::string& key) const {
  return groups.at(key);
}

void Params::list(std::ostream& out, std::string head) const {
  for (auto it=data.begin(); it!=data.end(); it++) {
    out << head << it->first;
    if (not it->second.empty()) {
      out << " = " << it->second;
    }
    out << std::endl;
  }
  for (auto it=groups.begin(); it!=groups.end(); it++) {
    out << head << it->first << " {" << std::endl;
    it->second.list(out,head+"  ");
    out << head << "}" << std::endl;
  }
}

template <>
void Params::add<std::string>(const std::string& key, const std::string& val) {
  data[key] = val;
}


template <>
std::string Params::get_<std::string>(const std::string& key) const {
  return get_(key);
}

template <>
bool Params::get_<bool>(const std::string& key) const {
  return exists(key);
}


void Config::load(uint argc, char** argv) {
  Params& params = getParams();
  params.load(argc,argv);
}

const std::string& Config::get(const std::string& key) {
  return getParams().get(key);
}

bool Config::groupExists(const std::string& key) {
  return getParams().groupExists(key);
}

bool Config::exists(const std::string& key) {
  return getParams().exists(key);
}

bool Config::hasValue(const std::string& key) {
  return getParams().hasValue(key);
}

const Params& Config::getGroup(const std::string& key) {
  return getParams().getGroup(key);
}


void Config::list(std::ostream& out, std::string head) {
  return getParams().list(out,head);
}


Params& Config::getParams() {
  static Params params;
  return params;
}


bool isSet(const std::string& key) {
  return Config::exists(key);
}
