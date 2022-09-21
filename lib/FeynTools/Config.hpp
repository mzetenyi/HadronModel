#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "require.h"
#include "halfint.hpp"

void remove_ws(std::string& str);

/**
   Templated conversion functions from string to numeric types.
*/
template <typename T>
T string_to(const std::string&);

template <> int string_to<int>(const std::string& s);
template <> long string_to<long>(const std::string& s);
template <> long long string_to<long long>(const std::string& s);
template <> unsigned int string_to<unsigned int>(const std::string& s);
template <> unsigned long string_to<unsigned long>(const std::string& s);
template <> unsigned long long string_to<unsigned long long>(const std::string& s);
template <> float string_to<float>(const std::string& s);
template <> double string_to<double>(const std::string& s);
template <> long double string_to<long double>(const std::string& s);
template <> halfint string_to<halfint>(const std::string& s);

class Params {
public:
  void add(const std::string&, const std::string&);

  void add(const std::string&);

  template <typename T> 
  void add(const std::string& key, const T& val) {
    add(key,std::to_string(val));
  }

  Params& addGroup(const std::string& key);

  void load(uint argc, char** argv);
  void load(std::istream& in);

  const std::string& get(const std::string&) const;

  template <typename T>
  T get(const std::string& key) const {
    size_t pos = key.find('.');
    if (pos==std::string::npos) {
      return get_<T>(key);
    }
    return getGroup(key.substr(0,pos)).get<T>(key.substr(pos+1));
  }

  bool groupExists(const std::string&) const;

  bool exists(const std::string&) const;

  bool hasValue(const std::string&) const;

  const Params& getGroup(const std::string& key) const;

  void list(std::ostream&, std::string head="") const;
private:
  std::map<std::string,std::string> data;
  std::map<std::string,Params> groups;

  void add_(const std::string&, const std::string&);
  const std::string& get_(const std::string&) const;

  template <typename T>
  T get_(const std::string& key) const {
    return string_to<T>(get_(key));
  }

};

template <>
std::string Params::get_<std::string>(const std::string& key) const;

template <>
bool Params::get_<bool>(const std::string& key) const;


class Config {
public:
  static void load(uint argc, char** argv);

  static const std::string& get(const std::string&);

  template <typename T>
  static T get(const std::string& key) {
    try {
      return getParams().get<T>(key);
    } catch(...) {
      std::cerr << "problem getting parameter " << key << std::endl;
      exit(0);
    }
  }

  template <typename T>
  static void add(const std::string& key, const T& val) {
    getParams().add(key,val);
  }

  static const Params& getGroup(const std::string& key);

  static bool groupExists(const std::string&);
  static bool exists(const std::string&);
  static bool hasValue(const std::string&);

  static void list(std::ostream&, std::string head="");

private:
  static Params& getParams();
};

bool isSet(const std::string&);

template <typename T>
T getParam(const std::string& key, T defaultVal) {
  if (isSet(key)) return Config::get<T>(key);
  Config::add(key,defaultVal);
  return defaultVal;
}

template <typename T>
T getParam(const std::string& key) {
  return Config::get<T>(key);
}

/// @brief Obtain a range of values to be used in range based for loop.
/// @param key a string used in Config to find min, max values and increment.
/// @param min default min value for range
/// @param max default max value for range
/// @param inc default increment for range
/// @return the generated range of values in a vector.
/// 
/// Searches for the strings "<key>min", "<key>max" and "d<key>" in Config parameters to obtain the 
/// min, max and increment values that are used to generate the values in the range according to
///    
///    for (double val(min), val<max, val+=inc) {
///         ret.push_back(val);
///    }
std::vector<double> getRange(const std::string& key, double min=0, double max=0, double inc=0);

#endif // CONFIG_HPP

