//----------------------------------------------------------------------
#include <iostream>
#include <map>
#include <fstream>
#include <stdlib.h>
#include "mpistream.h"
#include "parameter.h"
//----------------------------------------------------------------------
Parameter::Parameter(const char *filename) {
  valid = true;
  std::ifstream is(filename);
  if (is.fail()) {
    mout << "Could not open input file " << filename << std::endl;
    valid = false;
  }
  ReadFromStream(is);
}
//----------------------------------------------------------------------
Parameter::Parameter(std::istream &is) {
  valid = true;
  ReadFromStream(is);
}
//----------------------------------------------------------------------
void
Parameter::LoadFromFile(const char* filename) {
  std::ifstream is(filename);
  if (is.fail()) {
    mout << "Could not open input file " << filename << std::endl;
    valid = false;
  }
  ReadFromStream(is);
}
//----------------------------------------------------------------------
void
Parameter::ReadFromStream(std::istream &is) {
  std::string line;
  while (getline(is, line)) {
    size_t index = line.find("=");
    if (std::string::npos != index) {
      std::string key = line.substr(0, index);
      std::string value = line.substr(index + 1, line.length());
      params.insert(ptype::value_type(key, value));
    }
  }
}
//----------------------------------------------------------------------
bool
Parameter::Contains(std::string key) {
  if (params.find(key) == params.end()) {
    return false;
  } else {
    return true;
  }
}
//----------------------------------------------------------------------
std::string
Parameter::GetString(std::string key) {
  return GetStringDef(key, "");
}
//----------------------------------------------------------------------
std::string
Parameter::GetStringDef(std::string key, std::string value) {
  if (!Contains(key)) {
    return value;
  }
  return params[key];
}
//----------------------------------------------------------------------
double
Parameter::GetDouble(std::string key) {
  return GetDoubleDef(key, 0.0);
}
//----------------------------------------------------------------------
double
Parameter::GetDoubleDef(std::string key, double value) {
  if (!Contains(key)) {
    return value;
  }
  return atof(params[key].c_str());
}
//----------------------------------------------------------------------
int
Parameter::GetInteger(std::string key) {
  return GetIntegerDef(key, 0);
}
//----------------------------------------------------------------------
int
Parameter::GetIntegerDef(std::string key, int value) {
  if (!Contains(key)) {
    return value;
  }
  return atoi(params[key].c_str());
}
//----------------------------------------------------------------------
bool
Parameter::GetBoolean(std::string key) {
  return GetBooleanDef(key, false);
}
//----------------------------------------------------------------------
bool
Parameter::GetBooleanDef(std::string key, bool value) {
  if (!Contains(key)) {
    return value;
  }
  if ("yes" == params[key] || "Yes" == params[key]) {
    return true;
  } else {
    return false;
  }
}
//----------------------------------------------------------------------
void
Parameter::ShowAll(void) {
  ptype::iterator it = params.begin();
  while (it != params.end()) {
    mout << (*it).first << ":" << (*it).second << std::endl;
    ++it;
  }
}
//----------------------------------------------------------------------
