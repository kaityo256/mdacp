#ifndef parameter_h
#define parameter_h
//----------------------------------------------------------------------
#include <iostream>
#include <string>
#include <map>
//----------------------------------------------------------------------
typedef std::map<std::string, std::string> ptype;

class Parameter {
private:
  ptype params;
  bool valid;
  void ReadFromStream(std::istream &is);
public:
  Parameter(void) {valid = true;}
  Parameter(const char *filename);
  Parameter(std::istream &is);
  void LoadFromFile(const char *filename);
  bool IsValid(void) {return valid;}
  bool Contains(std::string key);
  int GetInteger(std::string key);
  int GetIntegerDef(std::string key, int value);
  double GetDouble(std::string key);
  double GetDoubleDef(std::string key, double value);
  std::string GetString(std::string key);
  std::string GetStringDef(std::string key, std::string value);
  bool GetBoolean(std::string key);
  bool GetBooleanDef(std::string key, bool value);
  void ShowAll(void);
};
//----------------------------------------------------------------------
#endif
