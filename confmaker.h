#ifndef confmaker_h
#define confmaker_h
//----------------------------------------------------------------------
#include <vector>
#include <algorithm>
#include "mdunit.h"
//----------------------------------------------------------------------
class ConfigurationMaker : public Executor {
private:
  Parameter *param;
public:
  ConfigurationMaker(Parameter *param_) {
    param = param_;
  };
  void Execute(MDUnit *mdu);
};
//----------------------------------------------------------------------
class SimpleConfigurationMaker : public Executor {
private:
  Parameter *param;
public:
  SimpleConfigurationMaker(Parameter *param_) {
    param = param_;
  };
  void Execute(MDUnit *mdu);
};
//----------------------------------------------------------------------
class ExactConfigurationMaker : public Executor {
private:
  Parameter *param;
  std::vector <unsigned long int> vec_delete;
  const double *L;
  int six, siy, siz;
public:
  ExactConfigurationMaker(Parameter *param_, const double *L_);
  void Execute(MDUnit *mdu);
};
//----------------------------------------------------------------------
class ExtractConfigurationMaker : public Executor {
private:
  Parameter *param;
  double density;
public:
  ExtractConfigurationMaker(Parameter *param_) {
    param = param_;
    density = param->GetDoubleDef("Density", 0.7);
  };
  void Execute(MDUnit *mdu);
};
//----------------------------------------------------------------------
#endif
