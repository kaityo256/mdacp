//----------------------------------------------------------------------
#ifndef benchmark_h
#define benchmark_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Benchmark : public Project {
private:
public:
  Benchmark(void) {
    ProjectManager::GetInstance().AddProject("Benchmark", this);
  };
  void Run(MDManager *mdm);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

