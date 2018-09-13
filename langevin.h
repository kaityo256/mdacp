//----------------------------------------------------------------------
#ifndef langevin_h
#define langevin_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Langevin : public Project {
private:
public:
  Langevin(void) {
    ProjectManager::GetInstance().AddProject("Langevin", this);
  };
  void Run(MDManager *mdm);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

