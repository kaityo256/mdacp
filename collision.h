//----------------------------------------------------------------------
#ifndef collision_h
#define collision_h
//----------------------------------------------------------------------
#include "projectmanager.h"
//----------------------------------------------------------------------
class Collision : public Project {
private:
public:
  Collision(void) {
    ProjectManager::GetInstance().AddProject("Collision", this);
  };
  void Run(MDManager *mdm);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

