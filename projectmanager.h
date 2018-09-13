#ifndef projectmanager_h
#define projectmanager_h
//----------------------------------------------------------------------
#include <map>
#include "mpistream.h"
#include "communicator.h"
#include "confmaker.h"
#include "mdmanager.h"
#include "parameter.h"
//----------------------------------------------------------------------
class Project {
public:
  virtual void Run(MDManager *) {}
};
//----------------------------------------------------------------------
typedef std::map<std::string, Project*> project_map;
//----------------------------------------------------------------------
class ProjectManager {
private:
  project_map projects;
  ProjectManager(void) {};
  Project * FindProject (const char *str);
public:
  static ProjectManager & GetInstance(void);
  void AddProject(const char *str, Project *p);
  void ExecuteProject(MDManager *mdm);
};
//----------------------------------------------------------------------
#endif

