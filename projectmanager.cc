#include <iostream>
#include "projectmanager.h"
#include "mpistream.h"
//----------------------------------------------------------------------
void
ProjectManager::AddProject(const char*str, Project *p) {
  if (FindProject(str) != NULL) {
    std::cout << "# Mode: " << str << " is already defined." << std::endl;
  }
  projects.insert(std::pair<std::string, Project*>(str, p));
}
//----------------------------------------------------------------------
Project *
ProjectManager::FindProject(const char *str) {
  project_map::iterator i;
  i = projects.find(str);
  if (i != projects.end()) {
    Project * p = i->second;
    return p;
  } else {
    return NULL;
  }
}
//----------------------------------------------------------------------
void
ProjectManager::ExecuteProject(MDManager *mdm) {
  Parameter *param = mdm->GetParameter();
  if (!param->Contains("Mode")) {
    mout << "Mode is not defined in inputfile." << std::endl;
    return;
  }

  std::string mode = param->GetString("Mode");
  mout << "# Mode = " << mode << std::endl;

  Project *p = FindProject(mode.c_str());
  if (NULL != p) {
    p->Run(mdm);
  } else {
    mout << "Mode " << mode << " is not found." << std::endl;
  }
}
//----------------------------------------------------------------------
ProjectManager & ProjectManager::GetInstance(void) {
  static ProjectManager pmanager;
  return pmanager;
}
//----------------------------------------------------------------------
