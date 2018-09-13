#include <iostream>
#include <omp.h>
#include <mpi.h>
#include "mpistream.h"
#include "mdmanager.h"
#include "projectmanager.h"
//----------------------------------------------------------------------
int
main(int argc, char **argv) {
  setvbuf(stdout, NULL, _IOLBF, 0);
  MDManager mdm(argc, argv);
  if (mdm.IsValid()) {
    ProjectManager::GetInstance().ExecuteProject(&mdm);
  } else {
    mout << "Program is aborted." << std::endl;
  }
}
//----------------------------------------------------------------------
