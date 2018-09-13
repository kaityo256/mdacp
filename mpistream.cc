#include <fstream>
#include <omp.h>
#include "mpistream.h"
//----------------------------------------------------------------------
MPIStream mout;
//----------------------------------------------------------------------
MPIStream & MPIStream::operator << (std::ostream & (*pf)(std::ostream &)) {
  os_backup << pf;
  if (0 == rank && omp_get_thread_num() == 0) {
    std::cout << oss.str() << pf;
  }
  oss.str("");
  oss.clear();
  return *this;
}
//----------------------------------------------------------------------
void
MPIStream::SaveToFile(std::string filename) {
  if (0 != rank || omp_get_thread_num() != 0) {
    return;
  }
  std::ofstream ofs(filename.c_str());
  ofs << os_backup.str() << std::endl;
}
//----------------------------------------------------------------------
