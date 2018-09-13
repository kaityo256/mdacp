//----------------------------------------------------------------------
#ifndef mpistream_h
#define mpistream_h
//----------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <string>
//----------------------------------------------------------------------
class MPIStream {
private:
  int rank;
  std::ostringstream oss;
  std::ostringstream os_backup;
public:
  MPIStream(void) {
    rank = 0;
  };
  void SetRank(int r) {rank = r;}
  template <class T>
  MPIStream & operator << (const T&a) {
    oss << a;
    os_backup << a;
    return *this;
  }
  MPIStream & operator << (std::ostream & (*pf)(std::ostream &));
  void SaveToFile(std::string filename);
};
//----------------------------------------------------------------------
extern MPIStream mout;
#endif
