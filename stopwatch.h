#ifndef stopwatch_h
#define stopwatch_h
//----------------------------------------------------------------------
#include <vector>
#include <fstream>
#include <mpi.h>
//----------------------------------------------------------------------
class StopWatch {
private:
  std::vector<double> data;
  double current_time;
  const char *basename;
  int id;
public:
  StopWatch(int rank, const char* bname) {
    basename = bname;
    id = rank;
  }
  ~StopWatch(void) {
    //if (id == 0)SaveToFile();
    //SaveToFile();
  }
  void Start(void) {
    current_time = Communicator::GetTime();
  }
  void Stop(void) {
    data.push_back(Communicator::GetTime() - current_time);
  }

  void SaveToFile(void) {
    std::vector<double> recvbuf;
    char filename[256];
    sprintf(filename, "%s%05d.dat", basename, id);
    std::ofstream ofs(filename);
    ofs.write((const char *)&data[0], sizeof(double)*data.size());
  }
};
//----------------------------------------------------------------------
#endif
