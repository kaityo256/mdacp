//----------------------------------------------------------------------
#ifndef parainfo_h
#define parainfo_h
//----------------------------------------------------------------------
#include <vector>
#include "mdconfig.h"
#include "parameter.h"
//----------------------------------------------------------------------
class ParaInfo {
private:
  std::vector<int> grid;
  std::vector<int> grid_x;
  std::vector<int> grid_y;
  std::vector<int> grid_z;
  int mpi_grid_size[D];
  int openmp_grid_size[D];
  const int num_procs;
  const int num_threads;
  void SetGrid(void);
  void ReadGlobalGrid(Parameter &param);
  void ReadLocalGrid(Parameter &param);
  bool valid;
  static const int diff[MAX_DIR][D];
public:
  ParaInfo(int np, int nt, Parameter &param);
  void ShowGrid(void);
  bool IsValid(void) {return valid;}
  void GetGridSize(int grid_size[D]);
  void GetGridPosition(int id, int grid_position[D]);
  int GetNeighborID(int id, int dir);
  int GetNeighborRank(int id, int dir);
  int Pos2ID(int pos[D]);
  bool IsOverBoundary(const int id, const int dir);
};
//----------------------------------------------------------------------
#endif
