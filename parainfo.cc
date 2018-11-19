#include <stdio.h>
#include <vector>
#include <mpi.h>
#include "mpistream.h"
#include "mdconfig.h"
#include "parainfo.h"

const int ParaInfo::diff[MAX_DIR][D] = {
  { -1, 0, 0},
  {1, 0, 0},
  {0, -1, 0},
  {0, 1, 0},
  {0, 0, -1},
  {0, 0, 1}
};
//----------------------------------------------------------------------
ParaInfo::ParaInfo(int np_, int nt_, Parameter &param):
  num_procs(np_), num_threads(nt_) {
  valid = true;
  grid.resize(num_procs * num_threads);
  grid_x.reserve(num_procs * num_threads);
  grid_y.reserve(num_procs * num_threads);
  grid_z.reserve(num_procs * num_threads);
  int g3[3] = {};
  MPI_Dims_create(num_procs, 3, g3);
  int gx = g3[X];
  int gy = g3[Y];
  int gz = g3[Z];
  int l3[3] = {};
  MPI_Dims_create(num_threads, 3, l3);
  // We reversed order so that gx*lx ~ gy*ly ~ gz*lz
  int lx = l3[Z];
  int ly = l3[Y];
  int lz = l3[X];
  mpi_grid_size[X] = gx;
  mpi_grid_size[Y] = gy;
  mpi_grid_size[Z] = gz;
  openmp_grid_size[X] = lx;
  openmp_grid_size[Y] = ly;
  openmp_grid_size[Z] = lz;
  ReadGlobalGrid(param);
  ReadLocalGrid(param);
  gx = mpi_grid_size[X];
  gy = mpi_grid_size[Y];
  gz = mpi_grid_size[Z];
  lx = openmp_grid_size[X];
  ly = openmp_grid_size[Y];
  lz = openmp_grid_size[Z];

  mout << "# Procs = (" << gx << "x" << gy << "x" << gz << ")" << std::endl;
  mout << "# Threads = (" << lx << "x" << ly << "x" << lz << ")" << std::endl;
  if (num_procs != gx * gy * gz) {
    show_error("Invalid Grid Numbers");
    valid = false;
    mout << "NumProcs = " << num_procs << std::endl;
    mout << "GridX=" << gx;
    mout << " GridY=" << gy;
    mout << " GridZ=" << gz;
    mout << " Total=" << gx*gy*gz << std::endl;
    return;
  }
  if (num_threads != lx * ly * lz) {
    show_error("Invalid LocalGrid Numbers");
    valid = false;
    mout << "NumThreads = " << num_threads << std::endl;
    mout << "LocalGridX=" << lx;
    mout << " LocalGridY=" << ly;
    mout << " LocalGridZ=" << lz;
    mout << " Total=" << lx*ly*lz << std::endl;
    return;
  }
  SetGrid();
}
//----------------------------------------------------------------------
void
ParaInfo::ReadGlobalGrid(Parameter &param) {
  if (!param.Contains("GridX")) {
    return;
  }
  if (!param.Contains("GridY")) {
    show_error("GridX is set, but GridY is not set");
    valid = false;
    return;
  }
  if (!param.Contains("GridZ")) {
    show_warning("GridX is set, but GridZ is not set");
    valid = false;
    return;
  }
  int gx = param.GetInteger("GridX");
  int gy = param.GetInteger("GridY");
  int gz = param.GetInteger("GridZ");
  mpi_grid_size[X] = gx;
  mpi_grid_size[Y] = gy;
  mpi_grid_size[Z] = gz;
}
//----------------------------------------------------------------------
void
ParaInfo::ReadLocalGrid(Parameter &param) {
  if (!param.Contains("LocalGridX")) {
    return;
  }
  if (!param.Contains("LocalGridY")) {
    show_error("LocalGridX is set, but LocalGridY is not set");
    valid = false;
    return;
  }
  if (!param.Contains("LocalGridZ")) {
    show_error("LocalGridX is set, but LocalGridZ is not set");
    valid = false;
    return;
  }

  const int lx = param.GetInteger("LocalGridX");
  const int ly = param.GetInteger("LocalGridY");
  const int lz = param.GetInteger("LocalGridZ");

  openmp_grid_size[X] = lx;
  openmp_grid_size[Y] = ly;
  openmp_grid_size[Z] = lz;
}
//----------------------------------------------------------------------
void
ParaInfo::SetGrid(void) {
  const int gx = mpi_grid_size[X];
  const int gy = mpi_grid_size[Y];
  const int lx = openmp_grid_size[X];
  const int ly = openmp_grid_size[Y];
  const int lz = openmp_grid_size[Z];
  const int glx = gx * lx;
  const int gly = gy * ly;
  for (int i = 0; i < num_procs; i++) {
    const int gix = i % gx;
    const int giy = (i / gx) % gy;
    const int giz = (i / gx / gy);
    for (int j = 0; j < num_threads; j++) {
      const int lix = j % lx;
      const int liy = (j / lx) % ly;
      const int liz = (j / lx / ly);
      const int ix = gix * lx + lix;
      const int iy = giy * ly + liy;
      const int iz = giz * lz + liz;
      int index = ix + iy * glx + iz * glx * gly;
      grid[index] = j + i * num_threads;
      grid_x.push_back(ix);
      grid_y.push_back(iy);
      grid_z.push_back(iz);
    }
  }
}
//----------------------------------------------------------------------
void
ParaInfo::GetGridSize(int grid_size[D]) {
  grid_size[X] = mpi_grid_size[X] * openmp_grid_size[X];
  grid_size[Y] = mpi_grid_size[Y] * openmp_grid_size[Y];
  grid_size[Z] = mpi_grid_size[Z] * openmp_grid_size[Z];
}
//----------------------------------------------------------------------
void
ParaInfo::GetGridPosition(int id, int grid_position[D]) {
  grid_position[X] = grid_x[id];
  grid_position[Y] = grid_y[id];
  grid_position[Z] = grid_z[id];
}
//----------------------------------------------------------------------
void
ParaInfo::ShowGrid(void) {
  printf("# %d %d %d\n", mpi_grid_size[X], mpi_grid_size[Y], mpi_grid_size[Z]);
  printf("# %d %d %d\n", openmp_grid_size[X], openmp_grid_size[Y], openmp_grid_size[Z]);
  const int glx = openmp_grid_size[X] * mpi_grid_size[X];
  for (unsigned int i = 0; i < grid.size(); i++) {
    if (i % glx == 0) printf("\n");
    printf("%d ", grid[i]);
  }
  printf("\n-----\n");
}
//----------------------------------------------------------------------
int
ParaInfo::Pos2ID(int pos[D]) {
  const int gx = mpi_grid_size[X];
  const int gy = mpi_grid_size[Y];
  const int lx = openmp_grid_size[X];
  const int ly = openmp_grid_size[Y];
  const int glx = gx * lx;
  const int gly = gy * ly;
  int index = pos[X] + pos[Y] * glx + pos[Z] * glx * gly;
  return grid[index];
}
//----------------------------------------------------------------------
int
ParaInfo::GetNeighborID(int id, int dir) {
  int pos[D];
  GetGridPosition(id, pos);
  const int gx = mpi_grid_size[X];
  const int gy = mpi_grid_size[Y];
  const int gz = mpi_grid_size[Z];
  const int lx = openmp_grid_size[X];
  const int ly = openmp_grid_size[Y];
  const int lz = openmp_grid_size[Z];
  const int glx = gx * lx;
  const int gly = gy * ly;
  const int glz = gz * lz;
  pos[X] += diff[dir][X];
  pos[Y] += diff[dir][Y];
  pos[Z] += diff[dir][Z];

  if (pos[X] < 0) {
    pos[X] += glx;
  } else if (pos[X] >= glx) {
    pos[X] -= glx;
  }
  if (pos[Y] < 0) {
    pos[Y] += gly;
  } else if (pos[Y] >= gly) {
    pos[Y] -= gly;
  }
  if (pos[Z] < 0) {
    pos[Z] += glz;
  } else if (pos[Z] >= glz) {
    pos[Z] -= glz;
  }
  return Pos2ID(pos);
}
//----------------------------------------------------------------------
int
ParaInfo::GetNeighborRank(int rank, int dir) {
  const int gx = mpi_grid_size[X];
  const int gy = mpi_grid_size[Y];
  const int gz = mpi_grid_size[Z];
  const int mx = rank % gx;
  const int my = (rank / gx) % gy;
  const int mz = (rank / gx / gy);
  int mx2 = mx + diff[dir][X];
  int my2 = my + diff[dir][Y];
  int mz2 = mz + diff[dir][Z];
  if (mx2 < 0) {
    mx2 += gx;
  } else if (mx2 >= gx) {
    mx2 -= gx;
  }
  if (my2 < 0) {
    my2 += gy;
  } else if (my2 >= gy) {
    my2 -= gy;
  }
  if (mz2 < 0) {
    mz2 += gz;
  } else if (mz2 >= gz) {
    mz2 -= gz;
  }
  const int rank2 = mx2 + my2 * gx + mz2 * gx * gy;
  return rank2;
}
//----------------------------------------------------------------------
bool
ParaInfo::IsOverBoundary(const int id, const int dir) {
  int pos[D];
  GetGridPosition(id, pos);
  const int gx = mpi_grid_size[X];
  const int gy = mpi_grid_size[Y];
  const int gz = mpi_grid_size[Z];
  const int lx = openmp_grid_size[X];
  const int ly = openmp_grid_size[Y];
  const int lz = openmp_grid_size[Z];
  const int glx = gx * lx;
  const int gly = gy * ly;
  const int glz = gz * lz;
  pos[X] += diff[dir][X];
  pos[Y] += diff[dir][Y];
  pos[Z] += diff[dir][Z];
  if (pos[X] < 0 || pos[X] >= glx) return true;
  if (pos[Y] < 0 || pos[Y] >= gly) return true;
  if (pos[Z] < 0 || pos[Z] >= glz) return true;
  return false;
}
//----------------------------------------------------------------------
