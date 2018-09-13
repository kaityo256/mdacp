//----------------------------------------------------------------------
#ifndef meshlist_h
#define meshlist_h
//----------------------------------------------------------------------
#include <vector>
#include "mdconfig.h"
#include "variables.h"
#include "mdrect.h"
//----------------------------------------------------------------------
class MeshList {
private:
  int number_of_pairs;
  double mesh_size_x;
  double mesh_size_y;
  double mesh_size_z;
  int mx, my, mz;
  int particle_position[N];
  int *key_particles;
  int *partner_particles;
  int * mesh_index;
  int * mesh_index2;
  int * mesh_particle_number;
  int sortbuf[N];
  int key_pointer[N];
  int key_pointer2[N];
  int number_of_mesh;
  int number_of_partners[N];
  int sorted_list[PAIRLIST_SIZE];
  int number_of_constructions;
  inline void RegisterPair(int index1, int index2);
  int sort_interval;

  void MakeListMesh(Variables *vars, SimulationInfo *sinfo, MDRect &myrect);
  void MakeMesh(Variables *vars, SimulationInfo *sinfo, MDRect &myrect);
  inline void index2pos(int index, int &ix, int &iy, int &iz);
  inline int pos2index(int ix, int iy, int iz);
  void SearchMesh(int index, Variables *vars, SimulationInfo *sinfo);
  void AppendList(int mx, int my, int mz, std::vector<int> &v);

public:
  MeshList(SimulationInfo *sinfo, MDRect &r);
  ~MeshList(void);

  void ChangeScale(SimulationInfo *sinfo, MDRect &myrect);

  int *GetKeyParticles(void) {return key_particles;}
  int *GetPartnerParticles(void) {return partner_particles;}
  int GetPairNumber(void) {return number_of_pairs;}
  int GetPartnerNumber(int i) {return number_of_partners[i];}
  int *GetSortedList(void) {return sorted_list;}
  int GetKeyPointer(int i) {return key_pointer[i];}
  int* GetKeyPointerP(void) {return key_pointer;}
  int* GetNumberOfPartners(void) {return number_of_partners;}
  void Sort(Variables *vars, SimulationInfo *sinfo, MDRect &myrect);
  int GetNumberOfConstructions(void) {return number_of_constructions;}
  void ClearNumberOfConstructions(void) {number_of_constructions = 0;}

  void MakeList(Variables *vars, SimulationInfo *sinfo, MDRect &myrect);
  void MakeListBruteforce(Variables *vars, SimulationInfo *sinfo, MDRect &myrect);
  void ShowPairs(void);
  void ShowSortedList(Variables *vars);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
