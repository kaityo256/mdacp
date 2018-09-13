//----------------------------------------------------------------------
#ifndef mdmanager_h
#define mdmanager_h
#include <vector>
#include "mdunit.h"
#include "parainfo.h"
#include "simulationinfo.h"
#include "parameter.h"
//----------------------------------------------------------------------
class MDManager {
private:
  int num_threads;
  int num_procs;
  int rank;
  ParaInfo *pinfo;
  SimulationInfo *sinfo;
  std::vector<MDUnit *> mdv;
  Parameter param;
  // For Grid Management
  std::vector<int> mpi_grid;
  std::vector<int> grid;
  std::vector<int> grid_x;
  std::vector<int> grid_y;
  std::vector<int> grid_z;
  int mpi_grid_size[D];
  int openmp_grid_size[D];
  bool IsMyUnit(int id);
  int GetLocalID(int id);
  bool IsPairListExpired(void);
  double s_time;
public:
  MDManager(int &argc, char** &argv);
  ~MDManager(void);
  bool IsValid(void);
  Parameter *GetParameter(void) {return &param;}
  //Getter, Setter
  MDUnit * GetMDUnit(int index) {return mdv[index];}
  int GetRank(void) {return rank;}
  int GetTotalThreads(void) {return num_threads;}
  int GetTotalUnits(void) {return num_threads * num_procs;}
  int GetTotalProcs(void) {return num_procs;}
  void GetGridSize(int g[D]) {pinfo->GetGridSize(g);}
  double * GetSystemSize(void) {return sinfo->L;}
  double GetTimeStep(void) {return sinfo->TimeStep;}


  //For MDUnit
  void Calculate(void);
  void CalculateForce(void);
  void CalculateNoseHoover(void);
  void SetSeed(const int seed);
  void CalculateLangevin(void);
  void SetInitialVelocity(double v0);
  void SaveConfiguration(void);
  void SaveAsCdviewSequential(void);
  void SaveAsCdview(const char *filename);
  void SaveAsDumpFile(const char *filename);
  void SaveAsDumpFileSequential(void);
  void SendParticlesSub(const int dir);
  void SendParticles(void);
  void MakePairList(void);
  void SendBorderParticles(void);
  void SendBorderParticlesSub(const int dir);
  void ExecuteAll(Executor *ex);
  void ExecuteAllSerial(Executor *ex);

  //For Observe
  double GetSimulationTime(void) {return s_time;}
  double ObserveDouble(DoubleObserver *obs);
  unsigned long int GetTotalParticleNumber(void);
  double KineticEnergy(void);
  double PotentialEnergy(void);
  double TotalEnergy(void);
  double Temperature(void);
  double ConfigurationTemperature(void);
  double Pressure(void);

  // Misc
  void SetControlTemperature(bool b) {sinfo->ControlTemperature = b;}
  void SetAimedTemperature(double t) {sinfo->AimedTemperature = t;}
  void ShowSystemInformation(void);
  void ChangeScale(double alpha);

};
//----------------------------------------------------------------------
#endif //mdmanager_h
