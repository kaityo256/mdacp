#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <stdlib.h>
#include "communicator.h"
#include "mpistream.h"
#include "mdmanager.h"
#include "observer.h"
#include "stopwatch.h"
//----------------------------------------------------------------------
MDManager::MDManager(int &argc, char ** &argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mout.SetRank(rank);
  mout << "# " << MDACP_VERSION << std::endl;
  std::string inputfile;
  if (argc > 1) {
    inputfile = argv[1];
  } else {
    mout << "# Input file is not specified. input.cfg is used." << std::endl;
    inputfile = "input.cfg";
  }
  param.LoadFromFile(inputfile.c_str());
  num_threads = omp_get_max_threads();
  mout << "# " << num_procs << " MPI Process(es), " << num_threads << " OpenMP Thread(s), Total " << num_procs * num_threads << " Unit(s)" << std::endl;

  pinfo = new ParaInfo(num_procs, num_threads, param);
  int grid_size[D];
  pinfo->GetGridSize(grid_size);
  sinfo = new SimulationInfo(param, grid_size);
  int tid;
  MDUnit *mdp;
  std::vector <MDUnit *> v;
  #pragma omp parallel shared(v) private(tid,mdp)
  {
    tid = omp_get_thread_num();
    mdp = new MDUnit(tid + rank * num_threads, sinfo, pinfo);
    #pragma omp critical
    v.push_back(mdp);
  }
  mdv.resize(num_threads);
  for (unsigned int i = 0; i < v.size(); i++) {
    const int local_id = GetLocalID(v[i]->GetID());
    mdv[local_id] = v[i];
  }
  s_time = 0.0;
}
//----------------------------------------------------------------------
MDManager::~MDManager(void) {
  for (unsigned int i = 0; i < mdv.size(); i++) {
    delete mdv[i];
  }
  delete pinfo;
  delete sinfo;
  MPI_Finalize();
}
//----------------------------------------------------------------------
bool
MDManager::IsValid(void) {
  return Communicator::AllReduceBoolean(pinfo->IsValid() && param.IsValid());
}
//----------------------------------------------------------------------
bool
MDManager::IsMyUnit(int id) {
  return (rank == (id / num_threads));
}
//----------------------------------------------------------------------
int
MDManager::GetLocalID(int id) {
  return (id % num_threads);
}
//----------------------------------------------------------------------
void
MDManager::SetInitialVelocity(double v0) {
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->SetInitialVelocity(v0);
  }
}
//----------------------------------------------------------------------
void
MDManager::SaveAsCdviewSequential(void) {
  static int index = 0;
  char filename[256];
  sprintf(filename, "conf%04d.cd", index);
  index++;
  SaveAsCdview(filename);
}
//----------------------------------------------------------------------
void
MDManager::SaveAsCdview(const char *filename) {
  if (0 == GetRank()) {
    std::ofstream ofs(filename);
    ofs.close();
  }
  for (int r = 0; r < GetTotalProcs(); r++) {
    Communicator::Barrier();
    if (r == GetRank()) {
      std::ofstream ofs(filename, std::ios::app);
      for (int i = 0; i < num_threads; i++) {
        mdv[i]->SaveAsCdview(ofs);
      }
      ofs.close();
    }
  }
}
//----------------------------------------------------------------------
void
MDManager::SaveConfiguration(void) {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->SaveConfiguration();
  }
}
//----------------------------------------------------------------------
void
MDManager::Calculate(void) {
  static StopWatch swAll(GetRank(), "all");
  static StopWatch swForce(GetRank(), "force");
  static StopWatch swComm(GetRank(), "comm");
  static StopWatch swPair(GetRank(), "pair");
  swAll.Start();
  if (IsPairListExpired()) {
    //mout << "# " << GetSimulationTime() << " # Expired!" << std::endl;
    swPair.Start();
    MakePairList();
    swPair.Stop();
  }
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->UpdatePositionHalf();
  }
  swComm.Start();
  SendBorderParticles();
  swComm.Stop();
  if (sinfo->ControlTemperature) {
    if (sinfo->HeatbathType == HT_NOSEHOOVER) {
      CalculateNoseHoover();
    } else {
      CalculateLangevin();
    }
  } else {
    swForce.Start();
    CalculateForce();
    swForce.Stop();
  }
  s_time += sinfo->TimeStep;
  swAll.Stop();
}
//----------------------------------------------------------------------
void
MDManager::CalculateForce(void) {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->CalculateForce();
    mdv[i]->UpdatePositionHalf();
  }
}
//----------------------------------------------------------------------
void
MDManager::CalculateNoseHoover(void) {
  double t = Temperature();
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->HeatbathZeta(t);
    mdv[i]->HeatbathMomenta();
    mdv[i]->CalculateForce();
    mdv[i]->HeatbathMomenta();
  }
  t = Temperature();
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->HeatbathZeta(t);
    mdv[i]->UpdatePositionHalf();
  }
}
//----------------------------------------------------------------------
void
MDManager::SetSeed(const int seed) {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->SetSeed(seed);
  }
}
//----------------------------------------------------------------------
void
MDManager::CalculateLangevin(void) {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->CalculateForce();
    mdv[i]->Langevin();
  }
}
//----------------------------------------------------------------------
void
MDManager::SendParticlesSub(const int dir) {
  const int o_dir = OppositeDir[dir];
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->MakeBufferForSendingParticle(dir);
  }
  std::vector<int> send_number;
  std::vector<int> recv_number;
  std::vector<ParticleInfo> send_buffer;
  std::vector<ParticleInfo> recv_buffer;
  for (int i = 0; i < num_threads; i++) {
    const int id = mdv[i]->GetID();
    const int id_dest = pinfo->GetNeighborID(id, dir);
    const int id_src = pinfo->GetNeighborID(id, o_dir);
    if (IsMyUnit(id_src)) {
      const int local_id_src  = GetLocalID(id_src);
      mdv[i]->ReceiveParticles(mdv[local_id_src]->send_buffer);
    }
    if (!IsMyUnit(id_dest)) {
      send_number.push_back(mdv[i]->send_buffer.size());
      send_buffer.insert(send_buffer.end(), mdv[i]->send_buffer.begin(), mdv[i]->send_buffer.end());
    }
  }
  const int dest_rank = pinfo->GetNeighborRank(rank, dir);
  const int src_rank = pinfo->GetNeighborRank(rank, o_dir);
  const int number = (int)send_number.size();
  Communicator::SendRecvVector(send_number, number, dest_rank, recv_number, number, src_rank);

  int send_sum = 0;
  int recv_sum = 0;
  for (int i = 0; i < number; i++) {
    send_sum += send_number[i];
    recv_sum += recv_number[i];
  }
  Communicator::SendRecvVector(send_buffer, send_sum, dest_rank, recv_buffer, recv_sum, src_rank);

  std::vector<ParticleInfo>::iterator it1 = recv_buffer.begin();
  std::vector<ParticleInfo>::iterator it2 = recv_buffer.begin();
  int index = 0;
  std::vector<ParticleInfo> temp_buffer;
  for (int i = 0; i < num_threads; i++) {
    const int id = mdv[i]->GetID();
    const int id_src = pinfo->GetNeighborID(id, o_dir);
    if (IsMyUnit(id_src)) {
      continue;
    }
    it1 = it2;
    it2 += recv_number[index];
    index++;
    temp_buffer.clear();
    temp_buffer.insert(temp_buffer.begin(), it1, it2);
    mdv[i]->ReceiveParticles(temp_buffer);
  }
}
//----------------------------------------------------------------------
void
MDManager::SendParticles(void) {
  for (int dir = 0; dir < MAX_DIR; dir++) {
    SendParticlesSub(dir);
  }
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->AdjustPeriodicBoundary();
  }
}
//----------------------------------------------------------------------
void
MDManager::MakePairList(void) {
  SendParticles();
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    const int pn = mdv[i]->GetParticleNumber();
    mdv[i]->SetTotalParticleNumber(pn);
  }
  for (int dir = 0; dir < MAX_DIR; dir++) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < num_threads; i++) {
      mdv[i]->FindBorderParticles(dir);
    }
    SendBorderParticlesSub(dir);
  }
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->MakePairList();
  }
}
//----------------------------------------------------------------------
void
MDManager::SendBorderParticles(void) {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    const int pn = mdv[i]->GetParticleNumber();
    mdv[i]->SetTotalParticleNumber(pn);
  }
  for (int dir = 0; dir < MAX_DIR; dir++) {
    SendBorderParticlesSub(dir);
  }
}
//----------------------------------------------------------------------
void
MDManager::SendBorderParticlesSub(const int dir) {
  const int o_dir = OppositeDir[dir];
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->MakeBufferForBorderParticles(dir);
  }
  std::vector<int> send_number;
  std::vector<int> recv_number;
  std::vector<ParticleInfo> send_buffer;
  std::vector<ParticleInfo> recv_buffer;
  for (int i = 0; i < num_threads; i++) {
    const int id = mdv[i]->GetID();
    const int id_dest = pinfo->GetNeighborID(id, dir);
    const int id_src = pinfo->GetNeighborID(id, o_dir);
    if (IsMyUnit(id_src)) {
      const int local_id_src  = GetLocalID(id_src);
      mdv[i]->ReceiveBorderParticles(mdv[local_id_src]->send_buffer);
    }
    if (!IsMyUnit(id_dest)) {
      send_number.push_back(mdv[i]->send_buffer.size());
      send_buffer.insert(send_buffer.end(), mdv[i]->send_buffer.begin(), mdv[i]->send_buffer.end());
    }
  }
  const int dest_rank = pinfo->GetNeighborRank(rank, dir);
  const int src_rank = pinfo->GetNeighborRank(rank, o_dir);
  const int number = (int)send_number.size();
  Communicator::SendRecvVector(send_number, number, dest_rank, recv_number, number, src_rank);

  int send_sum = 0;
  int recv_sum = 0;
  for (int i = 0; i < number; i++) {
    send_sum += send_number[i];
    recv_sum += recv_number[i];
  }
  Communicator::SendRecvVector(send_buffer, send_sum, dest_rank, recv_buffer, recv_sum, src_rank);

  std::vector<ParticleInfo>::iterator it1 = recv_buffer.begin();
  std::vector<ParticleInfo>::iterator it2 = recv_buffer.begin();
  int index = 0;
  std::vector<ParticleInfo> temp_buffer;
  for (int i = 0; i < num_threads; i++) {
    const int id = mdv[i]->GetID();
    const int id_src = pinfo->GetNeighborID(id, o_dir);
    if (IsMyUnit(id_src)) {
      continue;
    }
    it1 = it2;
    it2 += recv_number[index];
    temp_buffer.clear();
    temp_buffer.insert(temp_buffer.begin(), it1, it2);
    mdv[i]->ReceiveBorderParticles(temp_buffer);
    index++;
  }
}
//----------------------------------------------------------------------
void
MDManager::ShowSystemInformation(void) {
  const unsigned long int pn = GetTotalParticleNumber();
  sinfo->ShowAll(pn);
}
//----------------------------------------------------------------------
// For Observation
//----------------------------------------------------------------------
unsigned long int
MDManager::GetTotalParticleNumber(void) {
  unsigned long int pn = 0;
  for (int i = 0; i < num_threads; i++) {
    pn += static_cast<unsigned long int>(mdv[i]->GetParticleNumber());
  }
  pn = Communicator::AllReduceUnsignedLongInteger(pn);
  return pn;
}
//----------------------------------------------------------------------
bool
MDManager::IsPairListExpired(void) {
  bool expired = false;
  #pragma omp parallel for reduction(|:expired)
  for (int i = 0; i < num_threads; i++) {
    expired |= mdv[i]->IsPairListExpired();
  }
  expired = Communicator::AllReduceBoolean(expired);
  return expired;
}
//----------------------------------------------------------------------
double
MDManager::ObserveDouble(DoubleObserver *obs) {
  double e = 0.0;
  #pragma omp parallel for reduction(+:e)
  for (int i = 0; i < num_threads; i++) {
    e += mdv[i]->ObserveDouble(obs);
  }
  return Communicator::AllReduceDouble(e);
}
//----------------------------------------------------------------------
double
MDManager::Temperature(void) {
  return KineticEnergy() / 1.5;
}
//----------------------------------------------------------------------
double
MDManager::ConfigurationTemperature(void) {
  VirialObserver obs;
  const double pn = static_cast<double>(GetTotalParticleNumber());
  const double phi = ObserveDouble(&obs);
  return phi / pn;
}
//----------------------------------------------------------------------
double
MDManager::KineticEnergy(void) {
  KineticEnergyObserver obs;
  const double pn = static_cast<double>(GetTotalParticleNumber());
  return ObserveDouble(&obs) / pn;
}
//----------------------------------------------------------------------
double
MDManager::PotentialEnergy(void) {
  PotentialEnergyObserver obs(sinfo);
  const double pn = static_cast<double>(GetTotalParticleNumber());
  return ObserveDouble(&obs) / pn;
}
//----------------------------------------------------------------------
double
MDManager::TotalEnergy(void) {
  return KineticEnergy() + PotentialEnergy();
}
//----------------------------------------------------------------------
double
MDManager::Pressure(void) {
  VirialObserver obs;
  const double pn = static_cast<double>(GetTotalParticleNumber());
  const double phi = ObserveDouble(&obs) / pn;
  const double T = Temperature();
  const double V = sinfo->L[X] * sinfo->L[Y] * sinfo->L[Z];
  return (T - phi) * pn / V;
}
//----------------------------------------------------------------------
void
MDManager::ChangeScale(double alpha) {
  MakePairList();
  sinfo->L[X] = sinfo->L[X] * alpha;
  sinfo->L[Y] = sinfo->L[Y] * alpha;
  sinfo->L[Z] = sinfo->L[Z] * alpha;
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->ChangeScale(alpha);
  }
  MakePairList();
  mout << "# Change Scale " << alpha << std::endl;
}
//----------------------------------------------------------------------
void
MDManager::ExecuteAll(Executor *ex) {
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; i++) {
    mdv[i]->Execute(ex);
  }
}
//----------------------------------------------------------------------
void
MDManager::ExecuteAllSerial(Executor *ex) {
  for (int r = 0; r < GetTotalProcs(); r++) {
    Communicator::Barrier();
    if (r == GetRank()) {
      for (int i = 0; i < num_threads; i++) {
        mdv[i]->Execute(ex);
      }
    }
  }
}
//----------------------------------------------------------------------
void
MDManager::SaveAsDumpFileSequential(void) {
  static int index = 0;
  char filename[256];
  sprintf(filename, "conf%04d.dump", index);
  index++;
  SaveAsDumpFile(filename);
}
//----------------------------------------------------------------------
void
MDManager::SaveAsDumpFile(const char *filename) {
  class FileDumper : public Executor {
  private:
    const char *filename;
  public:
    FileDumper(const char *fn, const int rank) : filename(fn) {
      if (rank == 0) {
        std::ofstream ofs(filename, std::ios::out | std::ios::binary);
        ofs.close();
      }
    }
    void Execute(MDUnit *mdu) {
      std::ofstream ofs(filename, std::ios::out | std::ios::binary | std::ios::app);
      Variables *vars = mdu->GetVariables();
      const int pn = vars->GetParticleNumber();
      double (*q)[D] = vars->q;
      double (*p)[D] = vars->p;
      for (int i = 0; i < pn; i++) {
        ofs.write((char*) & (q[i][X]), sizeof(double));
        ofs.write((char*) & (q[i][Y]), sizeof(double));
        ofs.write((char*) & (q[i][Z]), sizeof(double));
        ofs.write((char*) & (p[i][X]), sizeof(double));
        ofs.write((char*) & (p[i][Y]), sizeof(double));
        ofs.write((char*) & (p[i][Z]), sizeof(double));
      }
      ofs.close();
    }
  };
  FileDumper fd(filename, GetRank());
  ExecuteAllSerial(&fd);
}
//----------------------------------------------------------------------
