//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string.h>
#include "mdunit.h"
#include "fcalculator.h"
//----------------------------------------------------------------------
MDUnit::MDUnit(int id_, SimulationInfo *si, ParaInfo *pi):
  id(id_) {
  vars = new Variables();
  plist = new PairList();
  sinfo = si;
  pinfo = pi;
  int grid_size[D];
  int grid_position[D];
  pinfo->GetGridSize(grid_size);
  pinfo->GetGridPosition(id, grid_position);
  double s[D], e[D];
  for (int d = 0; d < D; d++) {
    double ul = sinfo->L[d] / static_cast<double>(grid_size[d]);
    s[d] = ul * static_cast<double>(grid_position[d]);
    e[d] = s[d] + ul;
  }
  myrect = MDRect(s, e);
  mesh = new MeshList(sinfo, myrect);
  mt.seed(GetID());
}
//----------------------------------------------------------------------
MDUnit::~MDUnit(void) {
  delete vars;
  delete mesh;
  delete plist;
}
//----------------------------------------------------------------------
void
MDUnit::AddParticle(double x[D], int type) {
  double v[D];
  v[X] = 0.0;
  v[Y] = 0.0;
  v[Z] = 0.0;
  AddParticle(x, v, type);
}
//----------------------------------------------------------------------
void
MDUnit::AddParticle(double x[D], double v[D], int type) {
  if (myrect.IsInside(x)) {
    vars->AddParticle(x, v, type);
  }
}
//----------------------------------------------------------------------
void
MDUnit::SaveAsCdview(std::ofstream &ofs) {
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  int *type = vars->type;
  const int pn = vars->GetParticleNumber();
  for (int i = 0; i < pn; i++) {
    ofs << "0 ";
    ofs << type[i] << " ";
    ofs << q[i][X] << " ";
    ofs << q[i][Y] << " ";
    ofs << q[i][Z] << " ";
    ofs << p[i][X] << " ";
    ofs << p[i][Y] << " ";
    ofs << p[i][Z] << "\n";
  }
}
//----------------------------------------------------------------------
void
MDUnit::SaveConfiguration(void) {
  char filename[256];
  sprintf(filename, "conf%03d.cnf", id);
  std::ofstream fs(filename);
  vars->SaveConfiguration(fs);
}
//----------------------------------------------------------------------
void
MDUnit::MakeBufferForSendingParticle(int dir) {
  send_buffer.clear();
  int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  int *type = vars->type;
  int index = 0;
  for (int i = 0; i < pn; i++) {
    if (myrect.IsOverBoundary(dir, q[i])) {
      ParticleInfo pi(q[i][X], q[i][Y], q[i][Z], p[i][X], p[i][Y], p[i][Z], type[i]);
      send_buffer.push_back(pi);
    } else {
      q[index][X] = q[i][X];
      q[index][Y] = q[i][Y];
      q[index][Z] = q[i][Z];
      p[index][X] = p[i][X];
      p[index][Y] = p[i][Y];
      p[index][Z] = p[i][Z];
      type[index] = type[i];
      index++;
    }
  }
  vars->SetParticleNumber(index);
}
//----------------------------------------------------------------------
void
MDUnit::ReceiveParticles(std::vector<ParticleInfo> &recv_buffer) {
  //const unsigned int recv_number = recv_buffer.size() / 6;
  const unsigned int recv_number = recv_buffer.size();
  int index = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  int *type = vars->type;
  for (unsigned int i = 0; i < recv_number; i++) {
    q[index][X] = recv_buffer[i].q[X];
    q[index][Y] = recv_buffer[i].q[Y];
    q[index][Z] = recv_buffer[i].q[Z];
    p[index][X] = recv_buffer[i].p[X];
    p[index][Y] = recv_buffer[i].p[Y];
    p[index][Z] = recv_buffer[i].p[Z];
    type[index] = recv_buffer[i].type;
    index++;
  }
  vars->SetParticleNumber(index);
}
//----------------------------------------------------------------------
void
MDUnit::FindBorderParticles(const int dir) {
  border_particles[dir].clear();
  const int pn  = vars->GetTotalParticleNumber();
  double (*q)[D] = vars->q;
  for (int i = 0; i < pn; i++) {
    if (myrect.IsInsideEdge(dir, q[i], sinfo)) {
      border_particles[dir].push_back(i);
    }
  }
}
//----------------------------------------------------------------------
void
MDUnit::MakeBufferForBorderParticles(const int dir) {
  send_buffer.clear();
  double (*q)[D] = vars->q;
  double *L = GetSystemSize();
  double diff[D] = {0.0, 0.0, 0.0};
  if (pinfo->IsOverBoundary(GetID(), dir)) {
    switch (dir) {
    case D_LEFT:
      diff[X] = L[X];
      break;
    case D_RIGHT:
      diff[X] -= L[X];
      break;
    case D_BACK:
      diff[Y] = L[Y];
      break;
    case D_FORWARD:
      diff[Y] = -L[Y];
      break;
    case D_DOWN:
      diff[Z] = L[Z];
      break;
    case D_UP:
      diff[Z] = -L[Z];
      break;
    }
  }

  for (unsigned int i = 0; i < border_particles[dir].size(); i++) {
    const int j = border_particles[dir][i];
    double x[D] = {q[j][X] + diff[X], q[j][Y] + diff[Y], q[j][Z] + diff[Z]};
    ParticleInfo pi(x[X], x[Y], x[Z], 0, 0, 0, 0);
    send_buffer.push_back(pi);
    /*
    send_buffer.push_back(x[X]);
    send_buffer.push_back(x[Y]);
    send_buffer.push_back(x[Z]);
    */
  }
}
//----------------------------------------------------------------------
void
MDUnit::ReceiveBorderParticles(std::vector<ParticleInfo> &recv_buffer) {
  //const unsigned int recv_number = recv_buffer.size() / D;
  const unsigned int recv_number = recv_buffer.size();
  int index = vars->GetTotalParticleNumber();
  double (*q)[D] = vars->q;
  for (unsigned int i = 0; i < recv_number; i++) {
    q[i + index][X] = recv_buffer[i].q[X];
    q[i + index][Y] = recv_buffer[i].q[Y];
    q[i + index][Z] = recv_buffer[i].q[Z];
  }
  index = index + recv_number;
  vars->SetTotalParticleNumber(index);
}
//----------------------------------------------------------------------
void
MDUnit::MakePairList(void) {
  mesh->Sort(vars, sinfo, myrect);
  plist->Init(vars, sinfo);
  mesh->MakeList(vars, sinfo, myrect);
}
//----------------------------------------------------------------------
void
MDUnit::ChangeScale(double alpha) {
  myrect.ChangeScale(alpha);
  vars->ChangeScale(alpha);
  mesh->ChangeScale(sinfo, myrect);
}
//----------------------------------------------------------------------
