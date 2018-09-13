//----------------------------------------------------------------------
// MPI Communication Class
//----------------------------------------------------------------------
#include "communicator.h"
//----------------------------------------------------------------------
void
Communicator::Barrier(void) {
  MPI_Barrier(MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::SendRecvInteger(int &send_number, int dest_rank, int &recv_number, int src_rank) {
  MPI_Status st;
  MPI_Sendrecv(&send_number, 1, MPI_INT, dest_rank, 0, &recv_number, 1, MPI_INT, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
void
Communicator::SendRecvDouble(void *sendbuf, int send_number, int dest_rank,
                             void *recvbuf, int recv_number, int src_rank) {
  MPI_Status st;
  MPI_Sendrecv(sendbuf, send_number, MPI_DOUBLE, dest_rank, 0, recvbuf, recv_number, MPI_DOUBLE, src_rank, 0, MPI_COMM_WORLD, &st);
}
//----------------------------------------------------------------------
bool
Communicator::AllReduceBoolean(bool flag) {
  int n = 0;
  if (flag) {
    n = 1;
  }
  int sum = 0;
  MPI_Allreduce(&n, & sum, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  if (sum > 0) {
    return true;
  } else {
    return false;
  }
}
//----------------------------------------------------------------------
int
Communicator::AllReduceInteger(int value) {
  int sum  = 0;
  MPI_Allreduce(&value, & sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}
//----------------------------------------------------------------------
unsigned long int
Communicator::AllReduceUnsignedLongInteger(unsigned long int value) {
  unsigned long int  sum  = 0;
  MPI_Allreduce(&value, & sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}
//----------------------------------------------------------------------
double
Communicator::FindMaxDouble(double value) {
  double max = 0;
  MPI_Allreduce(&value, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return max;
}
//----------------------------------------------------------------------
double
Communicator::AllReduceDouble(double value) {
  double sum = 0;
  MPI_Allreduce(&value, & sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sum;
}
//----------------------------------------------------------------------
void
Communicator::AllReduceDoubleBuffer(double *sendbuf, int size, double *recvbuf) {
  MPI_Allreduce(sendbuf, recvbuf, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::AllGatherInteger(int *sendbuf, int number, int *recvbuf) {
  MPI_Allgather(sendbuf, number, MPI_INT, recvbuf, number, MPI_INT, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::AllGatherDouble(double *sendbuf, int number, double *recvbuf) {
  MPI_Allgather(sendbuf, number, MPI_DOUBLE, recvbuf, number, MPI_DOUBLE, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
void
Communicator::AllGatherIntegerVector(std::vector<int> &send_buffer, std::vector<int> &recv_buffer, int num_procs) {
  const int sendcount = static_cast<int>(send_buffer.size());
  recv_buffer.clear();
  recv_buffer.resize(sendcount * num_procs);
  MPI_Allgather(&send_buffer[0], sendcount, MPI_INT, &recv_buffer[0], sendcount, MPI_INT, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
double
Communicator::GetTime(void) {
  return MPI_Wtime();
}
//----------------------------------------------------------------------
void
Communicator::BroadcastInteger(int &value, int root) {
  MPI_Bcast(&value, 1, MPI_INT, root, MPI_COMM_WORLD);
}
//----------------------------------------------------------------------
