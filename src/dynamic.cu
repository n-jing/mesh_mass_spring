#include "mesh.h"
#include "remove_duplicate_vert.h"
#include "solver.h"
#include "writer.h"
#include "time_integral.h"
#include <iostream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <random>
#include <Eigen/Core>
#include <string>


using namespace std;
using namespace Eigen;
using namespace igl;

int main (int argc, char *argv[])
{
  remove_duplicate_vert(argv[1], "mesh.obj");
  MatrixXd V;
  MatrixXi F;
  readOBJ("mesh.obj", V, F);


  double *v;
  int *f;
  EdgeMesh *edge_mesh;
  cudaMallocManaged((void**)&edge_mesh, sizeof(EdgeMesh));
  cudaMallocManaged((void**)&v, V.size() * sizeof(double));
  cudaMallocManaged((void**)&f, F.size() * sizeof(int));
  cudaMallocManaged((void**)&edge_mesh->vert_, V.rows() * sizeof(EdgeMesh::Vert));
  cudaMallocManaged((void**)&edge_mesh->edge_, 3*F.rows()*sizeof(EdgeMesh::Edge));

  for (int i = 0; i < V.size(); ++i)
    v[i] = V[i];
  for (int i = 0; i < F.size(); ++i)
    f[i] = F[i];

  edge_mesh->init(v, V.rows(), V.cols(), f, F.rows(), F.cols());
  cerr << edge_mesh->get_vert_num() << " " << edge_mesh->get_edge_num() << endl;

  const int vert_num = edge_mesh->get_vert_num();
  const int var_num = 3 * (vert_num - 1);
  const double time = 4;
  const double delta_t = 1e-4;
  
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> vert_u(0, vert_num-1);

  edge_mesh->fixed_vert = vert_u(gen);

  double *var;
  double *speed;
  double *vert;
  cudaMallocManaged((void**)&var, var_num * sizeof(double));
  cudaMallocManaged((void**)&speed, var_num * sizeof(double));
  cudaMallocManaged((void**)&vert, 3*vert_num * sizeof(double));

  dim3 block_size(256);
  dim3 grid_size((edge_mesh->get_vert_num() + block_size.x - 1) / block_size.x);
  init_var_and_speed<<<grid_size, block_size>>>(var, speed, edge_mesh);

  get_vert(var, vert, edge_mesh);
  cudaDeviceSynchronize();
  write_mesh_to_vtk(&init_edge_vert[0], edge_mesh, "init_state.vtk");
  
  Integral integral = Integral::location_implicit;
  time_integral(var, speed, time, delta_t, integral, edge_mesh, vert);
  
  get_vert(var, vert, edge_mesh);
  cudaDeviceSynchronize();
  write_mesh_to_vtk(&balance_edge_vert[0], edge_mesh, "final_state.vtk");

  cudaFree(v);
  cudaFree(f);
  cudaFree(edge_mesh->vert_);
  cudaFree(edge_mesh->edge_);
  cudaFree(var);
  cudaFree(speed);
  cudaFree(vert);
  return 0;
}
