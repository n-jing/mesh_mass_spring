#include <iostream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <random>
#include <Eigen/Core>
#include <string>
#include "mesh.h"
#include "remove_duplicate_vert.h"
#include "solver.h"
#include "writer.h"


using namespace std;
using namespace Eigen;
using namespace igl;

extern EdgeMesh edge_mesh;

int main (int argc, char *argv[])
{
  remove_duplicate_vert(argv[1], "mesh.obj");
  MatrixXd V;
  MatrixXi F;
  readOBJ("mesh.obj", V, F);

  edge_mesh.init(V, F);
  cerr << edge_mesh.get_vert_num() << " " << edge_mesh.get_edge_num() << endl;

  const int max_iterator_num = 10000;
  const int bfgs_update_num = 10;
  const int hessinae_update_interval = 0;
  const bool with_hessian = false;
  const int vert_num = edge_mesh.get_vert_num();
  const int var_num = 3 * (vert_num - 1);

  
  random_device rd;
  mt19937 gen(rd());
  default_random_engine e;
  uniform_real_distribution<double> u(0, 0.1);
  uniform_int_distribution<int> vert_u(0, vert_num-1);

  const int fixed_vert = 0;
  edge_mesh.fixed_vert = fixed_vert;
  double var[var_num];
  for (int i = 0; i < vert_num; ++i)
  {
    if (i == fixed_vert)
      continue;
    const Vector3d &vert = edge_mesh.get_vert_coord(i);
    int var_id = i > fixed_vert ? i - 1 : i;
    for (int j = 0; j < 3; ++j)
      var[3*var_id + j] = vert[j];
  }

  vector<double> init_edge_vert = get_vert(var, edge_mesh);
  write_mesh_to_vtk(&init_edge_vert[0], edge_mesh, "init_state.vtk");

  construct_solver(var, bfgs_update_num, max_iterator_num, hessinae_update_interval, with_hessian);
  
  string out = "balance_state.vtk";
  vector<double> balance_edge_vert = get_vert(var, edge_mesh);
  write_mesh_to_vtk(&balance_edge_vert[0], edge_mesh, out.c_str());
  return 0;
}
