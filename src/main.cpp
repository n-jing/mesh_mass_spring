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
  const int bfgs_update_num = 7;
  const int hessinae_update_interval = 0;
  const bool with_hessian = false;
  const int vert_num = edge_mesh.get_vert_num();
  const int var_num = 3 * vert_num;

  default_random_engine e;
  uniform_real_distribution<double> u(0, 0.5);
  double var[var_num];
  for (int i = 0; i < vert_num; ++i)
  {
    const Vector3d &vert = edge_mesh.get_vert_coord(i);
    for (int j = 0; j < 3; ++j)
      var[3*i + j] = vert[j] + u(e) * (i != 0);
  }

  write_mesh_to_vtk(var, edge_mesh, "init_state.vtk");

  construct_solver(var, bfgs_update_num, max_iterator_num, hessinae_update_interval, with_hessian);
  string out = "balance_state.vtk";
  write_mesh_to_vtk(var, edge_mesh, out.c_str());
  return 0;
}
