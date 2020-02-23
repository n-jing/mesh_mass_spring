#include "writer.h"
#include "write_to_file.h"
#include <vector>
#include <iostream>
#include "mesh.h"


using namespace std;
using namespace Eigen;

int write_mesh_to_vtk(double *var, const EdgeMesh *const edge_mesh, const char *path)
{
  const size_t edge_num = edge_mesh->get_edge_num();
  vector<MatrixXd> vec_edge;
  for (size_t i = 0; i < edge_num; ++i)
  {
    const array<size_t, 2> &e = edge_mesh->get_edge(i);
    MatrixXd edge(3, 2);
    edge << var[3 * e[0]], var[3 * e[1]],
            var[3*e[0]+1], var[3*e[1]+1],
            var[3*e[0]+2], var[3*e[1]+2];
    vec_edge.push_back(edge);
  }

  write_to_vtk(vec_edge, path);
  return 0;
}
