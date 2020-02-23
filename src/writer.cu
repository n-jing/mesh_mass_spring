#include "writer.h"
#include "write_to_file.h"
#include <vector>
#include <iostream>
#include "mesh.h"


using namespace std;
using namespace Eigen;

__global__ void update_vert(double *var, double *vert,const EdgeMesh *const edge_mesh)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i == edge_mesh->fixed_vert)
  {
    const Vector3d &v = edge_mesh->get_vert_coord(i);
    for (int a = 0; a < 3; ++a)
      vert[3*i + a] = v[a];
    return ;
  } 
  int var_id = i > edge_mesh->fixed_vert ? i - 1 : i;
  for (int a = 0; a < 3; ++a)
    vert[3*i + a] = var[3*var_id + a];
}

int write_mesh_to_vtk(double *var, const EdgeMesh &edge_mesh, const char *path)
{
  const size_t edge_num = edge_mesh.get_edge_num();
  vector<MatrixXd> vec_edge;
  for (size_t i = 0; i < edge_num; ++i)
  {
    const array<size_t, 2> &e = edge_mesh.get_edge(i);
    MatrixXd edge(3, 2);
    edge << var[3 * e[0]], var[3 * e[1]],
            var[3*e[0]+1], var[3*e[1]+1],
            var[3*e[0]+2], var[3*e[1]+2];
    vec_edge.push_back(edge);
  }

  // write_to_vtk(vec_edge, path);
  return 0;
}
