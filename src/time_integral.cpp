#include "time_integral.h"
#include "mesh.h"
#include "writer.h"
#include "solver.h"
#include <iostream>
#include <vector>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;

extern const EdgeMesh edge_mesh;

void calc_force(const double *const var, double *g)
{
  const int vert_num = edge_mesh.get_vert_num();
  const int var_num = 3 * (vert_num - 1);
  for (int i = 0; i < var_num; ++i)
    g[i] = 0;

  const int edge_num = edge_mesh.get_edge_num();
  const int fixed_vert = edge_mesh.fixed_vert;
  for (int e = 0; e < edge_num; ++e)
  {
    array<size_t, 2> endpoint = edge_mesh.get_edge(e);
    double k_e = edge_mesh.get_edge_weight(e);
    array<Vector3d, 2> edge_vert = get_edge_vert(var, endpoint[1], endpoint[0]);
    double new_length = (edge_vert[1] - edge_vert[0]).norm();
    double origin_length = edge_mesh.get_edge_length(e);

    array<int, 2> dir = {-1, 1};
    for (int v = 0; v < 2; ++v)
    {
      if (endpoint[v] == fixed_vert)
        continue;
      int var_id = endpoint[v] > fixed_vert ? endpoint[v] - 1 : endpoint[v];
      const double coef = k_e * (new_length - origin_length);
      for (int i = 0; i < 3; ++i)
      {
        g[3*var_id + i] +=  coef *
          (edge_vert[1][i] - edge_vert[0][i]) * dir[v] / new_length;
      }
    }
  }

  for (int i = 0; i < vert_num; ++i)
  {
    if (i == fixed_vert)
      continue;
    double weight = edge_mesh.get_vert_weight(i);
    int var_id = i > fixed_vert ? i - 1 : i;
    for (int a = 0; a < 3; ++a)
      g[3*var_id + a] += weight * EdgeMesh::get_gravity()[a];
  }
}


void time_integral(double *const location, double *const speed, double time, double delta_t, Integral integral)
{
  const int vert_num = edge_mesh.get_vert_num();
  const int var_num = edge_mesh.get_vert_num() - 1;
  double force[3*var_num];
  const int fixed_vert = edge_mesh.fixed_vert;

  int count = 0;
  for (double t = 0; t < time; t += delta_t)
  {
    cerr << "t:" << t << endl;
    calc_force(location, force);
    for (int i = 0; i < vert_num; ++i)
    {
      if (i == fixed_vert)
        continue;
      double weight = edge_mesh.get_vert_weight(i);
      int var_id = i > fixed_vert ? i - 1 : i;
      for (int a = 0; a < 3; ++a)
      {
        if (integral == Integral::explicit_euler)
        {
          double acceleration = 1.0 / weight * force[3*var_id + a];
          speed[3*var_id + a] += acceleration * delta_t;
          location[3*var_id + a] += speed[3*var_id + a] * delta_t;
        }
        else if (integral == Integral::implicit_euler)
        {
          
        }
        else
        {
          cerr << "error in integral method" << endl;
          return ;
        }
      }
    }
    ++count;
    if (count % 1000 == 0)
    {
      string out = "state_" + to_string(count) + ".vtk";
      vector<double> vert = get_vert(location, edge_mesh);
      write_mesh_to_vtk(&vert[0], edge_mesh, out.c_str());
    }
  }
}
