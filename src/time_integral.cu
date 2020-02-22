#include "time_integral.h"
#include "mesh.h"
#include "writer.h"
#include "solver.h"
#include <iostream>
#include <vector>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;

void calc_force(const double *const var, double *g, const EdgeMesh *const edge_mesh)
{
  dim3 block_size(256);
  dim3 grid_size((edge_mesh->get_vert_num() + block_size.x - 1) / block_size.x);
  
  init_force<<<grid_size, block_size>>>(g, edge_mesh);

  const int edge_num = edge_mesh->get_edge_num();
  const int fixed_vert = edge_mesh->fixed_vert;
  for (int e = 0; e < edge_num; ++e)
  {
    array<size_t, 2> endpoint = edge_mesh->get_edge(e);
    double k_e = edge_mesh->get_edge_weight(e);
    array<Vector3d, 2> edge_vert = get_edge_vert(var, endpoint[1], endpoint[0]);
    double new_length = (edge_vert[1] - edge_vert[0]).norm();
    double origin_length = edge_mesh->get_edge_length(e);

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

}

__global__ void init_force(double *g, EdgeMesh *const edge_mesh)
{
  int vert_id = threadIdx.x + blockIdx.x * blockDim.x;
  
  if (vert_id == fixed_vert)
    return ;
  double weight = edge_mesh->get_vert_weight(vert_id);
  int var_id = vert_id > edge_mesh->fixed_vert ? vert_id - 1 : vert_id;
  for (int a = 0; a < 3; ++a)
    g[3*var_id + a] = weight * EdgeMesh::get_gravity()[a];
}


void time_integral(double *const location, double *const speed, double time, double delta_t, Integral integral, const EdgeMesh *const edge_mesh, double *vert)
{
  const int var_num = edge_mesh->get_vert_num() - 1;
  double *force;
  cudaMallocManaged((void**)&force, 3*var_num*sizeof(double));
  int count = 0;
  for (double t = 0; t < time; t += delta_t)
  {
    movement_integral(location, speed, force, delta_t, integral, edge_mesh);

    if (count++ % 100 == 0)
    {
      dim3 block_size(256);
      dim3 grid_size((edge_mesh->get_vert_num() + block_size.x-1)/block_size.x);
      update_vert<<<grid_size, block_size>>>(location, vert, edge_mesh);
      cudaDeviceSynchronize();

      string out = "state_" + to_string(count) + ".vtk";
      write_mesh_to_vtk(vert, edge_mesh, out.c_str());
    }
  }
  cudaFree(force);
}

void movement_integral(double *const location, double *const speed, double *const force, double delta_t, Integral integral, const EdgeMesh *const edge_mesh)
{
  dim3 block_size(256);
  dim3 grid_size((edge_mesh->get_vert_num() + block_size.x - 1) / block_size.x);

  if (integral == Integral::explicit_euler)
  {
    calc_force(location, force);
    explicit_integral<<<grid_size, block_size>>>
      (location, speed, force, delta_t, edge_mesh);
  }
  else if (integral == Integral::implicit_euler)
  {
    implicit_integral(location, speed, force, delta_t, edge_mesh);
  }
  else if (integral == Integral::location_implicit)
  {
    calc_force(location, force);
    location_semi_implicit_integral<<<grid_size, block_size>>>
      (location, speed, force, delta_t,edge_mesh);
  }
  else
  {
    speed_semi_implicit_integral(location, speed, force, delta_t, edge_mesh);
  }
  
}


__global__ void explicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh)
{
  int vert_id = threadIdx.x + blockIdx.x * blockDim.x;
  if (vert_id == edge_mesh->fixed_vert)
    return ;
  double weight = edge_mesh->get_vert_weight(vert_id);
  int var_id = vert_id > fixed_vert ? vert_id - 1 : vert_id;
  for (int a = 0; a < 3; ++a)
  {
    location[3*var_id + a] += speed[3*var_id + a] * delta_t;
    double acceleration = 1.0 / weight * force[3*var_id + a];
    speed[3*var_id + a] += acceleration * delta_t;
  }
}

__global__ void location_semi_implicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh)
{
  int vert_id = threadIdx.x + blockIdx.x * blockDim.x;
  if (vert_id == edge_mesh->fixed_vert)
    return ;
  double weight = edge_mesh->get_vert_weight(vert_id);
  int var_id = vert_id > fixed_vert ? vert_id - 1 : vert_id;
  for (int a = 0; a < 3; ++a)
  {
    double acceleration = 1.0 / weight * force[3*var_id + a];
    speed[3*var_id + a] += acceleration * delta_t;
    location[3*var_id + a] += speed[3*var_id + a] * delta_t;
  }
}

void speed_semi_implicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh)
{
  const int vert_num = edge_mesh->get_vert_num();
  const int fixed_vert = edge_mesh->fixed_vert;

  for (int i = 0; i < vert_num; ++i)
  {
    if (i == fixed_vert)
      continue;
    int var_id = i > fixed_vert ? i - 1 : i;
    for (int a = 0; a < 3; ++a)
      location[3*var_id + a] += speed[3*var_id + a] * delta_t;
  }

  calc_force(location, force);
  for (int i = 0; i < vert_num; ++i)
  {
    if (i == fixed_vert)
      continue;
    double weight = edge_mesh->get_vert_weight(i);
    int var_id = i > fixed_vert ? i - 1 : i;
    for (int a = 0; a < 3; ++a)
    {
      double acceleration = 1.0 / weight * force[3*var_id + a];
      speed[3*var_id + a] += acceleration * delta_t;
    }
  }
}

void implicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh)
{
  
}

__global__ void init_var_and_speed(double *var, double *speed, EdgeMesh *edge_mesh)
{
  int vert_id = threadIdx.x + blockIdx.x * blockDim.x;
  
  if (vert_id == edge_mesh->fixed_vert)
    return ;
  const Vector3d &vert = edge_mesh->get_vert_coord(vert_id);
  int var_id = vert_id > fixed_vert ? vert_id - 1 : vert_id;
  for (int j = 0; j < 3; ++j)
  {
    var[3*var_id + j] = vert[j];
    speed[3*var_id + j] = 0;
  } 
}
