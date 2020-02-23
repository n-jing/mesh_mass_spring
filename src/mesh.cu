#include "mesh.h"
#include <iostream>
#include <random>


using namespace std;
using namespace Eigen;


EdgeMesh edge_mesh;

__device__ __host__ const Eigen::Vector3d &EdgeMesh::get_gravity() const
{
  return gravity_;
}


__device__ __host__ EdgeMesh::EdgeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
  d_ = 0.1;
  gravity_ = {0, -9.8, 0};

  init(&V(0), V.rows(), V.cols(), &F(0), F.rows(), F.cols());
}


__device__ __host__ int EdgeMesh::init(
  const double *const V, int v_row, int v_col,
  const int *const F, int f_row, int f_col)
{
  for (size_t i = 0; i < v_row; ++i)
    add_vert(&V[i], i);

  int e_id = 0;
  double min_e = numeric_limits<double>::max();
  for (size_t i = 0; i < f_row; ++i)
  {
    for (size_t j = 0; j < f_col; ++j)
    {
      array<size_t, 2> e = {F[i, i + f_row * j],
                            F[i, i + f_row * (j+1) % f_col]};
      add_edge(e, e_id++);
      min_e = min((vert_[e[0]].v - vert_[e[1]].v).norm(), min_e);
    }
  }

  default_random_engine e;
  uniform_real_distribution<double> u(1, 2);
  for (size_t i = 0; i < v_row; ++i)
  {
#ifdef RANDOM
    set_vert_weight(i, u(e));
#else
    set_vert_weight(i, 1);
#endif
  }
  const size_t edge_num = sizeof(edge_) / sizeof(Edge);
  for (size_t i = 0; i < edge_num; ++i)
  {
    // m * g = k * delta x = k * d * x
    // where d is the coefficient
#ifdef RANDOM
    set_edge_stiffness(i, gravity_.norm() * u(e) / min_e / d_); 
#elif RIGID
    set_edge_stiffness(i, gravity_.norm() / min_e / d_ * 10000);
#else
    set_edge_stiffness(i, gravity_.norm() / min_e / d_);
#endif
  }
  return 0;
}



__device__ __host__ void EdgeMesh::set_gravity(const Eigen::Vector3d &g)
{
  gravity_ = g;
}

__device__ __host__ size_t EdgeMesh::get_vert_num() const
{
  return sizeof(vert_) / sizeof(Vert);
}

__device__ __host__ size_t EdgeMesh::get_edge_num() const
{
  return sizeof(edge_) / sizeof(Edge);
}

__device__ __host__ void EdgeMesh::set_vert_weight(size_t idx, double w)
{
  assert(idx < sizeof(vert_) / sizeof(Edge));
  vert_[idx].w = w;
}

__device__ __host__ void EdgeMesh::set_edge_stiffness(size_t idx, double k)
{
  assert(idx < sizeof(edge_) / sizeof(Edge));
  edge_[idx].k = k;
}

__device__ __host__ size_t EdgeMesh::add_vert(const double *const v, int v_id, double w)
{
  EdgeMesh::Vert vert(v, w);
  vert.id = v_id;

  vert_[v_id] = vert;
}

__device__ __host__ void EdgeMesh::add_edge(const std::array<size_t, 2> &e, int e_id, double k)
{
  EdgeMesh::Edge edge(e, k);
  edge.id = e_id;

  edge_[e_id] = edge;
}

__device__ __host__ const std::array<size_t, 2> &EdgeMesh::get_edge(size_t idx) const
{
  return edge_[idx].e;
}

__device__ __host__ double EdgeMesh::get_edge_length(size_t idx) const
{
  array<size_t, 2> endpoint = edge_[idx].e;
  const Vector3d &v1 = vert_[endpoint[0]].v;
  const Vector3d &v2 = vert_[endpoint[1]].v;
  return (v1 - v2).norm();
}


__device__ __host__ double EdgeMesh::get_edge_weight(size_t idx) const
{
  return edge_[idx].k;
}

__device__ __host__ const Eigen::Vector3d &EdgeMesh::get_vert_coord(size_t idx) const
{
  return vert_[idx].v;
}

__device__ __host__ double EdgeMesh::get_vert_weight(size_t idx) const
{
  return vert_[idx].w;
}
