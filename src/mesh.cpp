#include "mesh.h"
#include <iostream>
#include <random>


using namespace std;
using namespace Eigen;


EdgeMesh edge_mesh;

double EdgeMesh::d_ = 0.1;
Vector3d EdgeMesh::gravity_(0, -9.8, 0);

Eigen::Vector3d &EdgeMesh::get_gravity()
{
  return gravity_;
}


EdgeMesh::EdgeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
  init(V, F);
}

int EdgeMesh::init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
  const size_t vert_num = V.rows();
  const size_t face_num = F.rows();
  const size_t face_vert = F.cols();
  for (size_t i = 0; i < vert_num; ++i)
  {
    add_vert(V.row(i));
  }
  double min_e = numeric_limits<double>::max();
  for (size_t i = 0; i < face_num; ++i)
  {
    for (size_t j = 0; j < face_vert; ++j)
    {
      array<size_t, 2> e = {F(i, j), F(i, (j+1)%face_vert)};
      add_edge(e);
      min_e = min((vert_[e[0]].v - vert_[e[1]].v).norm(), min_e);
    }
  }

  default_random_engine e;
  uniform_real_distribution<double> u(1, 2);
  for (size_t i = 0; i < vert_num; ++i)
  {
#ifdef RANDOM
    set_vert_weight(i, u(e));
#else
    set_vert_weight(i, 1);
#endif
  }
  const size_t edge_num = edge_.size();
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

int EdgeMesh::init(const double *const V, int v_row, int v_col,
                   const double *const F, int f_row, int f_col)
{
  for (size_t i = 0; i < v_row; ++i)
    add_vert(V[i], i);

  int e_id = 0;
  double min_e = numeric_limits<double>::max();
  for (size_t i = 0; i < f_row; ++i)
  {
    for (size_t j = 0; j < f_col; ++j)
    {
      array<size_t, 2> e = {F(i, i + f_row * j),
                            F(i, i + f_row * (j+1)%face_vert)};
      add_edge(e, e_id++);
      min_e = min((vert_[e[0]].v - vert_[e[1]].v).norm(), min_e);
    }
  }

  default_random_engine e;
  uniform_real_distribution<double> u(1, 2);
  for (size_t i = 0; i < vert_num; ++i)
  {
#ifdef RANDOM
    set_vert_weight(i, u(e));
#else
    set_vert_weight(i, 1);
#endif
  }
  const size_t edge_num = edge_.size();
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



void EdgeMesh::set_gravity(const Eigen::Vector3d &g)
{
  gravity_ = g;
}

size_t EdgeMesh::get_vert_num() const
{
  return sizeof(vert_) / sizeof(Vert);
}

size_t EdgeMesh::get_edge_num() const
{
  return sizeof(edge_) / sizeof(Edge);
}

void EdgeMesh::set_vert_weight(size_t idx, double w)
{
  assert(idx < vert_.size());
  vert_[idx].w = w;
}

void EdgeMesh::set_edge_stiffness(size_t idx, double k)
{
  assert(idx < edge_.size());
  edge_[idx].k = k;
}

void EdgeMesh::add_vert(const double *const v, int v_id, double w)
{
  EdgeMesh::Vert vert(v, w);
  vert.id = v_id;

  vert_[v_id] = vert;
}

void EdgeMesh::add_edge(const std::array<size_t, 2> &e, int e_id, double k)
{
  EdgeMesh::Edge edge(e, k);
  edge.id = edge_.size();

  edge_[e_id] = edge;
}

const std::array<size_t, 2> &EdgeMesh::get_edge(size_t idx) const
{
  return edge_[idx].e;
}

double EdgeMesh::get_edge_length(size_t idx) const
{
  array<size_t, 2> endpoint = edge_[idx].e;
  const Vector3d &v1 = vert_[endpoint[0]].v;
  const Vector3d &v2 = vert_[endpoint[1]].v;
  return (v1 - v2).norm();
}


double EdgeMesh::get_edge_weight(size_t idx) const
{
  return edge_[idx].k;
}

const Eigen::Vector3d &EdgeMesh::get_vert_coord(size_t idx) const
{
  return vert_[idx].v;
}

double EdgeMesh::get_vert_weight(size_t idx) const
{
  return vert_[idx].w;
}
