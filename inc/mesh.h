#ifndef MESHH_JJ_H
#define MESHH_JJ_H

#include <vector>
#include <Eigen/Core>
#include <unordered_map>
#include <array>
#include <cuda.h>
#include <cuda_runtime.h>

class EdgeMesh
{
public:
  __device__ __host__ EdgeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
  __device__ __host__ EdgeMesh() { }
  struct Vert
  {
    Vert(const Eigen::Vector3d &vert, double weight) : v(vert), w(weight) { }
    size_t id;
    Eigen::Vector3d v;
    double w;
  };
  struct Edge
  {
    Edge(const std::array<size_t, 2> &edge, double ki) : e(edge), k(ki) { }
    size_t id;
    std::array<size_t, 2> e;
    double k;
  };
  __device__ __host__ void set_gravity(const Eigen::Vector3d &g);
  __device__ __host__ const Eigen::Vector3d &get_gravity() const;

  __device__ __host__ void set_displacement(double d);
  __device__ __host__ double get_displacement();
  int fixed_vert;
public:
  __device__ __host__ size_t get_vert_num() const;
  __device__ __host__ size_t get_edge_num() const;
  __device__ __host__ int init(const double *const V, int v_row, int v_col,
                               const int *const F, int f_row, int f_col);

public:
  __device__ __host__ size_t add_vert(const Eigen::Vector3d &vert, int v_id = 0, double w = 0);
  __device__ __host__ void add_edge(const std::array<size_t, 2> &e, int e_id = 0, double k = 0.1);
  __device__ __host__ const std::array<size_t, 2> &get_edge(size_t idx) const;
  __device__ __host__ double get_edge_length(size_t idx) const;
  __device__ __host__ double get_edge_weight(size_t idx) const;
  __device__ __host__ const Eigen::Vector3d &get_vert_coord(size_t idx) const;
  __device__ __host__ double get_vert_weight(size_t idx) const;

  __device__ __host__ void set_vert_weight(size_t idx, double w);
  __device__ __host__ void set_edge_stiffness(size_t idx, double k);
  
// private:
public:
  Vert *vert_; // vert
  Edge *edge_;
  double d_; // default displacement
  Eigen::Vector3d gravity_;  // gravity
  int vert_num;
  int edge_num;
};

#endif // MESHH_JJ_H
