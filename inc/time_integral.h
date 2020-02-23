#ifndef TIME_INTEGRAL_JJ_H
#define TIME_INTEGRAL_JJ_H


#include <array>
#include <Eigen/Core>

class EdgeMesh;

enum class Integral
{
  explicit_euler,
  implicit_euler,
  location_implicit,
  speed_implicit
};


void time_integral(double *const location, double *const speed, double time, double delta_t, Integral integral, const EdgeMesh *const edge_mesh, double *vert);

void calc_force(const double *const var, double *g, const EdgeMesh *const edge_mesh);

__global__ void explicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh);

__global__ void location_semi_implicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh);

void speed_semi_implicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh);

void implicit_integral(double *const location, double *const speed, double *const force, double delta_t, const EdgeMesh *const edge_mesh);


__global__ void init_force(double *g, const EdgeMesh *const edge_mesh);

void movement_integral(double *const location, double *const speed, double *const force, double delta_t, Integral integral, const EdgeMesh *const edge_mesh);

__global__ void init_var_and_speed(double *var, double *speed, EdgeMesh *edge_mesh);


std::array<Eigen::Vector3d, 2> get_edge_vert(const double* x, int vj, int vi, const EdgeMesh *const edge_mesh);

__global__ void update_vert(double *var,double *vert,const EdgeMesh *const edge_mesh);


#endif // TIME_INTEGRAL_JJ_H
