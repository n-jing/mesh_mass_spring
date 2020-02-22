#ifndef MESHH_JJ_H
#define MESHH_JJ_H

#include <vector>
#include <Eigen/Core>
#include <unordered_map>
#include <array>

class EdgeMesh
{
public:
  EdgeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
  EdgeMesh() { }
  struct Vert
  {
    Vert(const double *const vert, double weight)
      {
        v[0] = vert[0]; v[1] = vert[1]; v[2] = vert[2];
        w = weight;
      }
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
  static void set_gravity(const Eigen::Vector3d &g);
  static Eigen::Vector3d &get_gravity();

  static void set_displacement(double d);
  static double get_displacement();
  int fixed_vert;
public:
  size_t get_vert_num() const;
  size_t get_edge_num() const;
  int init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
  int init(const double *const V, int v_row, int v_col,
           const double *const F, int f_row, int f_col);

public:
  size_t add_vert(const double *const v, int v_id = 0, double w = 0);
  void add_edge(const std::array<size_t, 2> &e, int e_id = 0, double k = 0);
  const std::array<size_t, 2> &get_edge(size_t idx) const;
  double get_edge_length(size_t idx) const;
  double get_edge_weight(size_t idx) const;
  const Eigen::Vector3d &get_vert_coord(size_t idx) const;
  double get_vert_weight(size_t idx) const;

  void set_vert_weight(size_t idx, double w);
  void set_edge_stiffness(size_t idx, double k);
  
// private:
public:
  Vert *vert_; // vert
  Edge *edge_;
  static double d_; // default displacement
  static Eigen::Vector3d gravity_;  // gravity
};

#endif // MESHH_JJ_H
