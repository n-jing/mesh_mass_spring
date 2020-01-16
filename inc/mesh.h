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
  static void set_gravity(double g);
  static double get_gravity();
  static void set_displacement(double d);
  static double get_displacement();

public:
  size_t get_vert_num() const;
  size_t get_edge_num() const;
  int init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

public:
  size_t add_vert(const Eigen::Vector3d &v, double w = 0);
  size_t add_edge(const std::array<size_t, 2> &e, double k = 0);
  size_t add_edge(size_t v1, size_t v2, double k = 0);
  const std::array<size_t, 2> &get_edge(size_t idx) const;
  double get_edge_length(size_t idx) const;
  double get_edge_weight(size_t idx) const;
  const Eigen::Vector3d &get_vert_coord(size_t idx) const;
  double get_vert_weight(size_t idx) const;
  const std::vector<size_t> &get_neighbor_vert(size_t idx) const;

  void set_vert_weight(size_t idx, double w);
  void set_edge_stiffness(size_t idx, double k);
  
private:
  std::vector<Vert> vert_; // vert
  std::vector<Edge> edge_;
  std::unordered_map<size_t, std::vector<size_t>> neigh_vert_;
  static double g_; // gravity
  static double d_; // default displacement
};

#endif // MESHH_JJ_H
