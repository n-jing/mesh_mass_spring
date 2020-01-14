#ifndef MESHH_JJ_H
#define MESHH_JJ_H

#include <vector>
#include <Eigen/Core>
#include <unordered_map>

class EdgeMesh
{
public:
  EdgeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
  struct Vert
  {
    Vert(const Eigen::Vector3d &vert, double weight) : v(vert), w(weight) { }
    size_t id;
    Eigen::Vector3d v;
    double w;
  };
  struct Edge
  {
    Edge(const std::pair<size_t, size_t> &edge, double ki) : e(edge), k(ki) { }
    size_t id;
    std::pair<size_t, size_t> e;
    double k;
  };

public:
  size_t add_vert(const Eigen::Vector3d &v, double w);
  size_t add_edge(const std::pair<size_t, size_t> &e, double k);
  size_t add_edge(size_t v1, size_t v2, double k);
  const std::pair<size_t, size_t> &get_edge(size_t idx) const;
  double get_edge_weight(size_t idx) const;
  const Eigen::Vector3d &get_vert_coord(size_t idx) const;
  double get_vert_weight(size_t idx) const;
  const std::vector<size_t> &get_neighbor_vert(size_t idx) const;
  
private:
  std::vector<Vert> vert_; // vert
  std::vector<Edge> edge_;
  std::unordered_map<size_t, std::vector<size_t>> neigh_vert_;
};


#endif // MESHH_JJ_H
