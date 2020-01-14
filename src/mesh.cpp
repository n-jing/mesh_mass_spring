#include "mesh.h"


using namespace std;
using namespace Eigen;

EdgeMesh::EdgeMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
  
}



size_t EdgeMesh::add_vert(const Eigen::Vector3d &v, double w)
{
  EdgeMesh::Vert vert(v, w);
  vert.id = vert_.size();

  vert_.push_back(vert);
  return vert.id;
}

size_t EdgeMesh::add_edge(const std::pair<size_t, size_t> &e, double k)
{
  neigh_vert_[e.first].push_back(e.second);
  neigh_vert_[e.second].push_back(e.first);
  
  EdgeMesh::Edge edge(e, k);
  edge.id = edge_.size();

  edge_.push_back(edge);
  return edge.id;
}

size_t EdgeMesh::add_edge(size_t v1, size_t v2, double k)
{
  neigh_vert_[v1].push_back(v2);
  neigh_vert_[v2].push_back(v1);
  
  EdgeMesh::Edge edge({v1, v2}, k);
  edge.id = edge_.size();

  edge_.push_back(edge);
  return edge.id;
}

const std::pair<size_t, size_t> &EdgeMesh::get_edge(size_t idx) const
{
  return edge_.at(idx).e;
}

double EdgeMesh::get_edge_weight(size_t idx) const
{
  return edge_.at(idx).k;
}

const Eigen::Vector3d &EdgeMesh::get_vert_coord(size_t idx) const
{
  return vert_.at(idx).v;
}

double EdgeMesh::get_vert_weight(size_t idx) const
{
  return vert_.at(idx).w;
}

const std::vector<size_t> &EdgeMesh::get_neighbor_vert(size_t idx) const
{
  return neigh_vert_.at(idx);
}
