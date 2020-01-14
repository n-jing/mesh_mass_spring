#ifndef REMOVE_DUPLICATE_VERT
#define REMOVE_DUPLICATE_VERT


#include <string>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <unordered_map>

template <class T>
inline void hash_combine(std::size_t &seed, const T &v) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<typename T, size_t d>
std::size_t HashFunc(const T &key)
{
  size_t value = 0;
  for (size_t i = 0; i < d; ++i)
    hash_combine(value, key[i]);

  return value;
}

template<typename T, size_t d>
bool EqualKey(const T &lhs, const T &rhs)
{
  for (size_t i = 0; i < d; ++i)
  {
    if (lhs[i] != rhs[i])
      return false;
  }

  return true;
}

bool remove_duplicate_vert(const char *in, const char *out)
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;
  std::string ext(in);
  ext = ext.substr(ext.rfind(".") + 1);
  int is_in_valid = false;
  if (ext == "obj")
    is_in_valid = igl::readOBJ(in, V, F);
  else if (ext == "stl")
    is_in_valid = igl::readSTL(in, V, F, N);
  else if (ext == "off")
    is_in_valid = igl::readOFF(in, V, F);
  else
  {
    std::cout << "only obj, stl and off are supported" << std::endl;
    return false;
  }
  if (!is_in_valid)
    return false;
  std::function<size_t(const Eigen::Vector3d &)> HashFunc3 =
    std::bind(HashFunc<Eigen::Vector3d, 3>, std::placeholders::_1);
  std::function<bool(const Eigen::Vector3d &, const Eigen::Vector3d &)> EqualKey3 =
    std::bind(EqualKey<Eigen::Vector3d, 3>, std::placeholders::_1, std::placeholders::_2);

  std::unordered_map<
    Eigen::Vector3d, size_t, decltype(HashFunc3), decltype(EqualKey3)>
    map_vert(1, HashFunc3, EqualKey3);
  std::vector<Eigen::Vector3d> table_vert;
  std::unordered_map<size_t, size_t> map_old_to_new_vert_id;
  std::vector<Eigen::Array<size_t, 3, 1>> table_tri;

  const size_t vert_num = V.rows();
  for (size_t i = 0; i < vert_num; ++i)
  {
    const Eigen::Vector3d v = V.row(i);
    auto ret = map_vert.emplace(v, map_vert.size());
    if (ret.second == true)
    {
      table_vert.push_back(v);
      map_old_to_new_vert_id.emplace(i, ret.first->second);
    }
    else
    {
      map_old_to_new_vert_id.emplace(i, ret.first->second);
    }
  }
  const size_t num_face = F.rows();
  const size_t cols = F.cols();
  for (size_t c = 0; c < cols; ++c)
    for (size_t i = 0; i < num_face; ++i)
      F(i, c) = map_old_to_new_vert_id.at(F(i, c));

  Eigen::MatrixXd NewV(table_vert.size(), 3);
  for (size_t axis = 0; axis < 3; ++axis)
    for (size_t i = 0; i < table_vert.size(); ++i)
      NewV(i, axis) = table_vert[i][axis];

  bool is_out_valid = false;
  if (ext == "obj")
    is_out_valid = igl::writeOBJ(out, NewV, F);
  else if (ext == "stl")
    is_out_valid = igl::writeSTL(out, NewV, F);
  else if (ext == "off")
    is_out_valid = igl::writeOFF(out, NewV, F);

  return is_out_valid;
}

#endif // REMOVE_DUPLICATE_VERT
