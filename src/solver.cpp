#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif


#include "HLBFGS.h"
#include "Lite_Sparse_Matrix.h"
#include "solver.h"
#include "mesh.h"
#include <iostream>
#include <Eigen/Core>
#include "writer.h"

using namespace std;
using namespace Eigen;

extern const EdgeMesh edge_mesh;

int construct_solver(
  double *const var,
  const int bfgs_update_num,
  const int max_iterator_num,
  const int hessinae_update_interval,
  const bool with_hessian)
{
  const int var_num = 3 * (edge_mesh.get_vert_num() - 1);
  double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
  
  info[4] = max_iterator_num; // the max number of iterations
  info[6] = hessinae_update_interval; // the update interval of Hessian. (typical choices: 0-200)
  info[7] = with_hessian ? 1 : 0;

  parameter[5] = 1.0e-6; // stop criterion on relative error ||g|| /max(1, ||x||)
  // parameter[6] = 1.0e-6; // stop criterion on absolute error ||g||
  if (with_hessian)
  {
    cerr << "hessian is not supported" << endl;
    return 1;
  }
  else
  {
    //HLBFGS(N, M, INIT_X, EVALFUNC, EVALFUNC_H, HLBFGS_UPDATE_HESSIAN, NEWITERATION, PARAMETER, INFO); 
    HLBFGS(var_num, bfgs_update_num, var, eval_energy, 0, HLBFGS_UPDATE_Hessian, new_iteration, parameter, info);
    // HLBFGS(var_num, 0, var, eval_energy, 0, HLBFGS_UPDATE_Hessian, new_iteration, parameter, info);
  }
  return 0;
}



Lite_Sparse_Matrix* m_sparse_matrix = 0;

//////////////////////////////////////////////////////////////////////////
void eval_energy(int N, double* x, double *prev_x, double* f, double* g)
{
  const int vert_num = edge_mesh.get_vert_num();
  const int var_num = 3 * (vert_num - 1);
  for (int i = 0; i < var_num; ++i)
    g[i] = 0;
  const int edge_num = edge_mesh.get_edge_num();
  const int fixed_vert = edge_mesh.fixed_vert;
	*f = 0;
  for (int e = 0; e < edge_num; ++e)
  {
    array<size_t, 2> endpoint = edge_mesh.get_edge(e);
    double k_e = edge_mesh.get_edge_weight(e);
    array<Vector3d, 2> edge_vert = get_edge_vert(x, endpoint[1], endpoint[0]);
    double new_length = (edge_vert[1] - edge_vert[0]).norm();
    double origin_length = edge_mesh.get_edge_length(e);
    *f += 0.5 * k_e * (new_length - origin_length) * (new_length - origin_length);

    array<int, 2> dir = {1, -1};
    for (int v = 0; v < 2; ++v)
    {
      if (endpoint[v] == fixed_vert)
        continue;
      int var_id = endpoint[v] > fixed_vert ? endpoint[v] - 1 : endpoint[v];
      const double coef = 0.5 * k_e * 2 * (new_length - origin_length);
      for (int i = 0; i < 3; ++i)
      {
        g[3*var_id + i] +=  coef * 0.5 * 2 *
          (edge_vert[1][i] - edge_vert[0][i]) * dir[v] / new_length;
      }
    }
  }
  double ela = *f;
  // cerr << "ela:" << ela << endl;
  const Vector3d &horizontal = edge_mesh.get_vert_coord(fixed_vert);
  for (int i = 0; i < vert_num; ++i)
  {
    if (i == fixed_vert)
      continue;
    double weight = edge_mesh.get_vert_weight(i);
    int var_id = i > fixed_vert ? i - 1 : i;
    Vector3d v(x[3*var_id], x[3*var_id + 1], x[3*var_id + 2]);
    *f += weight * (horizontal - v).dot(EdgeMesh::get_gravity());

    for (int a = 0; a < 3; ++a)
      g[3*var_id + a] += weight * (-1) * EdgeMesh::get_gravity()[a];
  }
  
  // cerr << "gra:" << *f - ela << endl;

  // cerr << "F:" << *f << endl;
  // for (int i = 0; i < edge_mesh.get_vert_num() - 1; ++i)
  // {
  //   cerr << x[3*i] <<" " << x[3*i+1] <<" " << x[3*i+2] << endl;
  // }
  // cerr << "**********************" << endl;
  // for (int i = 0; i < edge_mesh.get_vert_num() - 1; ++i)
  // {
  //   cerr << g[3*i] <<" " << g[3*i+1] <<" " << g[3*i+2] << endl;
  // }
  // cerr << "========================" << endl;
  // getchar();
}

std::array<Vector3d, 2> get_edge_vert(const double* x, int vj, int vi)
{
  array<int, 2> vert_id = {vj, vi};
  array<Vector3d, 2> v;
  for (int i = 0; i < 2; ++i)
  {
    if (vert_id[i] == edge_mesh.fixed_vert)
      v[i] = edge_mesh.get_vert_coord(vert_id[i]);
    else
    {
      int var_id = vert_id[i] > edge_mesh.fixed_vert ? vert_id[i] - 1 : vert_id[i];
      v[i] << x[3 * var_id], x[3*var_id+1], x[3*var_id+2];
    }
  }
  return v;
}

//////////////////////////////////////////////////////////////////////////
void new_iteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
  static int id_file = 0;
  ++id_file;
  string str_f = "inter_f_" + to_string(id_file) + ".vtk";

  vector<double> init_edge_vert = get_vert(x, edge_mesh);
  write_mesh_to_vtk(&init_edge_vert[0], edge_mesh, str_f.c_str());

  // cerr << "F:" << *f << endl;
  // for (int i = 0; i < edge_mesh.get_vert_num() - 1; ++i)
  // {
  //   cerr << x[3*i] <<" " << x[3*i+1] <<" " << x[3*i+2] << endl;
  // }
  // cerr << "**********************" << endl;
  // for (int i = 0; i < edge_mesh.get_vert_num() - 1; ++i)
  // {
  //   cerr << g[3*i] <<" " << g[3*i+1] <<" " << g[3*i+2] << endl;
  // }
  // cerr << "========================" << endl;
  // getchar();
}
