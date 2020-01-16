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
  const int var_num = 3 * edge_mesh.get_vert_num();
  double parameter[20];
	int info[20];
	//initialize
	INIT_HLBFGS(parameter, info);
  
  
  info[4] = max_iterator_num; // the max number of iterations
  info[6] = hessinae_update_interval; // the update interval of Hessian. (typical choices: 0-200)
  info[7] = with_hessian ? 1 : 0;

  if (with_hessian)
  {
    cerr << "hessian is not supported" << endl;
    return 1;
  }
  else
  {
    //HLBFGS(N, M, INIT_X, EVALFUNC, EVALFUNC_H, HLBFGS_UPDATE_HESSIAN, NEWITERATION, PARAMETER, INFO); 
    HLBFGS(var_num, bfgs_update_num, var, eval_energy, 0, HLBFGS_UPDATE_Hessian, new_iteration, parameter, info);
  }
  return 0;
}



Lite_Sparse_Matrix* m_sparse_matrix = 0;

//////////////////////////////////////////////////////////////////////////
void eval_energy(int N, double* x, double *prev_x, double* f, double* g)
{
  for (int i = 0; i < edge_mesh.get_vert_num(); ++i)
  {
    g[3*i] = 0;
    g[3*i+1]  = 0;
    g[3*i+2] = 0;
  }

  const int edge_num = edge_mesh.get_edge_num();
	*f = 0;
  for (int e = 0; e < edge_num; ++e)
  {
    array<size_t, 2> endpoint = edge_mesh.get_edge(e);
    double k_e = edge_mesh.get_edge_weight(e);
    double new_length = get_edge_length(x, endpoint[1], endpoint[0]);
    double origin_length = edge_mesh.get_edge_length(e);
    *f += 0.5 * k_e * (new_length - origin_length) * (new_length - origin_length);
    array<int, 2> dir = {-1, 1};
    for (int v = 0; v < 2; ++v)
    {
      for (int i = 0; i < 3; ++i)
      {
        double coef = 0.5 * k_e * 2 * (new_length - origin_length);
        g[3*endpoint[v] + i] +=  coef * 0.5 * 2 * (x[3*endpoint[1] + i] - x[3*endpoint[0]+i]) * dir[v] / new_length;
      }
    }
  }
  g[0] = g[1] = g[2] = 0;
}

double get_edge_length(const double* x, int vj, int vi)
{
  Vector3d v1(x[3 * vj],
              x[3*vj+1],
              x[3*vj+2]);
  Vector3d v2(x[3 * vi],
              x[3*vi+1],
              x[3*vi+2]);
  return (v1 - v2).norm();
}


//////////////////////////////////////////////////////////////////////////
void new_iteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{
	std::cout << iter <<": " << call_iter <<" " << *f <<" " << *gnorm  << std::endl;
}
