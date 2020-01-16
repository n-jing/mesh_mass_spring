#ifndef SOLVER_JJ_H
#define SOLVER_JJ_H

class EdgeMesh;


int construct_solver(
  double *const init_var,
  const int bfgs_update_num = 7,
  const int max_iterator_num = 1000,
  const int hessinae_update_interval = 0,
  const bool with_hessian = false);

void eval_energy(int N, double* x, double *prev_x, double* f, double* g);

double get_edge_length(const double* x, int vj, int vi);

void new_iteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm);

int test();


#endif // SOLVER_JJ_H
