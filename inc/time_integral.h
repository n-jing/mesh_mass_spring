#ifndef TIME_INTEGRAL_JJ_H
#define TIME_INTEGRAL_JJ_H

class EdgeMesh;

enum class Integral
{
  explicit_euler,
  implicit_euler,
  semi_implicit
};


void time_integral(double *const location, double *const speed, double time, double delta_t, Integral integral);

void calc_force(const double *const var, double *g);

#endif // TIME_INTEGRAL_JJ_H
