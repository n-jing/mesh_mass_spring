#ifndef TIME_INTEGRAL_JJ_H
#define TIME_INTEGRAL_JJ_H

class EdgeMesh;

enum class Integral
{
  explicit_euler,
  implicit_euler,
  location_implicit,
  speed_implicit
};


void time_integral(double *const location, double *const speed, double time, double delta_t, Integral integral);

void calc_force(const double *const var, double *g);

void explicit_integral(double *const location, double *const speed, double *const force, double delta_t);

void location_semi_implicit_integral(double *const location, double *const speed, double *const force, double delta_t);

void speed_semi_implicit_integral(double *const location, double *const speed, double *const force, double delta_t);

void implicit_integral(double *const location, double *const speed, double *const force, double delta_t);

#endif // TIME_INTEGRAL_JJ_H
