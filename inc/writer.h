#ifndef WRITER
#define WRITER

#include <vector>

class EdgeMesh;


std::vector<double> get_vert(double *var, const EdgeMesh& edge_mesh);

int write_mesh_to_vtk(double *var, const EdgeMesh &edge_mesh, const char *path);

#endif // WRITER
