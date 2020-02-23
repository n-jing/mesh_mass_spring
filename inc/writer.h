#ifndef WRITER
#define WRITER

#include <vector>

class EdgeMesh;

int write_mesh_to_vtk(double *var, const EdgeMesh *const edge_mesh, const char *path);

#endif // WRITER
