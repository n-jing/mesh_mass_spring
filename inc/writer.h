#ifndef WRITER
#define WRITER

#include <vector>

class EdgeMesh;


__global__ void update_vert(double *var,double *vert,const EdgeMesh *const edge_mesh);

int write_mesh_to_vtk(double *var, const EdgeMesh &edge_mesh, const char *path);

#endif // WRITER
