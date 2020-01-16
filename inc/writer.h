#ifndef WRITER
#define WRITER

class EdgeMesh;
int write_mesh_to_vtk(double *var, const EdgeMesh &edge_mesh, const char *path);

#endif // WRITER
