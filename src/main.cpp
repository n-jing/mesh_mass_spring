#include <iostream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include "mesh.h"
#include "remove_duplicate_vert.h"


using namespace std;
using namespace Eigen;
using namespace igl;


int main (int argc, char *argv[])
{
  remove_duplicate_vert(argv[1], "mesh.obj");
  MatrixXd V;
  MatrixXi F;
  readOBJ("mesh.obj", V, F);

  EdgeMesh edge_mesh(V, F);
  
  return 0;
}
