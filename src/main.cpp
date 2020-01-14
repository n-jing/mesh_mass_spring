#include <iostream>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include "mesh.h"

using namespace std;
using namespace Eigen;
using namespace igl;


int main (int argc, char *argv[])
{
  MatrixXd V;
  MatrixXi F;
  readOBJ(argv[1], V, F);
  cerr << V << endl;
  getchar();
  cerr << F << endl;
  getchar();
  EdgeMesh edge_mesh(V, F);
  cerr << "Hello World!" << endl;
  return 0;
}
