#ifndef CP_CALC_H
#define CP_CALC_H

#include "StateMesh.h"
#include "TriangularMesh.h"
#include "Element.h"
#include "Face.h"
#include "Lagrange2DBasisFunctions.h"
#include "CurvedFace.h"
#include <Eigen/Dense>

void calcCP(const StateMesh& U, std::string filePath);

#endif