// Calculates Barth-Jespersen limiter
#pragma once
#include <vector>
#include <string>
#include </Eigen/Dense>
#include "/include/mesh/TriangularMesh.h"

Eigen::Matrix3d bj_limiter(const Eigen::Matrix3d& Lgrad, const TriangularMesh& triMesh)