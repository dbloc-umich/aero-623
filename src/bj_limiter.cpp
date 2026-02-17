#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include </Eigen/Dense>
#include "/include/mesh/TriangularMesh.h"

Eigen::Matrix3d bj_limiter(const Eigen::Matrix3d& Lgrad, const TriangularMesh& triMesh) {
    // std::array<std::array<double>> Li(triMesh.numElems(),3);
    // For each cell:
    Eigen::Matrix3d L_limit = Lgrad;
    for (i=0;i<triMesh.numElems();i++){
    // Calculate u at each node
        std::vector<std::vector<double>> uVals(4,4); //[u0, u1, u2, u3] at nodes; for each of the 4 states
        std::vector<std::vector<double>> alphaVals(4,3);
        Eigen::Vector2d centroid = triMesh.centroid(i);
        auto& elem = mesh.elem(i);
        for (k=0;k<3;k++){
            uVals[k][0] = triMesh.getState(i,k);//GET CURRENT STATE AT CELL i *******************
        }
        
        for (j=0;j<3;j++) { // iterating over each of the nodes
        // Compute ray from centroid to node
            Eigen::Vector2d nodePoint = triMesh.node(elem._pointID[j]);
            Eigen::Vector2d rayVec = nodePoint - centroid;
            for (k=0;k<3;k++){
                uVals[k][j+1] = uVals[k][0] + rayVec.dot(L0[i][k]); // SOMEHOW GET L0 FOR THE CELL AND WANT TO TREAT LIKE AN EIGEN::Vector2d *******************
            }
        }

        std::vector<double> alphaSet;
        for (int k=0;k<3;k++) { // iterate over each state
            double umin, umax, alpha;
            int idx;
            idx = std::min_element(uVals[k].begin(),uVals[k].end());
            umin = uVals[k][idx];
            idx = std::max_element(uVals[k].begin(),uVals[k].end());
            umax = uVals[k][idx];

            // find alpha
            for (int j=0;j<2;j++) { // iterate over 3 nodes
                if (uVals[k][j+1]-uVals[k][0] > 0) {
                    alphaVals[k][j] = std::min(1,(umax-uVals[k][0])/(uVals[k][j+1]-uVals[k][0]));
                }
                else if (uVals[k][j+1]-uVals[k][0] < 0) {
                    alphaVals[k][j] = std::min(1,(umin-uVals[k][0])/(uVals[k][j+1]-uVals[k][0]));
                }
                else {alphaVals[k][j] = 1;}
            }
            alpha = std::min(alphaVals[k][0],alphaVals[k][1],alphaVals[k][2]);
            L_limit[i,k] = Lgrad[i][k] * alpha;
        }
    }//repeat for all cells

    return L_limit;


}