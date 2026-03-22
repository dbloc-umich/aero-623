#include "FEAdvection.h"
#include "Element.h"
#include "Face.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "Lagrange2DBasisFunctions.h"
#include "FVFlux.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include <Eigen/LU>
#include <iostream>

// #define MONITOR
#ifdef MONITOR
#include "FreeStreamBC.h"
#include "InviscidWallBC.h"
#include "OutletBC.h"
#endif

static Eigen::Matrix<double,4,2> physicalFlux(const Eigen::Vector4d& U, double gamma)
{
    double gm1 = gamma - 1.0;
    double rho = U(0);
    double u = U(1)/rho;
    double v = U(2)/rho;
    double rhoE = U(3);
    double p = gm1*(rhoE - 0.5*rho*(u*u + v*v));
 
    Eigen::Matrix<double,4,2> F;

    F(0,0) = rho*u;
    F(1,0) = rho*u*u + p;
    F(2,0) = rho*u*v;
    F(3,0) = (rhoE + p)*u;

    F(0,1) = rho*v;
    F(1,1) = rho*u*v;
    F(2,1) = rho*v*v + p;
    F(3,1) = (rhoE + p)*v;
    return F;
}

Eigen::MatrixXd FEAdvection::computeResidual(const StateMesh& u) const{

    // // debug --------
    // std::cout << "Computing the residual..." << std::endl;
    // // --------------

    auto mesh = u.mesh();
    Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(u.stateCount(), u.cellCount()*u.Np());

    // // debug --------
    // std::cout << "Residual initialized to zero. Starting to compute contributions from each element..." << std::endl;
    // // --------------

    if (_flux){
        const ReferenceElement& ref = mesh->reference();
        std::size_t p = u.p();
        std::size_t Np = u.Np();
        Lagrange2DBasisFunctions PhiBasis(p);

        // //debug ---------
        // std::cout << "Reference element and basis functions initialized. Starting to loop over elements..." << std::endl;
        // // --------------

        const auto& intXi = ref.intXi(); // stored internal quad points location
        const auto& intW = ref.intW(); // stored internal quad weights
        const auto& intPhiXi = ref.intPhiXi(); // stored internal xi-derivative of phi at quad points
        const auto& intPhiEta = ref.intPhiEta(); // stored internal eta-derivative of phi ad quad points
        std::size_t intNq = intW.size();

        // //debug ---------
        // std::cout << "Internal quadrature points and weights retrieved. Starting to loop over elements..." << std::endl;
        // // --------------

        const auto& edgeW = ref.edgeW(); // stored edge quad weights
        std::size_t edgeNq = edgeW.size();
        
        for (std::size_t k = 0; k < mesh->numElems(); k++){
            Eigen::MatrixXd cell = u.cell(k); // block matrix with basis function weight
#ifdef MONITOR
            std::cout << "Element " << k << ":" << std::endl;
            std::cout << cell << std::endl;
#endif
            const Element& elem = mesh->elem(k);

            for (std::size_t i = 0; i < unsigned(Np); i++){
#ifdef MONITOR
                std::cout << "\tBasis function: phi" << i << std::endl;
#endif
                if (p > 0){
#ifdef MONITOR
                    std::cout << "\t\tComputing area integral." << std::endl;
#endif
                    Eigen::Vector4d integral = Eigen::Vector4d::Zero();
                    for (std::size_t nq = 0; nq < intNq; nq++){
                        double xi = intXi(0, nq);
                        double eta = intXi(1, nq);
                        Eigen::Vector4d u = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                        Eigen::Matrix<double,4,2> F = physicalFlux(u, _flux->gamma()); // get fluxes from states
#ifdef MONITOR
                        std::cout << "\t\t\tAt quadrature point [" << xi << ", " << eta << "], u = " << u.transpose() << ", Fx = " << F.col(0).transpose() << ", Fy = " << F.col(1).transpose() << std::endl;
#endif
                        double detJ = elem.detJacobian(nq);
                        Eigen::Vector2d gRef{intPhiXi(i,nq), intPhiEta(i,nq)}; // 2-by-1
                        Eigen::Matrix2d JT = elem.jacobian(nq).transpose(); // 2-by-2
                        Eigen::Vector2d gPhi = JT.lu().solve(gRef); // 2-by-1
                        integral += intW[nq] * detJ * F*gPhi; // (4-by-2)*(2-by-1) = (4-by-1)
                        // integral += intW[nq] * F*gRef;
                    }
#ifdef MONITOR
                    std::cout << "\t\t\tArea flux contribution is " << integral.transpose() << std::endl;
#endif
                    residual.col(k*Np+i) += integral;
                }

                // Line integral over each of the edges
                for (std::size_t edge = 0; edge < 3; edge++){
                    // std::cout << "Line #" << edge << " on element " << k << std::endl;
                    const auto& edgeXi = ref.edgeXi(edge); // stored edge phi derivatives at quad points
                    const auto& edgePhi = ref.edgePhi(edge); // stored edge phi values at quad points

                    auto faceID = elem.faceID(edge);
#ifdef MONITOR
                    std::cout << "\t\tComputing line integral on edge " << faceID << "." << std::endl;
#endif
                    const auto& face = mesh->face(faceID);
                    Eigen::Vector2d refNormal;
                    if (edge == 0) refNormal = {std::sqrt(2)/2, std::sqrt(2)/2};
                    else if (edge == 1) refNormal = {1, 0};
                    else refNormal = {0, -1};
                    if (face.isBoundaryFace()){
                        // Boundary edge, use boundary condition to compute flux
#ifdef MONITOR
                        std::cout << "\t\tThis is a boundary edge on " << face.title() << ". ";
#endif
                        std::shared_ptr<BoundaryCondition> bc;
                        if (face.title() == "Curve1") bc = u.bc(0);
                        else if (face.title() == "Curve3") bc = u.bc(1);
                        else if (face.title() == "Curve5") bc = u.bc(2);
                        else bc = u.bc(3);

                        // Transient cases, these lines do nothing if running in steady-state
                        if (auto inlet = dynamic_cast<InletBC*>(bc.get())){
                            double y0 = mesh->node(face.pointID(0))[1];
                            double y1 = mesh->node(face.pointID(1))[1];
                            inlet->setTransientRho((y1+y0)/2);
                        } else if (auto inlet = dynamic_cast<InletOutletBC*>(bc.get())){
                            double y0 = mesh->node(face.pointID(0))[1];
                            double y1 = mesh->node(face.pointID(1))[1];
                            inlet->setTransientRho((y1+y0)/2);                       
                        }
#ifdef MONITOR
                        if (dynamic_cast<FreeStreamBC*>(bc.get())) std::cout << "Freestream BC." << std::endl;
                        else if (dynamic_cast<InletBC*>(bc.get())) std::cout << "Inlet BC." << std::endl;
                        else if (dynamic_cast<InletOutletBC*>(bc.get())) std::cout << "Inlet-outlet BC." << std::endl;
                        else if (dynamic_cast<InviscidWallBC*>(bc.get())) std::cout << "Inviscid wall BC." << std::endl;
                        else if (dynamic_cast<OutletBC*>(bc.get())) std::cout << "Outlet BC." << std::endl;
                        else std::cout << "Unknown BC." << std::endl;
#endif
                        Eigen::Vector4d integral = Eigen::Vector4d::Zero();
                        for (std::size_t nq = 0; nq < edgeNq; nq++){
                            double xi = edgeXi(0, nq);
                            double eta = edgeXi(1, nq);
                            double phi = edgePhi(i, nq);
                            Eigen::Vector4d u = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                            Eigen::Vector2d normal = face.normal(nq);
                            Eigen::Vector4d F = bc->computeFlux(u, normal); // get numerical flux from boundary condition
                            // Eigen::Vector4d F = bc->computeFlux(u, refNormal);
#ifdef MONITOR
                            std::cout << "\t\t\tAt quadrature point [" << xi << ", " << eta << "], u = " << u.transpose() << ", F = " << F.transpose() << std::endl;
                            if (face.title() == "Curve3"){
                                if (u[1]*normal[0] + u[2]*normal[1] < 0) std::cout << "\t\t\tThis is an inlet condition." << std::endl;
                                else std::cout << "\t\t\tThis is an outlet condition." << std::endl;
                            }
#endif
                            integral -= edgeW[nq] * face.detJ(nq) * phi*F;
                            // integral -= edgeW[nq] * phi*F;
                        }
#ifdef MONITOR
                        std::cout << "\t\t\tBoundary flux contribution is " << integral.transpose() << std::endl;
#endif
                        residual.col(k*Np+i) += integral;
                    } else{
                        // Internal edge, use numerical flux function to compute flux
                        std::size_t kn;
                        Eigen::Vector2d normal = face.normal(); // linear face so normal does not vary along the edge
                        // std::cout << face.elemID().transpose() << std::endl;

                        // // ADD THIS CHECK:
                        // if (face.elemID(1) == -1){
                        //     throw std::runtime_error("ERROR: Face " + std::to_string(faceID) + 
                        //                             " has no right element but is not a boundary face.");
                        // }

                        // // Another sanity check:
                        // if (face.isBoundaryFace()){
                        //     throw std::runtime_error("ERROR: Face " + std::to_string(faceID) + 
                        //                             " is marked as a boundary face but has two neighboring elements.");
                        // }

                        // // Another sanity check:
                        // if (face.isPeriodicFace()){
                        //     throw std::runtime_error("ERROR: Face " + std::to_string(faceID) + 
                        //                             " is marked as a periodic face but has two neighboring elements. Periodic faces should be treated as boundary faces with special BCs.");
                        // }

                        if (face.elemID(0) == int(k)) kn = face.elemID(1);
                        else{
                            kn = face.elemID(0);
                            normal *= -1; // Element k is the R element, has to flip normal
                        }
#ifdef MONITOR
                        std::cout << "\t\tThis is an interior edge, facing neighbor element " << kn << ". Normal = " << normal.transpose() << std::endl;
#endif
                        // Find the local edge index of this edge on element kn
                        std::size_t edgeN;

                        // //debug --------
                        // std::cout << "\t\tFinding the local edge index of this face on the neighboring element..." << std::endl;
                        // // --------------

                        if (face.isPeriodicFace()) faceID = face.periodicFaceID();

                        // //debug --------
                        // std::cout << "\t\tThe periodic counterpart of this face is " << faceID << "." << std::endl;
                        // // --------------

                        if (mesh->elem(kn).faceID(0) == faceID) edgeN = 0;
                        else if (mesh->elem(kn).faceID(1) == faceID) edgeN = 1;
                        else edgeN = 2;
                        
                        // //debug --------
                        // std::cout << "\t\tThe local edge index on the neighboring element is " << edgeN << "." << std::endl;
                        // // --------------

                        Eigen::MatrixXd cellN = u.cell(kn);
                        Eigen::Matrix2Xd edgeXiN = ref.edgeXi(edgeN);
#ifdef MONITOR
                        std::cout << "\t\tThe local edge index on the neighboring element is " << edgeN << "." << std::endl;
                        // Eigen::Vector2d x0 = mesh->node(elem.pointID(0));
                        // Eigen::Vector2d x1 = mesh->node(elem.pointID(1));
                        // Eigen::Vector2d x2 = mesh->node(elem.pointID(2));

                        // const Element& elemN = mesh->elem(kn);
                        // Eigen::Vector2d x0n = mesh->node(elemN.pointID(0));
                        // Eigen::Vector2d x1n = mesh->node(elemN.pointID(1));
                        // Eigen::Vector2d x2n = mesh->node(elemN.pointID(2));
#endif

                        // std::cout << "\tThe actual integration" << std::endl;
                        Eigen::Vector4d integral = Eigen::Vector4d::Zero();
                        for (std::size_t nq = 0; nq < edgeNq; nq++){
                            // std::cout << "\t Quadrature " << nq << std::endl;
                            double xi = edgeXi(0, nq);
                            double eta = edgeXi(1, nq);
                            Eigen::Vector4d uL = PhiBasis.funcEval(xi, eta, cell); // left-state values at the quadrature point

                            double xiN = edgeXiN(0, edgeNq-1-nq);
                            double etaN = edgeXiN(1, edgeNq-1-nq);
                            Eigen::Vector4d uR = PhiBasis.funcEval(xiN, etaN, cellN); // right-state values at the quadrature point

#ifdef MONITOR
                            // std::cout << "\t\t\tThis quadrature point is " << (x0 + (x1-x0)*xi + (x2-x0)*eta).transpose() << " on elemL";
                            // std::cout << " and " << (x0n + (x1n-x0n)*xiN + (x2n-x0n)*etaN).transpose() << " on elemR. " << std::endl;
#endif

                            double phi = edgePhi(i, nq);
                            Eigen::Vector4d F = _flux->computeFlux(uL, uR, normal); // get numerical flux from boundary condition
#ifdef MONITOR
                            std::cout << "\t\t\tAt quadrature point [" << xi << ", " << eta << "], uL = " << uL.transpose() << ", uR = " << uR.transpose() << ", F = " << F.transpose() << std::endl;
#endif                   
                            integral -= edgeW[nq] * face.detJ(nq) * phi*F;
                            // integral -= edgeW[nq] * phi*F;
                        }
#ifdef MONITOR
                        std::cout << "\t\t\tInterior flux contribution is " << integral.transpose() << std::endl;
#endif
                        residual.col(k*Np+i) += integral;
                    }
                }
            }
#ifdef MONITOR
            std::cout << std::endl;
#endif
        }
        if (residual.array().isNaN().any()) throw std::runtime_error("ERROR: Unable to compute residual.");
    }
    return residual;
}