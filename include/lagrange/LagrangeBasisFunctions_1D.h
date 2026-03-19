#ifndef LAGRANGE_BASIS_FUNCTIONS__1D_H
#define LAGRANGE_BASIS_FUNCTIONS__1D_H

#include "Eigen/Dense"
class LagrangeBasisFunctions_1D{
    public:
    LagrangeBasisFunctions_1D(int p);

    int p() const noexcept{ return _p; }
    int Np() const noexcept{ return _Np; } 
    Eigen::MatrixXd phi() const noexcept{ return _phi; }
    Eigen::MatrixXd phix() const noexcept{ return _phix; }
    Eigen::Matrix2Xd nodes() const noexcept{ return _nodes; }

    // Evaluation of the basis functions and their derivatives
    Eigen::VectorXd evalPhi (double x) const noexcept; // _phi evaluated at x
    Eigen::VectorXd evalPhiX(double x) const noexcept; // _phix evaluated at x

    // Compute the Lagrange nodes
    Eigen::Matrix2Xd getLagrangeNodes() const noexcept;

    // Evaluation of a function (given by the associated weights) and its derivatives
    double funcEval (double x, const Eigen::VectorXd& coeff);
    double funcXEval(double x, const Eigen::VectorXd& coeff);

    // Multi-D overloads that returns an Ns-by-1 vector
    // coeff is of dimension Ns*Np, where Ns = number of states, Np = number of basis functions;
    Eigen::VectorXd funcEval (double x, const Eigen::MatrixXd& coeff);
    Eigen::VectorXd funcXEval(double x, const Eigen::MatrixXd& coeff);

    protected:
    int _p; // polynomial order
    int _Np; // number of basis polynomials
    Eigen::MatrixXd _phi; // matrix of basis function coefficients, each row is a basis with the monomial coefficients given in the columns
    Eigen::MatrixXd _phix; // matrix of x-derivative coefficients, each row is a basis with the monomial coefficients given in the columns
    Eigen::Matrix2Xd _nodes; // Lagrange nodes on a unit right triangle
};

#endif