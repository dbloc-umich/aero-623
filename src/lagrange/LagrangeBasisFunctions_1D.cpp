#include "Eigen/Dense"
#include "LagrangeBasisFunctions_1D.h"
#include <stdio.h>
#include <stdlib.h>

LagrangeBasisFunctions_1D::LagrangeBasisFunctions_1D(int p):
    _p(p),
    _Np(p+1),
    _phi(_Np, _Np),
    _phix(gradL1D_quad(getLagrangeNodes(), p))
{
    assert(p >= 0 && p <= 3);
    if (_p==0) _phi << 1;

    else if (_p==1){
        // variable order : 1, x
        _phi << 1, -1,
                0, 1;
    }

    else if (_p==2){
        // variable order : 1, x, x^2
        _phi << 1, -3, 2,
                0, 4, -4,
                0, -1, 2;
    }

    else if (_p==3){
        // variable order : 1, x, x^2, x^3
        _phi << 1, -11.0/2.0,       9.0,  -9.0/2.0,
                0,         9, -45.0/2.0,  27.0/2.0,
                0,  -9.0/2.0,      18.0, -27.0/2.0,
                0,         1,  -9.0/2.0,   9.0/2.0;
    }

    if (_p >= 1){
        _phix.col(0) = _phi.col(1);
    }

    if (_p >= 2){
        _phix.col(1) = 2.0 * _phi.col(2);
    }

    if (_p >= 3){
        _phix.col(2) = 3.0 * _phi.col(3);
    }
}

Eigen::VectorXd LagrangeBasisFunctions_1D::evalPhi(double x) const noexcept{
    if (_p == 0) return _phi;

    Eigen::MatrixXd eval = _phi;
    if (_p >= 1){
        eval.col(1) *= x;
    }

    if (_p >= 2){
        eval.col(2) *= x*x;
    }

    if (_p >= 3){
        eval.col(3) *= x*x*x;
    }

    return eval.rowwise().sum();

}

Eigen::VectorXd LagrangeBasisFunctions_1D::evalPhiX(double x) const noexcept{
    if (_p == 0) return Eigen::VectorXd::Zero(_Np);

    Eigen::MatrixXd eval = _phix;
    if (_p >= 2){
        eval.col(2) *= x*x;
    }

    if (_p >= 3){
        eval.col(3) *= x*x*x;
    }

    return eval.rowwise().sum();
}

// Have to fix here onwards
Eigen::Matrix2Xd LagrangeBasisFunctions_1D::getLagrangeNodes() const noexcept{
    Eigen::Matrix2Xd nodes(2, _Np);
    for (int i = 0; i < _Np; i++){
        nodes(0, i) = double(i)/double(_p); // x-coordinates of the nodes
        nodes(1, i) = 0.0; // y-coordinates of the nodes
    }
    return nodes;
}

double LagrangeBasisFunctions_1D::funcEval(double x, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _Np);
    if (_p == 0) return coeff[0];
    return evalPhi(x).dot(coeff);
}

Eigen::VectorXd LagrangeBasisFunctions_1D::funcEval(double x, const Eigen::MatrixXd& coeff){
    assert(coeff.cols() == _Np);
    if (_p == 0) return coeff.col(0);
    return coeff*evalPhi(x); // result is a Ns-by-1 vector
}

double LagrangeBasisFunctions_1D::funcXEval(double x, const Eigen::VectorXd& coeff){
    assert(coeff.size() == _Np);
    if (_p == 0) return 0;
    return evalPhiX(x).dot(coeff);
}

Eigen::VectorXd LagrangeBasisFunctions_1D::funcXEval(double x, const Eigen::MatrixXd& coeff){
    assert(coeff.rows() == _Np);
    if (_p == 0) return coeff.row(0).transpose();
    return coeff*evalPhiX(x); // result is a Ns-by-1 vector
}
// This is the end of the correct function everything beyond this is to be removed
/******************************************************************/


//   FUNCTION Definition: shapeL1D
int
shapeL1D(double sig, int q, double **pphi)
{
  /* Returns 1D Lagrange shape functions on [0,1]
     for order 1 <= q <= 3 at reference coordinate sig.
     pphi is reallocated to size (q+1). */
    int nshape;
    double *phi;
    nshape = q + 1;

    if ((*pphi) == NULL)
        (*pphi) = (double *) malloc(nshape*sizeof(double));
    else
        (*pphi) = (double *) realloc( (void *) (*pphi), nshape*sizeof(double));
    phi = *pphi;
    switch (q){

    // Shape functions for 1D Lagrange polynomials
    case 1:
        phi[0] = 1.0 - sig;
        phi[1] = sig;
        break;
        
    case 2:
        phi[0] = 2.0*sig*sig - 3.0*sig + 1.0;
        phi[1] = -4.0*sig*sig + 4.0*sig;
        phi[2] = 2.0*sig*sig - 1.0*sig;
        break;

    case 3:
        phi[0] = -9.0/2.0*sig*sig*sig + 9.0*sig*sig -11.0/2.0*sig +1.0;
        phi[1] = 27.0/2.0*sig*sig*sig - 45.0/2.0*sig*sig + 9.0*sig;
        phi[2] = -27.0/2.0*sig*sig*sig + 18.0*sig*sig - 9.0/2.0*sig;
        phi[3] = 9.0/2.0*sig*sig*sig - 9.0/2.0*sig*sig + sig;
        break;

    default:
        return -1;
    }
    return 0;
}

/******************************************************************/
//   FUNCTION Definition: gradientL1D
int 
gradientL1D(double sig, int q, double **qgphi)
{ 
  /* Returns dphi/dsig for 1D Lagrange shape functions on [0,1]
     for order 1 <= q <= 3 at reference coordinate sig.
     pgphi is reallocated to size (q+1). */

    int nshape;
    double *gphi;
    nshape = q+1;

    if ((*qgphi) == NULL)
        (*qgphi) = (double *) malloc(2*nshape*sizeof(double));
    else
        (*qgphi) = (double *) realloc( (void *) (*qgphi), 2*nshape*sizeof(double));

    gphi = *qgphi;
    
    switch (q){
        case 1:
            gphi[0] = -1.0;
            gphi[1] = 1.0;
            break;
                
        case 2:
            gphi[0] = 4.0*sig - 3.0;
            gphi[1] = -8.0*sig + 4.0;
            gphi[2] = 4.0*sig - 1.0;
            break;

        case 3:
            gphi[0] = -27.0/2.0*sig*sig + 18.0*sig - 11.0/2.0;
            gphi[1] = 81.0/2.0*sig*sig - 45.0*sig + 9.0;
            gphi[2] = -81.0/2.0*sig*sig + 36.0*sig - 9.0/2.0;
            gphi[3] = 27.0/2.0*sig*sig - 9.0*sig + 1.0;
            break;

        default:
            return -1;
    }
    return 0;
} // gradientL1D


/******************************************************************/
//   FUNCTION Definition: shapeL1D (Eigen wrapper)
Eigen::MatrixXd shapeL1D_quad(const Eigen::VectorXd& xiq, int p)
{
    int nQ = xiq.size();
    Eigen::MatrixXd phiq(nQ, p + 1);
    double *phi = NULL;
    for (int i = 0; i < nQ; i++) {
        shapeL1D(xiq(i), p, &phi);
        for (int j = 0; j < p + 1; j++)
            phiq(i, j) = phi[j];
    }
    free(phi);
    return phiq;
}

/******************************************************************/
//   FUNCTION Definition: gradL1D (Eigen wrapper)
Eigen::MatrixXd gradL1D_quad(const Eigen::VectorXd& xiq, int p)
{
    int nQ = xiq.size();
    Eigen::MatrixXd gphiq(nQ, p + 1);
    double *gphi = NULL;
    for (int i = 0; i < nQ; i++) {
        gradientL1D(xiq(i), p, &gphi);
        for (int j = 0; j < p + 1; j++)
            gphiq(i, j) = gphi[j];
    }
    free(gphi);
    return gphiq;
}