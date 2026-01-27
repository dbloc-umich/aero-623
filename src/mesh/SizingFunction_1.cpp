#include "mesh/SizingFunction.h"

#include <algorithm>
#include <cmath>

namespace mesh {

double SizingFunction::operator()(double d, double s) const
{
    d = std::max(0.0, d);

    // Surface Sizing
    const double gLE = 1 - std::exp(-std::pow((s - sLE_) / wLE_, 2));
    const double gTE = 1 - std::exp(-std::pow((s - sTE_) / wTE_, 2));
    
    // Surface Sizing variation
    const double hSurf = hMin_ + (hMid_ - hMin_) * gLE * gTE;

    // Growth Scaling
    const double gDeltaLE = 1 - std::exp(-std::pow((s - sLE_) / wDelta_, 2));
    const double gDeltaTE = 1 - std::exp(-std::pow((s - sTE_) / wDelta_, 2));

    // Growth Scaling along Surface
    const double delta = deltaMin_ + (deltaMax_ - deltaMin_) * gDeltaLE * gDeltaTE;

    // Final Sizing Function
    const double h = hSurf + (hInf_ - hSurf) * (1 - std::exp(-d / delta));

    // // wall transition
    // const double r  = d / std::max(1e-12, d0_);
    // const double rp = std::pow(r, std::max(1.0, pWall_));
    // const double t  = rp / (1.0 + rp);
    // double h = hMin_ + (hMax_ - hMin_) * t;

    // // LE/TE bump
    // const double sig = std::max(1e-12, xSigma_);
    // const double gLE = std::exp(-0.5 * std::pow((s - xLE_) / sig, 2));
    // const double gTE = std::exp(-0.5 * std::pow((s - xTE_) / sig, 2));
    // const double g   = std::max(gLE, gTE);

    // const double s = std::clamp(edgeStrength_, 0.0, 1.0);

    // // distance gate to prevent far-field LE/TE bias
    // const double dEdge = 2.0 * d0_;
    // const double q = d / std::max(1e-12, dEdge);
    // const double gNear = std::exp(-0.5 * q*q);

    // h *= (1.0 - s * g * gNear);

    return h;
}


} // namespace mesh

