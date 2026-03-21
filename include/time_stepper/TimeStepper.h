#ifndef TIME_STEPPER_H
#define TIME_STEPPER_H

#include "Eigen/Dense"
#include <memory>

class StateMesh;
class FVFlux;
class TimeStepper{
    public:
    TimeStepper(double CFL, double gamma, std::shared_ptr<FVFlux> flux): _CFL(CFL), _gamma(gamma), _flux(flux) {}
    virtual ~TimeStepper() = default;
    // u: state vector, s: wave-speed vector on the edges
    virtual Eigen::ArrayXd dt(const StateMesh& u) const noexcept = 0;

    double CFL() const noexcept{ return _CFL; }
    void setCFL(double CFL) noexcept{ _CFL = CFL; }

    protected:
    double _CFL;
    double _gamma;
    std::shared_ptr<FVFlux> _flux;
};

#endif