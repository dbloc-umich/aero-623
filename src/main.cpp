// src/main.cpp
#include <iostream>
#include <fstream>

#include "Constants.h"
#include "Element.h"
#include "Face.h"
#include "FEAdvection.h"
#include "FESteadySolver.h"
#include "FreeStreamBC.h"
#include "HLLEFlux.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "InviscidWallBC.h"
#include "LocalTimeStepper.h"
#include "OutletBC.h"
#include "RK4.h"
#include "RoeFlux.h"
#include "SSP_RK3.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include "Lagrange2DBasisFunctions.h"

int main() {
    std::shared_ptr<TriangularMesh> mesh;
    std::string meshName;
    do{
        std::cout << "Enter mesh name (\"test\", \"coarse\", \"fine\", \"finer\", or \"finest\"): ";
        std::cin >> meshName;
        std::transform(meshName.begin(), meshName.end(), meshName.begin(), [](unsigned char c){ return std::tolower(c); });
    } while (meshName != "test" && meshName != "coarse" && meshName != "fine" && meshName != "finer" && meshName != "finest");

    // Run first calculations with a low-order solver first
    if (meshName == "test") mesh = std::make_shared<TriangularMesh>("projects/Project-3/test2.gri", 0, 1, 2, false);
    else if (meshName == "coarse") mesh = std::make_shared<TriangularMesh>("projects/Project-2/mesh_refined_2394.gri", 0, 1, 2);
    else if (meshName == "fine") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined1.gri", 0, 1, 2);
    else if (meshName == "finer") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined2.gri", 0, 1, 2);
    else mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined3.gri", 0, 1, 2);

    // Inlet conditions
    double gamma = 1.4;
    double rho0 = 1;
    double a0 = 1;
    double p0 = rho0*a0*a0/gamma;
    double alpha = 50*mconst::pi/180;
    double pout = 0.7*p0;
    double M = 0.1;

    // Boundary conditions
    //std::shared_ptr<InletOutletBC> inlet = std::make_shared<InletOutletBC>(rho0, a0, alpha, pout, gamma);
    std::shared_ptr<InletBC> inlet = std::make_shared<InletBC>(rho0, a0, alpha, gamma);
    std::shared_ptr<BoundaryCondition> wall = std::make_shared<InviscidWallBC>(gamma);
    std::shared_ptr<BoundaryCondition> outlet = std::make_shared<OutletBC>(pout, gamma);
    std::shared_ptr<FreeStreamBC> freeStream = std::make_shared<FreeStreamBC>(gamma);
    std::vector<std::shared_ptr<BoundaryCondition>> bc{wall, inlet, wall, outlet};

    // Initialize to uniform inlet flow conditions
    StateMesh U(mesh, bc, 4, 0);
    U.state(0).fill(rho0);
    U.state(1).fill(rho0*M*a0*std::cos(alpha));
    U.state(2).fill(rho0*M*a0*std::sin(alpha));
    U.state(3).fill(p0/(gamma-1) + 0.5*rho0*M*M*a0*a0);

    // Solver
    std::shared_ptr<FVFlux> flux;
    std::string fluxName;
    do{
        std::cout << "Enter flux name (\"roe\" or \"hlle\"): ";
        std::cin >> fluxName;
        std::transform(fluxName.begin(), fluxName.end(), fluxName.begin(), [](unsigned char c){ return std::tolower(c); });
    } while (fluxName != "roe" && fluxName != "hlle");
    if (fluxName == "roe") flux = std::make_shared<RoeFlux>(gamma);
    else flux = std::make_shared<HLLEFlux>(gamma);

    //debug--------
    std::cout << "Setting up solver..." << std::endl;
    // ------------

    std::shared_ptr<Residual> residual = std::make_shared<FEAdvection>(flux);

    //debug--------
    std::cout << "Residual set up." << std::endl;
    // ------------

    std::shared_ptr<TimeIntegrator> integrator = std::make_shared<RK4>();

    //debug--------
    std::cout << "Time integrator set up." << std::endl;
    // ------------

    std::shared_ptr<TimeStepper> stepper = std::make_shared<LocalTimeStepper>(0.8, gamma, flux);

    //debug--------
    std::cout << "Time stepper set up." << std::endl;
    // ------------

    std::unique_ptr<Solver> solver = std::make_unique<FESteadySolver>(residual, integrator, stepper);

    // debug--------
    std::cout << "Solver set up." << std::endl;
    // ------------

    try{

        //debug--------
        std::cout << "Beginning solve..." << std::endl;
        // ------------

        solver->solve(U);

        // debug--------
        std::cout << "Solve complete." << std::endl;
        // ------------

        Eigen::MatrixXd results = solver->getResult().back(); // size = 1 if steady, more than 1 if unsteady

        // debug--------
        std::cout << "Results obtained from solver." << std::endl;
        // ------------

        std::vector<double> l1norm = solver->getNorm();

        //debug--------
        std::cout << "Solve complete. Writing results to file..." << std::endl;
        // ------------

        std::ofstream file;
        std::string resultFilePath = "projects/Project-3/results/";
        resultFilePath += meshName + "_mesh_steady_p0_q1_RK4_";
        resultFilePath += fluxName;

        //debug --------
        std::cout << "Result file path: " << resultFilePath << std::endl;
        // ------------
        
        file.open(resultFilePath + "_norm.txt");
        for (auto norm: l1norm) file << norm << "\n";
        file.close();

        // for (auto it: iter){
        //     if (steadyState == 0){
        //         std::string resultFilePathAtIter = resultFilePath + "_t_" + std::to_string(it*0.045*saveEveryNIterations);
        //         file.open(resultFilePathAtIter + ".txt");
        //     } else file.open(resultFilePath + ".txt");
        //     for (Eigen::Index e = 0; e < U.cellCount(); e++) file << results[it].col(e).transpose() << "\n";
        //     file.close();
        // }
        file.open(resultFilePath + ".txt");
        for (Eigen::Index i = 0; i < results.cols(); i++) file << results.col(i).transpose() << "\n";
        file.close();
    } catch (std::runtime_error& ex){
        std::cerr << ex.what() << std::endl;
    }

    std::size_t p, q;
    do{
        std::cout << "Enter the Lagrange polynomial order for solution approximation (p = 0, 1, 2, or 3): ";
        std::cin >> p;
    } while (p < 0 || p > 3);
    do{
        std::cout << "Enter the Lagrange polynomial order for geometry approximation (q = 1 or 3): ";
        std::cin >> q;
    } while (q != 1 && q != 3);

    if (p != 0 || q != 1){
        // Simulates a different 
        std::size_t r = 2*(p+q)+1;

        //debug
        std::cout << "Refining mesh and setting up high-order solver..." << std::endl;

        if (meshName == "test") mesh = std::make_shared<TriangularMesh>("projects/Project-3/test2.gri", p, q, r, false);
        else if (meshName == "coarse") mesh = std::make_shared<TriangularMesh>("projects/Project-2/mesh_refined_2394.gri", p, q, r);
        else if (meshName == "fine") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined1.gri", p, q, r);
        else if (meshName == "finer") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined2.gri", p, q, r);
        else mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined3.gri", p, q, r);

        Eigen::MatrixXd lowOrderSol = std::move(U.matrix());
        U = StateMesh(mesh, bc, 4, p);
        // Use the low-order solution as initial condition for the high-order solvers
        for (int k = 0; k < U.cellCount(); k++){
            auto cell = U.cell(k);
            for (int i = 0; i < int((p+1)*(p+2)/2); i++) cell.col(i) = lowOrderSol.col(k);
        }
        // std::cout << U.matrix() << std::endl;
        // std::cout << residual->computeResidual(U).lpNorm<1>() << std::endl;

        int timeOrder;
        do{
            std::cout << "Enter time integration order of accuracy (3 or 4): ";
            std::cin >> timeOrder;    
        } while (timeOrder != 3 && timeOrder != 4);
        if (timeOrder == 3) integrator = std::make_shared<SSP_RK3>();

        double cfl;
        do{
            std::cout << "Enter CFL number in (0, 1]: ";
            std::cin >> cfl;
        } while (cfl <= 0 || cfl > 1);

        //debug
        std::cout << "Setting CFL to " << cfl << "..." << std::endl;
        stepper->setCFL(cfl);

        //debug
        std::cout << "CFL set to " << stepper->CFL() << "." << std::endl;
        solver = std::make_unique<FESteadySolver>(residual, integrator, stepper);

        //debug
        std::cout << "Beginning high-order solve..." << std::endl;

        try{
            solver->solve(U);
        } catch (std::runtime_error& ex){
            std::cerr << ex.what() << std::endl;
        }
    }


    // int FVOrder;
    // do{
    //     std::cout << "Enter finite-volume order of accuracy (1 or 2): ";
    //     std::cin >> FVOrder;
    // } while (FVOrder != 1 && FVOrder != 2);
    // if (FVOrder == 1) residual = std::make_shared<FVAdvectionFirstOrder>(flux);
    // else{
    //     U.setGradientMethod(std::make_shared<HybridWalkPNGrad>());
    //     int useLimiter;
    //     do{
    //         std::cout << "Will a BJ limiter be used (0 = No, 1 = Yes)";
    //         std::cin >> useLimiter;
    //     } while (useLimiter != 0 && useLimiter != 1);
    //     residual = std::make_shared<FVAdvectionSecondOrder>(flux, useLimiter==1);
    // }

    // std::shared_ptr<TimeIntegrator> integrator;


    // double cfl;
    // do{
    //     std::cout << "Enter CFL number in (0, 1]: ";
    //     std::cin >> cfl;
    // } while (cfl <= 0 || cfl > 1);
    // std::shared_ptr<TimeStepper> stepper = std::make_shared<LocalTimeStepper>(cfl, gamma, flux);
    // std::unique_ptr<Solver> solver = std::make_unique<FESteadySolver>(residual, integrator, stepper);
    // int steadyState = 1;

    // int saveEveryNIterations, maxIterations;
    // do{
    //     std::cout << "Enter solver mode (0 = unsteady, 1 = steady): ";
    //     std::cin >> steadyState;
    // } while (steadyState != 0 && steadyState != 1);
    // if (steadyState == 1) solver = std::make_unique<FVSteadySolver>(residual, integrator, stepper);
    // else{
    //     inlet->setTransient(true);
    //     do{
    //         std::cout << "Enter the frequency (after every how many iterations) that data are saved: ";
    //         std::cin >> saveEveryNIterations;
    //     } while (saveEveryNIterations < 1);
    //     do{
    //         std::cout << "Enter the maximum number of iterations: ";
    //         std::cin >> maxIterations;
    //     } while (maxIterations < 1);
    //     solver = std::make_unique<FVUnSteadySolver>(residual, integrator, stepper, saveEveryNIterations, maxIterations);
    // }

    return 0;
}