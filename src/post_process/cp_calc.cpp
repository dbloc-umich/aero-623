#include <cp_calc.h>
#include <iostream>
#include <memory>
#include <cmath>
#include <fstream>

void calcCP(const StateMesh& U, std::string filePath) {
    std::string filePath_up = filePath + "_cpUp.csv";
    std::string filePath_low = filePath + "_cpLow.csv";
    std::string filePath_cxcy = filePath + "_cxcyBlade.csv";
    auto mesh = U.mesh();
    Eigen::ArrayXd sFace(mesh->numFaces());
    const ReferenceElement& ref = mesh->reference();
    const auto& edgeW = ref.edgeW(); // stored edge quad weights
    std::size_t nQ = edgeW.size();
    int p = U.p() /*polynomial order*/;
    int nP = U.Np(); /*number of basis functions*/
    Lagrange2DBasisFunctions PhiBasis(nP);

    int dim = 4; // numer of state variables
    double gamma = 1.4;
    double gmi = gamma-1;
    double P0 = 1/gamma;
    double Pout = P0*0.7;
    double c = 18.804;
    double Mout2 = (2/gmi)*(pow(P0/Pout,gmi/gamma)-1);
    double qout = 0.5*gamma*Pout*Mout2;
    std::vector<double> cpVals_low;
    std::vector<double> xVals_low;
    std::vector<double> cpVals_up;
    std::vector<double> xVals_up;
    double cx = 0;
    double cy = 0;

    // for (std::size_t i = 0; i < mesh->numFaces(); i++){
    for (std::size_t k = 0; k < mesh->numElems(); k++){ //loop over cells
        Eigen::MatrixXd cell = U.cell(k); // block matrix with basis function weight
        const Element& elem = mesh->elem(k);
        
        for (std::size_t edge = 0; edge < 3; edge++){ //loop over cell edges
            auto faceID = elem.faceID(edge);
            const auto& face = mesh->face(faceID);
            const std::string t = face.title();
            if (t == "Curve1" || t == "Curve5") { //check which curves are desired
                if (face.isCurvedFace()) {
                    const auto* curvedFace = dynamic_cast<const CurvedFace*>(&face);
                    const auto& edgeXi = ref.edgeXi(edge); // stored edge phi derivatives at quad points
                    // const auto& edgePhi = ref.edgePhi(edge); // stored edge phi values at quad points

                    std::vector<double> Pquad(nQ);
                    std::vector<double> cpVal(nQ);
                    std::vector<double> xTemp(nQ);

                    for (std::size_t q=0;q<nQ;q++) {
                        double xi = edgeXi(0, q);
                        double eta = edgeXi(1, q);
                        // double phi = edgePhi(k, q);
                        Eigen::Vector4d uquad = PhiBasis.funcEval(xi, eta, cell); // state values at the quadrature point
                        Eigen::Vector2d normal = face.normal(nQ);
                        double rho = uquad[0];
                        double rhou = uquad[1];
                        double rhov = uquad[2];
                        double rhoE = uquad[3];
                        Pquad[q] = (gamma-1)*(rhoE/rho - 0.5*rho*sqrt(pow(rhou,2) + pow(rhov,2))/rho);
                        cx += -Pquad[q]*normal[0]*edgeW[q]*face.detJ(q);
                        cy += -Pquad[q]*normal[1]*edgeW[q]*face.detJ(q);
                        cpVal[q] = (Pquad[q] - Pout)/qout/c;
                        xTemp[q] = curvedFace->xq(q).x();
                    }
                    // determine if upper or lower boundary
                    if (t== "Curve1") {
                        cpVals_up.insert(cpVals_up.end(), cpVal.begin(), cpVal.end());
                        xVals_up.insert(xVals_up.end(), xTemp.begin(), xTemp.end());
                    }
                    else if (t=="Curve5") {
                        cpVals_low.insert(cpVals_low.end(), cpVal.begin(), cpVal.end());
                        xVals_low.insert(xVals_low.end(), xTemp.begin(), xTemp.end());
                    }
            }
            else {
                std::cout << "One boundary face is not curved...error at cell " << k << " Edge " << edge << std::endl;
                return;
            }
            }
        }
    }
    std::cout << "Pressure Coefficients: Cx=" << cx << ", Cy=" << cy << std::endl;

    // Saving cp data to .csv to plot in python
    std::ofstream file1(filePath_low);
    std::vector<double>::iterator min_it_low = std::min_element(xVals_low.begin(), xVals_low.end());
    double min_val_low = *min_it_low; // shifting the x-coordinate to the leading edge of the airfoil
    file1 << "x, cp \n";
    for (std::size_t i = 0; i < xVals_low.size(); ++i)
        file1 << xVals_low[i]-min_val_low << "," << cpVals_low[i] << "\n";
    file1.close();

    std::ofstream file2(filePath_up);
    std::vector<double>::iterator min_it_up = std::min_element(xVals_up.begin(), xVals_up.end());
    double min_val_up = *min_it_up; // shifting the x-coordinate to the leading edge of the airfoil
    file2 << "x, cp \n";
    for (std::size_t i = 0; i < xVals_up.size(); ++i)
        file2 << xVals_up[i]-min_val_up << "," << cpVals_up[i] << "\n";
    file2.close();

    std::ofstream file3(filePath_cxcy);
    file3 << "cx, cy \n";
    file3 << cx << "," << cy << "\n";
    file3.close();
}
