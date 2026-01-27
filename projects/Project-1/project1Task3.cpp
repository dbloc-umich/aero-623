#include "../include/ReadGRI.h"
#include <fstream>
#include <sstream>
#include <iostream>

GRIData readGriFile(const std::string& filename) {

    std::cout << "Attempting to open: " << filename << std::endl;
    std::cout << "Current working directory: ";
    system("pwd");  // Linux/Mac
    // system("cd");  // Windows

    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Could not open file " + filename);

    GRIData data;

    file >> data.map.nNode >> data.map.nElemTot >> data.map.Dim;

    // Saves node coordinates
    data.map.nodeXYZ.reserve(data.map.nNode);
    for (int i=0;i<data.map.nNode;i++){
        double x, y, z = 0.0;
        file >> x >> y;
        if (data.map.Dim == 3) {
            file >> z;
        }
        data.map.nodeXYZ.emplace_back(std::vector<double>{x, y, z});
    }

    // Saving boundary data
    file >> data.map.nBGroup;
    data.boundaryGroup.nBFace.reserve(data.map.nBGroup);
    data.boundaryGroup.nf.reserve(data.map.nBGroup);
    data.boundaryGroup.Title.reserve(data.map.nBGroup);
    data.boundaryGroup.NB.reserve(data.map.nBGroup);
    for (int i=0;i<data.map.nBGroup;i++) {
        int nBFaceTemp, nfTemp;
        std::string title;
        file >> nBFaceTemp >> nfTemp >> title;

        std::vector<std::vector<int>> nbVec;
        nbVec.reserve(nBFaceTemp);
        for (int j=0;j<nBFaceTemp;j++) {
            std::vector<int> innerVec;
            innerVec.reserve(nfTemp);
            for (int k=0;k<nfTemp;k++) {
                int nbTemp;
                file >> nbTemp;
                innerVec.emplace_back(nbTemp);
            }
            nbVec.emplace_back(innerVec);
        }
        data.boundaryGroup.NB.emplace_back(std::move(nbVec));
        data.boundaryGroup.nBFace.emplace_back(nBFaceTemp);
        data.boundaryGroup.nf.emplace_back(nfTemp);
        data.boundaryGroup.Title.emplace_back(std::move(title));
    }

    // Saving element data
    // data.elementGroup.nElem.reserve(data.map.nElemGroup);
    // data.elementGroup.order.reserve(data.map.nElemGroup);
    // data.elementGroup.basis.reserve(data.map.nElemGroup);
    // for (int i=0;i<data.map.nElemGroup;i++) {
    int nElemCount = 0;
    while (nElemCount < data.map.nElemTot) {
        int nElemTemp, orderTemp;
        std::string basisTemp;
        file >> nElemTemp >> orderTemp >> basisTemp;
        std::vector<std::vector<int>> neVec;
        neVec.reserve(nElemTemp);
        nElemCount = nElemCount+nElemTemp;
        for (int j=0;j<nElemTemp;j++) {
            std::vector<int> innerVec2;
            int nn = 3; // CAUTION: This is hard-coded for triangular mesh right now!!!!!
            innerVec2.reserve(nn);
            for (int k=0;k<nn;k++) {
                int neTemp;
                file >> neTemp;
                innerVec2.emplace_back(neTemp);
            }
            neVec.emplace_back(innerVec2);
        }
        data.elementGroup.nElem.emplace_back(nElemTemp);
        data.elementGroup.order.emplace_back(orderTemp);
        data.elementGroup.basis.emplace_back(std::move(basisTemp));
        data.elementGroup.NE.emplace_back(std::move(neVec));
    };
    
    // Saving periodicity data
    std::string declarePeriodicity;
    file >> data.map.nPG >> declarePeriodicity;
    data.periodicGroup.nPGNode.reserve(data.map.nPG);
    data.periodicGroup.periodicity.reserve(data.map.nPG);
    data.periodicGroup.NP.reserve(data.map.nPG);
    for (int i=0;i<data.map.nPG;i++) {
        int npgNodeTemp;
        std::string periodTemp;
        std::vector<std::vector<int>> innerVec3;
        file >> npgNodeTemp >> periodTemp;
        for (int j=0;j<npgNodeTemp;j++) {
            int pair1, pair2;
            file >> pair1 >> pair2;
            innerVec3.emplace_back(std::vector<int>{pair1,pair2});
        }
        data.periodicGroup.NP.emplace_back(std::move(innerVec3));
        data.periodicGroup.nPGNode.emplace_back(npgNodeTemp);
        data.periodicGroup.periodicity.emplace_back(std::move(periodTemp));
    }

    return data;

}
