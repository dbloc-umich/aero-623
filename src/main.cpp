// #include "GaussKronrod.h"
#include "../include/ReadGRI.h"
#include <iostream>

int main(){
    try {
        GRIData gridData = readGriFile("/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/test.gri");

        // std::cout << "nNode: " << testData.map.nNode << std::endl;
        // std::cout << "nElemGroup: " <<testData.map.nElemGroup << std::endl;
        // std::cout << "Dim: " << testData.map.Dim << std::endl;

        // return 0;
        std::cout << "=== Grid Information ===" << std::endl;
        std::cout << "Nodes: " << gridData.map.nNode << std::endl;
        std::cout << "Element Groups: " << gridData.map.nElemTot << std::endl;
        std::cout << "Dimensions: " << gridData.map.Dim << std::endl;
        
        // Access node coordinates
        std::cout << "\n=== First 5 Node Coordinates ===" << std::endl;
        for (int i = 0; i < gridData.map.nNode; i++) {
            std::cout << "Node " << i << ": ";
            std::cout << "x=" << gridData.map.nodeXYZ[i][0] << ", ";
            std::cout << "y=" << gridData.map.nodeXYZ[i][1] << ", ";
            std::cout << "z=" << gridData.map.nodeXYZ[i][2] << std::endl;
        }
        
        // Access boundary group information
        std::cout << "\n=== Boundary Groups ===" << std::endl;
        for (int i = 0; i < gridData.map.nBGroup; i++) {
            std::cout << "Group " << i << ": " << gridData.boundaryGroup.Title[i];
            std::cout << " (Faces: " << gridData.boundaryGroup.nBFace[i] << ")" << std::endl;
        }
        
        // Access element group information
        std::cout << "\n=== Element Groups ===" << std::endl;
        int nGroupCount = 0;
        int iter = 0;
        // for (int i = 0; i < gridData.map.nElemTot; i++) {
        while (nGroupCount < gridData.map.nElemTot) {
            std::cout << "Group " << iter << ": ";
            std::cout << "Elements=" << gridData.elementGroup.nElem[iter] << ", ";
            std::cout << "Order=" << gridData.elementGroup.order[iter] << ", ";
            std::cout << "Basis=" << gridData.elementGroup.basis[iter] << std::endl;
            iter = iter+1;
            nGroupCount = nGroupCount+gridData.elementGroup.nElem[iter];
        }
        
        // Access periodic group information
        if (gridData.map.nPG > 0) {
            std::cout << "\n=== Periodic Groups ===" << std::endl;
            for (int i = 0; i < gridData.map.nPG; i++) {
                std::cout << "Group " << i << ": ";
                std::cout << gridData.periodicGroup.periodicity[i];
                std::cout << " (Node pairs: " << gridData.periodicGroup.nPGNode[i] << ")" << std::endl;
                for (int j=0;j<gridData.periodicGroup.nPGNode[i];j++) {
                    std::cout << "x: " << gridData.periodicGroup.NP[i][j][0] << "   y: " <<gridData.periodicGroup.NP[i][j][1] << std::endl;
                }
            }
        }
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}