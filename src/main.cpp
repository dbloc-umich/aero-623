// #include "GaussKronrod.h"
#include "../include/ReadGRI.h"
#include "../include/ReadConnData.h"
#include "../include/project1Task3.h"
#include <iostream>

int main(){
    // Testing GRI input
    /*try {
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
    } */


    // Mesh verification on test mesh
    std::string testGriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/test.gri";
    std::vector<std::string> testTxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testperiodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testI2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testIn.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testB2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testBn.txt"};
    std::cout << "------------Mesh Verification for Test Grid------------" <<std::endl;
    meshVerification(testGriFile, testTxtFiles);
    std::cout <<std::endl;
    
    // Mesh verification on coarse mesh
    std::string coarseGriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394.gri";
    std::vector<std::string> coarseTxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394periodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394I2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394In.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394B2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394Bn.txt"};
    std::cout << "------------Mesh Verification for Coarse Grid------------" <<std::endl;
    meshVerification(coarseGriFile, coarseTxtFiles);
    std::cout <<std::endl;

    std::string refinedGriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_9576.gri";
    std::vector<std::string> refinedtxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_9576periodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_9576I2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_9576In.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_9576B2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_9576Bn.txt"};
    std::cout << "------------Mesh Verification for Refined Grid------------" <<std::endl;
    meshVerification(refinedGriFile, refinedtxtFiles);
    std::cout <<std::endl;

}