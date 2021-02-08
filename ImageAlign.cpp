//
//  ImageAlign.cpp
//


#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <string>
#include "eigen/Eigen/Dense"

#define cimg_display 0
#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;


int main(int argc, const char * argv[]) {

    vector<int> xCoord, yCoord;
    
    for (int i = 1; i < 101; i++) {
        
        string input ("./lfw/");
        string output ("output");

        input += std::to_string(i);
        input += ".jpg";

        output += std::to_string(i);
        output += ".ppm";

        CImg<unsigned char> face(input.c_str());
        
        face.crop(0,0,200,200);
        face = face.get_RGBtoYCbCr().get_channel(0);
        face.save(output.c_str());

        
        cimg_forXY(face, x, y) {
            
            xCoord.push_back(x);
            yCoord.push_back(y);
        }
    }



    // Start PCA
    int N = xCoord.size();

    MatrixXd F(2,N);


    for(int i = 0; i < N; i++) {
        F(0,i) = xCoord[i];
        F(1,i) = yCoord[i];
    }

    Vector2d A = F.rowwise().mean();
    VectorXd ONES(N); ONES.setOnes();

    // N : nombre d'observations
    // F : 2 * N la matrice d'observations
    // A : 2 * 1 la matrice de moyenne
    // Q : la matrice de covariance.
    MatrixXd Q = (F - A * ONES.transpose()) * (F - A * ONES.transpose()).transpose() / (N - 1);

    // Calcule valeur propre, vecteur propre
    EigenSolver<MatrixXd> ES(Q);

    int maxIndex;

    // Conserve les reels, maxcoeff donne position de la plus grande valeur propre
    ES.eigenvalues().real().maxCoeff(&maxIndex);


    // stocke
    VectorXd eigen1((ES.eigenvectors().col(maxIndex)).real());

    return 0;
    }
