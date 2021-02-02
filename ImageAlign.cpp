//
//  ImageAlign.cpp
//

#include <iostream>
#include <vector>
#include <cmath>

#include "eigen/Eigen/Dense"

#define cimg_display 0
#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;

MatrixXd ImgToMatrixTresholder(CImg<unsigned char> image) {   


    int x = 0;
    int y = 0;

    vector<vector<int>> vec;
    CImg<unsigned char> img_to_grayscale = image.RGBtoYCbCr().channel(0);
    
    cimg_forXY(img_to_grayscale,x,y) { 

        
        float intensite = img_to_grayscale(x,y);

        if (intensite > 40) { 
            
            vector<int> coordonate;
            coordonate = {x,y};
            vec.push_back(coordonate);
            
        }
     }

     MatrixXd m(2,vec.size());

     for (int i = 0; i < vec.size(); i++)  {
         m(0,i) = vec[i][0];
         m(1,i) = vec[i][1];
     }

    return m;
}

int main(int argc, const char * argv[]) {
    
    vector<CIMG<unsigned chard>> faces;

    cv::glob("./lfw")


    
    CImg<unsigned char> image1("brain1.ppm");
    CImg<unsigned char> image2("brain2.ppm");

    MatrixXd image_filtered_1 = ImgToMatrixTresholder(image1);
    MatrixXd image_filtered_2 = ImgToMatrixTresholder(image2);

    VectorXd A = image_filtered_1.rowwise().mean();
    VectorXd ONES(image_filtered_1.size()); ONES.setOnes();    

    // N : nombre d'observations
    // F : 2 * N la matrice d'observations
    // A : 2 * 1 la matrice de moyenne
    // Q : la matrice de covariance.
    MatrixXd Q = (image_filtered_1 - A * ONES.transpose()) * (image_filtered_1 - A * ONES.transpose()).transpose() / (image_filtered_1.size() - 1);

    // Décomposition en valeur propre
    // Vp dans ES
    // Attention objet complexe dans notre ES
    EigenSolver<MatrixXd> ES(Q);

    int maxIndex;
    ES.eigenvalues().real().maxCoeff(&maxIndex);

    VectorXd eigen1((ES.eigenvectors().col(maxIndex)).real());
    VectorXd align1; align1 << 0,-1; align1.normalize();

    double prodscal(eigen1.dot(align1));


    image1.save("output.ppm");
    
    return 0;
}
