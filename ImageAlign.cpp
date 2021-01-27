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

        if (intensite > 0.4) { 
            
            vector<int> coordonate;
            coordonate = {x,y};
            vec.push_back(coordonate);
            
        }
     }
     MatrixXd m(2,vec.size());

    return m;
}

int main(int argc, const char * argv[]) {
    

    CImg<unsigned char> image1("brain1.ppm");
    CImg<unsigned char> image2("brain2.ppm");

    MatrixXd image_filtered_1 = ImgToMatrixTresholder(image1);
    MatrixXd image_filtered_2 = ImgToMatrixTresholder(image2);

    VectorXd image;

    // eigen sommation integree fct_mean

    // Mettre l'image 1 dans la même orientation que l'image 2 à l'aide de l'ACP

    image1.save("output.ppm");
    
    return 0;
}

//~!@#$%^&*()_+:"|"?><+_=-=-`]''\'/<
