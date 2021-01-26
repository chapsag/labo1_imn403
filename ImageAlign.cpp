//
//  ImageAlign.cpp
//

#include <iostream>
#include <vector>
#include <cmath>

#include "Eigen/Dense"

#define cimg_display 0
#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;

int main(int argc, const char * argv[]) {
    
    CImg<unsigned char> image1("brain1.ppm");
    CImg<unsigned char> image2("brain2.ppm");
    
    // Mettre l'image 1 dans la même orientation que l'image 2 à l'aide de l'ACP

    image1.save("output.ppm");
    
    return 0;
}
