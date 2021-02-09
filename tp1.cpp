/*
Alexandre Joanisse		joaa1801
Pierre-Emmanuel Goffi	gofp2701
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <string>
#include "eigen/Eigen/Dense"

#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;


bool file_exists(const string name) {
	ifstream f(name.c_str());
	return f.good();
}

CImg<unsigned char> VectorXd_to_CImg(VectorXd m,int image_width, int image_height)
{
	CImg<unsigned char> image(image_width, image_height, 1, 1, 0);
	for (int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < image_width; j++)
		{
			image(j, i) = m((i*image_width) + j, 0);
		}
	}
	return image;
}

void build_reduced_space(string data_base, string eigen_faces_directory, vector<CImg<unsigned char>> &Eigen_faces)
{
	string line, data_base_directory;
	int N = 0;
	int image_width, image_height;
	ifstream data(data_base);
	if (data.is_open())
	{
		getline(data, data_base_directory);

		getline(data, line);
		N = stoi(line);

		getline(data, line);
		image_width = stoi(line);

		getline(data, line);
		image_height = stoi(line);

		data.close();
	}
	else
	{
		cout << "Unable to open file";
		return;
	}

	MatrixXd F(image_width*image_height, N);

	for (int i = 0; i < N; i++) {


		string input = data_base_directory + to_string(i) + ".pgm";

		if (file_exists(input))
		{
			CImg<unsigned char> face(input.c_str());
			face.resize(image_width, image_height);

			cimg_forXY(face, x, y)
			{
				F((y*face.width()) + x, i) = face(x, y);
			}
		}
		else
		{
			cout << "Unable to open " + input;
			return;
		}
	}

	VectorXd ONES(N);
	ONES.setOnes();

	VectorXd m = F.rowwise().mean();

	CImg<unsigned char> mean = VectorXd_to_CImg(m, image_width, image_height);
	Eigen_faces.push_back(mean);
	mean.save((eigen_faces_directory + "/mean.pgm").c_str());

	// N : nombre d'observations
	// F : k^2 x N la matrice d'observations
	// m : 1 x N la matrice de moyenne
	// Q : la matrice de covariance.
	MatrixXd Q = ((F - (m * ONES.transpose())) * (F - (m * ONES.transpose())).transpose()) / (N - 1);

	// Calcule valeur propre, vecteur propre
	EigenSolver<MatrixXd> ES(Q);

	// Conserve les reels, maxcoeff donne position de la plus grande valeur propre
	int maxIndex;
	ES.eigenvalues().real().maxCoeff(&maxIndex);

	// stocke
	VectorXd eigen1((ES.eigenvectors().col(maxIndex)).real());

	VectorXd lambda = ES.eigenvalues().real();

	MatrixXd v = ES.eigenvectors().real();
	
	return;
}

void get_eigen_faces()
{
	//TODO
}

void image_to_reduce_space() 
{
	//TODO
}

void find_similar_image()
{
	//TODO
}

int main(int argc, const char * argv[]) {

	//vector of eigen faces, beginning by the mean.
	vector<CImg<unsigned char>> Eigen_faces;

	string b(argv[1]);
	string e(argv[2]);
	string s(argv[3]);
	int build = stoi(b);
	int get_ef = stoi(e);
	int similar = stoi(s);
	string data_base(argv[4]);
	string eigen_faces_directory(argv[5]);

	if (build)
			build_reduced_space(data_base, eigen_faces_directory, Eigen_faces);
	if (get_ef)
		get_eigen_faces();
	if (similar)
	{
		string similar_directory(argv[6]);
		find_similar_image();
	}
	
	return 0;
}
