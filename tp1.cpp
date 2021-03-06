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
#include <queue>
#include "eigen/Eigen/Dense"

#include "CImg.h"

using namespace Eigen;
using namespace cimg_library;
using namespace std;


bool file_exists(const string name) {
	ifstream f(name.c_str());
	return f.good();
}

//Convert eigen VectorXd to CImg and saves it.
CImg<unsigned char> VectorXd_to_CImg(VectorXd v,int image_width, int image_height, double min, double max)
{
	CImg<unsigned char> image(image_width, image_height, 1, 1, 0);
	if ((max - min) == 0)
		max = 255.0;
	for (int i = 0; i < image_height; i++)
	{
		for (int j = 0; j < image_width; j++)
		{
			double pix = v((i*image_width) + j, 0);
			int rescale = int(((pix - min) / (max - min))*255.0);
			image(j, i) = rescale;
		}
	}
	return image;
}

//Read image directory infos.
void read_txt_file(string data_base, string &data_base_directory, int &N, int &image_width, int &image_height)
{
	string line;
	
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
}

//Convert CImg to a eigen MatrixXd.
MatrixXd images_to_MatrixXd(int image_width, int image_height, int N, string directory)
{
	MatrixXd F(image_width*image_height, N);

	for (int i = 0; i < N; i++)
	{
		string input = directory + to_string(i) + ".pgm";

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
			return F;
		}
	}
	return F;
}

//Build eigen faces from N images from a image database.
void build_eigen_faces(string data_base, string eigen_faces_directory, vector<VectorXd> &Eigen_faces)
{
	string data_base_directory;
	int N = 0;
	int image_width, image_height;
	read_txt_file(data_base, data_base_directory, N, image_width, image_height);

	MatrixXd F = images_to_MatrixXd(image_width, image_height, N, data_base_directory);

	VectorXd ONES(N);
	ONES.setOnes();

	VectorXd m = F.rowwise().mean();

	Eigen_faces.push_back(m);
	CImg<unsigned char> mean = VectorXd_to_CImg(m, image_width, image_height, 0.0, 255.0);
	mean.save((eigen_faces_directory + "/mean.pgm").c_str());

	// Q : la matrice de covariance. Merci les mathématiques :)))
	MatrixXd A = F - (m * ONES.transpose());
	MatrixXd Q = A.transpose() * A;

	EigenSolver<MatrixXd> ES(Q);

	MatrixXd eigen_vectors = A * ES.eigenvectors().real();

	for (int n = 0; n < N; n++)
	{
		VectorXd eigen_vector = eigen_vectors.col(n);
		Eigen_faces.push_back(eigen_vector);

		double min = eigen_vector.minCoeff();
		double max = eigen_vector.maxCoeff();
		CImg<unsigned char> eigen_face = VectorXd_to_CImg(eigen_vector, image_width, image_height, min, max);
		eigen_face.save((eigen_faces_directory + "/eigen_face"+ to_string(n) +".pgm").c_str());
	}
	return;
}

//Get eigen faces from folder if already computed.
void get_eigen_faces(vector<VectorXd> &Eigen_faces, string eigen_faces_directory, int n_eigen_faces)
{
	if (file_exists(eigen_faces_directory + "/mean.pgm"))
	{
		CImg<unsigned char> mean((eigen_faces_directory + "/mean.pgm").c_str());
		VectorXd m(mean.width()*mean.height());
		cimg_forXY(mean, x, y)
		{
			m((y*mean.width()) + x) = mean(x, y);
		}
		Eigen_faces.push_back(m);
	}
	for (int n = 0; n < n_eigen_faces; n++)
	{
		string input = eigen_faces_directory + "/eigen_face" + to_string(n) + ".pgm";

		if (file_exists(input))
		{
			CImg<unsigned char> eigen_face(input.c_str());
			VectorXd v(eigen_face.width()*eigen_face.height());
			cimg_forXY(eigen_face, x, y)
			{
				v((y*eigen_face.width()) + x) = eigen_face(x, y);
			}
			Eigen_faces.push_back(v);
		}
		else
		{
			cout << "Unable to open " + input;
			return;
		}
	}
}

//Compute images eigen faces weigth.
MatrixXd image_weight(string images, vector<VectorXd> Eigen_faces, int n_eigen_faces)
{
	string directory;
	int N = 0;
	int image_width, image_height;
	read_txt_file(images, directory, N, image_width, image_height);
	
	MatrixXd I = images_to_MatrixXd(image_width, image_height, N, directory);
	MatrixXd weight(n_eigen_faces, I.cols());

	double w;
	VectorXd image_norm(I.rows());
	for (int i = 0; i < I.cols(); i++)
	{
		image_norm = (I.col(i) - Eigen_faces[0]);
		image_norm.normalize();
		for (int m = 1; m < n_eigen_faces + 1; m++)
		{
			w = Eigen_faces[m].transpose() * image_norm;
			weight(m - 1, i) = w;
		}
	}

	return weight;
}

//restore image from their weight.
void restore_image(vector<VectorXd> Eigen_faces, MatrixXd image_weight, string directory, int n_eigen_faces)
{
	VectorXd image(image_weight.rows());

	for (int i = 0; i < image_weight.cols(); i++)
	{
		image = Eigen_faces[0];
		for (int k = 0; k < n_eigen_faces; k++)
		{
			image += image_weight(k, i)*Eigen_faces[k+1];
		}

		int size = sqrt(image.rows());
		double min = image.minCoeff();
		double max = image.maxCoeff();
		CImg<unsigned char> eigen_face = VectorXd_to_CImg(image, size, size, min, max);
		eigen_face.save((directory + "restored_" + to_string(i) + ".pgm").c_str());
	}
	
}

//Finds k min from vector.
vector<int> k_min(vector<double> distances, int k)
{
	vector<int> index;
	vector<pair<double, int>> v_pair;

	for (int i = 0; i < distances.size(); i++)
	{
		pair<double, int> pair;
		pair.first = distances[i];
		pair.second = i;
		v_pair.push_back(pair);
	}
	sort(v_pair.begin(), v_pair.end());

	for (int i = 0; i < k; i++)
	{
		index.push_back(v_pair[i].second);
	}

	return index;
}

//Find k similar image from the database.
void find_similar_image(int k, string similar_directory, string data_base, vector<VectorXd> Eigen_faces, int n_eigen_faces, bool restore, string restore_dir)
{
	MatrixXd database_weight = image_weight(data_base, Eigen_faces, n_eigen_faces);
	MatrixXd compared_image_weight = image_weight(similar_directory, Eigen_faces, n_eigen_faces);

	if (restore)
	{
		restore_image(Eigen_faces, database_weight, restore_dir, n_eigen_faces);
	}

	VectorXd dif(compared_image_weight.rows());
	vector<double> distances;

	string dir;
	ifstream data(similar_directory);
	if (data.is_open())
	{
		getline(data, dir);
		data.close();
	}
	else
	{
		cout << "Unable to open file";
		return;
	}

	for (int i = 0; i < compared_image_weight.cols(); i++)
	{
		for (int d = 0; d < database_weight.cols(); d++)
		{
			double distance = 0.0;
			dif = compared_image_weight.col(i) - database_weight.col(d);
			for (int j = 0; j < dif.rows(); j++)
			{
				distance += dif(j)*dif(j);
			}
			distance = sqrt(distance);
			distances.push_back(distance);
		}

		vector<int> index = k_min(distances, k);

		ofstream myfile;
		myfile.open(dir + "similar_to_" + to_string(i+1) + ".pgm_list.txt");
		myfile << "Les " + to_string(k) + " images les plus similaire, de la plus similaise a la moins similaire, sont: \n";
		for (int i = 0; i < index.size(); i++)
		{
			myfile << to_string(index[i]) + ".pgm \n";
		}
	}
	
}

int main(int argc, const char * argv[]) {

	//vector of eigen faces, beginning by the mean.
	vector<VectorXd> Eigen_faces;
	vector<double> weigth;

	int build = atoi(argv[1]);
	int get_ef = atoi(argv[2]);
	int similar = atoi(argv[3]);
	string data_base(argv[4]);
	int n_eigen_faces = atoi(argv[5]);
	string eigen_faces_directory(argv[6]);
	string similar_directory(argv[7]);
	int k = atoi(argv[8]);
	int restore = atoi(argv[9]);
	string restore_directory(argv[10]);

	if (build)
		build_eigen_faces(data_base, eigen_faces_directory, Eigen_faces);
	else if (get_ef)
		get_eigen_faces(Eigen_faces, eigen_faces_directory, n_eigen_faces);
	else
	{
		cout << "Either build or get eigen_faces";
		return 0;
	}
	if (similar)
		find_similar_image(k, similar_directory, data_base, Eigen_faces, n_eigen_faces, restore, restore_directory);
	
	return 0;
}
