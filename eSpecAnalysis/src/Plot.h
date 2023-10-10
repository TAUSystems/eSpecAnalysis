#pragma once

#ifndef __Plot_h__
#define __Plot_h__

#include "matplotlibcpp.h"
#include "NumericalMethods.h"

namespace plt = matplotlibcpp;

void plot2D(std::vector<double> &x, std::vector<double> &y, std::string label) {
	int N = std::min((int)x.size(), (int)y.size());
	std::vector<double> xAxis;
	std::vector<double> yAxis;

	xAxis.resize(N, 0.0);
	yAxis.resize(N, 0.0);

	#pragma omp parallel 
	{
		for (int i = 0; i < N; i++) {
			xAxis[i] = x[i];
			yAxis[i] = y[i];
		}
	}

	std::map<std::string, std::string> keywords;
	keywords.insert(std::pair<std::string, std::string>("label", label));
	plt::plot(xAxis, yAxis, keywords);
}

void plot2DLog(std::vector<double>& x, std::vector<double>& y) {
	int N = std::min((int)x.size(), (int)y.size());
	std::vector<double> xAxis;
	std::vector<double> yAxis;

	xAxis.resize(N, 0.0);
	yAxis.resize(N, 0.0);

	#pragma omp parallel 
	{
		for (int i = 0; i < N; i++) {
			xAxis[i] = x[i];
			yAxis[i] = y[i];
		}
	}

	plt::semilogy(xAxis, yAxis);
}

void plot3D(std::vector<double>& x, std::vector<double>& y, imageBW& z, std::string label) {
	int Nx = std::min((int)x.size(), (int)z.sizeX());
	int Ny = std::min((int)y.size(), (int)z.sizeY());
	std::vector<std::vector<double>> xAxis;
	std::vector<std::vector<double>> yAxis;
	std::vector<std::vector<double>> zAxis;

	xAxis.resize(Nx);
	yAxis.resize(Nx);
	zAxis.resize(Nx);

	#pragma omp parallel for
		for (int i = 0; i < Nx; i++) {
			xAxis[i].resize(Ny, 0.0);
			yAxis[i].resize(Ny, 0.0);
			zAxis[i].resize(Ny, 0.0);
			for (int j = 0; j < Ny; j++) {
				xAxis[i][j] = x[i];
				yAxis[i][j] = y[j];
				zAxis[i][j] = z.value(i, j);
			}
		}

	std::map<std::string, std::string> keywords;
	keywords.insert(std::pair<std::string, std::string>("label", label));
	plt::contour(xAxis, yAxis, zAxis, keywords);
}

void plot3D(std::vector<double>& x, std::vector<double>& y, std::vector<std::vector<double>>& z, std::string label) {
	int Nx = std::min((int)x.size(), (int)z[0].size());
	int Ny = std::min((int)y.size(), (int)z[0].size());
	std::vector<std::vector<double>> xAxis;
	std::vector<std::vector<double>> yAxis;
	std::vector<std::vector<double>> zAxis;

	xAxis.resize(Nx);
	yAxis.resize(Nx);
	zAxis.resize(Nx);

	#pragma omp parallel for
		for (int i = 0; i < Nx; i++) {
			xAxis[i].resize(Ny, 0.0);
			yAxis[i].resize(Ny, 0.0);
			zAxis[i].resize(Ny, 0.0);
			for (int j = 0; j < Ny; j++) {
				xAxis[i][j] = x[i];
				yAxis[i][j] = y[j];
				zAxis[i][j] = z[i][j];
			}
		}

	std::map<std::string, std::string> keywords;
	keywords.insert(std::pair<std::string, std::string>("label", label));
	plt::plot_surface(xAxis, yAxis, zAxis, keywords);
}

void pltimshow(imageBW& image, bool invert, std::string label) {
	int Nx = (int)image.sizeX();
	int Ny = (int)image.sizeY();
	std::vector<float> buffer(Nx * Ny);

	if (invert) {
		#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				int Nj;
				for (int j = 0; j < Ny; j++) {
					Nj = Ny - 1 - j;
					buffer.at(Nx * j + i) = (float)image.value(i, Nj);
				}
			}
	}
	else {
		#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				int Nj;
				for (int j = 0; j < Ny; j++) {
					Nj = j;
					buffer.at(Nx * j + i) = (float)image.value(i, Nj);
				}
			}
	}

	const float* buffer_ptr = &(buffer[0]);
	const int colors = 1;

	std::map<std::string, std::string> keywords;
	keywords.insert(std::pair<std::string, std::string>("label", label));
	plt::imshow(buffer_ptr, Ny, Nx, colors, keywords);
}

#endif