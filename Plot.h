#pragma once

#ifndef __Plot_h__
#define __Plot_h__

#include "matplotlibcpp.h"
#include "NumericalMethods.h"
#include "fileIO.h"

namespace plt = matplotlibcpp;

void plot2D(std::vector<double> &x, std::vector<double> &y, std::string label) {
	int N = std::min((int)x.size(), (int)y.size());
	std::vector<double> xAxis;
	std::vector<double> yAxis;

	xAxis.resize(N, 0.0);
	yAxis.resize(N, 0.0);

	std::vector<int> indexT;
	indexT.resize(N);
	int countT = 0;
	std::generate(std::begin(indexT), std::end(indexT), [&] {
		return countT++;
	});

	concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
		xAxis[i] = x[i];
		yAxis[i] = y[i];
	});


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

	std::vector<int> indexT;
	indexT.resize(N);
	int countT = 0;
	std::generate(std::begin(indexT), std::end(indexT), [&] {
		return countT++;
		});

	concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
		xAxis[i] = x[i];
		yAxis[i] = y[i];
	});

	plt::semilogy(xAxis, yAxis);
}

void plot3D(std::vector<double>& x, std::vector<double>& y, imageBW& z, std::string label) {
	int Nx = std::min((int)x.size(), z.sizeX());
	int Ny = std::min((int)y.size(), z.sizeY());
	std::vector<std::vector<double>> xAxis;
	std::vector<std::vector<double>> yAxis;
	std::vector<std::vector<double>> zAxis;

	xAxis.resize(Nx);
	yAxis.resize(Nx);
	zAxis.resize(Nx);

	std::vector<int> indexT;
	indexT.resize(Nx);
	int countT = 0;
	std::generate(std::begin(indexT), std::end(indexT), [&] {
		return countT++;
	});

	concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
		xAxis[i].resize(Ny, 0.0);
		yAxis[i].resize(Ny, 0.0);
		zAxis[i].resize(Ny, 0.0);
		for (int j = 0; j < Ny; j++) {
			xAxis[i][j] = x[i];
			yAxis[i][j] = y[j];
			zAxis[i][j] = z.value(i,j);
		}
	});

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

	std::vector<int> indexT;
	indexT.resize(Nx);
	int countT = 0;
	std::generate(std::begin(indexT), std::end(indexT), [&] {
		return countT++;
		});

	concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
		xAxis[i].resize(Ny, 0.0);
		yAxis[i].resize(Ny, 0.0);
		zAxis[i].resize(Ny, 0.0);
		for (int j = 0; j < Ny; j++) {
			xAxis[i][j] = x[i];
			yAxis[i][j] = y[j];
			zAxis[i][j] = z[i][j];
		}
		});

	std::map<std::string, std::string> keywords;
	keywords.insert(std::pair<std::string, std::string>("label", label));
	plt::plot_surface(xAxis, yAxis, zAxis, keywords);
}

void imshow(imageBW& image, bool invert, std::string label) {
	int Nx = image.sizeX();
	int Ny = image.sizeY();
	std::vector<float> buffer(Nx * Ny);

	std::vector<int> indexT;
	indexT.resize(Nx);
	int countT = 0;
	std::generate(std::begin(indexT), std::end(indexT), [&] {
		return countT++;
	});

	if (invert) {
		concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
			int Nj;
			for (int j = 0; j < Ny; j++) {
				Nj = Ny - 1 - j;
				buffer.at(Nx * j + i) = (float)image.value(i, Nj);
			}
			});
	}
	else {
		concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
			int Nj;
			for (int j = 0; j < Ny; j++) {
				Nj = j;
				buffer.at(Nx * j + i) = (float)image.value(i, Nj);
			}
			});
	}

	const float* buffer_ptr = &(buffer[0]);
	const int colors = 1;

	std::map<std::string, std::string> keywords;
	keywords.insert(std::pair<std::string, std::string>("label", label));
	plt::imshow(buffer_ptr, Ny, Nx, colors, keywords);
}

#endif