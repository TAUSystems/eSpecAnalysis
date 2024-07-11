#pragma once

#ifndef __Analysis_h__
#define __Analysis_h__

#include "Calibration.h"

/*
void findSignalPeak(imageBW& image, std::vector<int>& peak, double& peakValue) {
	imageBW imSmooth = image;
	medianFilter(imSmooth, 2);
	int Nx = imSmooth.sizeX();
	int Ny = imSmooth.sizeY();
	peak.resize(2, 0);
	peakValue = imSmooth.value(0, 0);
	peak[0] = 0;
	peak[1] = 0;

	std::vector<int> linePeaksX,linePeaksY;
	std::vector<double> lineValuesX,lineValuesY;
	linePeaksX.resize(Nx, 0);
	lineValuesX.resize(Nx, 0.0);
	#pragma omp parallel for
	for (int i = 0; i < Nx; i++) {
		std::vector<double> line;
		line.resize(Ny, 0.0);
		double lineMin = 255.0;
		int linePI = linePeaksX[i];
		double linePV = imSmooth.value(i,linePI);
		double lineMean = 0.0;
		double lineSTD = 0.0;

		//Find minimum value
		for (int j = 0; j < Ny; j++) {
			line[j] = imSmooth.value(i, j);
			lineMean = lineMean + line[j];
			if (line[j] < lineMin) {
				lineMin = line[j];
			}
		}
		lineMean = lineMean / Ny;

		for (int j = 0; j < Ny; j++) {
			lineSTD = lineSTD + (line[j] - lineMean) * (line[j] - lineMean);
		}
		lineSTD = lineSTD / sqrt(Ny);
		lineMean = lineMean / Ny - lineMin + 1.0e-9;

		//Shift by minimum and find maximum location and value
		for (int j = 0; j < Ny; j++) {
			line[j] = line[j] - lineMin + 1.0e-9;
			if (line[j] > linePV) {
				if (line[j] < lineMean + 4.0 * lineSTD) {
					linePV = line[j];
					linePI = j;
				}
			}
		}

		//Fit to gaussian using quad regression
		double xValue, yValue, a1, a2, a0;
		int N, index;
		std::vector<double> sumX, sumY;
		cv::Mat quadRegMat = cv::Mat::zeros(3, 3, CV_64F);

		//Set fit window size
		N = (int)std::min(std::min((double)linePI, (double)(Ny - linePI)) - 1, 250.0);
		sumX.resize(5, 0.0);
		sumY.resize(3, 0.0);

		for (int i = 0; i < 2 * N + 1; i++) {
			xValue = (double)i + (double)linePI - (double)N;
			index = linePI + i - N;
			yValue = -log(line[index]/(linePV + 1.0e-9));
			
			sumX[0] = sumX[0] + 1;
			sumX[1] = sumX[1] + xValue;
			sumX[2] = sumX[2] + xValue * xValue;
			sumX[3] = sumX[3] + xValue * xValue * xValue;
			sumX[4] = sumX[4] + xValue * xValue * xValue * xValue;
			sumY[0] = sumY[0] + yValue;
			sumY[1] = sumY[1] + yValue * xValue;
			sumY[2] = sumY[2] + yValue * xValue * xValue;
		}

		quadRegMat = cv::Mat::zeros(3, 3, CV_64F);
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 0) = sumX[4];
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 1) = sumX[3];
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 2) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 0) = sumX[3];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 1) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 2) = sumX[1];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 0) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 1) = sumX[1];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 2) = sumX[0];

		quadRegMat = quadRegMat.inv();
		a2 = sumY[2] * quadRegMat.at<double>(0, 0) + sumY[1] * quadRegMat.at<double>(0, 1) + sumY[0] * quadRegMat.at<double>(0, 2);
		a1 = sumY[2] * quadRegMat.at<double>(1, 0) + sumY[1] * quadRegMat.at<double>(1, 1) + sumY[0] * quadRegMat.at<double>(1, 2);
		a0 = sumY[2] * quadRegMat.at<double>(2, 0) + sumY[1] * quadRegMat.at<double>(2, 1) + sumY[0] * quadRegMat.at<double>(2, 2);

		xValue = - a1 / (2.0 * a2);
		if (xValue >= 0.0 && xValue < Ny) {
			linePeaksX[i] = (int)xValue;
		}
		else {
			linePeaksX[i] = linePI;
		}
		lineValuesX[i] = a2 * xValue * xValue + a1 * xValue + a0;
		lineValuesX[i] = exp(-lineValuesX[i]) * (linePV + 1.0e-9) + lineMin - 1.0e-9;
	}

	for (int i = 0; i < Nx; i++) {
		if (lineValuesX[i] > peakValue) {
			peak[0] = i;
			peak[1] = linePeaksX[i];
			peakValue = lineValuesX[i];
		}
	}

	std::vector<double> line;
	line.resize(Nx, 0.0);
	double lineMin = 255.0;
	int linePI = 0;
	double linePV = imSmooth.value(linePI, peak[1]);

	//Find minimum value
	for (int j = 0; j < Nx; j++) {
		line[j] = imSmooth.value(j, peak[1]);
		if (line[j] < lineMin) {
			lineMin = line[j];
		}
	}

	//Shift by minimum and find maximum location and value
	for (int j = 0; j < Nx; j++) {
		line[j] = line[j] - lineMin + 1.0e-9;
		if (line[j] > linePV) {
			linePV = line[j];
			linePI = j;
		}
	}

	//Fit to gaussian using quad regression
	double xValue, yValue, a1, a2, a0;
	int N, index;
	std::vector<double> sumX, sumY;
	cv::Mat quadRegMat = cv::Mat::zeros(3, 3, CV_64F);

	//Set fit window size
	N = (int)std::min(std::min((double)linePI, (double)(Nx - linePI)) - 1, 250.0);
	sumX.resize(5, 0.0);
	sumY.resize(3, 0.0);

	for (int i = 0; i < 2 * N + 1; i++) {
		xValue = (double)i + (double)linePI - (double)N;
		index = linePI + i - N;
		yValue = -log(line[index] / (linePV + 1.0e-9));

		sumX[0] = sumX[0] + 1;
		sumX[1] = sumX[1] + xValue;
		sumX[2] = sumX[2] + xValue * xValue;
		sumX[3] = sumX[3] + xValue * xValue * xValue;
		sumX[4] = sumX[4] + xValue * xValue * xValue * xValue;
		sumY[0] = sumY[0] + yValue;
		sumY[1] = sumY[1] + yValue * xValue;
		sumY[2] = sumY[2] + yValue * xValue * xValue;
	}

	quadRegMat = cv::Mat::zeros(3, 3, CV_64F);
	quadRegMat.at<CvType<CV_64F>::type_t>(0, 0) = sumX[4];
	quadRegMat.at<CvType<CV_64F>::type_t>(0, 1) = sumX[3];
	quadRegMat.at<CvType<CV_64F>::type_t>(0, 2) = sumX[2];
	quadRegMat.at<CvType<CV_64F>::type_t>(1, 0) = sumX[3];
	quadRegMat.at<CvType<CV_64F>::type_t>(1, 1) = sumX[2];
	quadRegMat.at<CvType<CV_64F>::type_t>(1, 2) = sumX[1];
	quadRegMat.at<CvType<CV_64F>::type_t>(2, 0) = sumX[2];
	quadRegMat.at<CvType<CV_64F>::type_t>(2, 1) = sumX[1];
	quadRegMat.at<CvType<CV_64F>::type_t>(2, 2) = sumX[0];

	quadRegMat = quadRegMat.inv();
	a2 = sumY[2] * quadRegMat.at<double>(0, 0) + sumY[1] * quadRegMat.at<double>(0, 1) + sumY[0] * quadRegMat.at<double>(0, 2);
	a1 = sumY[2] * quadRegMat.at<double>(1, 0) + sumY[1] * quadRegMat.at<double>(1, 1) + sumY[0] * quadRegMat.at<double>(1, 2);
	a0 = sumY[2] * quadRegMat.at<double>(2, 0) + sumY[1] * quadRegMat.at<double>(2, 1) + sumY[0] * quadRegMat.at<double>(2, 2);

	xValue = -a1 / (2.0 * a2);
	if (xValue >= 0.0 && xValue < Nx) {
		peak[0] = (int)xValue;
	}
	else {
		peak[0] = linePI;
	}

	peakValue = imSmooth.value(peak[0], peak[1]);
}
*/

void findSignalPeak(imageBW& image, std::vector<int>& peak, double& peakValue) {
	imageBW imSmooth;
	image.copy(imSmooth);
	removeOutlier(imSmooth, 4.0);
	medianFilter(imSmooth, 12);
	int Nx = imSmooth.sizeX();
	int Ny = imSmooth.sizeY();
	peak.resize(2, 0);
	peakValue = imSmooth.value(0, 0);
	peak[0] = 0;
	peak[1] = 0;

	std::vector<double> lineValuesX, lineValuesY;
	lineValuesX.resize(Nx, 0.0);
	lineValuesY.resize(Ny, 0.0);
	
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			lineValuesX[i] = lineValuesX[i] + imSmooth.value(i, j) * imSmooth.value(i, j) / Nx;
			lineValuesY[j] = lineValuesY[j] + imSmooth.value(i, j) * imSmooth.value(i, j) / Ny;
		}
	}

	//Find X peak
	double lineMin = 255.0;
	double linePV = 0.0;
	int linePI = 0;

	//Find minimum value
	for (int i = 0; i < Nx; i++) {
		if (lineValuesX[i] < lineMin) {
			lineMin = lineValuesX[i];
		}
	}

	//Shift by minimum and find maximum location and value
	for (int i = 0; i < Nx; i++) {
		lineValuesX[i] = lineValuesX[i] - lineMin + 1.0e-9;
		if (lineValuesX[i] > linePV) {
			linePV = lineValuesX[i];
			linePI = i;
		}
	}

	double threshold = 10;
	if (linePV > threshold) {
		//Fit to gaussian using quad regression
		double xValue, yValue, a1, a2, a0;
		int N, index;
		std::vector<double> sumX, sumY;
		cv::Mat quadRegMat = cv::Mat::zeros(3, 3, CV_64F);

		//Set fit window size
		N = (int)std::min(std::min((double)linePI, (double)(Nx - linePI)) - 1, 250.0);
		sumX.resize(5, 0.0);
		sumY.resize(3, 0.0);

		for (int i = 0; i < 2 * N + 1; i++) {
			xValue = (double)i + (double)linePI - (double)N;
			index = linePI + i - N;
			yValue = -log(lineValuesX[index] / (linePV + 1.0e-9));

			sumX[0] = sumX[0] + 1;
			sumX[1] = sumX[1] + xValue;
			sumX[2] = sumX[2] + xValue * xValue;
			sumX[3] = sumX[3] + xValue * xValue * xValue;
			sumX[4] = sumX[4] + xValue * xValue * xValue * xValue;
			sumY[0] = sumY[0] + yValue;
			sumY[1] = sumY[1] + yValue * xValue;
			sumY[2] = sumY[2] + yValue * xValue * xValue;
		}

		quadRegMat = cv::Mat::zeros(3, 3, CV_64F);
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 0) = sumX[4];
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 1) = sumX[3];
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 2) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 0) = sumX[3];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 1) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 2) = sumX[1];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 0) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 1) = sumX[1];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 2) = sumX[0];

		quadRegMat = quadRegMat.inv();
		a2 = sumY[2] * quadRegMat.at<double>(0, 0) + sumY[1] * quadRegMat.at<double>(0, 1) + sumY[0] * quadRegMat.at<double>(0, 2);
		a1 = sumY[2] * quadRegMat.at<double>(1, 0) + sumY[1] * quadRegMat.at<double>(1, 1) + sumY[0] * quadRegMat.at<double>(1, 2);
		a0 = sumY[2] * quadRegMat.at<double>(2, 0) + sumY[1] * quadRegMat.at<double>(2, 1) + sumY[0] * quadRegMat.at<double>(2, 2);

		xValue = -a1 / (2.0 * a2);
		if (xValue >= 0.0 && xValue < Nx) {
			peak[0] = (int)xValue;
		}
		else {
			peak[0] = linePI;
		}
		quadRegMat.release();
		sumX.clear();
		sumY.clear();
	}
	else {
		peak[0] = linePI;
	}

	//Find Y peak
	lineMin = 255.0;
	linePV = 0.0;
	linePI = 0;

	//Find minimum value
	for (int i = 0; i < Ny; i++) {
		if (lineValuesY[i] < lineMin) {
			lineMin = lineValuesY[i];
		}
	}

	//Shift by minimum and find maximum location and value
	for (int i = 0; i < Ny; i++) {
		lineValuesY[i] = lineValuesY[i] - lineMin + 1.0e-9;
		if (lineValuesY[i] > linePV) {
			linePV = lineValuesY[i];
			linePI = i;
		}
	}

	if(linePI > threshold){
		//Fit to gaussian using quad regression
		double xValue, yValue, a1, a2, a0;
		int N, index;
		std::vector<double> sumX, sumY;
		cv::Mat quadRegMat = cv::Mat::zeros(3, 3, CV_64F);

		//Set fit window size
		N = (int)std::min(std::min((double)linePI, (double)(Ny - linePI)) - 1, 250.0);
		sumX.resize(5, 0.0);
		sumY.resize(3, 0.0);

		for (int i = 0; i < 2 * N + 1; i++) {
			xValue = (double)i + (double)linePI - (double)N;
			index = linePI + i - N;
			yValue = -log(lineValuesY[index] / (linePV + 1.0e-9));

			sumX[0] = sumX[0] + 1;
			sumX[1] = sumX[1] + xValue;
			sumX[2] = sumX[2] + xValue * xValue;
			sumX[3] = sumX[3] + xValue * xValue * xValue;
			sumX[4] = sumX[4] + xValue * xValue * xValue * xValue;
			sumY[0] = sumY[0] + yValue;
			sumY[1] = sumY[1] + yValue * xValue;
			sumY[2] = sumY[2] + yValue * xValue * xValue;
		}

		quadRegMat = cv::Mat::zeros(3, 3, CV_64F);
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 0) = sumX[4];
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 1) = sumX[3];
		quadRegMat.at<CvType<CV_64F>::type_t>(0, 2) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 0) = sumX[3];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 1) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(1, 2) = sumX[1];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 0) = sumX[2];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 1) = sumX[1];
		quadRegMat.at<CvType<CV_64F>::type_t>(2, 2) = sumX[0];

		quadRegMat = quadRegMat.inv();
		a2 = sumY[2] * quadRegMat.at<double>(0, 0) + sumY[1] * quadRegMat.at<double>(0, 1) + sumY[0] * quadRegMat.at<double>(0, 2);
		a1 = sumY[2] * quadRegMat.at<double>(1, 0) + sumY[1] * quadRegMat.at<double>(1, 1) + sumY[0] * quadRegMat.at<double>(1, 2);
		a0 = sumY[2] * quadRegMat.at<double>(2, 0) + sumY[1] * quadRegMat.at<double>(2, 1) + sumY[0] * quadRegMat.at<double>(2, 2);

		xValue = -a1 / (2.0 * a2);
		if (xValue >= 0.0 && xValue < Nx) {
			peak[1] = (int)xValue;
		}
		else {
			peak[1] = linePI;
		}
		quadRegMat.release();
		sumX.clear();
		sumY.clear();
	}
	else {
		peak[1] = linePI;
	}
	peakValue = imSmooth.value(peak[0], peak[1]);
	imSmooth.destroy();
	lineValuesX.clear();
	lineValuesY.clear();
	
}

void findPointing(spectrometer& eSpec, imageBW& image, std::vector<double>& rulerX, std::vector<double>& rulerY, std::vector<double>& pointing) {
	std::vector<int> peak;
	double peakValue;
	findSignalPeak(image, peak, peakValue);

	pointing.resize(2, 0.0);

	double x, y, z, phi, theta;
	phi = eSpec.phi(0);
	theta = eSpec.theta(0);
	
	
	z = eSpec.z(0) + rulerX[peak[0]] * std::cos(phi) + rulerY[peak[1]] * std::sin(theta);
	x = eSpec.x(0) + rulerX[peak[0]] * std::sin(phi);
	y = eSpec.y(0) + rulerY[peak[1]] *std::cos(theta);

	pointing[0] = std::atan2(z, x);
	pointing[1] = std::atan2(z, y);
	peak.clear();
}

void mRadAxis(spectrometer& eSpec, int& screen, std::vector<double>& rulerX, std::vector<double>& rulerY) {
	int Nx = (int)rulerX.size();
	int Ny = (int)rulerY.size();
	int N = std::max(Nx, Ny);

	double x, y, z, phi, theta;
	phi = eSpec.phi(0) / 180 * pi;
	theta = eSpec.theta(0) / 180 * pi;

	if (screen == 0) {
		for (int i = 0; i < N; i++) {
			if (i < Nx) {
				z = eSpec.z(0) + rulerX[i] * std::cos(phi);
				x = -(rulerX[i] * std::sin(phi) - eSpec.x(0));

				rulerX[i] = 1000.0 * std::atan2(x,z);
			}
			if (i < Ny) {
				z = eSpec.z(0) + rulerY[i] * std::sin(theta);
				y =  -(rulerY[i] * std::cos(theta) - eSpec.y(0));

				rulerY[i] = 1000.0 * std::atan2(y,z);
			}
		}
	}
	else {
		for (int i = 0; i < Ny; i++) {
			z = eSpec.z(0) + rulerY[i] * std::sin(phi) + rulerX[0] * std::cos(theta);
			y = rulerY[i] * std::cos(phi) - eSpec.x(0);

			rulerY[i] = 1000.0 * std::atan2(y, z);
		}
	}
}

void drawAxis(bool mode, double scale, spectrometer& eSpec, paramSpace & pSpace, int& screen, std::vector<double>& rulerX, std::vector<double>& rulerY, double pointing) {
	std::vector<double> pixelX, pixelY;
	int Nx = (int)rulerX.size();
	int Ny = (int)rulerY.size();
	pixelX.resize(Nx);
	pixelY.resize(Ny);

	double lw = 2.0 * scale;
	double txtSize = 11 * scale;

	int countT = 0;
	std::generate(std::begin(pixelX), std::end(pixelX), [&] {
		return countT++;
		});

	countT = 0;
	std::generate(std::begin(pixelY), std::end(pixelY), [&] {
		return countT++;
		});

	if (mode == 0) {
		std::vector<double> xAxis, yAxis;
		double value, eval, xZero, yZero;
		int xTickStart, yTickStart;
		eval = 0.0;
		FE1DInterp(rulerX, pixelX, eval, xZero);
		FE1DInterp(rulerY, pixelY, eval, yZero);

		xAxis.push_back(xZero);
		bool loop = 1;
		while (loop) {
			eval = eval - 5.0;
			FE1DInterp(rulerX, pixelX, eval, value);
			value = round(value * 10.0) / 10.0;
			if (value > 0 && value < Nx) {
				xAxis.push_back(value);
			}
			else {
				loop = 0;
			}
		}
		eval = 0.0;
		loop = 1;
		while (loop) {
			eval = eval + 5.0;
			FE1DInterp(rulerX, pixelX, eval, value);
			value = round(value * 10.0) / 10.0;
			if (value > 0 && value < Nx) {
				xAxis.insert(xAxis.begin(), value);
			}
			else {
				loop = 0;
			}
		}
		xTickStart = (int)(eval - 5.0);
		Nx = (int)xAxis.size();

		eval = 0.0;
		yAxis.push_back(yZero);
		loop = 1;
		while (loop) {
			eval = eval - 5.0;
			FE1DInterp(rulerY, pixelY, eval, value);
			value = round(value * 10.0) / 10.0;
			if (value > 0 && value < Ny) {
				yAxis.push_back(value);
			}
			else {
				loop = 0;
			}
		}
		eval = 0.0;
		loop = 1;
		while (loop) {
			eval = eval + 5.0;
			FE1DInterp(rulerY, pixelY, eval, value);
			value = round(value * 10.0) / 10.0;
			if (value > 0 && value < Ny) {
				yAxis.insert(yAxis.begin(), value);
			}
			else {
				loop = 0;
			}
		}
		yTickStart = (int)(eval - 5.0);
		Ny = (int)yAxis.size();

		std::vector<double> plotX, plotY;
		plotX.resize(2, 0.0);
		plotY.resize(2, 0.0);
		if (screen == 0) {
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
			for (int i = 0; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plotY[0] = yZero - 10;
				plotY[1] = yZero + 10;
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				if (xTickStart - 5 * i != 0) {
					//plt::text(xAxis[i] - 10, yZero + 25, std::to_string(xTickStart - 5 * i));
					if ((xTickStart - 5 * i) >= 0) {
						if ((xTickStart - 5 * i) > 0) {
							if ((xTickStart - 5 * i) < 10) {
								plt::text(xAxis[i] - 30, yZero + 50, " 0" + std::to_string(xTickStart - 5 * i));
							}
							else {
								plt::text(xAxis[i] - 30, yZero + 50, " " + std::to_string(xTickStart - 5 * i));
							}
						}
						else {
							plt::text(xAxis[i] - 30, yZero + 50, " 00");
						}

					}
					else {
						if ((xTickStart - 5 * i) <= -10) {
							plt::text(xAxis[i] - 30, yZero + 50, std::to_string(xTickStart - 5 * i));
						}
						else {
							plt::text(xAxis[i] - 30, yZero + 50, "-0" + std::to_string(abs(xTickStart - 5 * i)));
						}
					}
				}
			}
			for (int i = 0; i < Ny; i++) {
				plotX[0] = xZero - 10;
				plotX[1] = xZero + 10;
				plotY[0] = yAxis[i];
				plotY[1] = yAxis[i];
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				if (yTickStart - 5 * i != 0) {
					//plt::text(xZero + 25, yAxis[i] + 5, std::to_string(yTickStart - 5 * i));
					if ((yTickStart - 5 * i) >= 0) {
						if ((yTickStart - 5 * i) > 0) {
							if ((yTickStart - 5 * i) < 10) {
								plt::text(xZero + 25, yAxis[i] + 10, " 0" + std::to_string(yTickStart - 5 * i));
							}
							else {
								plt::text(xZero + 25, yAxis[i] + 10, " " + std::to_string(yTickStart - 5 * i));
							}
						}
						else {
							plt::text(xZero + 25, yAxis[i] + 10, " 00");
						}

					}
					else {
						if ((yTickStart - 5 * i) <= -10) {
							plt::text(xZero + 25, yAxis[i] + 10, std::to_string(yTickStart - 5 * i));
						}
						else {
							plt::text(xZero + 25, yAxis[i] + 10, "-0" + std::to_string(abs(yTickStart - 5 * i)));
						}
					}
				}
			}
			plotX[0] = 0;
			plotX[1] = (int)rulerX.size() - 1;
			plotY[0] = yZero;
			plotY[1] = yZero;
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
			plotX[0] = xZero;
			plotX[1] = xZero;
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
		}
		else {
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
			double locationX, locationY;
			for (int i = 1; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plotY[0] = yAxis[0];
				plotY[1] = 0;
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				locationY = 35;
				if ((xTickStart - 5 * i) % 10 == 0) {
					if ((xTickStart - 5 * i) == 0) {
						plt::text(xAxis[i] - 15, locationY, "00");
					}
					else {
						plt::text(xAxis[i] - 15, locationY, std::to_string(xTickStart - 5 * i));
					}
				}
			}
			plt::rcparams({ {"text.color", "k"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
			for (int i = 0; i < Ny; i++) {
				plotX[0] = xAxis[0];
				plotX[1] = 0;
				plotY[0] = yAxis[i];
				plotY[1] = yAxis[i];
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				locationX = -35;
				if ((yTickStart - 5 * i) >= 0) {
					if ((yTickStart - 5 * i) > 0) {
						if ((yTickStart - 5 * i) < 10) {
							plt::text(locationX, yAxis[i] + 5, " 0" + std::to_string(yTickStart - 5 * i));
						}
						else {
							plt::text(locationX, yAxis[i] + 5, " " + std::to_string(yTickStart - 5 * i));
						}
					}
					else {
						plt::text(locationX, yAxis[i] + 5, " 00");
					}
					
				}
				else {
					if ((yTickStart - 5 * i) <= -10) {
						plt::text(locationX, yAxis[i] + 5, std::to_string(yTickStart - 5 * i));
					}
					else {
						plt::text(locationX, yAxis[i] + 5, "-0" + std::to_string(abs(yTickStart - 5 * i)));
					}
				}
				
			}
			plotX[0] = 0;
			plotX[1] = (int)rulerX.size() - 1;
			plotY[0] = yAxis[0];
			plotY[1] = yAxis[0];
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
			plotX[0] = xAxis[0];
			plotX[1] = xAxis[0];
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
		}
	}
	else {
		int warning = 0;
		std::vector<double> mRadX = rulerX;
		std::vector<double> mRadY = rulerY;
		mRadAxis(eSpec, screen, mRadX, mRadY);

		std::vector<double> xAxis, yAxis;
		std::vector<int> xTick;
		double value, eval, xZero, yZero, ptMax, ptMin;
		int xTickStart, yTickStart;
		eval = 0.0;
		FE1DInterp(mRadX, pixelX, eval, xZero);
		FE1DInterp(mRadY, pixelY, eval, yZero);

		ptMax = 0;
		ptMin = 0;
		xAxis.push_back(xZero);
		bool loop = 1;
		if (screen == 0) {
			while (loop) {
				eval = eval - 5.0;
				FE1DInterp(mRadX, pixelX, eval, value);
				if (value > 0 && value < Nx) {
					xAxis.push_back(value);
				}
				else {
					loop = 0;
				}
			}
			eval = 0.0;
			loop = 1;
			while (loop) {
				eval = eval + 5.0;
				FE1DInterp(mRadX, pixelX, eval, value);
				if (value > 0 && value < Nx) {
					xAxis.insert(xAxis.begin(), value);
				}
				else {
					loop = 0;
				}
			}
			xTickStart = (int)(eval - 5.0);
		}
		else {
			int NE = (int)pSpace.energy(screen).size();
			int indexStart, indexEnd;
			std::vector<double> screenPos;
			std::vector<double> EnAxis = pSpace.energy(screen);
			std::vector<double> PtAxis = pSpace.pointing(screen);
			std::vector<std::vector<double>> pS = pSpace.parameterSpace(screen);
			ptMax = std::max(PtAxis.front(), PtAxis.back());
			ptMin = std::min(PtAxis.front(), PtAxis.back());
			double dE = abs(EnAxis[1] - EnAxis[0]);
			screenPos.resize(NE, 0.0);
			if (pointing > ptMax) {
				for (int i = 0; i < NE; i++) {
					if (ptMax == PtAxis.front()) {
						screenPos[i] = pS[i][0];
					}
					else {
						screenPos[i] = pS[i].back();
					}
					warning = 1;
				}
			}
			else {
				if (pointing < ptMin) {
					for (int i = 0; i < NE; i++) {
						if (ptMin == PtAxis.front()) {
							screenPos[i] = pS[i][0];
						}
						else {
							screenPos[i] = pS[i].back();
						}
						warning = -1;
					}
				}
				else {
					for (int i = 0; i < NE; i++) {
						FE2DInterp(EnAxis, PtAxis, pS, EnAxis[i], pointing, screenPos[i]);
					}
				}
			}

			double Elow, Ehigh;
			double screenLow = std::min(rulerX.front(), rulerX.back());
			double screenHigh = std::max(rulerX.front(), rulerX.back());
			indexStart = 0;
			indexEnd = NE - 1;
			if (screenPos.front() >= screenLow) {
				Elow = EnAxis.front();
				indexStart = 0;
			}
			else {
				FE1DInterp(screenPos, EnAxis, screenLow, Elow);
				Elow = ceil(Elow / 10.0) * 10.0;
				for (int i = 0; i < NE; i++) {
					if (abs(Elow - EnAxis[i]) < dE) {
						indexStart = i;
						Elow = EnAxis[indexStart];
						break;
					}
				}
			}
			if (screenPos.back() > screenHigh) {
				FE1DInterp(screenPos, EnAxis, screenHigh, Ehigh);
				Ehigh = floor(Ehigh / 10.0) * 10.0;
				for (int i = 0; i < NE; i++) {
					if (abs(Ehigh - EnAxis[NE - 1 - i]) < dE) {
						indexEnd = NE - 1 - i;
						Ehigh = EnAxis[indexEnd];
						break;
					}
				}
			}
			else {
				Ehigh = EnAxis.back();
				indexEnd = NE - 1;
			}
			if (indexEnd - indexStart != NE - 1) {
				std::vector<double> bufferE(EnAxis.begin() + indexStart,EnAxis.begin() + indexEnd + 1);
				std::vector<double> bufferSP(screenPos.begin() + indexStart, screenPos.begin() + indexEnd + 1);
				EnAxis.clear();
				EnAxis = bufferE;
				screenPos.clear();
				NE = (int)EnAxis.size();
				screenPos.resize(NE, 0.0);
				for (int i = 0; i < NE; i++) {
					FE1DInterp(rulerX, pixelX, bufferSP[i], screenPos[i]);
				}
				bufferE.clear();
				bufferSP.clear();
			}
			else {
				std::vector<double> bufferSP = screenPos;
				for (int i = 0; i < NE; i++) {
					FE1DInterp(rulerX, pixelX, bufferSP[i], screenPos[i]);
				}
				bufferSP.clear();
			}
			eval = Ehigh;

			xAxis.clear();
			xTick.clear();
			FE1DInterp(EnAxis, screenPos, eval, value);
			value = round(value * 10.0) / 10.0;
			xTick.push_back((int)eval);
			xAxis.push_back(value);
			loop = 1;
			while (loop) {
				eval = eval - 10;
				if (eval >= Elow) {
					FE1DInterp(EnAxis, screenPos, eval, value);
					value = round(value * 10.0) / 10.0;
					if (value > 0 && value < Nx) {
						if (abs(value - xAxis.back()) >= 15.0) {
							xTick.push_back((int)eval);
							xAxis.push_back(value);
						}
					}
					else {
						loop = 0;
					}
				}
				else {
					loop = 0;
				}
			}
			xTickStart = Ehigh;

			screenPos.clear();
			EnAxis.clear();
			PtAxis.clear();
			pS.clear();
		}
		Nx = (int)xAxis.size();

		eval = 0.0;
		yAxis.push_back(yZero);
		loop = 1;
		while (loop) {
			eval = eval - 5.0;
			FE1DInterp(mRadY, pixelY, eval, value);
			value = round(value * 10.0) / 10.0;
			if (value > 0 && value < Ny) {
				yAxis.push_back(value);
			}
			else {
				loop = 0;
			}
		}
		eval = 0.0;
		loop = 1;
		while (loop) {
			eval = eval + 5.0;
			FE1DInterp(mRadY, pixelY, eval, value);
			value = round(value * 10.0) / 10.0;
			if (value > 0 && value < Ny) {
				yAxis.insert(yAxis.begin(), value);
			}
			else {
				loop = 0;
			}
		}
		yTickStart = (int)(eval - 5.0);
		Ny = (int)yAxis.size();

		std::vector<double> plotX, plotY;
		plotX.resize(2, 0.0);
		plotY.resize(2, 0.0);
		if (screen == 0) {
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
			for (int i = 0; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plotY[0] = yZero - 25;
				plotY[1] = yZero + 25;
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				if (xTickStart - 5 * i != 0) {
					if ((xTickStart - 5 * i) >= 0) {
						if ((xTickStart - 5 * i) > 0) {
							if ((xTickStart - 5 * i) < 10) {
								plt::text(xAxis[i] - 50, yZero + 75, " 0" + std::to_string(xTickStart - 5 * i));
							}
							else {
								plt::text(xAxis[i] - 50, yZero + 75, " " + std::to_string(xTickStart - 5 * i));
							}
						}
						else {
							plt::text(xAxis[i] - 50, yZero + 75, " 00");
						}

					}
					else {
						if ((xTickStart - 5 * i) <= -10) {
							plt::text(xAxis[i] - 50, yZero + 75, std::to_string(xTickStart - 5 * i));
						}
						else {
							plt::text(xAxis[i] - 50, yZero + 75, "-0" + std::to_string(abs(xTickStart - 5 * i)));
						}
					}
				}
			}
			for (int i = 0; i < Ny; i++) {
				plotX[0] = xZero - 25;
				plotX[1] = xZero + 25;
				plotY[0] = yAxis[i];
				plotY[1] = yAxis[i];
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				if (yTickStart - 5 * i != 0) {
					if ((yTickStart - 5 * i) >= 0) {
						if ((yTickStart - 5 * i) > 0) {
							if ((yTickStart - 5 * i) < 10) {
								plt::text(xZero + 40, yAxis[i] + 15, " 0" + std::to_string(yTickStart - 5 * i));
							}
							else {
								plt::text(xZero + 40, yAxis[i] + 15, " " + std::to_string(yTickStart - 5 * i));
							}
						}
						else {
							plt::text(xZero + 40, yAxis[i] + 15, " 00");
						}

					}
					else {
						if ((yTickStart - 5 * i) <= -10) {
							plt::text(xZero + 40, yAxis[i] + 15, std::to_string(yTickStart - 5 * i));
						}
						else {
							plt::text(xZero + 40, yAxis[i] + 15, "-0" + std::to_string(abs(yTickStart - 5 * i)));
						}
					}
				}
			}
			plotX[0] = 0;
			plotX[1] = (int)rulerX.size() - 1;
			plotY[0] = yZero;
			plotY[1] = yZero;
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
			plotX[0] = xZero;
			plotX[1] = xZero;
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });

		}
		else {
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
			plotX[0] = 0;
			plotX[1] = (int)rulerX.size() - 1;
			double locationX, locationY;
			if (screen == 1) {
				plotY[0] = yAxis[0];
				plotY[1] = yAxis[0];
				locationY = 25;
			}
			else {
				plotY[0] = yAxis[0];
				plotY[1] = yAxis[0];
				locationY = 20;
			}
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
			if (screen == 1) {
				plotY[1] = yAxis[0] - 20;
			}
			else {
				plotY[1] = yAxis[0] - 20;
			}
			int lastLabel = 1;
			for (int i = 1; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				if (xTick[i] % 20 == 0) {
					if (i > 1){
						if (abs(xAxis[i] - xAxis[lastLabel]) > 60) {
							plt::text(xAxis[i] - 20, locationY, std::to_string(xTick[i]));
							lastLabel = i;
						}
					}
					else {
						plt::text(xAxis[i] - 20, locationY, std::to_string(xTick[i]));
					}
				}
				else {
					if (i == Nx - 1) {
						plt::text(xAxis[i] - 20, locationY, std::to_string(xTick[i]));
					}
				}
			}
			if (warning == 1) {
				std::string pValue = std::to_string(ptMax);
				pValue = pValue.substr(0, 4);
				plt::text((int)(0.775 * rulerX.size()), (int)(0.9 * rulerY.size()), std::string("Out of Range! Axis for ") + pValue + std::string(" mrad"));
			}
			else {
				if (warning == -1) {
					std::string pValue = std::to_string(ptMin);
					pValue = pValue.substr(0, 5);
					plt::text((int)(0.775 * rulerX.size()), (int)(0.9 * rulerY.size()), std::string("Out of Range! Axis for ") + pValue + std::string(" mrad"));
				}
			}
			plt::rcparams({ {"text.color", "k"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
			double dx = (int)rulerX.size();
			for (int i = 0; i < Nx; i++) {
				value = abs(xAxis[i] - 25);
				if (value < dx) {
					dx = value;
					plotX[0] = xAxis[i];
					plotX[1] = xAxis[i];
				}
			}
			if (dx > 60) {
				plotX[0] = 25;
				plotX[1] = 25;
			}
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
			plotX[1] = 0;
			for (int i = 0; i < Ny; i++) {
				plotY[0] = yAxis[i];
				plotY[1] = yAxis[i];
				if (yTickStart - 5 * i > -15 && yTickStart - 5 * i < 15) {
					plt::plot(plotX, plotY, { {"color","w"}, {"linewidth", std::to_string(lw)} });
				}
				locationX = -45;
				if ((yTickStart - 5 * i) >= 0) {
					if ((yTickStart - 5 * i) > 0) {
						if ((yTickStart - 5 * i) < 15) {
							if ((yTickStart - 5 * i) < 10) {
								plt::text(locationX, yAxis[i] + 5, " 0" + std::to_string(yTickStart - 5 * i));
							}
							else {
								plt::text(locationX, yAxis[i] + 5, " " + std::to_string(yTickStart - 5 * i));
							}
						}
					}
					else {
						plt::text(locationX, yAxis[i] + 5, " 00");
					}

				}
				else {
					if ((yTickStart - 5 * i) <= -10) {
						if ((yTickStart - 5 * i) > -15){
							plt::text(locationX, yAxis[i] + 5, std::to_string(yTickStart - 5 * i));
						}
					}
					else {
						plt::text(locationX, yAxis[i] + 5, "-0" + std::to_string(abs(yTickStart - 5 * i)));
					}
				}

			}
		}
		plotX.clear();
		plotY.clear();
		mRadX.clear();
		mRadY.clear();

		xAxis.clear();
		yAxis.clear();
		xTick.clear();
	}
	pixelX.clear();
	pixelY.clear();
}

void calMode(spectrometer& eSpec, std::string& pathCalibration, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler) {
	std::string listPath;
	std::vector<std::string> listFile;
	//listPath = "C:\\Users\\Xing\\Downloads\\Wakefield\\20230811\\pPointing";
	//std::vector<std::string> listPPointing = fileList(listPath);
	bool loadError = 0;

	std::vector<std::vector<int>> winRes;
	std::vector<std::vector<double>> zP;
	std::vector<int> buffer;
	buffer.resize(2, 0);
	H.clear();
	H.resize(3);
	xRuler.clear();
	xRuler.resize(3);
	yRuler.clear();
	yRuler.resize(3);
	winRes.clear();
	winRes.resize(3);
	zP.clear();
	zP.resize(3);

	imageBW imTeP, imTeSA, imTeSB;
	int mode = 1;
	screenCalibration calibration;
	int screen = 0;
	bool loop = 1;
	std::cout << "\nCalibration For Pointing Screen\n";
	while (loop) {
		listPath.clear();
		listFile.clear();
		listPath = eSpec.screenPath(screen);
		listDir(listPath, listFile);
		std::string refFile;
		findRef(listFile, refFile);
		if (refFile.find("-1NoRef") == std::string::npos) {
			calibration.loadCalibration(pathCalibration);
			screenCal(mode, screen, eSpec, calibration, H[0], xRuler[0], yRuler[0]);
			std::vector<double> viewRes = calibration.viewResolution(screen);
			imageBW image;
			getImage(refFile, image);
			perspectiveTransform(image, H[0], viewRes, imTeP);
			findZero(screen, imTeP, xRuler[0], yRuler[0], zP[0]);

			buffer[0] = imTeP.sizeX();
			buffer[1] = imTeP.sizeY();
			winRes[0] = buffer;
			loadError = 0;
			image.destroy();
			viewRes.clear();
		}
		else {
			std::cout << "\nNo Reference Found at: " << eSpec.screenPath(screen) << "\n";
			loadError = 1;
		}

		std::cout << "\nDo you want to rerun the calibration (Y / N) ?\n";
		std::string userInput;
		/*
		double ratio = (double)imTransform.sizeX() / (double)imTransform.sizeY();
		int resV = 720;
		int resH = (int)(ratio * resV);
		plt::figure_size(resH, resV);
		pltimshow(imTransform, 0, "");
		drawAxis(0, eSpec, screen, resXeP, resYeP);
		plt::title("ePointing Screen");
		plt::axis("off");
		plt::show();
		*/
		std::cin >> userInput;
		if (userInput.at(0) == 'y' || userInput.at(0) == 'Y' || userInput.at(0) == '1') {
			loop = 1;
		}
		else {
			loop = 0;
		}
	}
	screen = 1;
	loop = 1;
	if (loadError == 0) {
		std::cout << "\nCalibration For e Energy Screen A\n";
		while (loop) {
			listPath.clear();
			listFile.clear();
			listPath = eSpec.screenPath(screen);
			listDir(listPath, listFile);
			std::string refFile;
			findRef(listFile, refFile);
			if (refFile.find("-1NoRef") == std::string::npos) {
				calibration.loadCalibration(pathCalibration);
				screenCal(mode, screen, eSpec, calibration, H[1], xRuler[1], yRuler[1]);
				std::vector<double> viewRes = calibration.viewResolution(screen);
				imageBW image;
				getImage(refFile, image);
				perspectiveTransform(image, H[1], viewRes, imTeSA);
				findZero(screen, imTeSA, xRuler[1], yRuler[1], zP[1]);

				buffer[0] = imTeSA.sizeX();
				buffer[1] = imTeSA.sizeY();
				winRes[1] = buffer;
				loadError = 0;
				image.destroy();
				viewRes.clear();
			}
			else {
				std::cout << "\nNo Reference Found at: " << eSpec.screenPath(screen) << "\n";
				loadError = 1;
			}

			std::cout << "\nDo you want to rerun the calibration (Y / N) ?\n";
			std::string userInput;
			/*
			double ratio = (double)imTransform.sizeX() / (double)imTransform.sizeY();
			int resV = 720;
			int resH = (int)(ratio * resV);
			plt::figure_size(resH, resV);
			pltimshow(imTransform, 0, "");
			drawAxis(0, eSpec, screen, resXeSA, resYeSA);
			plt::title("eSpec Sceen A");
			plt::axis("off");
			plt::show();
			*/
			std::cin >> userInput;
			if (userInput.at(0) == 'y' || userInput.at(0) == 'Y' || userInput.at(0) == '1') {
				loop = 1;
			}
			else {
				loop = 0;
			}
		}
	}
	screen = 2;
	loop = 1;
	if (loadError == 0) {
		std::cout << "\nCalibration For e Energy Screen B\n";
		while (loop) {
			listPath.clear();
			listFile.clear();
			listPath = eSpec.screenPath(screen);
			listDir(listPath, listFile);
			std::string refFile;
			findRef(listFile, refFile);
			if (refFile.find("-1NoRef") == std::string::npos) {
				calibration.loadCalibration(pathCalibration);
				screenCal(mode, screen, eSpec, calibration, H[2], xRuler[2], yRuler[2]);
				std::vector<double> viewRes = calibration.viewResolution(screen);
				imageBW image;
				getImage(refFile, image);
				perspectiveTransform(image, H[2], viewRes, imTeSB);
				findZero(screen, imTeSB, xRuler[2], yRuler[2], zP[2]);

				buffer[0] = imTeSB.sizeX();
				buffer[1] = imTeSB.sizeY();
				winRes[2] = buffer;
				loadError = 0;
				image.destroy();
				viewRes.clear();
			}
			else {
				std::cout << "\nNo Reference Found at: " << eSpec.screenPath(screen) << "\n";
				loadError = 1;
			}

			std::cout << "\nDo you want to rerun the calibration (Y / N) ?\n";
			std::string userInput;
			/*
			double ratio = (double)imTransform.sizeX() / (double)imTransform.sizeY();
			int resV = 720;
			int resH = (int)(ratio * resV);
			plt::figure_size(resH, resV);
			pltimshow(imTransform, 0, "");
			drawAxis(0, eSpec, screen, resXeSB, resYeSB);
			plt::title("eSpec Sceen B");
			plt::axis("off");
			plt::show();
			*/
			std::cin >> userInput;
			if (userInput.at(0) == 'y' || userInput.at(0) == 'Y' || userInput.at(0) == '1') {
				loop = 1;
			}
			else {
				loop = 0;
			}
		}

		writeCalibration(pathCalibration, winRes, H, xRuler, yRuler, zP);

		screen = 0;
		pixelAxis(screen, winRes[0][0], winRes[0][1], xRuler[0], yRuler[0], zP[0]);
		screen = 1;
		pixelAxis(screen, winRes[1][0], winRes[1][1], xRuler[1], yRuler[1], zP[1]);
		screen = 2;
		pixelAxis(screen, winRes[2][0], winRes[2][1], xRuler[2], yRuler[2], zP[2]);
	}
	else {
		std::cout << "\nCalibration Failed Due To Missing Files.\n";
	}
	imTeP.destroy();
	imTeSA.destroy();
	imTeSB.destroy();
	listFile.clear();
	winRes.clear();
	zP.clear();
	buffer.clear();
}

bool loadFile(std::vector<std::string>& list, std::string& path, std::string& fileName, std::string& timeStamp, cv::Mat& H, std::vector<double>& viewRes, imageBW& output) {
	int fileLoopCount = 0;
	bool fileFound = 0;
	bool fileFindLoop = 1;
	int index = -1;
	double loopWait = 250;

	imageBW imBuffer;
	while (fileFindLoop) {
		listDir(path, list);
		fileFound = findFile(list, fileName, timeStamp, index);
		if (fileFound) {
			fileFindLoop = 0;
		}
		else {
			std::this_thread::sleep_for(std::chrono::milliseconds((long)loopWait));
		}
		if (fileLoopCount > 30) {
			fileFindLoop = 0;
			fileFound = 0;
			std::cout << "Can not find file.\n";
		}
		fileLoopCount = fileLoopCount + 1;
	}
	if (fileFound) {
		printf("Loading Image.\n");
		getImage(list[index], imBuffer);
		perspectiveTransform(imBuffer, H, viewRes, output);
		imBuffer.destroy();
		return 1;
	}
	else {
		imBuffer.destroy();
		return 0;
	}
}

void drawPointingAnalysis(spectrometer& eSpec, paramSpace& pSpace, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler, std::vector<double>& pxX, std::vector<double>& pxY, imageBW& imP, imageBW& imA, imageBW& imB, std::string outputName) {
	int screenA, screenB, screenP;
	screenA = 1;
	screenB = 2;
	screenP = 0;

	std::vector<double> sAEline, sBEline, APix, BPix;

	imageBW imASmooth, imBSmooth;

	imA.copy(imASmooth);
	imB.copy(imBSmooth);

	removeOutlier(imASmooth, 4.0);
	removeOutlier(imBSmooth, 4.0);
	medianFilter(imASmooth, 12.0);
	medianFilter(imBSmooth, 12.0);

	APix.resize(imA.sizeX(), 0.0);
	sAEline.resize(imA.sizeX(), 0.0);

	BPix.resize(imB.sizeX(), 0.0);
	sBEline.resize(imB.sizeX(), 0.0);

	#pragma omp parallel for
	for (int i = 0; i < imA.sizeX(); i++) {
		double buffer = 0.0;
		for (int j = 0; j < imA.sizeY(); j++) {
			buffer = buffer + imASmooth.value(i, j) / imA.sizeY();
		}
		sAEline[i] = buffer;
		APix[i] = (double)i;
	}
	#pragma omp parallel for
	for (int i = 0; i < imB.sizeX(); i++) {
		double buffer = 0.0;
		for (int j = 0; j < imB.sizeY(); j++) {
			buffer = buffer + imBSmooth.value(i, j) / imB.sizeY();
		}
		sBEline[i] = buffer;
		BPix[i] = (double)i;
	}
	//medianFilter(sAEline, 12.0);
	//medianFilter(sBEline, 12.0);

	double minBuffer = sAEline[0];
	double maxBuffer = sAEline[0];
	for (int i = 1; i < imA.sizeX(); i++) {
		if (sAEline[i] > maxBuffer) {
			maxBuffer = sAEline[i];
		}
		else {
			if (sAEline[i] < minBuffer) {
				minBuffer = sAEline[i];
			}
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < imA.sizeX(); i++) {
		sAEline[i] = (sAEline[i] - minBuffer) / (maxBuffer - minBuffer);
	}

	minBuffer = sBEline[0];
	maxBuffer = sBEline[0];
	for (int i = 1; i < imB.sizeX(); i++) {
		if (sBEline[i] > maxBuffer) {
			maxBuffer = sBEline[i];
		}
		else {
			if (sBEline[i] < minBuffer) {
				minBuffer = sBEline[i];
			}
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < imB.sizeX(); i++) {
		sBEline[i] = (sBEline[i] - minBuffer) / (maxBuffer - minBuffer);
	}

	int Apeak = 0;
	double maxSignal = 0.0;
	for (int i = 0; i < imA.sizeX(); i++) {
		if (sAEline[i] > maxSignal) {
			Apeak = i;
			maxSignal = sAEline[i];
		}
	}

	maxSignal = maxSignal / 2.0;

	int Alow, Ahigh;
	Alow = 0;
	for (int i = 0; i < Apeak; i++) {
		if (sAEline[i] >= maxSignal) {
			Alow = i;
			break;
		}
	}
	Ahigh = Apeak;
	for (int i = Apeak; i < imA.sizeX(); i++) {
		if (sAEline[i] <= maxSignal) {
			Ahigh = i;
			break;
		}
	}

	int Bpeak = 0;
	maxSignal = 0.0;
	for (int i = 0; i < imB.sizeX(); i++) {
		if (sBEline[i] > maxSignal) {
			Bpeak = i;
			maxSignal = sBEline[i];
		}
	}

	maxSignal = maxSignal / 2.0;

	int Blow, Bhigh;
	Blow = 0;
	for (int i = 0; i < Bpeak; i++) {
		if (sBEline[i] >= maxSignal) {
			Blow = i;
			break;
		}
	}
	Bhigh = Bpeak;
	for (int i = Bpeak; i < imB.sizeX(); i++) {
		if (sBEline[i] <= maxSignal) {
			Bhigh = i;
			break;
		}
	}


	std::vector<int> peak, peakBound;
	std::vector<double> pointX, pointY;
	pointX = xRuler[0];
	pointY = yRuler[0];
	mRadAxis(eSpec, screenP, pointX, pointY);
	double pointing, eval, dbuffer;
	std::vector<int> acceptanceBound;
	acceptanceBound.resize(4, 0);
	eval = -1.0 * eSpec.angleMax(0);
	FE1DInterp(pointX, pxX, eval, dbuffer);
	acceptanceBound[0] = (int)round(dbuffer);
	eval = eSpec.angleMax(0);
	FE1DInterp(pointX, pxX, eval, dbuffer);
	acceptanceBound[1] = (int)round(dbuffer);
	eval = -1.0 * eSpec.angleMax(1);
	FE1DInterp(pointY, pxY, eval, dbuffer);
	acceptanceBound[2] = (int)round(dbuffer);
	eval = eSpec.angleMax(1);
	FE1DInterp(pointY, pxY, eval, dbuffer);
	acceptanceBound[3] = (int)round(dbuffer);

	std::vector<double> boundboxX, boundboxY;
	boundboxX.resize(5, 0.0);
	boundboxY.resize(5, 0.0);
	boundboxX[0] = acceptanceBound[0];
	boundboxY[0] = acceptanceBound[2];
	boundboxX[1] = acceptanceBound[0];
	boundboxY[1] = acceptanceBound[3];
	boundboxX[2] = acceptanceBound[1];
	boundboxY[2] = acceptanceBound[3];
	boundboxX[3] = acceptanceBound[1];
	boundboxY[3] = acceptanceBound[2];
	boundboxX[4] = acceptanceBound[0];
	boundboxY[4] = acceptanceBound[2];

	
	printf("Loaded 3 Images.\n");
	imageBW imBufferS, imBufferL;
	imP.copy(imBufferL);
	int maxValue = 0;
	bool flagP = 0;
	bool flagB = 0;
	double peakValue, peakValueB, totalValue, acceptValue;
	removeOutlier(imBufferL, 4.0);
	medianFilter(imBufferL, 2);
	imBufferL.crop(acceptanceBound, imBufferS);
	findSignalPeak(imBufferS, peakBound, peakValueB);
	findSignalPeak(imBufferL, peak, peakValue);
	peakBound[0] = peakBound[0] + acceptanceBound[0];
	peakBound[1] = peakBound[1] + acceptanceBound[2];
	//eval = (double)((double)imBufferL.sizeY() - 1 - peakBound[1]);
	eval = (double)(peakBound[1]);
	FE1DInterp(pxY, pointY, eval, pointing);
	maxValue = (int)round(peakValue * 10000);
	Sum(imP, totalValue);
	Sum(imBufferS, acceptValue);
	totalValue = round(totalValue / ((double)(imP.sizeX() * imP.sizeY())) * 10000);
	acceptValue = round(acceptValue / ((double)(imBufferS.sizeX() * imBufferS.sizeY())) * 10000);
	imBufferS.destroy();
	imBufferL.destroy();
	printf("Found Electron Pointing.\n");

	std::vector<double> screenPos;
	std::vector<double> EnAxis = pSpace.energy(1);
	std::vector<double> PtAxis = pSpace.pointing(1);
	std::vector<std::vector<double>> pS = pSpace.parameterSpace(1);
	double ptMax = std::max(PtAxis.front(), PtAxis.back());
	double ptMin = std::min(PtAxis.front(), PtAxis.back());
	int NE = EnAxis.size();
	screenPos.resize(NE, 0.0);
	if (pointing > ptMax) {
		for (int i = 0; i < NE; i++) {
			if (ptMax == PtAxis.front()) {
				screenPos[i] = pS[i][0];
			}
			else {
				screenPos[i] = pS[i].back();
			}
		}
	}
	else {
		if (pointing < ptMin) {
			for (int i = 0; i < NE; i++) {
				if (ptMin == PtAxis.front()) {
					screenPos[i] = pS[i][0];
				}
				else {
					screenPos[i] = pS[i].back();
				}
			}
		}
		else {
			for (int i = 0; i < NE; i++) {
				FE2DInterp(EnAxis, PtAxis, pS, EnAxis[i], pointing, screenPos[i]);
			}
		}
	}

	double ACenterEn, ASpreadEn, BCenterEn, BSpreadEn, inputBuffer, outputBuffer;
	inputBuffer = (double)Apeak;
	FE1DInterp(APix, xRuler[1], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, ACenterEn);

	inputBuffer = (double)Alow;
	FE1DInterp(APix, xRuler[1], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, ASpreadEn);

	inputBuffer = (double)Ahigh;
	FE1DInterp(APix, xRuler[1], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, outputBuffer);
	ASpreadEn = ASpreadEn - outputBuffer;


	EnAxis = pSpace.energy(2);
	PtAxis = pSpace.pointing(2);
	pS = pSpace.parameterSpace(2);
	ptMax = std::max(PtAxis.front(), PtAxis.back());
	ptMin = std::min(PtAxis.front(), PtAxis.back());
	NE = EnAxis.size();
	screenPos.resize(NE, 0.0);
	if (pointing > ptMax) {
		for (int i = 0; i < NE; i++) {
			if (ptMax == PtAxis.front()) {
				screenPos[i] = pS[i][0];
			}
			else {
				screenPos[i] = pS[i].back();
			}
		}
	}
	else {
		if (pointing < ptMin) {
			for (int i = 0; i < NE; i++) {
				if (ptMin == PtAxis.front()) {
					screenPos[i] = pS[i][0];
				}
				else {
					screenPos[i] = pS[i].back();
				}
			}
		}
		else {
			for (int i = 0; i < NE; i++) {
				FE2DInterp(EnAxis, PtAxis, pS, EnAxis[i], pointing, screenPos[i]);
			}
		}
	}

	inputBuffer = (double)Bpeak;
	FE1DInterp(BPix, xRuler[2], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, BCenterEn);

	inputBuffer = (double)Blow;
	FE1DInterp(BPix, xRuler[2], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, BSpreadEn);

	inputBuffer = (double)Bhigh;
	FE1DInterp(BPix, xRuler[2], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, outputBuffer);
	BSpreadEn = BSpreadEn - outputBuffer;


	#pragma omp parallel for
	for (int i = 0; i < imA.sizeX(); i++) {
		sAEline[i] = imA.sizeY() * (1.0 - sAEline[i] / 3.0);
	}

	#pragma omp parallel for
	for (int i = 0; i < imB.sizeX(); i++) {
		sBEline[i] = imB.sizeY() * (1.0 - sBEline[i] / 3.0);
	}

	double scale = 2;
	double lw = 2.0 * scale;
	double txtSize = 11 * scale;

	
	size_t resV, resH;
	double ratio;
	resH = 2224 * scale;
	double spX, spY;
	spY = ((double)imA.sizeY() + (double)imB.sizeY());
	spX = (std::max((double)imP.sizeX(), (double)imB.sizeX()));
	ratio = spY / (spX + spY);
	ratio = ratio * resH - 7;
	resV = (size_t)round(ratio);

	std::vector<double> drawLineX, drawLineY;
	drawLineX.resize(2, 0.0);
	drawLineY.resize(2, 0.0);

	printf("Drawing Image.\n");
	plt::figure_size(resH, resV);
	

	printf("Drawing Pointing Image.\n");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 0, 2, 1);
	pltimshow(imP, 0, "");

	plt::plot(boundboxX, boundboxY, { {"color","b"}, {"linewidth", std::to_string(lw)} });

	/*
	drawLineX[0] = 0;
	drawLineX[1] = (int)imP.sizeX() - 1;
	drawLineY[0] = (int)imP.sizeY() - 1 - peak[1];
	//drawLineY[0] = peak[1];
	drawLineY[1] = drawLineY[0];
	plt::plot(drawLineX, drawLineY, { {"color","k"} });

	drawLineX[0] = peak[0];
	drawLineX[1] = drawLineX[0];
	drawLineY[0] = 0;
	drawLineY[1] = (int)imP.sizeY() - 1;
	plt::plot(drawLineX, drawLineY, { {"color","k"} });
	*/

	
	drawLineX[0] = 0;
	drawLineX[1] = (int)imP.sizeX() - 1;
	//drawLineY[0] = (int)imP.sizeY() - 1 - peakBound[1];
	drawLineY[0] = peakBound[1];
	drawLineY[1] = drawLineY[0];
	plt::plot(drawLineX, drawLineY, { {"color","r"}, {"linewidth", std::to_string(lw)} });

	drawLineX[0] = peakBound[0];
	drawLineX[1] = drawLineX[0];
	drawLineY[0] = 0;
	drawLineY[1] = (int)imP.sizeY() - 1;
	plt::plot(drawLineX, drawLineY, { {"color","r"}, {"linewidth", std::to_string(lw)} });
	

	plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
	std::string pValue = std::to_string(pointing);
	std::string mValue = std::to_string(maxValue);
	std::string tValue = std::to_string((int)totalValue);
	std::string aValue = std::to_string((int)acceptValue);
	pValue = pValue.substr(0, 4);
	plt::text((int)(0.775 * imP.sizeX()), (int)(0.975 * imP.sizeY()), pValue + std::string(" mrad"));
	plt::text((int)(0.025 * imP.sizeX()), (int)(0.975 * imP.sizeY()), std::string("Max Px: ") + mValue + std::string("/10000"));
	plt::text((int)(0.025 * imP.sizeX()), (int)(0.100 * imP.sizeY()), std::string("AVG ToT Signal: ") + tValue);
	plt::text((int)(0.025 * imP.sizeX()), (int)(0.050 * imP.sizeY()), std::string("AVG Red Signal: ") + aValue);
	if (flagP) {
		if (flagB) {
			plt::text((int)(0.800 * imP.sizeX()), (int)(0.925 * imP.sizeY()), std::string("Error P,B"));
		}
		else {
			plt::text((int)(0.800 * imP.sizeX()), (int)(0.925 * imP.sizeY()), std::string("Error P"));
		}
	}
	else {
		if (flagB) {
			plt::text((int)(0.800 * imP.sizeX()), (int)(0.925 * imP.sizeY()), std::string("Error B"));
		}
	}
	printf("Drawing Pointing Axis: %f.\n", pointing);
	drawAxis(1, scale, eSpec, pSpace, screenP, xRuler[screenP], yRuler[screenP], pointing);

	printf("Drawing Energy A Image.\n");
	plt::axis("off");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 1, 1, (int)(spX / spY));
	//removeOutlier(imA, 4.0);
	//medianFilter(imA, 2);
	pltimshow(imASmooth, 0, "");
	plt::plot(APix, sAEline, { {"color","w"}, {"linewidth", std::to_string(lw)} });

	printf("Drawing Energy A Axis.\n");
	drawAxis(1, scale, eSpec, pSpace, screenA, xRuler[screenA], yRuler[screenA], pointing);
	plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
	std::string AEValue = std::to_string((int)std::round(ACenterEn));
	plt::text((int)(0.825 * imA.sizeX()), (int)(0.250 * imA.sizeY()), std::string("Centroid Energy: ") + AEValue + std::string(" MeV"));
	std::string AESValue = std::to_string((int)std::round(ASpreadEn));
	plt::text((int)(0.825 * imA.sizeX()), (int)(0.325 * imA.sizeY()), std::string("Energy Spread: ") + AESValue + std::string(" MeV"));
	plt::axis("off");

	printf("Drawing Energy B Image.\n");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 1, 1, 1, (int)(spX / spY));
	//removeOutlier(imB, 4.0);
	//medianFilter(imB, 2);
	pltimshow(imBSmooth, 0, "");
	plt::plot(BPix, sBEline, { {"color","w"}, {"linewidth", std::to_string(lw)} });

	printf("Drawing Energy B Axis.\n");
	drawAxis(1, scale, eSpec, pSpace, screenB, xRuler[screenB], yRuler[screenB], pointing);
	plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
	std::string BEValue = std::to_string((int)std::round(BCenterEn));
	plt::text((int)(0.825 * imB.sizeX()), (int)(0.250 * imB.sizeY()), std::string("Centroid Energy: ") + BEValue + std::string(" MeV"));
	std::string BESValue = std::to_string((int)std::round(BSpreadEn));
	plt::text((int)(0.825 * imB.sizeX()), (int)(0.325 * imB.sizeY()), std::string("Energy Spread: ") + BESValue + std::string(" MeV"));
	plt::axis("off");

	printf("Drawing Canvas.\n");
	plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.075}, {"hspace",0.0} });
	plt::draw();

	printf("Saving Analysis.\n");
	std::string outputHR = eSpec.analysisPath() + "/" + outputName + "-HR.png";
	std::string outputLR = eSpec.analysisPath() + "/" + outputName + ".png";
	plt::save(outputHR);
	plt::close();
	//double scaling = 0.5;
	//resizeImage(scaling, outputHR, outputLR);
	printf("Analysis Saved.\n");

	sAEline.clear();
	sBEline.clear();
	APix.clear();
	BPix.clear();
	peak.clear();
	peakBound.clear();
	pointX.clear();
	pointY.clear();
	acceptanceBound.clear();
	screenPos.clear();
	EnAxis.clear();
	PtAxis.clear();
	pS.clear();
	drawLineX.clear();
	drawLineY.clear();
	imASmooth.destroy();
	imBSmooth.destroy();
}

void drawPointingAnalysisManual(spectrometer& eSpec, paramSpace& pSpace, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler, std::vector<double>& pxX, std::vector<double>& pxY, imageBW& imP, imageBW& imA, imageBW& imB, std::string outputName) {
	int screenA, screenB, screenP;
	screenA = 1;
	screenB = 2;
	screenP = 0;

	imageBW imASmooth, imBSmooth;

	imA.copy(imASmooth);
	imB.copy(imBSmooth);

	removeOutlier(imASmooth, 4.0);
	removeOutlier(imBSmooth, 4.0);
	medianFilter(imASmooth, 12.0);
	medianFilter(imBSmooth, 12.0);

	std::vector<double> sAEline, sBEline, APix, BPix, AEn, BEn;
	std::vector<double> PXPix, PYPix, PXSum, PYSum;
	PXSum.resize(imP.sizeX(), 0.0);
	PXPix.resize(imP.sizeX(), 0.0);

	PYSum.resize(imP.sizeY(), 0.0);
	PYPix.resize(imP.sizeY(), 0.0);

	APix.resize(imA.sizeX(), 0.0);
	AEn.resize(imA.sizeX(), 0.0);
	sAEline.resize(imA.sizeX(), 0.0);

	BPix.resize(imB.sizeX(), 0.0);
	BEn.resize(imB.sizeX(), 0.0);
	sBEline.resize(imB.sizeX(), 0.0);
	#pragma omp parallel for
	for (int i = 0; i < imP.sizeX(); i++) {
		double buffer = 0.0;
		for (int j = 0; j < imP.sizeY(); j++) {
			buffer = buffer + imP.value(i, j)/imP.sizeY();
		}
		PXSum[i] = buffer;
		PXPix[i] = (double)i;
	}
	#pragma omp parallel for
	for (int i = 0; i < imP.sizeY(); i++) {
		double buffer = 0.0;
		for (int j = 0; j < imP.sizeX(); j++) {
			buffer = buffer + imP.value(j, imP.sizeY() - i) / imP.sizeX();
		}
		PYSum[i] = buffer;
		PYPix[i] = (double)i;
	}
	#pragma omp parallel for
	for (int i = 0; i < imA.sizeX(); i++) {
		double buffer = 0.0;
		for (int j = 0; j < imA.sizeY(); j++) {
			buffer = buffer + imASmooth.value(i, j)/imA.sizeY();
		}
		sAEline[i] = buffer;
		APix[i] = (double)i;
	}
	#pragma omp parallel for
	for (int i = 0; i < imB.sizeX(); i++) {
		double buffer = 0.0;
		for (int j = 0; j < imB.sizeY(); j++) {
			buffer = buffer + imBSmooth.value(i, j)/imB.sizeY();
		}
		sBEline[i] = buffer;
		BPix[i] = (double)i;
	}
	//medianFilter(sAEline, 16.0);
	//medianFilter(sBEline, 16.0);

	double minBuffer = sAEline[0];
	double maxBuffer = sAEline[0];
	for (int i = 1; i < imA.sizeX(); i++) {
		if (sAEline[i] > maxBuffer) {
			maxBuffer = sAEline[i];
		}
		else {
			if (sAEline[i] < minBuffer) {
				minBuffer = sAEline[i];
			}
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < imA.sizeX(); i++) {
		sAEline[i] = (sAEline[i] - minBuffer)/(maxBuffer - minBuffer);
	}

	minBuffer = sBEline[0];
	maxBuffer = sBEline[0];
	for (int i = 1; i < imB.sizeX(); i++) {
		if (sBEline[i] > maxBuffer) {
			maxBuffer = sBEline[i];
		}
		else {
			if (sBEline[i] < minBuffer) {
				minBuffer = sBEline[i];
			}
		}
	}
	#pragma omp parallel for
	for (int i = 0; i < imB.sizeX(); i++) {
		sBEline[i] = (sBEline[i] - minBuffer) / (maxBuffer - minBuffer);
	}

	std::vector<int> peak, peakBound;
	std::vector<double> pointX, pointY;
	pointX = xRuler[0];
	pointY = yRuler[0];
	mRadAxis(eSpec, screenP, pointX, pointY);
	double pointing, eval, dbuffer;
	std::vector<int> acceptanceBound;
	acceptanceBound.resize(4, 0);
	eval = -1.0 * eSpec.angleMax(0);
	FE1DInterp(pointX, pxX, eval, dbuffer);
	acceptanceBound[0] = (int)round(dbuffer);
	eval = eSpec.angleMax(0);
	FE1DInterp(pointX, pxX, eval, dbuffer);
	acceptanceBound[1] = (int)round(dbuffer);
	eval = -1.0 * eSpec.angleMax(1);
	FE1DInterp(pointY, pxY, eval, dbuffer);
	acceptanceBound[2] = (int)round(dbuffer);
	eval = eSpec.angleMax(1);
	FE1DInterp(pointY, pxY, eval, dbuffer);
	acceptanceBound[3] = (int)round(dbuffer);

	std::vector<double> boundboxX, boundboxY;
	boundboxX.resize(5, 0.0);
	boundboxY.resize(5, 0.0);
	boundboxX[0] = acceptanceBound[0];
	boundboxY[0] = acceptanceBound[2];
	boundboxX[1] = acceptanceBound[0];
	boundboxY[1] = acceptanceBound[3];
	boundboxX[2] = acceptanceBound[1];
	boundboxY[2] = acceptanceBound[3];
	boundboxX[3] = acceptanceBound[1];
	boundboxY[3] = acceptanceBound[2];
	boundboxX[4] = acceptanceBound[0];
	boundboxY[4] = acceptanceBound[2];


	printf("Loaded 3 Images.\n");
	imageBW imBufferS;
	int maxValue = 0;
	bool flagP = 0;
	bool flagB = 0;
	double peakValue, peakValueB, totalValue, acceptValue;
	imP.crop(acceptanceBound, imBufferS);
	int peakBuffer;
	peak.resize(2, 0);
	peakBound.resize(2, 0);

	printf("Enter Horizontal Pixel Value for Peak: ");
	plt::figure();
	plt::plot(PXPix,PXSum);
	plt::show();
	std::cin >> peakBuffer;
	peak[0] = peakBuffer;

	printf("Enter Vertical Pixel Value for Peak: ");
	plt::figure();
	plt::plot(PYPix, PYSum);
	plt::show();
	std::cin >> peakBuffer;
	peak[1] = peakBuffer;

	peakValue = imP.value(peak[0], peak[1]);
	if (peak[0] < acceptanceBound[1]) {
		if (peak[0] > acceptanceBound[0]) {
			peakBound[0] = peak[0];
		}
		else {
			peakBound[0] = acceptanceBound[0];
		}
	}
	else {
		peakBound[0] = acceptanceBound[1];
	}
	if (peak[1] < acceptanceBound[3]) {
		if (peak[1] > acceptanceBound[2]) {
			peakBound[1] = peak[1];
		}
		else {
			peakBound[1] = acceptanceBound[2];
		}
	}
	else {
		peakBound[1] = acceptanceBound[3];
	}
	peakValueB = imP.value(peakBound[0], peakBound[1]);
	eval = (double)((double)imP.sizeY() - 1 - peak[1]);
	//eval = peak[1];
	FE1DInterp(pxY, pointY, eval, pointing);
	maxValue = (int)round(peakValue * 10000);
	Sum(imP, totalValue);
	Sum(imBufferS, acceptValue);
	totalValue = round(totalValue / ((double)(imP.sizeX() * imP.sizeY())) * 10000);
	acceptValue = round(acceptValue / ((double)(imBufferS.sizeX() * imBufferS.sizeY())) * 10000);
	imBufferS.destroy();
	printf("Found Electron Pointing.\n");
	
	int Apeak, Alow, Ahigh;
	printf("Enter Screen A Pixel Value for Peak: ");
	plt::figure();
	plt::plot(APix,sAEline);
	plt::show();
	std::cin >> Apeak;

	printf("Enter Screen A Pixel Value for FWHM Lower Bound: ");
	plt::figure();
	plt::plot(APix, sAEline);
	plt::show();
	std::cin >> Alow;

	printf("Enter Screen A Pixel Value for FWHM Upper Bound: ");
	plt::figure();
	plt::plot(APix, sAEline);
	plt::show();
	std::cin >> Ahigh;

	int Bpeak, Blow, Bhigh;
	printf("Enter Screen B Pixel Value for Peak: ");
	plt::figure();
	plt::plot(BPix, sBEline);
	plt::show();
	std::cin >> Bpeak;

	printf("Enter Screen B Pixel Value for FWHM Lower Bound: ");
	plt::figure();
	plt::plot(BPix, sBEline);
	plt::show();
	std::cin >> Blow;

	printf("Enter Screen B Pixel Value for FWHM Upper Bound: ");
	plt::figure();
	plt::plot(BPix, sBEline);
	plt::show();
	std::cin >> Bhigh;
	
	std::vector<double> screenPos;
	std::vector<double> EnAxis = pSpace.energy(1);
	std::vector<double> PtAxis = pSpace.pointing(1);
	std::vector<std::vector<double>> pS = pSpace.parameterSpace(1);
	double ptMax = std::max(PtAxis.front(), PtAxis.back());
	double ptMin = std::min(PtAxis.front(), PtAxis.back());
	int NE = EnAxis.size();
	screenPos.resize(NE, 0.0);
	if (pointing > ptMax) {
		for (int i = 0; i < NE; i++) {
			if (ptMax == PtAxis.front()) {
				screenPos[i] = pS[i][0];
			}
			else {
				screenPos[i] = pS[i].back();
			}
		}
	}
	else {
		if (pointing < ptMin) {
			for (int i = 0; i < NE; i++) {
				if (ptMin == PtAxis.front()) {
					screenPos[i] = pS[i][0];
				}
				else {
					screenPos[i] = pS[i].back();
				}
			}
		}
		else {
			for (int i = 0; i < NE; i++) {
				FE2DInterp(EnAxis, PtAxis, pS, EnAxis[i], pointing, screenPos[i]);
			}
		}
	}
	
	double ACenterEn, ASpreadEn, BCenterEn, BSpreadEn, inputBuffer, outputBuffer;
	inputBuffer = (double)Apeak;
	FE1DInterp(APix, xRuler[1], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, ACenterEn);

	inputBuffer = (double)Alow;
	FE1DInterp(APix, xRuler[1], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, ASpreadEn);

	inputBuffer = (double)Ahigh;
	FE1DInterp(APix, xRuler[1], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, outputBuffer);
	ASpreadEn = ASpreadEn - outputBuffer;

	printf("Screen A Centeroid Energy: %0.2e MeV\n", ACenterEn);
	printf("Screen A Energy Spread: %0.2e MeV\n", ASpreadEn);
	

	EnAxis = pSpace.energy(2);
	PtAxis = pSpace.pointing(2);
	pS = pSpace.parameterSpace(2);
	ptMax = std::max(PtAxis.front(), PtAxis.back());
	ptMin = std::min(PtAxis.front(), PtAxis.back());
	NE = EnAxis.size();
	screenPos.resize(NE, 0.0);
	if (pointing > ptMax) {
		for (int i = 0; i < NE; i++) {
			if (ptMax == PtAxis.front()) {
				screenPos[i] = pS[i][0];
			}
			else {
				screenPos[i] = pS[i].back();
			}
		}
	}
	else {
		if (pointing < ptMin) {
			for (int i = 0; i < NE; i++) {
				if (ptMin == PtAxis.front()) {
					screenPos[i] = pS[i][0];
				}
				else {
					screenPos[i] = pS[i].back();
				}
			}
		}
		else {
			for (int i = 0; i < NE; i++) {
				FE2DInterp(EnAxis, PtAxis, pS, EnAxis[i], pointing, screenPos[i]);
			}
		}
	}

	inputBuffer = (double)Bpeak;
	FE1DInterp(BPix, xRuler[2], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, BCenterEn);

	inputBuffer = (double)Blow;
	FE1DInterp(BPix, xRuler[2], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, BSpreadEn);

	inputBuffer = (double)Bhigh;
	FE1DInterp(BPix, xRuler[2], inputBuffer, outputBuffer);
	inputBuffer = outputBuffer;
	FE1DInterp(screenPos, EnAxis, inputBuffer, outputBuffer);
	BSpreadEn = BSpreadEn - outputBuffer;

	printf("Screen B Centeroid Energy: %0.2e MeV\n", BCenterEn);
	printf("Screen B Energy Spread: %0.2e MeV\n", BSpreadEn);
	
	#pragma omp parallel for
	for (int i = 0; i < imA.sizeX(); i++) {
		sAEline[i] = imA.sizeY() * (1.0 - sAEline[i] / 3.0);
	}

	#pragma omp parallel for
	for (int i = 0; i < imB.sizeX(); i++) {
		sBEline[i] = imB.sizeY() * (1.0 - sBEline[i] / 3.0);
	}

	double scale = 2;
	double lw = 2.0 * scale;
	double txtSize = 11 * scale;

	size_t resV, resH;
	double ratio;
	resH = 2224 * scale;
	double spX, spY;
	spY = ((double)imA.sizeY() + (double)imB.sizeY());
	spX = (std::max((double)imP.sizeX(), (double)imB.sizeX()));
	ratio = spY / (spX + spY);
	ratio = ratio * resH - 7;
	resV = (size_t)round(ratio);

	std::vector<double> drawLineX, drawLineY;
	drawLineX.resize(2, 0.0);
	drawLineY.resize(2, 0.0);

	printf("Drawing Image.\n");
	plt::figure_size(resH, resV);

	printf("Drawing Pointing Image.\n");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 0, 2, 1);
	pltimshow(imP, 0, "");


	plt::plot(boundboxX, boundboxY, { {"color","b"}, {"linewidth", std::to_string(lw)} });

	/*
	drawLineX[0] = 0;
	drawLineX[1] = (int)imP.sizeX() - 1;
	drawLineY[0] = (int)imP.sizeY() - 1 - peak[1];
	drawLineY[1] = drawLineY[0];
	plt::plot(drawLineX, drawLineY, { {"color","k"} });

	drawLineX[0] = peak[0];
	drawLineX[1] = drawLineX[0];
	drawLineY[0] = 0;
	drawLineY[1] = (int)imP.sizeY() - 1;
	plt::plot(drawLineX, drawLineY, { {"color","k"} });
	*/

	
	drawLineX[0] = 0;
	drawLineX[1] = (int)imP.sizeX() - 1;
	drawLineY[0] = (int)imP.sizeY() - 1 - peakBound[1];
	drawLineY[1] = drawLineY[0];
	plt::plot(drawLineX, drawLineY, { {"color","r"}, {"linewidth", std::to_string(lw)} });

	drawLineX[0] = peakBound[0];
	drawLineX[1] = drawLineX[0];
	drawLineY[0] = 0;
	drawLineY[1] = (int)imP.sizeY() - 1;
	plt::plot(drawLineX, drawLineY, { {"color","r"}, {"linewidth", std::to_string(lw)} });
	

	plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
	std::string pValue = std::to_string(pointing);
	std::string mValue = std::to_string(maxValue);
	std::string tValue = std::to_string((int)totalValue);
	std::string aValue = std::to_string((int)acceptValue);
	pValue = pValue.substr(0, 4);
	plt::text((int)(0.775 * imP.sizeX()), (int)(0.975 * imP.sizeY()), pValue + std::string(" mrad"));
	plt::text((int)(0.025 * imP.sizeX()), (int)(0.975 * imP.sizeY()), std::string("Max Px: ") + mValue + std::string("/10000"));
	plt::text((int)(0.025 * imP.sizeX()), (int)(0.100 * imP.sizeY()), std::string("AVG ToT Signal: ") + tValue);
	plt::text((int)(0.025 * imP.sizeX()), (int)(0.050 * imP.sizeY()), std::string("AVG Red Signal: ") + aValue);
	if (flagP) {
		if (flagB) {
			plt::text((int)(0.800 * imP.sizeX()), (int)(0.925 * imP.sizeY()), std::string("Error P,B"));
		}
		else {
			plt::text((int)(0.800 * imP.sizeX()), (int)(0.925 * imP.sizeY()), std::string("Error P"));
		}
	}
	else {
		if (flagB) {
			plt::text((int)(0.800 * imP.sizeX()), (int)(0.925 * imP.sizeY()), std::string("Error B"));
		}
	}
	printf("Drawing Pointing Axis: %f.\n", pointing);
	drawAxis(1, scale, eSpec, pSpace, screenP, xRuler[screenP], yRuler[screenP], pointing);
	plt::axis("off");

	printf("Drawing Energy A Image.\n");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 1, 1, (int)(spX / spY));
	//removeOutlier(imA, 4.0);
	//medianFilter(imA, 2);
	pltimshow(imASmooth, 0, "");
	plt::plot(APix, sAEline, { {"color","w"}, {"linewidth", std::to_string(lw)} });
	printf("Drawing Energy A Axis.\n");
	drawAxis(1, scale, eSpec, pSpace, screenA, xRuler[screenA], yRuler[screenA], pointing);
	plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
	std::string AEValue = std::to_string((int)std::round(ACenterEn));
	plt::text((int)(0.825 * imA.sizeX()), (int)(0.25 * imA.sizeY()), std::string("Centroid Energy: ") + AEValue + std::string(" MeV"));
	std::string AESValue = std::to_string((int)std::round(ASpreadEn));
	plt::text((int)(0.825 * imA.sizeX()), (int)(0.325 * imA.sizeY()), std::string("Energy Spread: ") + AESValue + std::string(" MeV"));
	plt::axis("off");

	printf("Drawing Energy B Image.\n");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 1, 1, 1, (int)(spX / spY));
	//removeOutlier(imB, 4.0);
	//medianFilter(imB, 2);
	pltimshow(imBSmooth, 0, "");
	plt::plot(BPix, sBEline, { {"color","w"}, {"linewidth", std::to_string(lw)} });

	printf("Drawing Energy B Axis.\n");
	drawAxis(1, scale, eSpec, pSpace, screenB, xRuler[screenB], yRuler[screenB], pointing);
	plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"}, {"font.size", std::to_string(txtSize)} });
	std::string BEValue = std::to_string((int)std::round(BCenterEn));
	plt::text((int)(0.825 * imB.sizeX()), (int)(0.25 * imB.sizeY()), std::string("Centroid Energy: ") + BEValue + std::string(" MeV"));
	std::string BESValue = std::to_string((int)std::round(BSpreadEn));
	plt::text((int)(0.825 * imB.sizeX()), (int)(0.325 * imB.sizeY()), std::string("Energy Spread: ") + BESValue + std::string(" MeV"));
	plt::axis("off");

	printf("Drawing Canvas.\n");
	plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.075}, {"hspace",0.0} });
	plt::draw();

	printf("Saving Analysis.\n");
	std::string outputHR = eSpec.analysisPath() + "/" + outputName + "-HR.png";
	std::string outputLR = eSpec.analysisPath() + "/" + outputName + ".png";
	plt::save(outputHR);
	plt::close();
	//double scaling = 0.5;
	//resizeImage(scaling, outputHR, outputLR);
	printf("Analysis Saved.\n");

	sAEline.clear();
	sBEline.clear();
	APix.clear();
	BPix.clear();
	peak.clear();
	peakBound.clear();
	pointX.clear();
	pointY.clear();
	acceptanceBound.clear();
	screenPos.clear();
	EnAxis.clear();
	PtAxis.clear();
	pS.clear();
	drawLineX.clear();
	drawLineY.clear();
	imASmooth.destroy();
	imBSmooth.destroy();
}

void pointingMode(double& rate, double& timeout, spectrometer& eSpec, screenCalibration& calibration, paramSpace & pSpace, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler) {
	int screenA, screenB, screenP;
	screenA = 1;
	screenB = 2;
	screenP = 0;
	std::string pathA, pathB, pathP, fileName, outputName, timeStamp;
	pathA = eSpec.screenPath(screenA);
	pathB = eSpec.screenPath(screenB);
	pathP = eSpec.screenPath(screenP);
	std::vector<double> viewResA, viewResB, viewResP, lineBuffer, pxX, pxY;
	viewResA = calibration.viewResolution(screenA);
	viewResB = calibration.viewResolution(screenB);
	viewResP = calibration.viewResolution(screenP);
	int pathLength = (int)pathA.length();
	std::vector<std::string> listRef, listUpdate, listB, listP;
	listDir(pathA, listRef);
	
	bool loop = 1;
	bool fileFound = 0;
	imageBW imA, imB, imP;

	int N = (int)std::max(viewResP[0], viewResP[1]);
	for (int i = 0; i < N; i++) {
		if (i < viewResP[0]) {
			pxX.push_back((double)i);
		}
		if (i < viewResP[1]) {
			pxY.push_back((double)i);
		}
	}
	
	int i = 0;
	std::vector<int> updateStatus;
	while (loop) {
		scanNewFile(pathA, listRef, updateStatus);
		for (int i = 0; i < (int)updateStatus.size(); i++) {
			if (updateStatus[i] == 1) {
				clearCMD();
				printf("Found New File.\n");
				uint fileCount = 0;
				fileName = listRef[i].substr(pathLength + 1, listRef[i].length() - pathLength - 1);
				outputName = fileName.substr(0, fileName.length() - 5);\
				int strStart = fileName.find("-");
				fileName = fileName.substr(strStart + 1, fileName.length() - strStart - 1);
				timeStamp = fileName;
				strStart = fileName.find("-");
				timeStamp = timeStamp.substr(0, strStart - 3);
				fileName = fileName.substr(strStart + 1, fileName.length() - strStart - 1);
				
				bool fileLoad;
				fileLoad = loadFile(listB, pathB, fileName, timeStamp, H[screenB], viewResB, imB);
				if (fileLoad) {
					fileCount = fileCount + 5;
				}

				fileLoad = loadFile(listP, pathP, fileName, timeStamp, H[screenP], viewResP, imP);
				if (fileLoad) {
					fileCount = fileCount + 3;
				}

				fileLoad = loadFile(listRef, pathA, fileName, timeStamp, H[screenA], viewResA, imA);
				if (fileLoad) {
					fileCount = fileCount + 1;
				}

				/*
				printf("Loading eScreen A Image.\n");
				getImage(listRef[i], imBuffer);
				perspectiveTransform(imBuffer, H[screenA], viewResA, imA);
				removeOutlier(imA, 4.0);
				medianFilter(imA, 2);
				fileCount = fileCount + 1;
				*/
				if (fileCount == 9) {
					drawPointingAnalysis(eSpec, pSpace, xRuler, yRuler, pxX, pxY, imP, imA, imB, outputName);
					imA.destroy();
					imB.destroy();
					imP.destroy();
					
				}
				else {
					switch (fileCount) {
					case 1:
						printf("Could Not Find Pointing and eScreen B.\n");
						break;
					case 3:
						printf("Could Not Find eScreen A and eScreen B.\n");
						break;
					case 4:
						printf("could Not Find eScreen B.\n");
						break;
					case 5:
						printf("Could Not Find Pointing and eScreen A.\n");
						break;
					case 6:
						printf("Could Not Find Pointing.\n");
						break;
					case 8:
						printf("Could Not Find eScreen A.\n");
						break;
					}
				}
			}
		}
		//plt::show(false);
		//plt::pause(rate);
		std::this_thread::sleep_for(std::chrono::milliseconds((long)(1000/rate)));
		if (i > (int)round(timeout / rate)) {
			loop = 0;
		}
		i++;
	}
	viewResA.clear();
	viewResB.clear();
	viewResP.clear();
	lineBuffer.clear();
	pxX.clear();
	pxY.clear();
	listRef.clear();
	listUpdate.clear();
	listB.clear();
	listP.clear();
	updateStatus.clear();
}

void pointingModeManual(double& rate, double& timeout, spectrometer& eSpec, screenCalibration& calibration, paramSpace& pSpace, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler) {
	int screenA, screenB, screenP;
	screenA = 1;
	screenB = 2;
	screenP = 0;
	std::string pathA, pathB, pathP, fileName, outputName, timeStamp, userInput;
	pathA = eSpec.screenPath(screenA);
	pathB = eSpec.screenPath(screenB);
	pathP = eSpec.screenPath(screenP);
	std::vector<double> viewResA, viewResB, viewResP, lineBuffer, pxX, pxY;
	viewResA = calibration.viewResolution(screenA);
	viewResB = calibration.viewResolution(screenB);
	viewResP = calibration.viewResolution(screenP);
	int pathLength = (int)pathP.length();
	std::vector<std::string> listA, listB, listP;

	bool loop = 1;
	bool fileFound = 0;
	imageBW imA, imB, imP;

	int N = (int)std::max(viewResP[0], viewResP[1]);
	for (int i = 0; i < N; i++) {
		if (i < viewResP[0]) {
			pxX.push_back((double)i);
		}
		if (i < viewResP[1]) {
			pxY.push_back((double)i);
		}
	}

	while (loop) {
		printf("Enter File Name (use Pointing Screen): ");
		std::cin >> fileName;
		uint fileCount = 0;
		outputName = fileName.substr(0, fileName.length() - 5); \
			int strStart = fileName.find("-");
		fileName = fileName.substr(strStart + 1, fileName.length() - strStart - 1);
		timeStamp = fileName;
		strStart = fileName.find("-");
		timeStamp = timeStamp.substr(0, strStart - 3);
		fileName = fileName.substr(strStart + 1, fileName.length() - strStart - 1);

		bool fileLoad;
		fileLoad = loadFile(listP, pathP, fileName, timeStamp, H[screenP], viewResP, imP);
		if (fileLoad) {
			fileCount = fileCount + 3;
		}

		fileLoad = loadFile(listA, pathA, fileName, timeStamp, H[screenA], viewResA, imA);
		if (fileLoad) {
			fileCount = fileCount + 1;
		}

		fileLoad = loadFile(listB, pathB, fileName, timeStamp, H[screenB], viewResB, imB);
		if (fileLoad) {
			fileCount = fileCount + 5;
		}

		if (fileCount == 9) {
			drawPointingAnalysisManual(eSpec, pSpace, xRuler, yRuler, pxX, pxY, imP, imA, imB, outputName);
			imA.destroy();
			imB.destroy();
			imP.destroy();
		}
		else {
			switch (fileCount) {
			case 1:
				printf("Could Not Find Pointing and eScreen B.\n");
				break;
			case 3:
				printf("Could Not Find eScreen A and eScreen B.\n");
				break;
			case 4:
				printf("could Not Find eScreen B.\n");
				break;
			case 5:
				printf("Could Not Find Pointing and eScreen A.\n");
				break;
			case 6:
				printf("Could Not Find Pointing.\n");
				break;
			case 8:
				printf("Could Not Find eScreen A.\n");
				break;
			}
		}

		printf("Analyze Another Image (Y/N): ");
		std::cin >> userInput;
		if (strcmp(userInput.c_str(), "N") == 0 || strcmp(userInput.c_str(), "n") == 0 || strcmp(userInput.c_str(), "0") == 0) {
			loop = 0;
		}
	}
	viewResA.clear();
	viewResB.clear();
	viewResP.clear();
	lineBuffer.clear();
	pxX.clear();
	pxY.clear();
	listA.clear();
	listB.clear();
	listP.clear();
}

void viewMode(double& rate, double& timeout, spectrometer& eSpec, screenCalibration& calibration, int& screen, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler) {
    std::string listPath = eSpec.screenPath(screen);
	std::vector<double> viewRes = calibration.viewResolution(screen);
	std::vector<std::string> listRef, listUpdate;
	paramSpace null;

    bool loop = 1;
    imageBW imBuffer, im;
    std::vector<double> line;
	int i = 0;
	std::vector<int> updateStatus;
	while (loop) {
		scanNewFile(listPath, listRef, updateStatus);
		for (int i = 0; i < (int)updateStatus.size(); i++) {
			if (updateStatus[i]) {
				plt::close();

				getImage(listRef[i], imBuffer);
				perspectiveTransform(imBuffer, H[screen], viewRes, im);
				imBuffer.destroy();
				removeOutlier(im, 4.0);
				medianFilter(im, 8);

    
				size_t resV, resH;
				double ratio;
				double scale = 2.0;
				if (screen == 0) {
					resV = 720;
					ratio = (double)im.sizeX() / (double)im.sizeY();
					ratio = ratio * (resV + 7);
					resH = (size_t)round(ratio);
				}
				else {
					resH = 1900;
					ratio = (double)im.sizeY() / (double)im.sizeX();
					ratio = ratio * resH;
					resV = (size_t)round(ratio);
				}


                std::vector<double> axis;
                double limL, limU;
				double pointing = 0.0;

                plt::figure_size(resH, resV);
                plt::subplot2grid(8, 8, 0, 0, 7, 7);
                pltimshow(im, 1, "");
				drawAxis(0, scale, eSpec, null, screen, xRuler[screen], yRuler[screen],pointing);
                plt::axis("off");
                plt::subplot2grid(8, 8, 0, 7, 7, 1);
                lineOut(1, 1, im, 0, line);
                axis.resize((int)line.size(), 0.0);
                for (int j = 0; j < (int)line.size(); j++) {
                    axis[j] = (double)j;
                }
                plt::plot(line, axis);
                limL = 0;
                limU = (double)line.size() - 1;
                plt::xlim(0.0, 1.0);
                plt::ylim(limL, limU);
                plt::tick_params({ {"labelsize","0"},{"direction","in"} });
                plt::subplot2grid(8, 8, 7, 0, 1, 7);
                lineOut(1, 1, im, 1, line);
                axis.resize((int)line.size(), 0.0);
                for (int j = 0; j < (int)line.size(); j++) {
                    axis[j] = (double)j;
                }
                plt::plot(axis, line);
                limL = 0;
                limU = (double)line.size() - 1;
                plt::xlim(limL, limU);
                plt::ylim(0.0, 1.0);
                plt::tick_params({ {"labelsize","0"},{"direction","in"} });
				plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.0}, {"hspace",0.0} });
                plt::draw();

				axis.clear();
				im.destroy();

            }
        }
        plt::show(false);
		plt::pause(rate);
        //std::this_thread::sleep_for(std::chrono::milliseconds((long)rate));
		if (i > (int)round(timeout / rate)) {
			loop = 0;
		}
		i++;
    }
    plt::close();

	viewRes.clear();
	listRef.clear();
	listUpdate.clear();
	line.clear();
	updateStatus.clear();
}

#endif