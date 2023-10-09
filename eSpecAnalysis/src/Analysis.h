#pragma once

#ifndef __Analysis_h__
#define __Analysis_h__

#include "Calibration.h"

void findSignalPeak(imageBW& image, std::vector<int>& peak, double& peakValue) {
	int Nx = image.sizeX();
	int Ny = image.sizeY();

	imageBW imSmooth = image;
	medianFilter(imSmooth, 3);
	peak.resize(2, 0);
	peakValue = imSmooth.value(0, 0);
	peak[0] = 0;
	peak[1] = 0;
	double C;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++){
			C = imSmooth.value(i, j);
			if (C > peakValue) {
				peak[0] = i;
				peak[1] = j;
				peakValue = C;
			}
		}
	}

	std::vector<double> lineY, dLineY;
	std::vector<int> countY;
	double meanY = 0.0;
	lineY.resize(Ny, 0.0);
	dLineY.resize(Ny, 0.0);
	countY.resize(Ny, 0);
	for (int i = 0; i < Ny; i++) {
		if (i > 0 && i < Ny - 1) {
			dLineY[i] = (imSmooth.value(peak[0], i + 1) - imSmooth.value(peak[0], i - 1)) / peakValue;
		}
		else {
			dLineY[i] = 0.0;
		}
		lineY[i] = imSmooth.value(peak[0], i);
		meanY = meanY + lineY[i];
	}
	meanY = meanY / Ny;

	int counter;
	bool loop = 1;
	int i = 0;
	while(loop) {
		counter = 0;
		for (int j = 1; j < Ny - i - 4; j++) {
			if (dLineY[i + j] < 0.0) {
				counter = counter + 1;
			}
			else {
				if (dLineY[i + j] > 0.0) {
					counter = counter + 1;
				}
				else {
					if (dLineY[i + j + 1] > 0.0 || dLineY[i + j + 2] > 0.0 || dLineY[i + j + 3] > 0.0) {
						counter = counter + 1;
					}
					else {
						break;
					}
				}
			}
		}
		if (counter > 0) {
			for (int j = 0; j < counter; j++) {
				countY[i + j + 1] = counter;
			}
			i = i + counter;
		}
		else {
			i++;
		}
		if (i > Ny - 2) {
			loop = 0;
		}
	}
	
	loop = 1;
	i = 1;
	while (loop) {
		if (dLineY[i] != 0.0 && dLineY[i + 1] != 0.0 && dLineY[i + 2] != 0.0) {
			for (int j = 0; j < countY[i] + 1; j++) {
				lineY[i + j] = ceil((lineY[i - 1] + lineY[i + countY[i]]) / 2.0);
			}
			i = i + countY[i] + 1;
		}
		else {
			i++;
		}
		
		if (i > Ny - 1) {
			loop = 0;
		}
	}


	std::vector<double> lineX, dLineX;
	std::vector<int> countX;
	double meanX = 0.0;
	lineX.resize(Nx, 0.0);
	dLineX.resize(Nx, 0.0);
	countX.resize(Nx, 0);
	for (int i = 0; i < Nx; i++) {
		if (i > 0 && i < Nx - 1) {
			dLineX[i] = (imSmooth.value(i + 1, peak[1]) - imSmooth.value(i - 1, peak[1])) / peakValue;
		}
		else {
			dLineX[i] = 0.0;
		}
		lineX[i] = imSmooth.value(i, peak[0]);
		meanX = meanX + lineX[i];
	}
	meanX = meanX / Nx;

	counter;
	loop = 1;
	i = 0;
	while (loop) {
		counter = 0;
		for (int j = 1; j < Nx - i - 4; j++) {
			if (dLineX[i + j] < 0.0) {
				counter = counter + 1;
			}
			else {
				if (dLineX[i + j] > 0.0) {
					counter = counter + 1;
				}
				else {
					if (dLineX[i + j + 1] > 0.0 || dLineX[i + j + 2] > 0.0 || dLineX[i + j + 3] > 0.0) {
						counter = counter + 1;
					}
					else {
						break;
					}
				}
			}
		}
		if (counter > 0) {
			for (int j = 0; j < counter; j++) {
				countX[i + j + 1] = counter;
			}
			i = i + counter;
		}
		else {
			i++;
		}
		if (i > Nx - 2) {
			loop = 0;
		}
	}

	loop = 1;
	i = 1;
	while (loop) {
		if (dLineX[i] != 0.0 && dLineX[i + 1] != 0.0 && dLineX[i + 2] != 0.0) {
			for (int j = 0; j < countX[i] + 1; j++) {
				lineX[i + j] = ceil((lineX[i - 1] + lineX[i + countX[i]]) / 2.0);
			}
			i = i + countX[i] + 1;
		}
		else {
			i++;
		}
		if (i > Nx - 1) {
			loop = 0;
		}
	}

	double xValue, yValue, a1, a2, a0;
	int N;
	std::vector<double> sumX, sumY;
	cv::Mat quadRegMat = cv::Mat::zeros(3, 3, CV_64F);

	N = (int)std::min(std::min((double)peak[0], (double)(Nx - peak[0])) - 1, 250.0);
	sumX.resize(5, 0.0);
	sumY.resize(3, 0.0);

	for (int i = 0; i < 2 * N + 1; i++) {
		xValue = (double)i + (double)peak[0] - (double)N;
		yValue = lineX[peak[0] + i - N];
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
	if (xValue >= 0 && xValue < (int)lineX.size()) {
		peak[0] = xValue;
	}

	/*
	std::vector<double> fX;
	fX.resize(Nx, 0.0);
	for (int i = 0; i < Nx; i++) {
		fX[i] = a0 + a1 * (double)i + a2 * (double)i * (double)i;
	}
	*/

	N = (int)std::min(std::min((double)peak[1], (double)(Ny - peak[1])) - 1, 250.0);
	sumX.resize(5, 0.0);
	sumY.resize(3, 0.0);

	for (int i = 0; i < 2 * N + 1; i++) {
		xValue = (double)i + (double)peak[1] - (double)N;
		yValue = lineY[peak[1] + i - N];
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

	yValue = -a1 / (2.0 * a2);
	if (yValue >= 0 && yValue < (int)lineY.size()) {
		peak[1] = yValue;
	}

	/*
	std::vector<double> fY;
	fY.resize(Ny, 0.0);
	for (int i = 0; i < Ny; i++) {
		fY[i] = a0 + a1 * (double)i + a2 * (double)i * (double)i;
	}


	std::cout << peak[0] << "\t" << peak[1] << "\n";

	plt::figure();
	plt::subplot(1, 2, 1);
	plt::plot(lineX);
	plt::plot(fX);
	plt::subplot(1, 2, 2);
	plt::plot(lineY);
	plt::plot(fY);
	plt::show();
	*/

}

void findPointing(spectrometer& eSpec, imageBW& image, std::vector<double>& rulerX, std::vector<double>&rulerY, std::vector<double>& pointing) {
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

void drawAxis(bool mode, spectrometer& eSpec, paramSpace & pSpace, int& screen, std::vector<double>& rulerX, std::vector<double>& rulerY, double pointing) {
	std::vector<double> pixelX, pixelY;
	int Nx = (int)rulerX.size();
	int Ny = (int)rulerY.size();
	pixelX.resize(Nx);
	pixelY.resize(Ny);

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
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"} });
			for (int i = 0; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plotY[0] = yZero - 10;
				plotY[1] = yZero + 10;
				plt::plot(plotX, plotY, { {"color","w"} });
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
				plt::plot(plotX, plotY, { {"color","w"} });
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
			plt::plot(plotX, plotY, { {"color","w"} });
			plotX[0] = xZero;
			plotX[1] = xZero;
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"} });
		}
		else {
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"} });
			double locationX, locationY;
			for (int i = 1; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plotY[0] = yAxis[0];
				plotY[1] = 0;
				plt::plot(plotX, plotY, { {"color","w"} });
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
			plt::rcparams({ {"text.color", "k"}, {"font.weight", "bold"} });
			for (int i = 0; i < Ny; i++) {
				plotX[0] = xAxis[0];
				plotX[1] = 0;
				plotY[0] = yAxis[i];
				plotY[1] = yAxis[i];
				plt::plot(plotX, plotY, { {"color","w"} });
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
			plt::plot(plotX, plotY, { {"color","w"} });
			plotX[0] = xAxis[0];
			plotX[1] = xAxis[0];
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"} });
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
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"} });
			for (int i = 0; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plotY[0] = yZero - 25;
				plotY[1] = yZero + 25;
				plt::plot(plotX, plotY, { {"color","w"} });
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
				plt::plot(plotX, plotY, { {"color","w"} });
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
			plt::plot(plotX, plotY, { {"color","w"} });
			plotX[0] = xZero;
			plotX[1] = xZero;
			plotY[0] = 0;
			plotY[1] = (int)rulerY.size() - 1;
			plt::plot(plotX, plotY, { {"color","w"} });
		}
		else {
			plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"} });
			plotX[0] = 0;
			plotX[1] = (int)rulerX.size() - 1;
			double locationX, locationY;
			if (screen == 1) {
				plotY[0] = yAxis[1];
				plotY[1] = yAxis[1];
				locationY = 25;
			}
			else {
				plotY[0] = yAxis[0];
				plotY[1] = yAxis[0];
				locationY = 20;
			}
			plt::plot(plotX, plotY, { {"color","w"} });
			if (screen == 1) {
				plotY[1] = yAxis[1] - 20;
			}
			else {
				plotY[1] = yAxis[0] - 20;
			}
			int lastLabel = 1;
			for (int i = 1; i < Nx; i++) {
				plotX[0] = xAxis[i];
				plotX[1] = xAxis[i];
				plt::plot(plotX, plotY, { {"color","w"} });
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
			plt::rcparams({ {"text.color", "k"}, {"font.weight", "bold"} });
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
			plt::plot(plotX, plotY, { {"color","w"} });
			plotX[1] = 0;
			for (int i = 0; i < Ny; i++) {
				plotY[0] = yAxis[i];
				plotY[1] = yAxis[i];
				if (yTickStart - 5 * i > -15 && yTickStart - 5 * i < 15) {
					plt::plot(plotX, plotY, { {"color","w"} });
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
	}
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
}

void pointingMode(double& rate, double& timeout, spectrometer& eSpec, screenCalibration& calibration, paramSpace & pSpace, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler) {
	int screenA, screenB, screenP;
	screenA = 1;
	screenB = 2;
	screenP = 0;
	std::string pathA, pathB, pathP, fileName, outputName;
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
	
	double loopWait = 250;
	bool fileFindLoop = 1;
	int fileLoopCount = 0;
	bool loop = 1;
	bool flagP = 0;
	bool flagB = 0;
	bool fileFound = 0;
	imageBW imBuffer, imA, imB, imP;

	int N = (int)std::max(viewResP[0], viewResP[1]);
	for (int i = 0; i < N; i++) {
		if (i < viewResP[0]) {
			pxX.push_back((double)i);
		}
		if (i < viewResP[1]) {
			pxY.push_back((double)i);
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
	
	int i = 0;
	std::vector<int> updateStatus;
	while (loop) {
		scanNewFile(pathA, listRef, updateStatus);
		for (int i = 0; i < (int)updateStatus.size(); i++) {
			if (updateStatus[i]) {
				plt::close();
				uint fileCount = 0;
				fileName = listRef[i].substr(pathLength + 1, listRef[i].length() - pathLength - 1);
				outputName = fileName.substr(0, fileName.length() - 5);
				int strStart = fileName.find("-");
				fileName = fileName.substr(strStart + 1, fileName.length() - strStart - 1);
				strStart = fileName.find("-");
				fileName = fileName.substr(strStart + 1, fileName.length() - strStart - 1);

				fileLoopCount = 0;
				fileFound = 0;
				fileFindLoop = 1;
				int index = -1;
				while (fileFindLoop) {
					listDir(pathB, listB);
					fileFound = findFile(listB, fileName, index);
					if (fileFound) {
						fileFindLoop = 0;
					}
					if (fileLoopCount > 30) {
						fileFindLoop = 0;
						fileFound = 0;
						std::cout << "Can not find eScreen B file.\n";
					}
					fileLoopCount = fileLoopCount + 1;
					std::this_thread::sleep_for(std::chrono::milliseconds((long)loopWait));
				}
				if (fileFound) {
					getImage(listB[index], imBuffer);
					perspectiveTransform(imBuffer, H[screenB], viewResB, imB);
					removeOutlier(imB, 4.0);
					medianFilter(imB, 2);
					fileCount = fileCount + 5;
				}

				fileLoopCount = 0;
				fileFound = 0;
				fileFindLoop = 1;
				index = -1;
				while (fileFindLoop) {
					listDir(pathP, listP);
					fileFound = findFile(listP, fileName, index);
					if (fileFound) {
						fileFindLoop = 0;
					}
					if (fileLoopCount > 30) {
						fileFindLoop = 0;
						fileFound = 0;
						std::cout << "Can not find Pointing file.\n";
					}
					fileLoopCount = fileLoopCount + 1;
					std::this_thread::sleep_for(std::chrono::milliseconds((long)loopWait));
				}
				if (fileFound) {
					getImage(listP[index], imBuffer);
					perspectiveTransform(imBuffer, H[screenP], viewResP, imP);
					removeOutlier(imP, 4.0);
					medianFilter(imP, 2);
					fileCount = fileCount + 3;
				}

				getImage(listRef[i], imBuffer);
				perspectiveTransform(imBuffer, H[screenA], viewResA, imA);
				removeOutlier(imA, 4.0);
				medianFilter(imA, 2);
				fileCount = fileCount + 1;

				if (fileCount == 9) {
					double peakValue, peakValueB;
					imP.crop(acceptanceBound, imBuffer);
					findSignalPeak(imBuffer, peakBound, peakValueB);
					findSignalPeak(imP, peak, peakValue);
					peakBound[0] = peakBound[0] + acceptanceBound[0];
					peakBound[1] = peakBound[1] + acceptanceBound[2];
					eval = (double)((double)imP.sizeY() - 1 - peakBound[1]);
					FE1DInterp(pxY, pointY, eval, pointing);


					size_t resV, resH;
					double ratio;
					resH = 1900;
					double spX, spY;
					spY = ((double)imA.sizeY() + (double)imB.sizeY());
					spX = (std::max((double)imP.sizeX(), (double)imB.sizeX()));
					ratio = spY / (spX + spY);
					ratio = ratio * resH - 7;
					resV = (size_t)round(ratio);

					std::vector<double> drawLineX, drawLineY;
					drawLineX.resize(2, 0.0);
					drawLineY.resize(2, 0.0);


					plt::figure_size(resH, resV);
					plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 0, 2, 1);
					pltimshow(imP, 1, "");


					plt::plot(boundboxX, boundboxY, { {"color","r"} });

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

					drawLineX[0] = 0;
					drawLineX[1] = (int)imP.sizeX() - 1;
					drawLineY[0] = (int)imP.sizeY() - 1 - peakBound[1];
					drawLineY[1] = drawLineY[0];
					plt::plot(drawLineX, drawLineY, { {"color","b"} });

					drawLineX[0] = peakBound[0];
					drawLineX[1] = drawLineX[0];
					drawLineY[0] = 0;
					drawLineY[1] = (int)imP.sizeY() - 1;
					plt::plot(drawLineX, drawLineY, { {"color","b"} });

					plt::rcparams({ {"text.color", "w"}, {"font.weight", "bold"} });
					std::string pValue = std::to_string(pointing);
					pValue = pValue.substr(0, 4);
					plt::text((int)(0.775 * imP.sizeX()), (int)(0.975 * imP.sizeY()), pValue + std::string(" mrad"));
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

					drawAxis(1, eSpec, pSpace, screenP, xRuler[screenP], yRuler[screenP], pointing);

					plt::axis("off");
					plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 1, 1, (int)(spX / spY));
					pltimshow(imA, 1, "");
					drawAxis(1, eSpec, pSpace, screenA, xRuler[screenA], yRuler[screenA], pointing);
					plt::axis("off");
					plt::subplot2grid(2, (int)((spY + spX) / spY), 1, 1, 1, (int)(spX / spY));
					pltimshow(imB, 1, "");
					drawAxis(1, eSpec, pSpace, screenB, xRuler[screenB], yRuler[screenB], pointing);
					plt::axis("off");
					plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.075}, {"hspace",0.0} });
					plt::draw();

					outputName = eSpec.analysisPath() + "/" + outputName + ".png";
					plt::save(outputName);
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
		plt::show(false);
		plt::pause(rate);
		if (i > (int)round(timeout / rate)) {
			loop = 0;
		}
		i++;
	}
	plt::close();
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
				removeOutlier(im, 4.0);
				medianFilter(im, 2);

    
				size_t resV, resH;
				double ratio;
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
				drawAxis(0, eSpec, null, screen, xRuler[screen], yRuler[screen],pointing);
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
}

#endif