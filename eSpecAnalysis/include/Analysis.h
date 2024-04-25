#pragma once

#ifndef __Analysis_h__
#define __Analysis_h__

#include "Calibration.h"

/**
 *
 * @param image The image to analyze.
 * @param peak output (x, y)
 * @param peakValue output (smoothed) value at peak
 */
void findSignalPeak(imageBW& image, std::vector<int>& peak, double& peakValue) {
	imageBW imSmooth = image;
	medianFilter(imSmooth, 3);
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


/**
 * @brief Convert angle axes to mrad
 * 
 * For axes representing transverse angles, calculate their values in mrad. For 
 * the pointing screen, that's both x and y axes. For LowEnergy and HighEnergy
 * that's just the y axis.
 * 
 * @param eSpec The spectrometer object.
 * @param screen The screen to find the axis/axes of.
 * @param rulerX The x-axis that is converted inplace to angle in mrad
 * @param rulerY The y-axis that is converted inplace to angle in mrad
 */
void mRadAxis(spectrometer& eSpec, const ScreenName& screen, std::vector<double>& rulerX, std::vector<double>& rulerY) {
	int Nx = (int)rulerX.size();
	int Ny = (int)rulerY.size();
	int N = std::max(Nx, Ny);

	double x, y, z, phi, theta;
	phi = eSpec.phi(0) / 180 * pi;
	theta = eSpec.theta(0) / 180 * pi;

	if (screen == Pointing) {
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

void drawAxis(bool mode, spectrometer& eSpec, paramSpace & pSpace, const ScreenName& screen, std::vector<double>& rulerX, std::vector<double>& rulerY, double pointing) {
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
		if (screen == Pointing) {
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
		else { // screen is LowEnergy or HighEnergy
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
	else {  // mode != 0 
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
		if (screen == Pointing) {
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
		else {  // screen is LowEnergy or HighEnergy
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
		if (screen == Pointing) {
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

void loadFile_old(std::vector<std::string>& list, std::string& path, std::string& filename, std::string& timeStamp, cv::Mat& H, std::vector<double>& viewRes, uint& fileCount, imageBW& output) {
	int fileLoopCount = 0;
	bool fileFound = 0;
	bool fileFindLoop = 1;
	int index = -1;
	double loopWait = 250;

	imageBW imBuffer;
	while (fileFindLoop) {
		listDir(path, list);
		fileFound = findFile(list, filename, timeStamp, index);
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
		removeOutlier(output, 4.0);
		medianFilter(output, 2);
		fileCount = fileCount + 5;
	}
}

void loadFile(std::string& filepath, cv::Mat& H, std::vector<double>& viewRes, imageBW& output) {
	
	imageBW imBuffer;

	getImage(filepath, imBuffer);
	perspectiveTransform(imBuffer, H, viewRes, output);
	removeOutlier(output, 4.0);
	medianFilter(output, 2);

}



/**
 * @brief Computes the spectrum and draws it.
 *
 * @param eSpec The spectrometer object.
 * @param pSpace The paramSpace object.
 * @param xRuler x-axis in millimeters.
 * @param yRuler y-axis in millimeters
 * @param pxX 0..Nx-1
 * @param pxY 0..Ny-1
 * @param imP Pointing image
 * @param imA low energy image
 * @param imB high energy image
 * @param filepath_spectrum The filepath for the spectrum & pointing png file
 */
void drawPointingAnalysis(spectrometer& eSpec, paramSpace& pSpace, 
						  std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler, std::vector<double>& pxX, std::vector<double>& pxY, 
						  imageBW& imP, imageBW& imA, imageBW& imB, std::string filepath_spectrum
						 ) {

	// peak corresponds to the whole image, and peakBound only to that within
	// the desired max transverse angle
	std::vector<int> peak, peakBound;
	std::vector<double> pointX, pointY;
	pointX = xRuler[0];
	pointY = yRuler[0];
	mRadAxis(eSpec, Pointing, pointX, pointY);
	double pointing, eval, dbuffer;
	
	// acceptanceBound is the pixel values on the pointing screen corresponding 
	// to the desired maximum transverse angles, in mrad, as [xmin, xmax, ymin, ymax]
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
	
	// analyze Pointing image
	imageBW imBuffer;
	int maxValue = 0;
	bool flagP = 0;
	bool flagB = 0;
	// peakValue and totalValue correspond to the whole image, and peakValueB 
	// and acceptValue only to that within the desired max transverse angle
	double peakValue, peakValueB, totalValue, acceptValue;
	imP.crop(acceptanceBound, imBuffer);
	findSignalPeak(imBuffer, peakBound, peakValueB);
	findSignalPeak(imP, peak, peakValue);
	// put peakBound back into the whole image coordinates
	peakBound[0] = peakBound[0] + acceptanceBound[0];
	peakBound[1] = peakBound[1] + acceptanceBound[2];
	eval = (double)((double)imP.sizeY() - 1 - peakBound[1]);
	FE1DInterp(pxY, pointY, eval, pointing);
	maxValue = (int)round(peakValue * 10000);
	Sum(imP, totalValue);
	Sum(imBuffer, acceptValue);
	// totalValue and acceptValue are actualy mean pixel values
	totalValue = round(totalValue / ((double)(imP.sizeX() * imP.sizeY())) * 10000);
	acceptValue = round(acceptValue / ((double)(imBuffer.sizeX() * imBuffer.sizeY())) * 10000);
	printf("Found Pointing.\n");

	size_t resV, resH;
	double ratio;
	resH = 2224;
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

	drawAxis(1, eSpec, pSpace, Pointing, xRuler[Pointing], yRuler[Pointing], pointing);

	plt::axis("off");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 1, 1, (int)(spX / spY));
	pltimshow(imA, 1, "");
	drawAxis(1, eSpec, pSpace, LowEnergy, xRuler[LowEnergy], yRuler[LowEnergy], pointing);
	plt::axis("off");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 1, 1, 1, (int)(spX / spY));
	pltimshow(imB, 1, "");
	drawAxis(1, eSpec, pSpace, HighEnergy, xRuler[HighEnergy], yRuler[HighEnergy], pointing);
	plt::axis("off");
	plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.075}, {"hspace",0.0} });
	plt::draw();

	printf("Saving Analysis.\n");
	plt::save(filepath_spectrum);
	double scaling = 0.5;
	resizeImage(scaling, filepath_spectrum, filepath_spectrum);
	printf("Analysis Saved.\n");
}



/**
 *
 * @param filepath_eScreenA The file path of the low energy image.
 * @param filepath_eScreenB The file path of the high energy image.
 * @param filepath_ePointing The file path of the pointing image.
 * @param filepath_spectrum The file path to save the resulting spectrum image.
 * @param rate Not used.
 * @param timeout Not used.
 * @param eSpec The spectrometer object.
 * @param calibration The screen calibration object.
 * @param pSpace The parameter space object.
 * @param H Homography matrices, one for each screen
 * @param xRuler x-axis in millimeters, one for each screen
 * @param yRuler y-axis in millimeters, one for each screen
 */
void pointingMode(std::string filepath_eScreenA, std::string filepath_eScreenB, std::string filepath_ePointing, std::string filepath_spectrum,
				  double& rate, double& timeout, 
				  spectrometer& eSpec, screenCalibration& calibration, paramSpace & pSpace, 
				  std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler
				 ) {

	std::vector<double> viewResA, viewResB, viewResP, lineBuffer, pxX, pxY;
	// (winSize_x, winSize_y) for each screen
	viewResA = calibration.viewResolution(LowEnergy);
	viewResB = calibration.viewResolution(HighEnergy);
	viewResP = calibration.viewResolution(Pointing);
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
	
	plt::close();
	imA.destroy();
	imB.destroy(); 
	imP.destroy();
	uint fileCount = 0;

	loadFile(filepath_eScreenA,  H[LowEnergy], viewResA, imA);
	loadFile(filepath_eScreenB,  H[HighEnergy], viewResB, imB);
	loadFile(filepath_ePointing, H[Pointing], viewResP, imP);

	// historical. Used to check if the files were found.
	fileCount = 15;

	if (fileCount == 15) {
		drawPointingAnalysis(eSpec, pSpace, xRuler, yRuler, pxX, pxY, imP, imA, imB, filepath_spectrum);
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

#endif