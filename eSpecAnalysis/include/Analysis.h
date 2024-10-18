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
	imageBW imSmooth;
	image.copy(imSmooth);
	removeOutlier(imSmooth, 4.0);
	medianFilter(imSmooth, 4);
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
	double lineSum = 0.0;
	double linePV = 0.0;
	double lineCV = 0.0;
	int linePI = 0;
	int indexBuffer = 0;

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

	#pragma omp parallel for reduction(+:linePV,lineSum)
	for (int i = 0; i < Nx; i++) {
		lineCV = lineCV + lineValuesX[i] * (double)i;
		lineSum = lineSum + lineValuesX[i];
	}
	indexBuffer = (int)round(lineCV / lineSum);
	if (lineValuesX[indexBuffer] > linePV) {
		linePI = indexBuffer;
		linePV = lineValuesX[linePI];
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
	lineSum = 0.0;
	linePV = 0.0;
	lineCV = 0.0;
	linePI = 0;
	indexBuffer = 0;

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

	#pragma omp parallel for reduction(+:linePV,lineSum)
	for (int i = 0; i < Ny; i++) {
		lineCV = lineCV + lineValuesY[i] * (double)i;
		lineSum = lineSum + lineValuesY[i];
	}
	indexBuffer = (int)round(lineCV / lineSum);
	if (lineValuesY[indexBuffer] > linePV) {
		linePI = indexBuffer;
		linePV = lineValuesX[linePI];
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


/**
 * @brief Convert angle axes to mrad
 * 
 * For axes representing transverse angles, calculate their values in mrad. For 
 * the pointing screen, that's both x and y axes. For LowEnergy and HighEnergy
 * that's just the y axis.
 * 
 * @param eSpec The spectrometer object.
 * @param screen The screen to find the axis/axes of.
 * @param xAxis The x-axis that is converted inplace to angle in mrad
 * @param yAxis The y-axis that is converted inplace to angle in mrad
 */
void mRadAxis(spectrometer& eSpec, const ScreenName& screen, std::vector<double>& xAxis, std::vector<double>& yAxis) {
	int Nx = (int)xAxis.size();
	int Ny = (int)yAxis.size();
	int N = std::max(Nx, Ny);

	double x, y, z, phi, theta;
	phi = eSpec.phi(0) / 180 * pi;
	theta = eSpec.theta(0) / 180 * pi;

	if (screen == Pointing) {
		for (int i = 0; i < N; i++) {
			if (i < Nx) {
				z = eSpec.z(0) + xAxis[i] * std::cos(phi);
				x = -(xAxis[i] * std::sin(phi) - eSpec.x(0));

				xAxis[i] = 1000.0 * std::atan2(x,z);
			}
			if (i < Ny) {
				z = eSpec.z(0) + yAxis[i] * std::sin(theta);
				y =  -(yAxis[i] * std::cos(theta) - eSpec.y(0));

				yAxis[i] = 1000.0 * std::atan2(y,z);
			}
		}
	}
	else {
		for (int i = 0; i < Ny; i++) {
			z = eSpec.z(0) + yAxis[i] * std::sin(phi) + xAxis[0] * std::cos(theta);
			y = yAxis[i] * std::cos(phi) - eSpec.x(0);

			yAxis[i] = 1000.0 * std::atan2(y, z);
		}
	}
}


void getPointing(imageBW& pointingImage, 
				  std::vector<double>& xAxis, std::vector<double>& yAxis,
				  std::vector<double>& pxX, std::vector<double>& pxY,
				  spectrometer& eSpec, 
				  double verticalPointingAngle
				 ) {
	// peak corresponds to the whole image, and peakBound only to that within
	// the desired max transverse angle
	std::vector<int> peak, peakBound;

	std::vector<double> xAxisPointing, yAxisPointing;
	xAxisPointing = xAxis;
	yAxisPointing = yAxis;

	mRadAxis(eSpec, Pointing, xAxisPointing, yAxisPointing);
	// the transverse angle in mrad along the Pointing screen y-axis at the peak
	double eval, dbuffer;
	
	// acceptanceBound is the pixel values on the pointing screen corresponding 
	// to the desired maximum transverse angles, in mrad, as [xmin, xmax, ymin, ymax]
	std::vector<int> acceptanceBound;
	acceptanceBound.resize(4, 0);
	eval = -1.0 * eSpec.angleMax(0);
	FE1DInterp(xAxisPointing, pxX, eval, dbuffer);
	acceptanceBound[0] = (int)round(dbuffer);
	eval = eSpec.angleMax(0);
	FE1DInterp(xAxisPointing, pxX, eval, dbuffer);
	acceptanceBound[1] = (int)round(dbuffer);
	eval = -1.0 * eSpec.angleMax(1);
	FE1DInterp(yAxisPointing, pxY, eval, dbuffer);
	acceptanceBound[2] = (int)round(dbuffer);
	eval = eSpec.angleMax(1);
	FE1DInterp(yAxisPointing, pxY, eval, dbuffer);
	acceptanceBound[3] = (int)round(dbuffer);

	// analyze Pointing image
	imageBW imBuffer;
	int maxValue = 0;
	bool flagP = 0;
	bool flagB = 0;
	// peakValue and totalValue correspond to the whole image, and peakValueB 
	// and acceptValue only to that within the desired max transverse angle
	double peakValue, peakValueB, totalValue, acceptValue;
	pointingImage.crop(acceptanceBound, imBuffer);
	findSignalPeak(imBuffer, peakBound, peakValueB);
	findSignalPeak(pointingImage, peak, peakValue);
	// put peakBound back into the whole image coordinates
	peakBound[0] = peakBound[0] + acceptanceBound[0];
	peakBound[1] = peakBound[1] + acceptanceBound[2];

	// get peak y location in mrad
	eval = (double)((double)pointingImage.sizeY() - 1 - peakBound[1]);
	FE1DInterp(pxY, yAxisPointing, eval, verticalPointingAngle);

}



/**
 * @brief Converts axes of a screen from spatial (mm) into angle or energy
 *
 * For Pointing, that's x-axis and y-axis to angle in mrad
 * For LowEnergy and HighEnergy, that's x-axis to energy in MeV and y-axis to angle in mrad
 * 
 * @param xAxis The x-axis in millimeter.
 * @param yAxis The y-axis in millimeter.
 * @param screen Which screen this is.
 * @param pSpace The trajectory endpoint surfaces.
 * @param eSpec The spectrometer.
 */
void transformAxes(std::vector<double>& xAxis, std::vector<double>& yAxis, const ScreenName& screen, trajectoryEndpointSurfaces& pSpace, double verticalPointingAngle, spectrometer& eSpec) {
	std::vector<double> pixelX, pixelY;
	int Nx = (int)xAxis.size();
	int Ny = (int)yAxis.size();
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


	// convert axes representing angles into angle values in mrad
	int warning = 0;
	mRadAxis(eSpec, screen, xAxis, yAxis);

	// convert axes representing energies into energy values in MeV
	if (screen == LowEnergy || screen == HighEnergy)  {

		int NE;
		int indexStart, indexEnd;
		// screenPos is the trajectory endpoint values for angle = verticalPointingAngle
		std::vector<double> screenPos;

		std::vector<double> EnAxis = pSpace.getEnergyAxis(screen);
		std::vector<double> PtAxis = pSpace.getPointingAxis(screen);
		std::vector<std::vector<double>> pS = pSpace.getTrajectoryEndpointSurface(screen);


		double ptMax, ptMin, enMin, enMax;
		double ptEval, xEval;

		ptMax = std::max(PtAxis.front(), PtAxis.back());
		ptMin = std::min(PtAxis.front(), PtAxis.back());
		enMax = std::max(EnAxis.front(), EnAxis.back());
		enMin = std::min(EnAxis.front(), EnAxis.back());
		double dE = abs(EnAxis[1] - EnAxis[0]);

		ptEval = std::clamp(verticalPointingAngle, ptMin, ptMax);

		// get screen x position in mm for each energy value in pS for given 
		// vertical pointing angle
		NE = (int)pSpace.getEnergyAxis(screen).size();
		screenPos.resize(NE, 0.0);
		for (int i = 0; i < NE; i++) {
			FE1DInterp(PtAxis, pS[i], ptEval, screenPos[i]);
		}

		// get energy corresponding to each xAxis value 
		NE = (int)xAxis.size();
		for (int i = 0; i < NE; i++) {
			xEval = std::clamp(xAxis[i], enMin, enMax);
			FE1DInterp(screenPos, EnAxis, xEval, xAxis[i]);
		}
	}
}

void loadFile(std::string& filepath, cv::Mat& H, std::vector<double>& viewRes, imageBW& output) {
	
	imageBW imBuffer;

	getImage(filepath, imBuffer);
	perspectiveTransform(imBuffer, H, viewRes, output);
	removeOutlier(output, 4.0);
	medianFilter(output, 2);
	imBuffer.destroy();

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
 * @param xAxes x-axis in millimeters, one for each screen
 * @param yAxes y-axis in millimeters, one for each screen
 */
void pointingMode(std::string filepath_eScreenA, std::string filepath_eScreenB, std::string filepath_ePointing, std::string filepath_spectrum,
				  double& rate, double& timeout, 
				  spectrometer& eSpec, screenCalibration& calibration, trajectoryEndpointSurfaces & pSpace, 
				  std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xAxes, std::vector<std::vector<double>>& yAxes
				 ) {

	std::vector<double> viewResA, viewResB, viewResP, pxX, pxY;
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

	uint fileCount = 0;

	loadFile(filepath_eScreenA,  H[LowEnergy], viewResA, imA);
	loadFile(filepath_eScreenB,  H[HighEnergy], viewResB, imB);
	loadFile(filepath_ePointing, H[Pointing], viewResP, imP);

	// historical. Used to check if the files were found.
	fileCount = 15;

	if (fileCount == 15) {

		double verticalPointingAngle;
		getPointing(imP, xAxes[Pointing], yAxes[Pointing], pxX, pxY, eSpec, verticalPointingAngle);

		transformAxes(xAxes[LowEnergy], yAxes[LowEnergy], LowEnergy, pSpace, verticalPointingAngle, eSpec);
		transformAxes(xAxes[HighEnergy], yAxes[HighEnergy], HighEnergy, pSpace, verticalPointingAngle, eSpec);
		transformAxes(xAxes[Pointing], yAxes[Pointing], Pointing, pSpace, verticalPointingAngle, eSpec);

		saveCroppedTransformedImage(imA, xAxes[LowEnergy], yAxes[LowEnergy], filepath_eScreenA);
		saveCroppedTransformedImage(imB, xAxes[HighEnergy], yAxes[HighEnergy], filepath_eScreenB);
		saveCroppedTransformedImage(imP, xAxes[Pointing], yAxes[Pointing], filepath_ePointing);
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
	
	imA.destroy();
	imB.destroy(); 
	imP.destroy();

	viewResA.clear();
	viewResB.clear();
	viewResP.clear();
	pxX.clear();
	pxY.clear();
}

#endif
