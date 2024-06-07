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
 * Draws the axis on the screen.
 *
 * @param convertToEnergyAndAngle Flag indicating whether axes need to be converted
 * @param eSpec The spectrometer object.
 * @param pSpace Trajectory endpoint surfaces
 * @param screen Which screen's axis to draw.
 * @param xAxis The X-axis values in millimeter.
 * @param yAxis The Y-axis values in millimeter.
 * @param verticalPointingAngle The angle in mrad of the peak's y value on the pointing screen
 */
void drawAxis(bool convertToEnergyAndAngle, spectrometer& eSpec, trajectoryEndpointSurfaces & pSpace, const ScreenName& screen, std::vector<double>& xAxis, std::vector<double>& yAxis, double verticalPointingAngle) {
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
 * @param pSpace The trajectoryEndpointSurfaces object.
 * @param xAxes x-axis in millimeters, one for each screen
 * @param yAxes y-axis in millimeters, one for each screen
 * @param pxX 0..Nx-1
 * @param pxY 0..Ny-1
 * @param imP Pointing image
 * @param imA low energy image
 * @param imB high energy image
 * @param filepath_spectrum The filepath for the spectrum & pointing png file
 */
void drawPointingAnalysis(spectrometer& eSpec, trajectoryEndpointSurfaces& pSpace, 
						  std::vector<std::vector<double>>& xAxes, std::vector<std::vector<double>>& yAxes, std::vector<double>& pxX, std::vector<double>& pxY, 
						  imageBW& imP, imageBW& imA, imageBW& imB, std::string filepath_spectrum
						 ) {

	// peak corresponds to the whole image, and peakBound only to that within
	// the desired max transverse angle
	std::vector<int> peak, peakBound;
	std::vector<double> xAxisPointing, yAxisPointing;
	xAxisPointing = xAxes[Pointing];
	yAxisPointing = yAxes[Pointing];
	mRadAxis(eSpec, Pointing, xAxisPointing, yAxisPointing);
	// the transverse angle in mrad along the Pointing screen y-axis at the peak
	double verticalPointingAngle;
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
	FE1DInterp(pxY, yAxisPointing, eval, verticalPointingAngle);
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
	std::string pValue = std::to_string(verticalPointingAngle);
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

	drawAxis(1, eSpec, pSpace, Pointing, xAxes[Pointing], yAxes[Pointing], verticalPointingAngle);

	plt::axis("off");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 0, 1, 1, (int)(spX / spY));
	pltimshow(imA, 1, "");
	drawAxis(1, eSpec, pSpace, LowEnergy, xAxes[LowEnergy], yAxes[LowEnergy], verticalPointingAngle);
	plt::axis("off");
	plt::subplot2grid(2, (int)((spY + spX) / spY), 1, 1, 1, (int)(spX / spY));
	pltimshow(imB, 1, "");
	drawAxis(1, eSpec, pSpace, HighEnergy, xAxes[HighEnergy], yAxes[HighEnergy], verticalPointingAngle);
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
 * @param xAxes x-axis in millimeters, one for each screen
 * @param yAxes y-axis in millimeters, one for each screen
 */
void pointingMode(std::string filepath_eScreenA, std::string filepath_eScreenB, std::string filepath_ePointing, std::string filepath_spectrum,
				  double& rate, double& timeout, 
				  spectrometer& eSpec, screenCalibration& calibration, trajectoryEndpointSurfaces & pSpace, 
				  std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xAxes, std::vector<std::vector<double>>& yAxes
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
		drawPointingAnalysis(eSpec, pSpace, xAxes, yAxes, pxX, pxY, imP, imA, imB, filepath_spectrum);

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
}

#endif