#pragma once

#ifndef __Calibration_h__
#define __Calibration_h__

#include "Plot.h"

void homographyMat(std::vector<cv::Point2d> points, cv::Mat& matrixH) {
    std::vector<cv::Point2d> sourcePoints;
    std::vector<cv::Point2d> targetPoints;

    targetPoints.push_back(points[0]);
    targetPoints.push_back(points[1]);
    targetPoints.push_back(points[2]);
    targetPoints.push_back(points[3]);



    sourcePoints.push_back(points[4]);
    sourcePoints.push_back(points[5]);
    sourcePoints.push_back(points[6]);
    sourcePoints.push_back(points[7]);

    matrixH = cv::findHomography(sourcePoints, targetPoints);
}

void transformPixel(cv::Mat& matrixH, cv::Point3d& pixel) {
    double x, y, z;
    x = matrixH.at<double>(0, 0) * pixel.x + matrixH.at<double>(0, 1) * pixel.y + matrixH.at<double>(0, 2) * pixel.z;
    y = matrixH.at<double>(1, 0) * pixel.x + matrixH.at<double>(1, 1) * pixel.y + matrixH.at<double>(1, 2) * pixel.z;
    z = matrixH.at<double>(2, 0) * pixel.x + matrixH.at<double>(2, 1) * pixel.y + matrixH.at<double>(2, 2) * pixel.z;

    pixel.x = x;
    pixel.y = y;
    pixel.z = z;
}

void perspectiveTransform(imageBW& input, cv::Mat matrixH, std::vector<double>& windowSize, imageBW& output) {
    size_t Nx = (size_t)windowSize[0];
    size_t Ny = (size_t)windowSize[1];

    std::vector<size_t> bounds;
    bounds.resize(4, 0);
    bounds[0] = 0;
    bounds[1] = Nx;
    bounds[2] = 0;
    bounds[3] = Ny;
    input.crop(bounds, output);

    #pragma omp parallel for
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                cv::Point3d pixel((double)i, (double)j, 1.0);
                transformPixel(matrixH, pixel);
                int iPrime = (int)round(pixel.x / pixel.z);
                int jPrime = (int)round(pixel.y / pixel.z);
                jPrime = input.sizeY() - 1 - jPrime;
                if (iPrime < 0) {
                    iPrime = 0;
                }
                if (jPrime < 0) {
                    jPrime = 0;
                }
                if (iPrime > (input.sizeX() - 1)) {
                    iPrime = (input.sizeX() - 1);
                }
                if (jPrime > (input.sizeY() - 1)) {
                    jPrime = (input.sizeY() - 1);
                }
                double value = input.value(iPrime, jPrime);
                //imageInterp(input, iPrime, jPrime, value);
                output.definePixel(i, j, value);
            }
        }
}

void edgeFind1D(int mode, double threshold, std::vector<double>& input) {
    int N = (int)input.size();

    std::vector<double> edge;
    edge.resize(N, 0.0);
    #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            if (i == 0 || i == N - 1) {
                edge[i] = 0;
            }
            else {
                edge[i] = input[i - 1] - 2 * input[i] + input[i + 1];
            }
        }
    medianFilter(edge, 1);
    /*
    if (mode == 1) {
        plt::plot(edge);
        plt::show();
    }
    */
    #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            if (edge[i] > threshold) {
                edge[i] = 1;
            }
            else {
                edge[i] = 0;
            }
        }
    std::vector<double> dedge;
    dedge.resize(N, 0.0);
    #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            dedge[i] = abs(edge[i + 1] - edge[i]);
        }
    #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            if (dedge[i] < 0.1) {
                edge[i] = edge[i];
            }
            else {
                edge[i] = 0;
            }
        }
    input = edge;
}

void findRuler(int& mode, int& screen, imageBW& image, std::vector<double>& threshold, std::vector<cv::Point2d>& screenQuad, std::vector<double>& windowSize, std::vector<double>& rulerX, std::vector<double>& rulerY) {
    std::vector<double> line;
    imageBW imCrop;
    std::vector<int> bounds;
    bounds.resize(4, 0);

    //Horizontal Axis
    if (screen == 0) {
        bounds[0] = (int)round(screenQuad[4].x - windowSize[0] / 2 + 10);
        bounds[1] = (int)round(screenQuad[4].x + windowSize[0] / 2 - 10);
        bounds[2] = (int)screenQuad[5].y - 10;
        bounds[3] = (int)screenQuad[5].y + 10;
    }
    else {
        bounds[0] = (int)screenQuad[4].x + 10;
        bounds[1] = (int)screenQuad[6].x - 10;
        bounds[2] = (int)screenQuad[5].y - 20;
        bounds[3] = (int)screenQuad[5].y - 10;
    }
    image.crop(bounds, imCrop);
    lineOut(1, 1.0, imCrop, 1, line);
    medianFilter(line, 1);
    Contrast(threshold[0], line);
    if (mode == 1) {
        plt::figure;
        plt::title("Horizontal Threshold");
    }
    edgeFind1D(mode, threshold[1], line);



    int Nim = image.sizeX();
    int N = (int)line.size();
    std::vector<int> peak;
    for (int i = 0; i < N; i++) {
        if (line[i] > 0.9) {
            peak.push_back(i);
        }
    }
    line.clear();

    double averageSpacing = 0.0;
    rulerX.clear();
    std::vector<double> spacing;
    double value;
    N = (int)peak.size();
    for (int i = 0; i < N - 1; i++) {
        value = (double)peak[i + 1] - (double)peak[i];
        if (value >= 2) {
            if (value == 1 || value == 2) {
                if (rulerX.size() > 0) {
                    rulerX.back() = rulerX.back() + 1;
                    spacing.back() = spacing.back() + 1;
                    averageSpacing = averageSpacing + 1;
                }
            }
            else {
                rulerX.push_back(peak[i + 1] + bounds[0]);
                spacing.push_back(value);
                averageSpacing = averageSpacing + value;
            }
        }
    }
    N = (int)spacing.size();
    averageSpacing = averageSpacing / N;
    int loop = 1;
    while (loop > 0) {
        for (int i = 0; i < N; i++) {
            if (spacing[i] < 0.75 * averageSpacing) {
                rulerX.erase(rulerX.begin() + i + 1);
                break;
            }
            if (i == N - 1) {
                loop = 0;
            }
        }
        N = (int)rulerX.size();
        spacing.clear();
        averageSpacing = 0.0;
        for (int i = 0; i < N - 1; i++) {
            value = (double)rulerX[i + 1] - (double)rulerX[i];
            if (value > 2) {
                spacing.push_back(value);
                averageSpacing = averageSpacing + value;
            }
        }
        N = (int)spacing.size();
        averageSpacing = averageSpacing / N;
        if (loop > 64) {
            loop = 0;
        }
        else {
            if (loop == 0) {
                loop = 0;
            }
            else {
                loop = loop + 1;
            }
        }
    }
    loop = 1;
    while (loop > 0) {
        for (int i = 0; i < N; i++) {
            if (spacing[i] > 1.75 * averageSpacing) {
                int Nfill = round(spacing[i] / averageSpacing) - 1;
                for (int j = 0; j < Nfill; j++) {
                    rulerX.insert(rulerX.begin() + i + j + 1, rulerX[i + j] + round(averageSpacing));
                }
                break;
            }
            if (i == N - 1) {
                loop = 0;
            }
        }
        N = (int)rulerX.size();
        spacing.clear();
        averageSpacing = 0.0;
        for (int i = 0; i < N - 1; i++) {
            value = (double)rulerX[i + 1] - (double)rulerX[i];
            if (value > 2) {
                spacing.push_back(value);
                averageSpacing = averageSpacing + value;
            }
        }
        N = (int)spacing.size();
        averageSpacing = averageSpacing / N;
        if (loop > 64) {
            loop = 0;
        }
        else {
            if (loop == 0) {
                loop = 0;
            }
            else {
                loop = loop + 1;
            }
        }
    }
    double stdev, dx, dy;
    N = (int)rulerX.size();
    stdev = 0.0;
    for (int i = 0; i < N - 1; i++) {
        stdev = stdev + (rulerX[i + 1] - rulerX[i] - averageSpacing) * (rulerX[i + 1] - rulerX[i] - averageSpacing);
    }
    stdev = sqrt(stdev/(N-1));
    for (int i = 0; i < N - 1; i++) {
        dx = abs(rulerX[i + 1] - rulerX[i] - averageSpacing);
        if (dx > stdev) {
            rulerX[i + 1] = rulerX[i] + averageSpacing;
        }
    }
    /*
    if (screen == 0) {
        for (int i = 0; i < N; i++) {
            rulerX[i] = rulerX[i] - averageSpacing / 2.0;
        }
    }
    */
    loop = 1;
    while (loop) {
        value = round(rulerX[0] - averageSpacing);
        if (value > 0) {
            rulerX.insert(rulerX.begin(), value);
        }
        else {
            loop = 0;
        }
    }
    loop = 1;
    while (loop) {
        value = round(rulerX.back() + averageSpacing);
        if (value < Nim) {
            rulerX.push_back(value);
        }
        else {
            loop = 0;
        }
    }
    peak.clear();
    spacing.clear();


    //Vertical Axis
    if (screen == 0) {
        bounds[0] = (int)screenQuad[4].x - 10;
        bounds[1] = (int)screenQuad[4].x + 10;
        bounds[2] = (int)round(screenQuad[5].y - windowSize[1] / 2 + 10);
        bounds[3] = (int)round(screenQuad[5].y + windowSize[1] / 2 - 10);
    }
    else {
        bounds[0] = (int)screenQuad[4].x - 20;
        bounds[1] = (int)screenQuad[4].x - 10;
        bounds[2] = (int)screenQuad[4].y + 10;
        bounds[3] = (int)screenQuad[6].y - 10;
    }
    image.crop(bounds, imCrop);
    lineOut(1, 1.0, imCrop, 0, line);
    medianFilter(line, 1);
    Contrast(threshold[2], line);
    if (mode == 1) {
        plt::figure;
        plt::title("Vertical Threshold");
    }
    edgeFind1D(mode, threshold[3], line);


    Nim = image.sizeY();
    N = (int)line.size();
    peak.clear();
    for (int i = 0; i < N; i++) {
        if (line[i] > 0.9) {
            peak.push_back(i);
        }
    }
    line.clear();

    averageSpacing = 0.0;
    rulerY.clear();
    spacing.clear();
    value;
    N = (int)peak.size();
    for (int i = 0; i < N - 1; i++) {
        value = (double)peak[i + 1] - (double)peak[i];
        if (value >= 2) {
            if (value == 1 || value == 2) {
                if (rulerY.size() > 0) {
                    rulerY.back() = rulerY.back() + 1;
                    spacing.back() = spacing.back() + 1;
                    averageSpacing = averageSpacing + 1;
                }
            }
            else {
                rulerY.push_back(peak[i + 1] + bounds[2]);
                spacing.push_back(value);
                averageSpacing = averageSpacing + value;
            }
        }
    }
    N = (int)spacing.size();
    averageSpacing = averageSpacing / N;
    loop = 1;
    while (loop > 0) {
        for (int i = 0; i < N; i++) {
            if (spacing[i] < 0.75 * averageSpacing) {
                rulerY.erase(rulerY.begin() + i + 1);
                break;
            }
            if (i == N - 1) {
                loop = 0;
            }
        }
        N = (int)rulerY.size();
        spacing.clear();
        averageSpacing = 0.0;
        for (int i = 0; i < N - 1; i++) {
            value = (double)rulerY[i + 1] - (double)rulerY[i];
            if (value > 2) {
                spacing.push_back(value);
                averageSpacing = averageSpacing + value;
            }
        }
        N = (int)spacing.size();
        averageSpacing = averageSpacing / N;
        if (loop > 64) {
            loop = 0;
        }
        else {
            if (loop == 0) {
                loop = 0;
            }
            else {
                loop = loop + 1;
            }
        }
    }
    loop = 1;
    while (loop > 0) {
        for (int i = 0; i < N; i++) {
            if (spacing[i] > 1.75 * averageSpacing) {
                int Nfill = round(spacing[i] / averageSpacing) - 1;
                for (int j = 0; j < Nfill; j++) {
                    rulerY.insert(rulerY.begin() + i + 1 + j, rulerY[i + j] + round(averageSpacing));
                }
                break;
            }
            if (i == N - 1) {
                loop = 0;
            }
        }
        N = (int)rulerY.size();
        spacing.clear();
        averageSpacing = 0.0;
        for (int i = 0; i < N - 1; i++) {
            value = (double)rulerY[i + 1] - (double)rulerY[i];
            if (value > 2) {
                spacing.push_back(value);
                averageSpacing = averageSpacing + value;
            }
        }
        N = (int)spacing.size();
        averageSpacing = averageSpacing / N;
        if (loop > 64) {
            loop = 0;
        }
        else {
            if (loop == 0) {
                loop = 0;
            }
            else {
                loop = loop + 1;
            }
        }
    }
    N = (int)rulerY.size();
    stdev = 0.0;
    for (int i = 0; i < N - 1; i++) {
        stdev = stdev + (rulerY[i + 1] - rulerY[i] - averageSpacing) * (rulerY[i + 1] - rulerY[i] - averageSpacing);
    }
    stdev = sqrt(stdev / (N-1));
    for (int i = 0; i < N - 1; i++) {
        dy = abs(rulerY[i + 1] - rulerY[i] - averageSpacing);
        if (dy > stdev) {
            rulerY[i + 1] = rulerY[i] + averageSpacing;
        }
    }
    /*
    if (screen == 0) {
        for (int i = 0; i < N; i++) {
            rulerY[i] = rulerY[i] - averageSpacing / 2.0;
        }
    }
    */
    loop = 1;
    while (loop) {
        value = round(rulerY[0] - averageSpacing);
        if (value > 0) {
            rulerY.insert(rulerY.begin(), value);
        }
        else {
            loop = 0;
        }
    }
    loop = 1;
    while (loop) {
        value = round(rulerY.back() + averageSpacing);
        if (value < Nim) {
            rulerY.push_back(value);
        }
        else {
            loop = 0;
        }
    }
    peak.clear();
    spacing.clear();
}

void findZero(int& screen, imageBW& image, std::vector<double>& rulerX, std::vector<double>& rulerY, std::vector<double>& zeroPoint) {
    double Nx = (double)image.sizeX();
    double Ny = (double)image.sizeY();
    int nRx = (int)rulerX.size();
    int nRy = (int)rulerY.size();


    double dv, value;
    int pxYCenter = 0;
    int pxXCenter = 0;
    dv = Ny;
    for (int i = 0; i < nRy; i++) {
        value = abs(Ny / 2.0 - rulerY[i]);
        if (value < dv) {
            dv = value;
            pxYCenter = rulerY[i];
        }
    }

    if (screen == 0) {
        dv = Nx;
        for (int i = 0; i < nRx; i++) {
            value = abs(Nx / 2.0 - rulerX[i]);
            if (value < dv) {
                dv = value;
                pxXCenter = rulerX[i];
            }
        }
    }
    else {
        imageBW imCrop;
        std::vector<int> bounds;
        std::vector<double> line;
        bounds.resize(4, 0.0);
        bounds[0] = (int)round(3.0 / 4.0 * Nx);
        bounds[1] = Nx - 1;
        bounds[2] = (int)round(1.0 / 1.1 * (double)pxYCenter);
        bounds[3] = (int)round(1.1 * (double)pxYCenter);
        image.crop(bounds, imCrop);
        lineOut(1, 1.0, imCrop, 1, line);
        medianFilter(line, 10);
        edgeFind1D(0, 0.005, line);
        for (int i = 0; i < (int)line.size(); i++) {
            if (line[i] > 0.9) {
                pxXCenter = i + bounds[0];
            }
        }
    }
    zeroPoint.resize(2, 0.0);
    zeroPoint[0] = pxXCenter;
    zeroPoint[1] = pxYCenter;
}

void pixelAxis(int& screen, int Nx, int Ny, std::vector<double>& rulerX, std::vector<double>& rulerY, std::vector<double>& zeroPoint) {
    int nRx = (int)rulerX.size();
    int nRy = (int)rulerY.size();

    int n = std::max(nRx, nRy);
    std::vector<double> mmX, mmY;
    mmX.resize(nRx, 0.0);
    mmY.resize(nRy, 0.0);
    if (screen == 0) {
        for (int i = 0; i < n; i++) {
            if (i < nRx) {
                mmX[i] = i;
            }
            if (i < nRy) {
                mmY[i] = i;
            }
        }
    }
    else {
        for (int i = 0; i < n; i++) {
            if (i < nRx) {
                mmX[i] = i;
            }
            if (i < nRy) {
                mmY[i] = i;
            }
        }
    }

    std::vector<double> resX, resY, pxX, pxY;
    resX = rulerX;
    resY = rulerY;
    rulerX.clear();
    rulerY.clear();
    int N = std::max(Nx, Ny);
    for (int i = 0; i < N; i++) {
        if (i < Nx) {
            pxX.push_back(i);
        }
        if (i < Ny) {
            pxY.push_back(i);
        }
    }

    double pxIndex, value, valueMin, valueMax, dv;

    
    bool loop = 1;
    valueMin = std::min(resX.front(), resX.back());
    valueMax = std::max(resX.front(), resX.back());
    rulerX.push_back(0.0);
    pxIndex = valueMin + 1.0;
    while (loop) {
        if (pxIndex >= valueMin && pxIndex <= valueMax) {
            FE1DInterp(resX, mmX, pxIndex, value);
            value = round(value * 1000.0) / 1000.0;
            rulerX.push_back(value);
            pxIndex = pxIndex + 1.0;
        }
        else {
            loop = 0;
        }
    }
    loop = 1;
    dv = rulerX[1] - rulerX[0];
    pxIndex = valueMin - 1.0;
    while (loop) {
        if (pxIndex >= 0 && pxIndex < Nx) {
            value = rulerX.front() - dv;
            rulerX.insert(rulerX.begin(), value);
            pxIndex = pxIndex - 1;
        }
        else {
            loop = 0;
        }
    }
    loop = 1;
    dv = rulerX[rulerX.size() - 2] - rulerX[rulerX.size() - 3];
    pxIndex = rulerX.size();
    while (loop) {
        if (pxIndex >= 0 && pxIndex < Nx) {
            value = rulerX.back() + dv;
            rulerX.push_back(value);
            pxIndex = pxIndex + 1;
        }
        else {
            loop = 0;
        }
    }
    

    loop = 1;
    valueMin = std::min(resY.front(), resY.back());
    valueMax = std::max(resY.front(), resY.back());
    rulerY.push_back(0.0);
    pxIndex = valueMin + 1.0;
    while (loop) {
        if (pxIndex >= valueMin && pxIndex <= valueMax) {
            FE1DInterp(resY, mmY, pxIndex, value);
            value = round(value * 1000.0) / 1000.0;
            rulerY.push_back(value);
            pxIndex = pxIndex + 1.0;
        }
        else {
            loop = 0;
        }
    }
    loop = 1;
    dv = rulerY[1] - rulerY[0];
    pxIndex = valueMin - 1.0;
    while (loop) {
        if (pxIndex >= 0 && pxIndex < Ny) {
            value = rulerY.front() - dv;
            rulerY.insert(rulerY.begin(), value);
            pxIndex = pxIndex - 1;
        }
        else {
            loop = 0;
        }
    }
    loop = 1;
    dv = rulerY[rulerY.size() - 2] - rulerY[rulerY.size() - 3];
    pxIndex = rulerY.size();
    while (loop) {
        if (pxIndex >= 0 && pxIndex < Ny) {
            value = rulerY.back() + dv;
            rulerY.push_back(value);
            pxIndex = pxIndex + 1;
        }
        else {
            loop = 0;
        }
    }

    linReg(pxX, rulerX);
    linReg(pxY, rulerY);

    double zPX, zPY;
    FE1DInterp(pxX, rulerX, zeroPoint[0], zPX);
    FE1DInterp(pxY, rulerY, zeroPoint[1], zPY);
    zPX = round(zPX * 1000.0) / 1000.0;
    zPY = round(zPY * 1000.0) / 1000.0;
    if (screen == 0) {
        for (int i = 0; i < N; i++) {
            if (i < Nx) {
                rulerX[i] = (zPX - rulerX[i]);
            }
            if (i < Ny) {
                rulerY[i] = (zPY - rulerY[i]);
            }
        }
    }
    else {
        for (int i = 0; i < N; i++) {
            if (i < Nx) {
                rulerX[i] = (zPX - rulerX[i]);
            }
            if (i < Ny) {
                rulerY[i] = (zPY - rulerY[i]);
            }
        }
    }
}

void screenCal(int mode, int& screen, spectrometer& eSpec, screenCalibration& calibration, cv::Mat& transformMatH, std::vector<double>& rulerX, std::vector<double>& rulerY) {
    std::string screenPath = eSpec.screenPath(screen);
    std::vector<cv::Point2d> viewQuad = calibration.viewQuad(screen);
    std::vector<double> viewRes = calibration.viewResolution(screen);
    std::vector<double> viewThreshold = calibration.viewThreshold(screen);
    std::vector<std::string> folderList; 
    listDir(screenPath, folderList);
    std::string refFile;
    findRef(folderList, refFile);
    imageBW image, imTransform;

    homographyMat(viewQuad, transformMatH);

    getImage(refFile, image);
    perspectiveTransform(image, transformMatH, viewRes, imTransform);
    findRuler(mode, screen, imTransform, viewThreshold, viewQuad, viewRes, rulerX, rulerY);

    if (mode == 1) {
        imageBW imPlot;
        std::vector<int> bounds;
        bounds.resize(4, 0);
        if (screen == 0) {
            bounds[0] = (viewQuad[0].x + viewQuad[2].x - viewRes[0]) / 2;
            bounds[1] = (viewQuad[0].x + viewQuad[2].x + viewRes[0]) / 2;
            bounds[2] = (viewQuad[1].y + viewQuad[1].y - viewRes[1]) / 2;
            bounds[3] = (viewQuad[1].y + viewQuad[1].y + viewRes[1]) / 2;

            image.crop(bounds, imPlot);

            size_t resV = 720;
            double ratio = 2.0 * (double)imPlot.sizeX() / (double)imPlot.sizeY();
            ratio = ratio * (resV + resV / 100.0);
            size_t resH = (size_t)round(ratio);

            std::vector<double> plotX, plotY;
            plotX.resize(5, 0.0);
            plotY.resize(5, 0.0);

            plt::figure_size(resH, resV);
            plt::subplot(1, 2, 1);

            plotX[0] = viewQuad[0].x - bounds[0];
            plotX[1] = viewQuad[1].x - bounds[0];
            plotX[2] = viewQuad[2].x - bounds[0];
            plotX[3] = viewQuad[3].x - bounds[0];
            plotX[4] = viewQuad[0].x - bounds[0];
            plotY[0] = image.sizeY() - 1 - viewQuad[0].y - bounds[2];
            plotY[1] = image.sizeY() - 1 - viewQuad[1].y - bounds[2];
            plotY[2] = image.sizeY() - 1 - viewQuad[2].y - bounds[2];
            plotY[3] = image.sizeY() - 1 - viewQuad[3].y - bounds[2];
            plotY[4] = image.sizeY() - 1 - viewQuad[0].y - bounds[2];

            pltimshow(imPlot, 0, "");
            plt::plot(plotX, plotY, { {"color","w"} });
            plt::axis("off");

            int Nx = (int)rulerX.size();
            int Ny = (int)rulerY.size();
            plotX.resize(2, 0.0);
            plotY.resize(2, 0.0);
            plt::subplot(1, 2, 2);
            pltimshow(imTransform, 0, "");
            for (int i = 0; i < Nx; i++) {
                //if (i % 10 == 2) {
                plotX[0] = rulerX[i];
                plotX[1] = plotX[0];
                plotY[0] = viewQuad[5].y - 10;
                plotY[1] = viewQuad[5].y + 10;
                plt::plot(plotX, plotY, { {"color","w"} });
                //}
            }
            for (int i = 0; i < Ny; i++) {
                //if (i % 10 == 2) {
                plotY[0] = rulerY[i];
                plotY[1] = plotY[0];
                plotX[0] = viewQuad[4].x - 10;
                plotX[1] = viewQuad[4].x + 10;
                plt::plot(plotX, plotY, { {"color","w"} });
                //}
            }
            plt::axis("off");
            plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.0}, {"hspace",0.0} });
            plt::show();
        }
        else {
            bounds[0] = (int)std::max(std::min(std::min(viewQuad[0].x, viewQuad[1].x), std::min(viewQuad[2].x, viewQuad[3].x)) - 100, 0.0);
            bounds[1] = (int)std::min(std::max(std::max(viewQuad[0].x, viewQuad[1].x), std::max(viewQuad[2].x, viewQuad[3].x)) + 300, (double)(image.sizeX() - 1));
            bounds[2] = (int)std::max(std::min(std::min(viewQuad[0].y, viewQuad[1].y), std::min(viewQuad[2].y, viewQuad[3].y)) - 100, 0.0);
            bounds[3] = (int)std::min(std::max(std::max(viewQuad[0].y, viewQuad[1].y), std::max(viewQuad[2].y, viewQuad[3].y)) + 100, (double)(image.sizeY() - 1));

            image.crop(bounds, imPlot);

            size_t resV = 720;
            double ratio = 0.5 * (double)imPlot.sizeX() / (double)imPlot.sizeY();
            ratio = ratio * (resV + resV / 100.0);
            size_t resH = (size_t)round(ratio);

            std::vector<double> plotX, plotY;
            plotX.resize(5, 0.0);
            plotY.resize(5, 0.0);

            plt::figure_size(resH, resV);
            plt::subplot(2, 1, 1);

            plotX[0] = viewQuad[0].x - bounds[0];
            plotX[1] = viewQuad[1].x - bounds[0];
            plotX[2] = viewQuad[2].x - bounds[0];
            plotX[3] = viewQuad[3].x - bounds[0];
            plotX[4] = viewQuad[0].x - bounds[0];
            plotY[0] = image.sizeY() - 1 - viewQuad[0].y - bounds[2];
            plotY[1] = image.sizeY() - 1 - viewQuad[1].y - bounds[2];
            plotY[2] = image.sizeY() - 1 - viewQuad[2].y - bounds[2];
            plotY[3] = image.sizeY() - 1 - viewQuad[3].y - bounds[2];
            plotY[4] = image.sizeY() - 1 - viewQuad[0].y - bounds[2];

            pltimshow(imPlot, 0, "");
            plt::plot(plotX, plotY, { {"color","w"} });
            plt::axis("off");

            int Nx = (int)rulerX.size();
            int Ny = (int)rulerY.size();
            plotX.resize(2, 0.0);
            plotY.resize(2, 0.0);
            plt::subplot(2, 1, 2);
            pltimshow(imTransform, 0, "");
            for (int i = 0; i < Nx; i++) {
                //if (i % 10 == 2) {
                plotX[0] = rulerX[i];
                plotX[1] = plotX[0];
                plotY[0] = viewQuad[5].y - 10;
                plotY[1] = viewQuad[5].y + 10;
                plt::plot(plotX, plotY, { {"color","w"} });
                //}
            }
            for (int i = 0; i < Ny; i++) {
                //if (i % 10 == 2) {
                plotY[0] = rulerY[i];
                plotY[1] = plotY[0];
                plotX[0] = viewQuad[4].x - 10;
                plotX[1] = viewQuad[4].x + 10;
                plt::plot(plotX, plotY, { {"color","w"} });
                //}
            }
            plt::axis("off");
            plt::subplots_adjust({ {"left",0.05},{"right",0.95},{"top", 0.95},{"bottom",0.04}, {"wspace", 0.0}, {"hspace",0.0} });
            plt::show();
        }
    }
}

void writeCalibration(std::string calPath, std::vector<std::vector<int>>& winRes, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler, std::vector<std::vector<double>>& zP) {
    std::fstream output;
    std::string filePath = calPath + "/perspective.cache";
    output.open(filePath, std::fstream::out);
    int N = (int)H.size();
    for (int n = 0; n < N; n++) {
        output << winRes[n][0] << ",\t" << winRes[n][1] << ",\t" << zP[n][0] << ",\t" << zP[n][1] << "\n";
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i == 2 && j == 2) {
                    output << H[n].at<double>(i, j);
                }
                else {
                    output << H[n].at<double>(i, j) << ",\t";
                }
            }
        }
        output << "\n";
        int Nx = (int)xRuler[n].size();
        for (int i = 0; i < Nx; i++) {
            if (i == Nx - 1) {
                output << xRuler[n][i];
            }
            else {
                output << xRuler[n][i] << ",\t";
            }
        }
        output << "\n";
        int Ny = (int)yRuler[n].size();
        for (int i = 0; i < Ny; i++) {
            if (i == Ny - 1) {
                output << yRuler[n][i];
            }
            else {
                output << yRuler[n][i] << ",\t";
            }
        }
        output << "\n";
    }
    output.close();
}

void readCalibration(std::string calPath, std::vector<cv::Mat>& H, std::vector<std::vector<double>>& xRuler, std::vector<std::vector<double>>& yRuler) {
    std::string filePath = calPath + "/perspective.cache";

    std::vector<std::vector<int>> winRes;
    std::vector<std::vector<double>> zP;
    winRes.clear();
    winRes.resize(3);
    H.clear();
    H.resize(3);
    zP.clear();
    zP.resize(3);
    xRuler.clear();
    xRuler.resize(3);
    yRuler.clear();
    yRuler.resize(3);

    size_t index;
    std::string strValue;
    std::vector<std::string> buffer;
    readFile(filePath, buffer);
    for (int n = 0; n < 3; n++) {
        std::vector<int> bufferWR(2, 0);
        std::vector<double> bufferZP(2, 0.0);
        index = buffer[4 * n].find(",");
        strValue = buffer[4 * n].substr(0, index);
        bufferWR[0] = (int)std::stod(strValue.c_str());
        index = index + 1;
        buffer[4 * n] = buffer[4 * n].substr(index, buffer[4 * n].length());
        index = buffer[4 * n].find(",");
        strValue = buffer[4 * n].substr(0, index);
        bufferWR[1] = (int)std::stod(strValue.c_str());
        index = index + 1;
        buffer[4 * n] = buffer[4 * n].substr(index, buffer[4 * n].length());
        index = buffer[4 * n].find(",");
        strValue = buffer[4 * n].substr(0, index);
        bufferZP[0] = std::stod(strValue.c_str());
        index = index + 1;
        buffer[4 * n] = buffer[4 * n].substr(index, buffer[4 * n].length());
        index = buffer[4 * n].find(",");
        strValue = buffer[4 * n].substr(0, index);
        bufferZP[1] = std::stod(strValue.c_str());

        winRes[n] = bufferWR;
        zP[n] = bufferZP;

        cv::Mat bufferH = cv::Mat::zeros(3, 3, CV_64F);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                index = buffer[4 * n + 1].find(",");
                strValue = buffer[4 * n + 1].substr(0, index);
                bufferH.at<CvType<CV_64F>::type_t>(i, j) = std::stod(strValue.c_str());
                index = index + 1;
                buffer[4 * n + 1] = buffer[4 * n + 1].substr(index, buffer[4 * n + 1].length());
            }
        }
        H[n] = bufferH;

        bool loop = 1;
        std::vector<double> ruler;
        while (loop) {
            if (buffer[4 * n + 2].find(",") != std::string::npos) {
                index = buffer[4 * n + 2].find(",");
                strValue = buffer[4 * n + 2].substr(0, index);
                ruler.push_back(std::stod(strValue.c_str()));
                index = index + 1;
                buffer[4 * n + 2] = buffer[4 * n + 2].substr(index, buffer[4 * n + 2].length());
            }
            else {
                ruler.push_back(std::stod(buffer[4 * n + 2].c_str()));
                loop = 0;
            }
        }
        xRuler[n] = ruler;

        loop = 1;
        ruler.clear();
        while (loop) {
            if (buffer[4 * n + 3].find(",") != std::string::npos) {
                index = buffer[4 * n + 3].find(",");
                strValue = buffer[4 * n + 3].substr(0, index);
                ruler.push_back(std::stod(strValue.c_str()));
                index = index + 1;
                buffer[4 * n + 3] = buffer[4 * n + 3].substr(index, buffer[4 * n + 3].length());
            }
            else {
                ruler.push_back(std::stod(buffer[4 * n + 3].c_str()));
                loop = 0;
            }
        }
        yRuler[n] = ruler;
    }
    int screen = 0;
    pixelAxis(screen, winRes[0][0], winRes[0][1], xRuler[0], yRuler[0], zP[0]);
    screen = 1;
    pixelAxis(screen, winRes[1][0], winRes[1][1], xRuler[1], yRuler[1], zP[1]);
    screen = 2;
    pixelAxis(screen, winRes[2][0], winRes[2][1], xRuler[2], yRuler[2], zP[2]);
}
#endif