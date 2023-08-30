#pragma once

#ifndef __NumericalMethods_h__
#define __NumericalMethods_h__

#include <vector>
#include "fileIO.h"

//Composite Simpson's 3/8 Rule
void simpsonInt(double& dx, std::vector<double>& y, double& output) {
    output = 0.0;
    size_t N = y.size()-1;
    int remainder = N % 3;
    size_t Nend = floor(N / 3);
    if (remainder == 0) {
        Nend++;
    }
    else if (remainder == 2) {
        Nend++;
    }
    for (int i = 1; i < Nend; i++) {
        output = output + y[3 * i - 3] + 3.0 * y[3 * i - 2] + 3.0 * y[3 * i - 1] + y[3 * i];
    }
    output = output * dx * 3.0 / 8.0;
    if (remainder > 0) {
        output = output + dx / 3.0 * (y[3 * (Nend - 1)] + 4.0 * y[3 * (Nend - 1) + 1] + y[3 * (Nend - 1) + 2]);
        if (remainder == 1) {
            output = output + dx / 3.0 * (y[3 * (Nend - 1) + 2] + 4.0 * y[3 * (Nend - 1) + 3] + y[3 * (Nend - 1) + 4]);
        }
    }
}

//Finite Element Interpolation
void FE1DInterp(std::vector<double> &x, std::vector<double> &y, double &xEval, double &output) {
    output = 0.0;
    int Lx = (int)x.size();
    int NxE = (int)floor((Lx - 2) / 2) + 1;
    int eleX;
    double mX0, mX;

    std::vector<double> psi;
    psi.resize(3, 0.0);
    if (abs(2 * NxE - Lx) > pow(10.0, -3.0)) {
        eleX = NxE;
        if (x.front() < x.back()) {
            for (int i = 0; i < NxE; i++) {
                if (xEval <= x[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NxE; i++) {
                if (xEval >= x[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }

        if (eleX < NxE) {
            mX0 = (2 * x[2 * eleX + 1] - (x[2 * eleX + 2] + x[2 * eleX])) / (x[2 * eleX + 2] - x[2 * eleX]);
            mX = (2 * xEval - (x[2 * eleX + 2] + x[2 * eleX])) / (x[2 * eleX + 2] - x[2 * eleX]);

            psi[0] = (mX - mX0) * (mX - 1.0) / (2 * (1 + mX0));
            psi[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
            psi[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
        }
        else {
            mX = (2 * xEval - (x[2 * eleX + 1] + x[2 * eleX])) / (x[2 * eleX + 1] - x[2 * eleX]);

            psi[0] = -(mX - 1.0) / 2.0;
            psi[1] = (mX + 1.0) / 2.0;
        }
    }
    else {
        eleX = NxE - 1;
        for (int i = 0; i < NxE; i++) {
            if (abs(xEval - x[2 * (i + 1)]) <= abs(x[2 * (i + 1)] - x[2 * i])) {
                eleX = i;
                break;
            }
        }
        mX0 = (2 * x[2 * eleX + 1] - (x[2 * eleX + 2] + x[2 * eleX])) / (x[2 * eleX + 2] - x[2 * eleX]);
        mX = (2 * xEval - (x[2 * eleX + 2] + x[2 * eleX])) / (x[2 * eleX + 2] - x[2 * eleX]);

        psi[0] = (mX - mX0) * (mX - 1.0) / (2 * (1.0 + mX0));
        psi[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
        psi[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
    }
    if (eleX > NxE - 1) {
        output = y[2 * eleX] * psi[0] + y[2 * eleX + 1] * psi[1];
    }
    else{
        output = y[2 * eleX] * psi[0] + y[2 * eleX + 1] * psi[1] + y[2 * eleX + 2] * psi[2];
    }
}

void FE2DInterp(std::vector<double>& xAxis, std::vector<double>& yAxis, std::vector<std::vector<double>>& surface, double& xEval, double& yEval, double& output) {
    output = 0.0;
    int Lx = (int)surface.size();
    int Ly = (int)surface[0].size();
    int NxE = (int)floor((Lx - 2) / 2) + 1;
    int NyE = (int)floor((Ly - 2) / 2) + 1;
    int eleX, eleY;
    double mX0, mX, mY0, mY;

    std::vector<double> psiX, psiY;
    psiX.resize(3, 0.0);
    if (abs(2 * NxE - Lx) > pow(10.0, -3.0)) {
        eleX = NxE;
        if (xAxis.front() < xAxis.back()) {
            for (int i = 0; i < NxE; i++) {
                if (xEval <= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NxE; i++) {
                if (xEval >= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }

        if (eleX < NxE) {
            mX0 = (2 * xAxis[2 * eleX + 1] - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);
            mX = (2 * xEval - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);

            psiX[0] = (mX - mX0) * (mX - 1.0) / (2 * (1 + mX0));
            psiX[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
            psiX[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
        }
        else {
            mX = (2 * xEval - (xAxis[2 * eleX + 1] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 1] - xAxis[2 * eleX]);

            psiX[0] = -(mX - 1.0) / 2.0;
            psiX[1] = (mX + 1.0) / 2.0;
        }
    }
    else {
        eleX = NxE - 1;
        if (xAxis.front() < xAxis.back()) {
            for (int i = 0; i < NxE; i++) {
                if (xEval <= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NxE; i++) {
                if (xEval >= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }

        mX0 = (2 * xAxis[2 * eleX + 1] - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);
        mX = (2 * xEval - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);

        psiX[0] = (mX - mX0) * (mX - 1.0) / (2 * (1.0 + mX0));
        psiX[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
        psiX[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
    }

    psiY.resize(3, 0.0);
    if (abs(2 * NyE - Ly) > pow(10.0, -3.0)) {
        eleY = NyE;
        if (yAxis.front() < yAxis.back()) {
            for (int i = 0; i < NyE; i++) {
                if (yEval <= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NyE; i++) {
                if (yEval >= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }

        if (eleY < NyE) {
            mY0 = (2 * yAxis[2 * eleY + 1] - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);
            mY = (2 * yEval - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);

            psiY[0] = (mY - mY0) * (mY - 1.0) / (2 * (1 + mY0));
            psiY[1] = (mY + 1.0) * (mY - 1.0) / ((mY0 + 1.0) * (mY0 - 1.0));
            psiY[2] = (mY + 1.0) * (mY - mY0) / (2 * (1.0 - mY0));
        }
        else {
            mY = (2 * yEval - (yAxis[2 * eleY + 1] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 1] - yAxis[2 * eleY]);

            psiY[0] = -(mY - 1.0) / 2.0;
            psiY[1] = (mY + 1.0) / 2.0;
        }
    }
    else {
        eleY = NyE - 1;
        if (yAxis.front() < yAxis.back()) {
            for (int i = 0; i < NyE; i++) {
                if (yEval <= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NyE; i++) {
                if (yEval >= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }

        mY0 = (2 * yAxis[2 * eleY + 1] - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);
        mY = (2 * yEval - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);

        psiY[0] = (mY - mY0) * (mY - 1.0) / (2 * (1.0 + mY0));
        psiY[1] = (mY + 1.0) * (mY - 1.0) / ((mY0 + 1.0) * (mY0 - 1.0));
        psiY[2] = (mY + 1.0) * (mY - mY0) / (2 * (1.0 - mY0));
    }

    if (eleX > NxE - 1) {
        if (eleY > NyE - 1) {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1];
        }
        else {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1] \
                + surface[2 * eleX][2 * eleY + 2] * psiX[0] * psiY[2] \
                + surface[2 * eleX + 1][2 * eleY + 2] * psiX[1] * psiY[2];
        }
    }
    else {
        if (eleY > NyE - 1) {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX + 2][2 * eleY] * psiX[2] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1] \
                + surface[2 * eleX + 2][2 * eleY + 1] * psiX[2] * psiY[1];
        }
        else {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX + 2][2 * eleY] * psiX[2] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1] \
                + surface[2 * eleX + 2][2 * eleY + 1] * psiX[2] * psiY[1] \
                + surface[2 * eleX][2 * eleY + 2] * psiX[0] * psiY[2] \
                + surface[2 * eleX + 1][2 * eleY + 2] * psiX[1] * psiY[2] \
                + surface[2 * eleX + 2][2 * eleY + 2] * psiX[2] * psiY[2];
        }
    }
}

void imageInterp(imageBW & image, double& xEval, double& yEval, double& output) {
    output = 0.0;
    int Lx = (int)image.sizeX();
    int Ly = (int)image.sizeY();
    int NxE = (int)floor((Lx - 2) / 2) + 1;
    int NyE = (int)floor((Ly - 2) / 2) + 1;
    int eleX, eleY;
    double mX0, mX, mY0, mY;

    std::vector<double> psiX, psiY;
    psiX.resize(3, 0.0);
    if (abs(2 * NxE - Lx) > pow(10.0, -3.0)) {
        eleX = NxE;
        for (int i = 0; i < NxE; i++) {
            if (xEval <= (double)(2 * i + 2)) {
                eleX = i;
                break;
            }
        }

        if (eleX < NxE - 1) {
            mX0 = (2 * (double)(2 * eleX + 1) - ((double)(2 * eleX + 2) + (double)(2 * eleX))) / ((double)(2 * eleX + 2) - (double)(2 * eleX));
            mX = (2 * xEval - ((double)(2 * eleX + 2) + (double)(2 * eleX))) / ((double)(2 * eleX + 2) - (double)(2 * eleX));

            psiX[0] = (mX - mX0) * (mX - 1.0) / (2 * (1 + mX0));
            psiX[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
            psiX[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
        }
        else {
            mX = (2 * xEval - ((double)(2 * eleX + 1) + (double)(2 * eleX))) / ((double)(2 * eleX + 1) - (double)(2 * eleX));

            psiX[0] = -(mX - 1.0) / 2.0;
            psiX[1] = (mX + 1.0) / 2.0;
        }
    }
    else {
        eleX = NxE - 1;
        for (int i = 0; i < NxE; i++) {
            if (xEval <= (double)(2 * i + 2)) {
                eleX = i;
                break;
            }
        }

        mX0 = (2 * (double)(2 * eleX + 1) - ((double)(2 * eleX + 2) + (double)(2 * eleX))) / ((double)(2 * eleX + 2) - (double)(2 * eleX));
        mX = (2 * xEval - ((double)(2 * eleX + 2) + (double)(2 * eleX))) / ((double)(2 * eleX + 2) - (double)(2 * eleX));

        psiX[0] = (mX - mX0) * (mX - 1.0) / (2 * (1.0 + mX0));
        psiX[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
        psiX[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
    }

    psiY.resize(3, 0.0);
    if (abs(2 * NyE - Ly) > pow(10.0, -3.0)) {
        eleY = NyE;
        for (int i = 0; i < NyE; i++) {
            if (yEval <= (double)(2 * i + 2)) {
                eleY = i;
                break;
            }
        }

        if (eleY < NyE - 1) {
            mY0 = (2 * (double)(2 * eleY + 1) - ((double)(2 * eleY + 2) + (double)(2 * eleY))) / ((double)(2 * eleY + 2) - (double)(2 * eleY));
            mY = (2 * yEval - ((double)(2 * eleY + 2) + (double)(2 * eleY))) / ((double)(2 * eleY + 2) - (double)(2 * eleY));

            psiY[0] = (mY - mY0) * (mY - 1.0) / (2 * (1 + mY0));
            psiY[1] = (mY + 1.0) * (mY - 1.0) / ((mY0 + 1.0) * (mY0 - 1.0));
            psiY[2] = (mY + 1.0) * (mY - mY0) / (2 * (1.0 - mY0));
        }
        else {
            mY = (2 * yEval - ((double)(2 * eleY + 1) + (double)(2 * eleY))) / ((double)(2 * eleY + 1) - (double)(2 * eleY));

            psiY[0] = -(mY - 1.0) / 2.0;
            psiY[1] = (mY + 1.0) / 2.0;
        }
    }
    else {
        eleY = NyE - 1;
        for (int i = 0; i < NyE; i++) {
            if (yEval <= (double)(2 * i + 2)) {
                eleY = i;
                break;
            }
        }

        mY0 = (2 * (double)(2 * eleY + 1) - ((double)(2 * eleY + 2) + (double)(2 * eleY))) / ((double)(2 * eleY + 2) - (double)(2 * eleY));
        mY = (2 * yEval - ((double)(2 * eleY + 2) + (double)(2 * eleY))) / ((double)(2 * eleY + 2) - (double)(2 * eleY));

        psiY[0] = (mY - mY0) * (mY - 1.0) / (2 * (1.0 + mY0));
        psiY[1] = (mY + 1.0) * (mY - 1.0) / ((mY0 + 1.0) * (mY0 - 1.0));
        psiY[2] = (mY + 1.0) * (mY - mY0) / (2 * (1.0 - mY0));
    }

    if (eleX > NxE - 1) {
        if (eleY > NyE - 1) {
            output = image.value(2 * eleX, 2 * eleY) * psiX[0] * psiY[0] \
                + image.value(2 * eleX + 1, 2 * eleY) * psiX[1] * psiY[0] \
                + image.value(2 * eleX, 2 * eleY + 1) * psiX[0] * psiY[1] \
                + image.value(2 * eleX + 1, 2 * eleY + 1) * psiX[1] * psiY[1];
        }
        else {
            output = image.value(2 * eleX, 2 * eleY) * psiX[0] * psiY[0] \
                + image.value(2 * eleX + 1, 2 * eleY) * psiX[1] * psiY[0] \
                + image.value(2 * eleX, 2 * eleY + 1) * psiX[0] * psiY[1] \
                + image.value(2 * eleX + 1, 2 * eleY + 1) * psiX[1] * psiY[1] \
                + image.value(2 * eleX, 2 * eleY + 2) * psiX[0] * psiY[2] \
                + image.value(2 * eleX + 1, 2 * eleY + 2) * psiX[1] * psiY[2];
        }
    }
    else {
        if (eleY > NyE - 1) {
            output = image.value(2 * eleX, 2 * eleY) * psiX[0] * psiY[0] \
                + image.value(2 * eleX + 1, 2 * eleY) * psiX[1] * psiY[0] \
                + image.value(2 * eleX + 2, 2 * eleY) * psiX[2] * psiY[0] \
                + image.value(2 * eleX, 2 * eleY + 1) * psiX[0] * psiY[1] \
                + image.value(2 * eleX + 1, 2 * eleY + 1) * psiX[1] * psiY[1] \
                + image.value(2 * eleX + 2, 2 * eleY + 1) * psiX[2] * psiY[1];
        }
        else {
            output = image.value(2 * eleX, 2 * eleY) * psiX[0] * psiY[0] \
                + image.value(2 * eleX + 1, 2 * eleY) * psiX[1] * psiY[0] \
                + image.value(2 * eleX + 2, 2 * eleY) * psiX[2] * psiY[0] \
                + image.value(2 * eleX, 2 * eleY + 1) * psiX[0] * psiY[1] \
                + image.value(2 * eleX + 1, 2 * eleY + 1) * psiX[1] * psiY[1] \
                + image.value(2 * eleX + 2, 2 * eleY + 1) * psiX[2] * psiY[1] \
                + image.value(2 * eleX, 2 * eleY + 2) * psiX[0] * psiY[2] \
                + image.value(2 * eleX + 1, 2 * eleY + 2) * psiX[1] * psiY[2] \
                + image.value(2 * eleX + 2, 2 * eleY + 2) * psiX[2] * psiY[2];
        }
    }
}

void imageInterp(std::vector<double>& xAxis, std::vector<double>& yAxis, std::vector<std::vector<double>>& surface, double& xEval, double& yEval, double& output) {
    output = 0.0;
    int Lx = (int)surface.size();
    int Ly = (int)surface[0].size();
    int NxE = (int)floor((Lx - 2) / 2) + 1;
    int NyE = (int)floor((Ly - 2) / 2) + 1;
    int eleX, eleY;
    double mX0, mX, mY0, mY;

    std::vector<double> psiX, psiY;
    psiX.resize(3, 0.0);
    if (abs(2 * NxE - Lx) > pow(10.0, -3.0)) {
        eleX = NxE;
        if (xAxis.front() < xAxis.back()) {
            for (int i = 0; i < NxE; i++) {
                if (xEval <= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NxE; i++) {
                if (xEval >= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }

        if (eleX < NxE) {
            mX0 = (2 * xAxis[2 * eleX + 1] - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);
            mX = (2 * xEval - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);

            psiX[0] = (mX - mX0) * (mX - 1.0) / (2 * (1 + mX0));
            psiX[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
            psiX[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
        }
        else {
            mX = (2 * xEval - (xAxis[2 * eleX + 1] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 1] - xAxis[2 * eleX]);

            psiX[0] = -(mX - 1.0) / 2.0;
            psiX[1] = (mX + 1.0) / 2.0;
        }
    }
    else {
        eleX = NxE - 1;
        if (xAxis.front() < xAxis.back()) {
            for (int i = 0; i < NxE; i++) {
                if (xEval <= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NxE; i++) {
                if (xEval >= xAxis[2 * i + 2]) {
                    eleX = i;
                    break;
                }
            }
        }

        mX0 = (2 * xAxis[2 * eleX + 1] - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);
        mX = (2 * xEval - (xAxis[2 * eleX + 2] + xAxis[2 * eleX])) / (xAxis[2 * eleX + 2] - xAxis[2 * eleX]);

        psiX[0] = (mX - mX0) * (mX - 1.0) / (2 * (1.0 + mX0));
        psiX[1] = (mX + 1.0) * (mX - 1.0) / ((mX0 + 1.0) * (mX0 - 1.0));
        psiX[2] = (mX + 1.0) * (mX - mX0) / (2 * (1.0 - mX0));
    }

    psiY.resize(3, 0.0);
    if (abs(2 * NyE - Ly) > pow(10.0, -3.0)) {
        eleY = NyE;
        if (yAxis.front() < yAxis.back()) {
            for (int i = 0; i < NyE; i++) {
                if (yEval <= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NyE; i++) {
                if (yEval >= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }

        if (eleY < NyE) {
            mY0 = (2 * yAxis[2 * eleY + 1] - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);
            mY = (2 * yEval - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);

            psiY[0] = (mY - mY0) * (mY - 1.0) / (2 * (1 + mY0));
            psiY[1] = (mY + 1.0) * (mY - 1.0) / ((mY0 + 1.0) * (mY0 - 1.0));
            psiY[2] = (mY + 1.0) * (mY - mY0) / (2 * (1.0 - mY0));
        }
        else {
            mY = (2 * yEval - (yAxis[2 * eleY + 1] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 1] - yAxis[2 * eleY]);

            psiY[0] = -(mY - 1.0) / 2.0;
            psiY[1] = (mY + 1.0) / 2.0;
        }
    }
    else {
        eleY = NyE - 1;
        if (yAxis.front() < yAxis.back()) {
            for (int i = 0; i < NyE; i++) {
                if (yEval <= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < NyE; i++) {
                if (yEval >= yAxis[2 * i + 2]) {
                    eleY = i;
                    break;
                }
            }
        }

        mY0 = (2 * yAxis[2 * eleY + 1] - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);
        mY = (2 * yEval - (yAxis[2 * eleY + 2] + yAxis[2 * eleY])) / (yAxis[2 * eleY + 2] - yAxis[2 * eleY]);

        psiY[0] = (mY - mY0) * (mY - 1.0) / (2 * (1.0 + mY0));
        psiY[1] = (mY + 1.0) * (mY - 1.0) / ((mY0 + 1.0) * (mY0 - 1.0));
        psiY[2] = (mY + 1.0) * (mY - mY0) / (2 * (1.0 - mY0));
    }

    if (eleX > NxE - 1) {
        if (eleY > NyE - 1) {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1];
        }
        else {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1] \
                + surface[2 * eleX][2 * eleY + 2] * psiX[0] * psiY[2] \
                + surface[2 * eleX + 1][2 * eleY + 2] * psiX[1] * psiY[2];
        }
    }
    else {
        if (eleY > NyE - 1) {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX + 2][2 * eleY] * psiX[2] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1] \
                + surface[2 * eleX + 2][2 * eleY + 1] * psiX[2] * psiY[1];
        }
        else {
            output = surface[2 * eleX][2 * eleY] * psiX[0] * psiY[0] \
                + surface[2 * eleX + 1][2 * eleY] * psiX[1] * psiY[0] \
                + surface[2 * eleX + 2][2 * eleY] * psiX[2] * psiY[0] \
                + surface[2 * eleX][2 * eleY + 1] * psiX[0] * psiY[1] \
                + surface[2 * eleX + 1][2 * eleY + 1] * psiX[1] * psiY[1] \
                + surface[2 * eleX + 2][2 * eleY + 1] * psiX[2] * psiY[1] \
                + surface[2 * eleX][2 * eleY + 2] * psiX[0] * psiY[2] \
                + surface[2 * eleX + 1][2 * eleY + 2] * psiX[1] * psiY[2] \
                + surface[2 * eleX + 2][2 * eleY + 2] * psiX[2] * psiY[2];
        }
    }
}

void transpose(std::vector<std::vector<double>> &matrix) {
    int Nx = (int)matrix.size();
    int Ny = (int)matrix[0].size();

    std::vector<std::vector<double>> output;
    output.resize(Ny);
    for (int i = 0; i < Ny; i++) {
        output[i].resize(Nx);
        for (int j = 0; j < Nx; j++) {
            output[i][j] = matrix[j][i];
        }
    }
    matrix = output;
}

void lineOut(bool normalize, double power, imageBW& Image, int AxisSum, std::vector<double>& output) {
    int Ni, Nj;
    if (AxisSum == 1) {
        Ni = Image.sizeX();
        Nj = Image.sizeY();

        std::vector<int> indexT;
        indexT.resize(Ni);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        output.resize(Ni,0.0);
        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            for (int j = 0; j < Nj; j++) {
                output[i] = output[i] + std::pow(Image.value(i, j), power);
            }
        });

        if (normalize) {
            double maxValue = -pow(2.0, 32.0);
            double minValue = pow(2.0, 32.0);
            for (int i = 0; i < Ni; i++) {
                if (output[i] < minValue) {
                    minValue = output[i];
                }
                if (output[i] > maxValue) {
                    maxValue = output[i];
                }
            }
            concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
                output[i] = (output[i] - minValue) / (maxValue - minValue);
            });
        }
    }
    else {
        Ni = Image.sizeY();
        Nj = Image.sizeX();

        std::vector<int> indexT;
        indexT.resize(Ni);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        output.resize(Ni, 0.0);
        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            for (int j = 0; j < Nj; j++) {
                output[i] = output[i] + std::pow(Image.value(j, i),power);
            }
        });

        if (normalize) {
            double maxValue = -pow(2.0, 32.0);
            double minValue = pow(2.0, 32.0);
            for (int i = 0; i < Ni; i++) {
                if (output[i] < minValue) {
                    minValue = output[i];
                }
                if (output[i] > maxValue) {
                    maxValue = output[i];
                }
            }
            concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
                output[i] = (output[i] - minValue) / (maxValue - minValue);
            });
        }
    }
}

void medianFilter(std::vector<double>& data, int windowRadius) {
    int N = (int)data.size();

    std::vector<int> indexT;
    indexT.resize(N);
    int countT = 0;
    std::generate(std::begin(indexT), std::end(indexT), [&] {
        return countT++;
    });

    std::vector<double> filtered;
    filtered.resize(N, 0.0);
    concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
        if (i >= windowRadius) {
            if (i < N - windowRadius) {
                std::vector<double> window;
                window.resize(2 * windowRadius + 1, 0.0);
                for (int j = 0; j < 2 * windowRadius + 1; j++) {
                    window[j] = data[i + j - windowRadius];
                }
                std::sort(window.begin(), window.end());
                filtered[i] = window[windowRadius + 1];
            }
            else {
                std::vector<double> window;
                window.resize(2 * windowRadius + 1, 0.0);
                int j0 = N - 1 - i - 2 * windowRadius;
                for (int j = 0; j < 2 * windowRadius + 1; j++) {
                    window[j] = data[i + j + j0];
                }
                std::sort(window.begin(), window.end());
                filtered[i] = window[windowRadius + 1];
            }
        }
        else {
            std::vector<double> window;
            window.resize(2 * windowRadius + 1, 0.0);
            for (int j = 0; j < 2 * windowRadius + 1; j++) {
                window[j] = data[j];
            }
            std::sort(window.begin(), window.end());
            filtered[i] = window[windowRadius + 1];
        }
    });

    data = filtered;
}

void medianFilter(imageBW& data, int windowRadius) {
    int Nx = data.sizeX();
    int Ny = data.sizeY();
    int Nw = 2 * windowRadius + 1;

    std::vector<int> indexT;
    indexT.resize(Nx);
    int countT = 0;
    std::generate(std::begin(indexT), std::end(indexT), [&] {
        return countT++;
        });

    imageBW filtered;
    filtered.resize(Nx, Ny);
    int indexW = (int)((Nw * Nw + 1) / 2);
    concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
        if (i >= windowRadius) {
            if (i < Nx - windowRadius) {
                for (int j = 0; j < Ny; j++) {
                    if (j >= windowRadius) {
                        if (j < Ny - windowRadius) {
                            std::vector<double> window(Nw * Nw, 0.0);
                            for (int k = 0; k < Nw; k++) {
                                for (int l = 0; l < Nw; l++) {
                                    window[k * Nw + l] = data.value(i + k - windowRadius, j + l - windowRadius);
                                }
                            }
                            std::sort(window.begin(), window.end());
                            filtered.definePixel(i, j, window[indexW]);
                            window.clear();
                        }
                        else {
                            std::vector<double> window(Nw * Nw, 0.0);
                            int l0 = Ny - 1 - j - 2 * windowRadius;
                            for (int k = 0; k < Nw; k++) {
                                for (int l = 0; l < Nw; l++) {
                                    window[k * Nw + l] = data.value(i + k - windowRadius, j + l + l0);
                                }
                            }
                            std::sort(window.begin(), window.end());
                            filtered.definePixel(i, j, window[indexW]);
                            window.clear();
                        }
                    }
                    else {
                        std::vector<double> window(Nw * Nw, 0.0);
                        for (int k = 0; k < Nw; k++) {
                            for (int l = 0; l < Nw; l++) {
                                window[k * Nw + l] = data.value(i + k - windowRadius, l);
                            }
                        }
                        std::sort(window.begin(), window.end());
                        filtered.definePixel(i, j, window[indexW]);
                        window.clear();
                    }
                    
                }
            }
            else {
                int k0 = Nx - 1 - i - 2 * windowRadius;
                for (int j = 0; j < Ny; j++) {
                    if (j >= windowRadius) {
                        if (j < Ny - windowRadius) {
                            std::vector<double> window(Nw * Nw, 0.0);
                            for (int k = 0; k < Nw; k++) {
                                for (int l = 0; l < Nw; l++) {
                                    window[k * Nw + l] = data.value(i + k + k0, j + l - windowRadius);
                                }
                            }
                            std::sort(window.begin(), window.end());
                            filtered.definePixel(i, j, window[indexW]);
                            window.clear();
                        }
                        else {
                            std::vector<double> window(Nw * Nw, 0.0);
                            int l0 = Ny - 1 - j - 2 * windowRadius;
                            for (int k = 0; k < Nw; k++) {
                                for (int l = 0; l < Nw; l++) {
                                    window[k * Nw + l] = data.value(i + k + k0, j + l + l0);
                                }
                            }
                            std::sort(window.begin(), window.end());
                            filtered.definePixel(i, j, window[indexW]);
                            window.clear();
                        }
                    }
                    else {
                        std::vector<double> window(Nw * Nw, 0.0);
                        for (int k = 0; k < Nw; k++) {
                            for (int l = 0; l < Nw; l++) {
                                window[k * Nw + l] = data.value(i + k + k0, l);
                            }
                        }
                        std::sort(window.begin(), window.end());
                        filtered.definePixel(i, j, window[indexW]);
                        window.clear();
                    }

                }
            }
        }
        else {
            for (int j = 0; j < Ny; j++) {
                if (j >= windowRadius) {
                    if (j < Ny - windowRadius) {
                        std::vector<double> window(Nw* Nw, 0.0);
                        for (int k = 0; k < Nw; k++) {
                            for (int l = 0; l < Nw; l++) {
                                window[k * Nw + l] = data.value(k, j + l - windowRadius);
                            }
                        }
                        std::sort(window.begin(), window.end());
                        filtered.definePixel(i, j, window[indexW]);
                        window.clear();
                    }
                    else {
                        std::vector<double> window(Nw* Nw, 0.0);
                        int l0 = Ny - 1 - j - 2 * windowRadius;
                        for (int k = 0; k < Nw; k++) {
                            for (int l = 0; l < Nw; l++) {
                                window[k * Nw + l] = data.value(k, j + l + l0);
                            }
                        }
                        std::sort(window.begin(), window.end());
                        filtered.definePixel(i, j, window[indexW]);
                        window.clear();
                    }
                }
                else {
                    std::vector<double> window(Nw* Nw, 0.0);
                    for (int k = 0; k < Nw; k++) {
                        for (int l = 0; l < Nw; l++) {
                            window[k * Nw + l] = data.value(k, l);
                        }
                    }
                    std::sort(window.begin(), window.end());
                    filtered.definePixel(i, j, window[indexW]);
                    window.clear();
                }

            }
        }
    });

    data = filtered;
}

void removeOutlier(imageBW& data, double sigmaOrder) {
    int Nx = data.sizeX();
    int Ny = data.sizeY();

    double meanValue = 0.0;

    std::vector<int> indexT;
    indexT.resize(Nx);
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            meanValue = meanValue + data.value(i, j);
        }
        indexT[i] = i;
    }
    meanValue = meanValue / (Nx * Ny);
    double stdValue = 0.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            stdValue = stdValue + (data.value(i, j) - meanValue) * (data.value(i, j) - meanValue);
        }
    }
    stdValue = sqrt(stdValue / (Nx * Ny));

    double limit = meanValue + sigmaOrder * stdValue;
    concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
        for (int j = 0; j < Ny; j++) {
            if (data.value(i, j) > limit) {
                data.definePixel(i, j, limit);
            }
        }
    });
}

/*
void findCenter(imageBW& Image, std::vector<int>& center) {
    int Nx = Image.sizeX();
    int Ny = Image.sizeY();

    std::vector<double> xlinePointing, ylinePointing;
    lineOut(1, 1.0, Image, 0, xlinePointing);
    lineOut(1, 1.0, Image, 1, ylinePointing);

    int windowRadius = 2;
    medianFilter(xlinePointing, windowRadius);
    medianFilter(ylinePointing, windowRadius);

    double minXL, minXH, minYL, minYH, pxValue;
    minXL = xlinePointing[Nx];
    minXH = minXL;
    minYL = ylinePointing[Ny];
    minYH = minYL;
    int indexXL, indexXH, indexYL, indexYH;
    indexXL = 0;
    indexXH = Nx - 1;
    indexYL = 0;
    indexYH = Ny - 1;

    int N = std::max(Nx, Ny);
    for (int i = 0; i < N; i++) {
        if (i < Nx) {
            pxValue = xlinePointing[i];
            if (pxValue < minXL) {
                minXL = pxValue;
                indexXL = i;
            }
            pxValue = xlinePointing[Nx - 1 - i];
            if (pxValue < minXH) {
                minXH = pxValue;
                indexXH = Nx - 1 - i;
            }
        }
        if (i < Ny) {
            pxValue = ylinePointing[i];
            if (pxValue < minYL) {
                minYL = pxValue;
                indexYL = i;
            }
            pxValue = ylinePointing[Ny - 1 - i];
            if (pxValue < minYH) {
                minYH = pxValue;
                indexYH = Ny - 1 - i;
            }
        }
    }

    center.resize(2);
    center[0] = round((indexXL + indexXH) / 2);
    center[1] = round((indexYL + indexYH) / 2);
}
*/
/*
void designerMatrix(std::vector<std::vector<int>>& input, std::vector<std::vector<double>>& matrixM) {
    //input should be 4 pairwise points (x,y);
    std::vector<std::vector<int>> uv;
    uv.resize(4);
    int xH, xL, yH, yL;
    xH = std::max(std::max(input[0][0], input[1][0]), std::max(input[2][0], input[3][0]));
    xL = std::min(std::min(input[0][0], input[1][0]), std::min(input[2][0], input[3][0]));
    yH = std::max(std::max(input[0][1], input[1][1]), std::max(input[2][1], input[3][1]));
    yL = std::min(std::min(input[0][1], input[1][1]), std::min(input[2][1], input[3][1]));

    uv[0].resize(2,0);
    uv[0][0] = xL;
    uv[0][1] = yL;

    uv[1].resize(2, 0);
    uv[1][0] = xL;
    uv[1][1] = yH;

    uv[2].resize(2, 0);
    uv[2][0] = xH;
    uv[2][1] = yH;

    uv[3].resize(2, 0);
    uv[3][0] = xH;
    uv[3][1] = yL;

    std::vector<std::vector<double>> matrixA;
    matrixA.resize(8);
    
    matrixA[0].resize(9, 0.0);
    matrixA[1].resize(9, 0.0);
    matrixA[2].resize(9, 0.0);
    matrixA[3].resize(9, 0.0);
    matrixA[4].resize(9, 0.0);
    matrixA[5].resize(9, 0.0);
    matrixA[6].resize(9, 0.0);
    matrixA[7].resize(9, 0.0);

    matrixA[0][0] = (double)input[0][0] * (double)uv[0][1];
    matrixA[0][1] = (double)input[0][1] * (double)uv[0][1];
    matrixA[0][2] = 1.0 * (double)uv[0][1];

    matrixA[0][3] = -1.0 * (double)input[0][0] * ((double)uv[0][0] + 1.0);
    matrixA[0][4] = -1.0 * (double)input[0][1] * ((double)uv[0][0] + 1.0);
    matrixA[0][5] = -1.0 * ((double)uv[0][0] + 1.0);

    matrixA[0][6] = (double)input[0][0] * (double)uv[0][1];
    matrixA[0][7] = (double)input[0][1] * (double)uv[0][1];
    matrixA[0][8] = 1.0 * (double)uv[0][1];

    matrixA[1][0] = (double)input[0][0] * ((double)uv[0][1] + 1.0);
    matrixA[1][1] = (double)input[0][1] * ((double)uv[0][1] + 1.0);
    matrixA[1][2] = 1.0 * ((double)uv[0][1] + 1.0);

    matrixA[1][3] = -1.0 * (double)input[0][0] * (double)uv[0][0];
    matrixA[1][4] = -1.0 * (double)input[0][1] * (double)uv[0][0];
    matrixA[1][5] = -1.0 * (double)uv[0][0];

    matrixA[1][6] = -1.0 * (double)input[0][0] * (double)uv[0][0];
    matrixA[1][7] = -1.0 * (double)input[0][1] * (double)uv[0][0];
    matrixA[1][8] = -1.0 * (double)uv[0][0];


    matrixA[2][0] = (double)input[1][0] * (double)uv[1][1];
    matrixA[2][1] = (double)input[1][1] * (double)uv[1][1];
    matrixA[2][2] = 1.0 * (double)uv[1][1];

    matrixA[2][3] = -1.0 * (double)input[1][0] * ((double)uv[1][0] + 1.0);
    matrixA[2][4] = -1.0 * (double)input[1][1] * ((double)uv[1][0] + 1.0);
    matrixA[2][5] = -1.0 * ((double)uv[1][0] + 1.0);

    matrixA[2][6] = (double)input[1][0] * (double)uv[1][1];
    matrixA[2][7] = (double)input[1][1] * (double)uv[1][1];
    matrixA[2][8] = 1.0 * (double)uv[1][1];

    matrixA[3][0] = (double)input[1][0] * ((double)uv[1][1] + 1.0);
    matrixA[3][1] = (double)input[1][1] * ((double)uv[1][1] + 1.0);
    matrixA[3][2] = 1.0 * ((double)uv[1][1] + 1.0);

    matrixA[3][3] = -1.0 * (double)input[1][0] * (double)uv[1][0];
    matrixA[3][4] = -1.0 * (double)input[1][1] * (double)uv[1][0];
    matrixA[3][5] = -1.0 * (double)uv[1][0];

    matrixA[3][6] = -1.0 * (double)input[1][0] * (double)uv[1][0];
    matrixA[3][7] = -1.0 * (double)input[1][1] * (double)uv[1][0];
    matrixA[3][8] = -1.0 * (double)uv[1][0];


    matrixA[4][0] = (double)input[2][0] * (double)uv[2][1];
    matrixA[4][1] = (double)input[2][1] * (double)uv[2][1];
    matrixA[4][2] = 1.0 * (double)uv[2][1];

    matrixA[4][3] = -1.0 * (double)input[2][0] * ((double)uv[2][0] + 1.0);
    matrixA[4][4] = -1.0 * (double)input[2][1] * ((double)uv[2][0] + 1.0);
    matrixA[4][5] = -1.0 * ((double)uv[2][0] + 1.0);

    matrixA[4][6] = (double)input[2][0] * (double)uv[2][1];
    matrixA[4][7] = (double)input[2][1] * (double)uv[2][1];
    matrixA[4][8] = 1.0 * (double)uv[2][1];

    matrixA[5][0] = (double)input[2][0] * ((double)uv[2][1] + 1.0);
    matrixA[5][1] = (double)input[2][1] * ((double)uv[2][1] + 1.0);
    matrixA[5][2] = 1.0 * ((double)uv[2][1] + 1.0);

    matrixA[5][3] = -1.0 * (double)input[2][0] * (double)uv[2][0];
    matrixA[5][4] = -1.0 * (double)input[2][1] * (double)uv[2][0];
    matrixA[5][5] = -1.0 * (double)uv[2][0];

    matrixA[5][6] = -1.0 * (double)input[2][0] * (double)uv[2][0];
    matrixA[5][7] = -1.0 * (double)input[2][1] * (double)uv[2][0];
    matrixA[5][8] = -1.0 * (double)uv[2][0];


    matrixA[6][0] = (double)input[3][0] * (double)uv[3][1];
    matrixA[6][1] = (double)input[3][1] * (double)uv[3][1];
    matrixA[6][2] = 1.0 * (double)uv[3][1];

    matrixA[6][3] = -1.0 * (double)input[3][0] * ((double)uv[3][0] + 1.0);
    matrixA[6][4] = -1.0 * (double)input[3][1] * ((double)uv[3][0] + 1.0);
    matrixA[6][5] = -1.0 * ((double)uv[3][0] + 1.0);

    matrixA[6][6] = (double)input[3][0] * (double)uv[3][1];
    matrixA[6][7] = (double)input[3][1] * (double)uv[3][1];
    matrixA[6][8] = 1.0 * (double)uv[3][1];

    matrixA[7][0] = (double)input[3][0] * ((double)uv[3][1] + 1.0);
    matrixA[7][1] = (double)input[3][1] * ((double)uv[3][1] + 1.0);
    matrixA[7][2] = 1.0 * ((double)uv[3][1] + 1.0);

    matrixA[7][3] = -1.0 * (double)input[3][0] * (double)uv[3][0];
    matrixA[7][4] = -1.0 * (double)input[3][1] * (double)uv[3][0];
    matrixA[7][5] = -1.0 * (double)uv[3][0];

    matrixA[7][6] = -1.0 * (double)input[3][0] * (double)uv[3][0];
    matrixA[7][7] = -1.0 * (double)input[3][1] * (double)uv[3][0];
    matrixA[7][8] = -1.0 * (double)uv[3][0];

    matrixM.resize(9);
    for (int i = 0; i < 9; i++) {
        matrixM[i].resize(9, 0.0);
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 8; k++) {
                matrixM[i][j] = matrixM[i][j] + matrixA[k][i] * matrixA[k][j];
            }
        }
    }
}
*/

void Contrast(double scale, std::vector<double>& input) {
    int N = (int)input.size();

    std::vector<int> indexT;
    indexT.resize(N - 1);
    int countT = 1;
    std::generate(std::begin(indexT), std::end(indexT), [&] {
        return countT++;
    });

    double avg = 0.0;
    for (int i = 0; i < N; i++) {
        avg = avg + input[i];
    }
    avg = avg / N;

    concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
        input[i] = (scale * (input[i] - avg) + avg);
    });
    input[0] = input[1];
    input[N - 1] = input[N - 2];
}

void Average(std::vector<double>& input, double& avg) {
    int N = (int)input.size();
    avg = 0.0;
    for (int i = 0; i < N; i++) {
        avg = avg + input[i];
    }
    avg = avg / N;
}

void Difference(std::vector<double>& input, std::vector<double>& output) {
    int N = (int)input.size() - 1;
    output.resize(N, 0.0);
    
    std::vector<int> indexT;
    indexT.resize(N);
    int countT = 0;
    std::generate(std::begin(indexT), std::end(indexT), [&] {
        return countT++;
    });
    concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
        output[i] = input[i + 1] - input[i];
    });
}

void Deviation(std::vector<double>& input, double& mean, double& stdev) {
    int N = (int)input.size();
    stdev = 0.0;
    for (int i = 0; i < N; i++) {
        stdev = stdev + (input[i] - mean) * (input[i] - mean);
    }
    stdev = sqrt(stdev / N);
}

/*
double pxScale(std::vector<double>& input) {
    int N = (int)input.size();

    std::vector<int> peak;
    for (int i = 0; i < N; i++) {
        if (input[i] > 0.9) {
            peak.push_back(i);
        }
    }

    bool loop = 1;
    double averageSpacing = 0.0;
    std::vector<double> spacing;
    double value;
    N = (int)peak.size();
    for (int i = 0; i < N - 1; i++) {
        value = (double)peak[i + 1] - (double)peak[i];
        if (value > 1) {
            spacing.push_back(value);
            averageSpacing = averageSpacing + value;
        }
    }
    N = (int)spacing.size();
    averageSpacing = averageSpacing / (N-1);
    int iStart = 0;
    while (loop) {
        for (int i = 0; i < N; i++) {
            if (spacing[i] > 1.75 * averageSpacing) {
                spacing.erase(spacing.begin() + i);
                break;
            }
            if (i == N - 1) {
                loop = 0;
            }
        }
        N = (int)spacing.size();
        averageSpacing = 0.0;
        for (int i = 0; i < N; i++) {
            averageSpacing = averageSpacing + spacing[i];
        }
        averageSpacing = averageSpacing / (N-1);
    }
    return averageSpacing;
}
*/

void linReg(std::vector<double>& x, std::vector<double>& y) {
    int N = std::min((int)x.size(), (int)y.size());

    double sumX, sumY, sumXY, sumX2;
    sumX = 0;
    sumY = 0;
    sumXY = 0;
    sumX2 = 0;
    for (int i = 0; i < N ; i++) {
        sumX = sumX + x[i];
        sumY = sumY + y[i];
        sumXY = sumXY + x[i] * y[i];
        sumX2 = sumX2 + x[i] * x[i];
    }

    double a0, a1;

    a1 = (N * sumXY - sumX * sumY) / (N * sumX2 - sumX * sumX);
    a0 = (sumX2 * sumY - sumX * sumXY) / (N * sumX2 - sumX * sumX);

    for (int i = 0; i < N; i++) {
        y[i] = a0 + a1 * x[i];
    }
}

void findMax(std::vector<double>& input, int& indexMax) {
    int N = (int)input.size();

    indexMax = 0;
    double valueMax = 0.0;
    for (int i = 0; i < N; i++) {
        if (input[i] > valueMax) {
            indexMax = i;
            valueMax = input[i];
        }
    }
}

#endif