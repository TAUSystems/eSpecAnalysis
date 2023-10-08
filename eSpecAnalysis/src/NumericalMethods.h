#pragma once

#ifndef __NumericalMethods_h__
#define __NumericalMethods_h__

#include "fileIO.h"

#define pi 3.14159265358979323846

//Composite Simpson's 3/8 Rule
void simpsonInt(double& dx, std::vector<double>& y, double& output) {
    output = 0.0;
    size_t N = y.size()-1;
    int remainder = N % 3;
    size_t Nend = (size_t)floor(N / 3);
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
        Ni = (int)Image.sizeX();
        Nj = (int)Image.sizeY();

        
        output.resize(Ni,0.0);
        #pragma omp parallel 
        {
            for (int i = 0; i < Ni; i++) {
                for (int j = 0; j < Nj; j++) {
                    output[i] = output[i] + std::pow(Image.value(i, j), power);
                }
            }
        }

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
            #pragma omp parallel 
            {
                for (int i = 0; i < Ni; i++) {
                    output[i] = (output[i] - minValue) / (maxValue - minValue);
                }
            }
        }
    }
    else {
        Ni = (int)Image.sizeY();
        Nj = (int)Image.sizeX();

        std::vector<int> indexT;
        indexT.resize(Ni);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        output.resize(Ni, 0.0);
        #pragma omp parallel 
        {
            for (int i = 0; i < Ni; i++) {
                for (int j = 0; j < Nj; j++) {
                    output[i] = output[i] + std::pow(Image.value(j, i), power);
                }
            }
        }

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
            #pragma omp parallel 
            {
                for (int i = 0; i < Ni; i++) {
                    output[i] = (output[i] - minValue) / (maxValue - minValue);
                }
            }
        }
    }
}

void medianFilter(std::vector<double>& data, int windowRadius) {
    int N = (int)data.size();

    std::vector<double> filtered;
    filtered.resize(N, 0.0);
    #pragma omp parallel 
    {
        for (int i = 0; i < N; i++) {
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
        }
    }

    data = filtered;
}

void medianFilter(imageBW& data, int windowRadius) {
    int Nx = (int)data.sizeX();
    int Ny = (int)data.sizeY();
    int Nw = 2 * windowRadius + 1;

    
    imageBW filtered;
    filtered.resize(Nx, Ny);
    int indexW = (int)((Nw * Nw + 1) / 2);
    #pragma omp parallel 
    {
        for (int i = 0; i < Nx; i++) {
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
                            std::vector<double> window(Nw * Nw, 0.0);
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
                            std::vector<double> window(Nw * Nw, 0.0);
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
                        std::vector<double> window(Nw * Nw, 0.0);
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
        }
    }

    data = filtered;
}

void removeOutlier(imageBW& data, double sigmaOrder) {
    int Nx = (int)data.sizeX();
    int Ny = (int)data.sizeY();

    double meanValue = 0.0;
    
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            meanValue = meanValue + data.value(i, j);
        }
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
    #pragma omp parallel 
    {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                if (data.value(i, j) > limit) {
                    data.definePixel(i, j, limit);
                }
            }
        }
    }
}

void Contrast(double scale, std::vector<double>& input) {
    int N = (int)input.size();

    double avg = 0.0;
    for (int i = 0; i < N; i++) {
        avg = avg + input[i];
    }
    avg = avg / N;

    #pragma omp parallel 
    {
        for (int i = 1; i < N - 1; i++) {
            input[i] = (scale * (input[i] - avg) + avg);
        }
    }
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
    
    #pragma omp parallel 
    {
        for (int i = 0; i < N; i++) {
            output[i] = input[i + 1] - input[i];
        }
    }
}

void Deviation(std::vector<double>& input, double& mean, double& stdev) {
    int N = (int)input.size();
    stdev = 0.0;
    for (int i = 0; i < N; i++) {
        stdev = stdev + (input[i] - mean) * (input[i] - mean);
    }
    stdev = sqrt(stdev / N);
}

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