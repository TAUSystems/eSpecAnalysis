// eSpecAnalysis.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "Plot.h"

#include "Calibration.h"
#include "Analysis.h"

//#include <chrono>
//#include <thread>

int main() {
    std::string pathRoot = pathGet();

    std::string pathSetting = pathRoot;
    pathSetting = pathSetting + "\\settings.cfg";
    std::vector<std::string> settingsParameters = readFile(pathSetting);
    int opMode;
    double rate, timeout;
    getOpParameters(settingsParameters, opMode, rate, timeout);
    rate = 1.0 / (2.0 * rate);

    spectrometer eSpec;
    eSpec.generate(settingsParameters);

    std::string pathCalibration = pathRoot;
    pathCalibration = pathCalibration + "\\Calibration";
    screenCalibration calibration;
    calibration.loadCalibration(pathCalibration);
    
    std::vector<cv::Mat> H, Hread;
    std::vector<std::vector<double>> xRuler, yRuler;

    struct stat sb;

    /*
    int screen = 1;
    readCalibration(pathCalibration, H, xRuler, yRuler);
    paramSpace pSpace;
    pSpace.loadMap(pathCalibration);
    std::vector<double> EnAxis = pSpace.energy(screen);
    std::vector<double> PtAxis = pSpace.pointing(screen);
    std::vector<std::vector<double>> pS = pSpace.parameterSpace(screen);
    int NE = (int)EnAxis.size();
    int NP = (int)PtAxis.size();

    int pIndex = NP-1;
    double p0 = PtAxis[pIndex];
    std::vector<double> pSI, pST, dpS;
    pSI.resize(NE, 0.0);
    pST.resize(NE, 0.0);
    dpS.resize(NE, 0.0);
    for (int i = 0; i < NE; i++) {
        FE2DInterp(EnAxis, PtAxis, pS, EnAxis[i], p0, pSI[i]);
        pST[i] = pS[i][pIndex];
        dpS[i] = abs((pSI[i] - pST[i]) / pST[i]);
        std::cout << std::setprecision(3) << std::scientific << i << "\t" << EnAxis[i] << "\t" << pSI[i] << "\t" << pST[i] << "\t" << dpS[i] << "\n";
    }
    plt::figure();
    plt::subplot(1, 2, 1);
    plt::plot(EnAxis, pSI);
    plt::plot(EnAxis, pST);
    plt::subplot(1, 2, 2);
    plt::plot(EnAxis, dpS);
    plt::show();
    */

    
    int screen;
    switch(opMode) {
        case 0:
            std::cout << "\nCalibration Mode\n";
            calMode(eSpec, pathCalibration, H, xRuler, yRuler);
            break;
        case 1:
            screen = 1;
            std::cout << "\nAnalysis with Pointing Screen\n";
            if (stat(eSpec.screenPath(screen).c_str(), &sb) == 0) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                paramSpace pSpace;
                pSpace.loadMap(pathCalibration);
                pointingMode(rate, timeout, eSpec, calibration, pSpace, H, xRuler, yRuler);
            }
            else {
                std::cout << "\nData Path Not Found\n";
            }
            break;
        case 2:
            std::cout << "\nCurrently Not Implemented\n";
            //readCalibration(pathCalibration, H, xRuler, yRuler);
            break;
        case 3:
            std::cout << "\nCurrently Not Implemented\n";
            //readCalibration(pathCalibration, H, xRuler, yRuler);
            break;

        case 4:
            screen = 0;
            std::cout << "\nPointing View Auto Update View\n";
            if (stat(eSpec.screenPath(screen).c_str(), &sb) == 0) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                viewMode(rate, timeout, eSpec, calibration, screen, H, xRuler, yRuler);
            }
            else {
                std::cout << "\nData Path Not Found\n";
            }
            break;
        case 5:
            screen = 1;
            std::cout << "\neScreen A View Auto Update View\n";
            if (stat(eSpec.screenPath(screen).c_str(), &sb) == 0) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                viewMode(rate, timeout, eSpec, calibration, screen, H, xRuler, yRuler);
            }
            else {
                std::cout << "\nData Path Not Found\n";
            }
            break;
        case 6:
            screen = 2;
            std::cout << "\neScreen B View Auto Update View\n";
            if (stat(eSpec.screenPath(screen).c_str(), &sb) == 0) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                viewMode(rate, timeout, eSpec, calibration, screen, H, xRuler, yRuler);
            }
            else {
                std::cout << "\nData Path Not Found\n";
            }
            break;
    }
    system("PAUSE");
}
