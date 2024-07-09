
#include "Analysis.h"

int main() {
    int threadCount = omp_get_max_threads();
    omp_set_num_threads(threadCount);

	std::string pathRoot = getPath();
	std::string pathSettings = pathRoot + "/settings.cfg";
	std::string pathCalibration = pathRoot + "/Calibration";
	std::vector<std::string> settingsFile;
	std::vector<cv::Mat> H;
	std::vector<std::vector<double>> xRuler, yRuler;
	spectrometer eSpec;
	screenCalibration calibration;
	int opMode;
	double rate, timeout;
	struct stat sb;

	readFile(pathSettings, settingsFile);
	getOpParameters(settingsFile, opMode, rate, timeout);
	rate = 1.0 / (2.0 * rate);

    eSpec.generate(settingsFile);
	calibration.loadCalibration(pathCalibration);

    int screen;
    switch (opMode) {
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
        screen = 1;
        std::cout << "\nManual Analysis with Pointing Screen\n";
        std::cout << eSpec.screenPath(screen).c_str() << "\n";
        if (stat(eSpec.screenPath(screen).c_str(), &sb) == 0) {
            readCalibration(pathCalibration, H, xRuler, yRuler);
            paramSpace pSpace;
            pSpace.loadMap(pathCalibration);
            pointingModeManual(rate, timeout, eSpec, calibration, pSpace, H, xRuler, yRuler);
        }
        else {
            std::cout << "\nData Path Not Found\n";
        }
        break;
    case 3:
        std::cout << "\nCurrently Not Implemented\n";
        readCalibration(pathCalibration, H, xRuler, yRuler);
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
	return 0;
}
