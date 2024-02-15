
#include "main.h"
#include "Analysis.h"


void analyze_espec_images(std::string filename) {
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
        case 1:
            screen = 1;
            // std::cout << "\nAnalysis with Pointing Screen\n";
            if (stat(eSpec.screenPath(screen).c_str(), &sb) == 0) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                paramSpace pSpace;
                pSpace.loadMap(pathCalibration);
                pointingMode(filename, rate, timeout, eSpec, calibration, pSpace, H, xRuler, yRuler);
            }
            else {
                std::cout << "\nData Path Not Found\n";
            }
            break;

        default:
            std::cout << "\nOnly opMode 1 (pointingMode) is implemented.\n";
            break;
    }

	return 0;
}
