
#include "main.h"
#include "Analysis.h"


int analyze_espec_images(std::string filepath_eScreenA, std::string filepath_eScreenB, std::string filepath_ePointing, std::string filepath_spectrum) {

    printf("Starting analyze_espec_images()\n");

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
            printf("In case 1\n");
            if (true || (stat(eSpec.screenPath(screen).c_str(), &sb) == 0)) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                paramSpace pSpace;
                pSpace.loadMap(pathCalibration);
                printf("calling pointingMode\n");
                pointingMode(filepath_eScreenA, filepath_eScreenB, filepath_ePointing, filepath_spectrum, 
                             rate, timeout, eSpec, calibration, pSpace, H, xRuler, yRuler
                            );
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


int main( void ) {
    analyze_espec_images("/home/reinier/data/eSpecAnalysis_test_data/shot-10192023102619-0-eScreenA.tiff",
                         "/home/reinier/data/eSpecAnalysis_test_data/shot-10192023102619-0-eScreenB.tiff",
                         "/home/reinier/data/eSpecAnalysis_test_data/shot-10192023102619-0-ePointing.tiff",
                         "/home/reinier/data/eSpecAnalysis_test_data/shot-10192023102619-0-spectrum.png"
    );
}