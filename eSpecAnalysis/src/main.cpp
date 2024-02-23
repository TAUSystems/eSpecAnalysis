
#include "main.h"
#include "Analysis.h"


int analyze_espec_images(const char* filepath_eScreenA, const char* filepath_eScreenB, const char* filepath_ePointing, const char* filepath_spectrum) {
    printf("%s\n", filepath_eScreenA);

    int threadCount = omp_get_max_threads();
    omp_set_num_threads(threadCount);

    const std::filesystem::path configPath = getConfigPath();
    std::string pathSettings = (configPath / "settings.cfg").string();
    std::string pathCalibration = (configPath / "Calibration").string();
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
            if (true || (stat(eSpec.screenPath(screen).c_str(), &sb) == 0)) {
                readCalibration(pathCalibration, H, xRuler, yRuler);
                paramSpace pSpace;
                pSpace.loadMap(pathCalibration);
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


int main( int argc, char* argv[] ) {
    if (argc < 5) { 
        printf("Not enough arguments. Usage: \n  eSpecAnalysis low-energy-image-filename high-energy-image-filename pointing-image-filename spectrum-output-filename\n");
        return 1;
    }

    return analyze_espec_images(argv[1], argv[2], argv[3], argv[4]);
}