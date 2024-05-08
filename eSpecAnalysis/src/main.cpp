
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
	// homography matrices for each of the three screens
    std::vector<cv::Mat> homographyMatrices;
    // x-axis and y-axis rulers for each of the three screens
	std::vector<std::vector<double>> xRuler, yRuler;
	spectrometer eSpec;
	screenCalibration calibration;
	int operatingMode;
	double rate, timeout;
	struct stat sb;

	readFile(pathSettings, settingsFile);
	getOpParameters(settingsFile, operatingMode, rate, timeout);
	rate = 1.0 / (2.0 * rate);

    eSpec.generate(settingsFile);
	// load calibration data, used in operatingMode 1 to get window size
    calibration.loadCalibration(pathCalibration);

    int screen;
    switch (operatingMode) {
        case 1:
            screen = 1;
            // std::cout << "\nAnalysis with Pointing Screen\n";
            if (true || (stat(eSpec.screenPath(screen).c_str(), &sb) == 0)) {
                // reads perspective.cache
                readCalibration(pathCalibration, homographyMatrices, xRuler, yRuler);
                trajectoryEndpointSurfaces pSpace;
                pSpace.loadMap(pathCalibration);
                pointingMode(filepath_eScreenA, filepath_eScreenB, filepath_ePointing, filepath_spectrum, 
                             rate, timeout, eSpec, calibration, pSpace, homographyMatrices, xRuler, yRuler
                            );
            }
            else {
                std::cout << "\nData Path Not Found\n";
            }
            break;

        default:
            std::cout << "\nOnly operatingMode 1 (pointingMode) is implemented.\n";
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