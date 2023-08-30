#pragma once

#pragma once

#ifndef __fileIO_h__
#define __fileIO_h__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ppl.h>
#include <ctime>
#include <vector>

#include <windows.h>
#include <filesystem>
#include <sys/stat.h>
#include "tiffio.h"

#include "opencv2/core/core.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"

template<int I>
struct CvType {};
template<>
struct CvType<CV_64F> { typedef double type_t; };

std::string getDate() {
    std::time_t timeRaw;
    std::tm* timeInfo;
    char buffer[80];
    std::time(&timeRaw);
    timeInfo = std::localtime(&timeRaw);
    std::strftime(buffer, 80, "%Y%m%d", timeInfo);

    std::string dateString(buffer);
    return dateString;
}

class imageBW {
    double** data;
    int* size;
public:
    imageBW() {
        data = (double**)malloc(sizeof(double*));
        data[0] = (double*)malloc(sizeof(double));
        data[0][0] = 0.0;
        size = (int*)malloc(sizeof(int) * 2);
        size[0] = 1;
        size[1] = 1;
    }

    void resize(int& Nx, int& Ny) {
        size[0] = Nx;
        size[1] = Ny;
        data = (double**)malloc(sizeof(double*) * size[0]);

        std::vector<int> indexT;
	    indexT.resize(size[0]);
	    int countT = 0;
	    std::generate(std::begin(indexT), std::end(indexT), [&] {
		    return countT++;
	    });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            data[i] = (double*)malloc(sizeof(double) * size[1]);
            for (int j = 0; j < size[1]; j++) {
                data[i][j] = 0.0;
            }
        });
    }

    void resize(uint32& Nx, uint32& Ny) {
        size[0] = (int)Nx;
        size[1] = (int)Ny;
        data = (double**)malloc(sizeof(double*) * size[0]);

        std::vector<int> indexT;
        indexT.resize(size[0]);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            data[i] = (double*)malloc(sizeof(double) * size[1]);
            for (int j = 0; j < size[1]; j++) {
                data[i][j] = 0.0;
            }
        });
    }

    void destroy() {
        for (int i = 0; i < size[0]; i++) {
            free(data[i]);
        }
        free(data);
        free(size);
    }

    void definePixel(int indexX, int indexY, double value) {
        data[indexX][indexY] = value;
    }

    void definePixel(uint32 indexX, uint32 indexY, double value) {
        data[indexX][indexY] = value;
    }

    void definePixel(int indexX, uint32 indexY, double value) {
        data[indexX][indexY] = value;
    }

    void definePixel(uint32 indexX, int indexY, double value) {
        data[indexX][indexY] = value;
    }

    double value(int indexX, int indexY) {
        return data[indexX][indexY];
    }

    int sizeX() {
        return size[0];
    }

    int sizeY() {
        return size[1];
    }
    void crop(std::vector<int> bounds, imageBW& output) {
        int NxS, NxE, NyS, NyE;
        
        if (bounds[1] > size[0]) {
            NxE = size[0];
        }
        else {
            NxE = bounds[1];
        }
        if (bounds[0] < 0) {
            NxS = 0;
        }
        else {
            NxS = bounds[0];
        }
        if (bounds[3] > size[1]) {
            NyE = size[0];
        }
        else {
            NyE = bounds[3];
        }
        if (bounds[2] < 0) {
            NyS = 0;
        }
        else {
            NyS = bounds[2];
        }
        int Nx = NxE - NxS;
        int Ny = NyE - NyS;
        output.resize(Nx, Ny);

        std::vector<int> indexT;
        indexT.resize(Nx);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            int Ni, Nj;
            for (int j = 0; j < Ny; j++) {
                Ni = NxS + i;
                Nj = NyS + j;
                output.definePixel(i,j,data[Ni][Nj]);
            }
        });
    }
};


class imageColor {
    double*** data;
    int* size;
public:
    imageColor() {
        data = (double***)malloc(sizeof(double**));
        data[0] = (double**)malloc(sizeof(double*));
        data[0][0] = (double*)malloc(sizeof(double) * 4);
        data[0][0][0] = 0.0;
        data[0][0][1] = 0.0;
        data[0][0][2] = 0.0;
        data[0][0][3] = 0.0;
        size = (int*)malloc(sizeof(int) * 2);
        size[0] = 1;
        size[1] = 1;
    }

    void resize(int& Nx, int& Ny) {
        size[0] = Nx;
        size[1] = Ny;
        data = (double***)malloc(sizeof(double**) * size[0]);

        std::vector<int> indexT;
        indexT.resize(size[0]);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            data[i] = (double**)malloc(sizeof(double*) * size[1]);
            for (int j = 0; j < size[1]; j++) {
                data[i][j] = (double*)malloc(sizeof(double) * 4);
                data[i][j][0] = 0.0;
                data[i][j][1] = 0.0;
                data[i][j][2] = 0.0;
                data[i][j][3] = 0.0;
            }
        });
    }

    void resize(uint32& Nx, uint32& Ny) {
        size[0] = (int)Nx;
        size[1] = (int)Ny;
        data = (double***)malloc(sizeof(double**) * size[0]);

        std::vector<int> indexT;
        indexT.resize(size[0]);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            data[i] = (double**)malloc(sizeof(double*) * size[1]);
            for (int j = 0; j < size[1]; j++) {
                data[i][j] = (double*)malloc(sizeof(double) * 4);
                data[i][j][0] = 0.0;
                data[i][j][1] = 0.0;
                data[i][j][2] = 0.0;
                data[i][j][3] = 0.0;
            }
        });
    }

    void destroy() {
        for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
                free(data[i][j]);
            }
            free(data[i]);
        }
        free(data);
        free(size);
    }

    void definePixel(int& indexX, int& indexY, double r, double g, double b, double a) {
        data[indexX][indexY][0] = r;
        data[indexX][indexY][1] = g;
        data[indexX][indexY][2] = b;
        data[indexX][indexY][3] = a;
    }

    void definePixel(uint32& indexX, uint32& indexY, double r, double g, double b, double a) {
        data[indexX][indexY][0] = r;
        data[indexX][indexY][1] = g;
        data[indexX][indexY][2] = b;
        data[indexX][indexY][3] = a;
    }

    void definePixel(int& indexX, uint32& indexY, double r, double g, double b, double a) {
        data[indexX][indexY][0] = r;
        data[indexX][indexY][1] = g;
        data[indexX][indexY][2] = b;
        data[indexX][indexY][3] = a;
    }

    void definePixel(uint32& indexX, int& indexY, double r, double g, double b, double a) {
        data[indexX][indexY][0] = r;
        data[indexX][indexY][1] = g;
        data[indexX][indexY][2] = b;
        data[indexX][indexY][3] = a;
    }

    double r(int& indexX, int& indexY) {
        return data[indexX][indexY][0];
    }

    double g(int& indexX, int& indexY) {
        return data[indexX][indexY][1];
    }

    double b(int& indexX, int& indexY) {
        return data[indexX][indexY][2];
    }

    double a(int& indexX, int& indexY) {
        return data[indexX][indexY][3];
    }

    imageBW color2bw() {
        imageBW output;
        output.resize(size[0], size[1]);
        double pxValue;

        std::vector<int> indexT;
        indexT.resize(size[0]);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            for (int j = 0; j < size[1]; j++) {
                if (data[i][j] != 0) {
                    pxValue = 0.2990 * data[i][j][0] + 0.5870 * data[i][j][1] + 0.1140 * data[i][j][2];
                    output.definePixel(i, j, pxValue);
                }
            }
        });

        return output;
    }
    int sizeX() {
        return size[0];
    }
    int sizeY() {
        return size[1];
    }
};


class spectrometer {
    std::vector<std::string> path;
    double** screen;
    double* angle;
    std::string analysis;

public:
    spectrometer() {
        path.resize(3);
        screen = (double**)malloc(sizeof(double*) * 3);
        screen[0] = (double*)malloc(sizeof(double) * 5);
        screen[1] = (double*)malloc(sizeof(double) * 5);
        screen[2] = (double*)malloc(sizeof(double) * 5);

        screen[0][0] = 0.0;
        screen[0][1] = 0.0;
        screen[0][2] = 0.0;
        screen[0][3] = 0.0;
        screen[0][4] = 0.0;

        screen[1][0] = 0.0;
        screen[1][1] = 0.0;
        screen[1][2] = 0.0;
        screen[1][3] = 0.0;
        screen[1][4] = 0.0;

        screen[2][0] = 0.0;
        screen[2][1] = 0.0;
        screen[2][2] = 0.0;
        screen[2][3] = 0.0;
        screen[2][4] = 0.0;

        angle = (double*)malloc(sizeof(double) * 2);
        angle[0] = 0.0;
        angle[1] = 0.0;
    }
    void generate(std::vector<std::string>& settingsParameters) {
        int N = (int)settingsParameters.size();

        std::vector<int> indexT;
        indexT.resize(N);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        std::string rootPath;
        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            if (settingsParameters[i].find("dirRoot") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                rootPath = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                indexStart = rootPath.find("\"");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = rootPath.find("\"", indexStart + 1);
                    rootPath = rootPath.substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                indexStart = rootPath.find("\'");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = rootPath.find("\'", indexStart + 1);
                    rootPath = rootPath.substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                if (rootPath.find(" ") == 0) {
                    rootPath = rootPath.substr(1, rootPath.length() - 1);
                }
            }
            if (settingsParameters[i].find("dirPointing") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                path[0] = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                indexStart = path[0].find("\"");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = path[0].find("\"", indexStart + 1);
                    path[0] = path[0].substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                indexStart = path[0].find("\'");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = path[0].find("\'", indexStart + 1);
                    path[0] = path[0].substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                if (path[0].find(" ") == 0) {
                    path[0] = path[0].substr(1, path[0].length() - 1);
                }
            }

            if (settingsParameters[i].find("dirEScreenA") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                path[1] = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                indexStart = path[1].find("\"");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = path[1].find("\"", indexStart + 1);
                    path[1] = path[1].substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                indexStart = path[1].find("\'");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = path[1].find("\'", indexStart + 1);
                    path[1] = path[1].substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                if (path[1].find(" ") == 0) {
                    path[1] = path[1].substr(1, path[1].length() - 1);
                }
            }

            if (settingsParameters[i].find("dirEScreenB") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                path[2] = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                indexStart = path[2].find("\"");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = path[2].find("\"", indexStart + 1);
                    path[2] = path[2].substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                indexStart = path[2].find("\'");
                if (indexStart != std::string::npos) {
                    size_t indexEnd = path[2].find("\'", indexStart + 1);
                    path[2] = path[2].substr(indexStart + 1, indexEnd - indexStart - 1);
                }
                if (path[2].find(" ") == 0) {
                    path[2] = path[2].substr(1, path[2].length() - 1);
                }
            }

            if (settingsParameters[i].find("xPointing") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[0][0] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("yPointing") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[0][1] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("zPointing") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[0][2] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("phiPointing") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[0][3] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("thetaPointing") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[0][4] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("xMaxAngle") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                angle[0] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("yMaxAngle") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                angle[1] = std::stod(strValue.c_str());
            }

            if (settingsParameters[i].find("xEScreenA") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[1][0] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("yEScreenA") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[1][1] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("zEScreenA") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[1][2] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("phiEScreenA") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[1][3] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("thetaEScreenA") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[1][4] = std::stod(strValue.c_str());
            }

            if (settingsParameters[i].find("xEScreenB") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[2][0] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("yEScreenB") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[2][1] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("zEScreenB") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[2][2] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("phiEScreenB") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[2][3] = std::stod(strValue.c_str());
            }
            if (settingsParameters[i].find("thetaEScreenB") != std::string::npos) {
                size_t indexStart = settingsParameters[i].find("=") + 1;
                std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
                screen[2][4] = std::stod(strValue.c_str());
            }
        });
        std::string dateString = getDate();
        path[0] = rootPath + "\\" + dateString + "\\" + path[0];
        path[1] = rootPath + "\\" + dateString + "\\" + path[1];
        path[2] = rootPath + "\\" + dateString + "\\" + path[2];
        analysis = rootPath + "\\" + dateString;

        struct stat sb;
        if (stat(analysis.c_str(), &sb) == 0) {
            analysis = analysis + "\\Analysis";
            std::filesystem::create_directory(analysis);
        }
        else {
            std::cout << "\nData Path Not Found\n";
        }
    }

    std::string screenPath(int indexScreen) {
        return path[indexScreen];
    }

    double x(int indexScreen) {
        return screen[indexScreen][0];
    }

    double y(int indexScreen) {
        return screen[indexScreen][1];
    }

    double z(int indexScreen) {
        return screen[indexScreen][2];
    }
    double phi(int indexScreen) {
        return screen[indexScreen][3];
    }

    double theta(int indexScreen) {
        return screen[indexScreen][4];
    }

    double angleMax(int dir) {
        return angle[dir];
    }

    std::string analysisPath() {
        return analysis;
    }
};


std::string pathGet() {
    TCHAR buffer[MAX_PATH] = { 0 };
    GetModuleFileName(NULL, buffer, MAX_PATH);
    std::wstring::size_type pos = std::wstring(buffer).find_last_of(L"\\/");
    std::wstring path = std::wstring(buffer).substr(0, pos);
    std::string output(path.begin(), path.end());
    return output;
}

std::vector<std::string> fileList(std::string& path) {
    std::vector<std::string> fileName;
    std::string file;
    bool imgFile;
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        imgFile = 0;
        file = entry.path().generic_string();
        if (file.find(".tiff") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".tif") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".Tiff") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".Tif") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".TIFF") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".TIF") != std::string::npos) {
            imgFile = 1;
        }
        /*
        if (file.find(".png") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".Png") != std::string::npos) {
            imgFile = 1;
        }
        if (file.find(".PNG") != std::string::npos) {
            imgFile = 1;
        }
        */
        if (imgFile) {
            fileName.push_back(file);
        }
    }
    return fileName;
}

std::string findRef(std::vector<std::string>& list) {
    std::string output = "-1NoRef";
    int N = (int)list.size();

    for (int i = 0; i < N; i++) {
        if (list[i].find("RefImage") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("refimage") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("REFIMAGE") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("refIMAGE") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("REFimage") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("REFImage") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("refimg") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("RefImg") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("REFIMG") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("refIMG") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("REFimg") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
        if (list[i].find("REFImg") != std::string::npos) {
            if (list[i].find(".tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tiff") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".Tif") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIFF") != std::string::npos) {
                output = list[i];
                break;
            }
            if (list[i].find(".TIF") != std::string::npos) {
                output = list[i];
                break;
            }
            break;
        }
    }
    return output;
}

std::vector<std::string> readFile(std::string& pathFile) {
    std::ifstream file;
    file.open(pathFile);
    std::vector<std::string> fileOutput;
    std::string tempLine;
    if (file.is_open()) {
        while (file) {
            std::getline(file, tempLine);
            fileOutput.push_back(tempLine);
        }
    }
    size_t N = fileOutput.size() - 1;
    fileOutput.resize(N);
    return fileOutput;
}

bool readImage(std::string& pathFile, imageColor& Image) {
    const char* fileName = pathFile.c_str();
    TIFF* input = TIFFOpen(fileName, "r");
    if (input) {
        uint32 width, height;
        size_t npixels;
        uint32* raster;

        TIFFGetField(input, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(input, TIFFTAG_IMAGELENGTH, &height);
        npixels = (size_t)(width);
        npixels = npixels * (size_t)(height);
        raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));
        Image.resize(width, height);
        if (raster != NULL) {
            if (TIFFReadRGBAImage(input, width, height, raster, 0)) {
                std::vector<int> indexT;
                indexT.resize((int)width);
                int countT = 0;
                std::generate(std::begin(indexT), std::end(indexT), [&] {
                    return countT++;
                });

                concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
                    double pixelValueR, pixelValueG, pixelValueB, pixelValueA;
                    for (uint32 j = 0; j < height; j++) {
                        pixelValueR = (double)TIFFGetR(raster[j * width + (uint32)i]);
                        pixelValueG = (double)TIFFGetG(raster[j * width + (uint32)i]);
                        pixelValueB = (double)TIFFGetB(raster[j * width + (uint32)i]);
                        pixelValueA = (double)TIFFGetA(raster[j * width + (uint32)i]);
                        Image.definePixel(i, j, pixelValueR, pixelValueG, pixelValueB, pixelValueA);
                    }
                });
            }
            _TIFFfree(raster);
        }
        TIFFClose(input);
        return 0;
    }
    else {
        return 1;
    }
}

bool readImage(std::string& pathFile, imageBW& Image) {
    const char* fileName = pathFile.c_str();
    TIFF* input = TIFFOpen(fileName, "r");
    if (input) {
        uint32 width, height;
        size_t npixels;
        uint32* raster;

        TIFFGetField(input, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(input, TIFFTAG_IMAGELENGTH, &height);
        npixels = (size_t)width;
        npixels = npixels * (size_t)height;
        raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32));
        Image.resize(width, height);
        if (raster != NULL) {
            if (TIFFReadRGBAImage(input, width, height, raster, 0)) {
                std::vector<int> indexT;
                indexT.resize((int)width);
                int countT = 0;
                std::generate(std::begin(indexT), std::end(indexT), [&] {
                    return countT++;
                });

                concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
                    double pixelValue;
                    for (uint32 j = 0; j < height; j++) {
                        if (TIFFGetA(raster[j * width + (uint32)i]) != 0) {
                            pixelValue = 0.2990 * (double)TIFFGetR(raster[j * width + (uint32)i]) + 0.5870 * (double)TIFFGetG(raster[j * width + (uint32)i]) + 0.1140 * (double)TIFFGetB(raster[j * width + (uint32)i]);
                            Image.definePixel(i, j, pixelValue);
                        }
                    }
                });
            }
            _TIFFfree(raster);
        }
        TIFFClose(input);
        return 0;
    }
    else {
        return 1;
    }
}

class screenCalibration {
    std::vector<std::vector<cv::Point2d>> calPoints;
    std::vector<std::vector<double>> winSize;
    std::vector<std::vector<double>> threshold;
public:
    screenCalibration() {
        calPoints.resize(4);
        calPoints[0].resize(8, cv::Point2d(0.0, 0.0));
        calPoints[1].resize(8, cv::Point2d(0.0, 0.0));
        calPoints[2].resize(8, cv::Point2d(0.0, 0.0));
        calPoints[3].resize(8, cv::Point2d(0.0, 0.0));

        winSize.resize(4);
        winSize[0].resize(2, 0.0);
        winSize[1].resize(2, 0.0);
        winSize[2].resize(2, 0.0);
        winSize[3].resize(2, 0.0);

        threshold.resize(4);
        threshold[0].resize(4, 0.0);
        threshold[1].resize(4, 0.0);
        threshold[2].resize(4, 0.0);
        threshold[3].resize(4, 0.0);
    }
    void loadCalibration(std::string& pathCalibration) {
        std::string fileCalibration = pathCalibration + "\\perspective.cal";
        std::vector<std::string> calibration = readFile(fileCalibration);


        int N = (int)calibration.size();
        std::vector<int> indexT;
        indexT.resize(N);
        int countT = 0;
        std::generate(std::begin(indexT), std::end(indexT), [&] {
            return countT++;
        });

        concurrency::parallel_for_each(std::begin(indexT), std::end(indexT), [&](int i) {
            if (calibration[i].find("xTEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][0].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][0].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xLEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][3].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yLEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][3].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][2].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][2].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xREPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][1].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yREPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[0][1].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xPadEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[0][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yPadEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[0][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xContrastEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[0][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[0][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yContrastEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[0][2] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdEPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[0][3] = std::stod(strValue.c_str());
            }

            if (calibration[i].find("xTLEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][0].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTLEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][0].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBLEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][3].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBLEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][3].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBREScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][2].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBREScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][2].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xTREScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][1].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTREScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[1][1].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xPadEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[1][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yPadEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[1][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xContrastEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[1][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[1][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yContrastEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[1][2] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdEScreenA") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[1][3] = std::stod(strValue.c_str());
            }

            if (calibration[i].find("xTLEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][0].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTLEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][0].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBLEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][3].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBLEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][3].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBREScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][2].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBREScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][2].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xTREScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][1].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTREScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[2][1].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xPadEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[2][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yPadEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[2][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xContrastEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[2][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[2][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yContrastEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[2][2] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdEScreenB") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[2][3] = std::stod(strValue.c_str());
            }

            if (calibration[i].find("xTLPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][0].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTLPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][0].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBLPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][3].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBLPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][3].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xBRPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][2].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yBRPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][2].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xTRPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][1].x = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yTRPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                calPoints[3][1].y = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xPadPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[3][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yPadPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                winSize[3][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xContrastPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[3][0] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[3][1] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("yContrastPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[3][2] = std::stod(strValue.c_str());
            }
            if (calibration[i].find("xThresholdPPointing") != std::string::npos) {
                size_t indexStart = calibration[i].find("=") + 1;
                std::string strValue = calibration[i].substr(indexStart, calibration[i].length() - indexStart);
                threshold[3][3] = std::stod(strValue.c_str());
            }
        });

        for (int i = 0; i < 4; i++) {
            if (i == 0) {
                double xL, xH, yL, yH, xC, yC;
                xH = std::max(std::max(calPoints[i][0].x, calPoints[i][1].x), std::max(calPoints[i][2].x, calPoints[i][3].x));
                xL = std::min(std::min(calPoints[i][0].x, calPoints[i][1].x), std::min(calPoints[i][2].x, calPoints[i][3].x));
                yH = std::max(std::max(calPoints[i][0].y, calPoints[i][1].y), std::max(calPoints[i][2].y, calPoints[i][3].y));
                yL = std::min(std::min(calPoints[i][0].y, calPoints[i][1].y), std::min(calPoints[i][2].y, calPoints[i][3].y));
                xC = (xH - xL)/2 + winSize[i][0];
                yC = (yH - yL)/2 + winSize[i][1];

                xH = xH - xL + winSize[i][0];
                xL = winSize[i][0];
                yH = yH - yL + winSize[i][1];
                yL = winSize[i][1];

                calPoints[i][4].x = xC;
                calPoints[i][4].y = yL;
                calPoints[i][7].x = xL;
                calPoints[i][7].y = yC;
                calPoints[i][6].x = xC;
                calPoints[i][6].y = yH;
                calPoints[i][5].x = xH;
                calPoints[i][5].y = yC;

                winSize[i][0] = xH + winSize[i][0];
                winSize[i][1] = yH + winSize[i][1];
            }
            else {
                double xL, xH, yL, yH;
                xH = std::max(std::max(calPoints[i][0].x, calPoints[i][1].x), std::max(calPoints[i][2].x, calPoints[i][3].x));
                xL = std::min(std::min(calPoints[i][0].x, calPoints[i][1].x), std::min(calPoints[i][2].x, calPoints[i][3].x));
                yH = std::max(std::max(calPoints[i][0].y, calPoints[i][1].y), std::max(calPoints[i][2].y, calPoints[i][3].y));
                yL = std::min(std::min(calPoints[i][0].y, calPoints[i][1].y), std::min(calPoints[i][2].y, calPoints[i][3].y));

                xH = xH - xL + winSize[i][1];
                xL = winSize[i][1];
                yH = yH - yL + winSize[i][1];
                yL = winSize[i][1];

                calPoints[i][4].x = xL;
                calPoints[i][4].y = yL;
                calPoints[i][7].x = xL;
                calPoints[i][7].y = yH;
                calPoints[i][6].x = xH;
                calPoints[i][6].y = yH;
                calPoints[i][5].x = xH;
                calPoints[i][5].y = yL;

                winSize[i][0] = xH + winSize[i][0];
                winSize[i][1] = yH + winSize[i][1];
            }
        }
    }

    std::vector<cv::Point2d> viewQuad(int screen) {
        return calPoints[screen];
    }
    std::vector<double> viewResolution(int screen) {
        return winSize[screen];
    }
    std::vector<double> viewThreshold(int screen) {
        return threshold[screen];
    }
};

class paramSpace {
    std::vector<std::vector<std::vector<double>>> ps;
    std::vector<std::vector<double>> energyAxis;
    std::vector<std::vector<double>> pointingAxis;
public:
    paramSpace() {
        ps.resize(2);
        energyAxis.resize(2);
        pointingAxis.resize(2);
    }
    void loadMap(std::string calPath) {
        std::string filePath;
        double dE, Estart, dp, pStart;
        size_t index;
        std::string strValue;
        std::vector<std::string> buffer;
        std::vector<double> bufferM;
        bool loop;
        int m, n, N;


        filePath = calPath + "\\eScreenA.map";
        buffer = readFile(filePath);
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        Estart = std::stod(strValue.c_str());
        index = index + 1;
        buffer[0] = buffer[0].substr(index, buffer[0].length());
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        dE = std::stod(strValue.c_str());
        index = index + 1;
        buffer[0] = buffer[0].substr(index, buffer[0].length());
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        pStart = std::stod(strValue.c_str());
        index = index + 1;
        buffer[0] = buffer[0].substr(index, buffer[0].length());
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        dp = std::stod(strValue.c_str());

        m = (int)buffer.size();
        ps[0].clear();
        for (int i = 1; i < m; i++) {
            loop = 1;
            bufferM.clear();
            while (loop) {
                if (buffer[i].find(",") != std::string::npos) {
                    index = buffer[i].find(",");
                    strValue = buffer[i].substr(0, index);
                    bufferM.push_back(std::stod(strValue.c_str()));
                    index = index + 1;
                    buffer[i] = buffer[i].substr(index, buffer[i].length());
                }
                else {
                    bufferM.push_back(std::stod(buffer[i].c_str()));
                    loop = 0;
                }
            }
            ps[0].push_back(bufferM);
        }
        m = m - 1;
        n = (int)ps[0][0].size();
        N = std::max(n, m);
        for (int i = 0; i < N; i++) {
            if (i < m) {
                energyAxis[0].push_back(Estart + dE * (double)i);
            }
            if (i < n) {
                pointingAxis[0].push_back(pStart + dp * (double)i);
            }
        }

        filePath = calPath + "\\eScreenB.map";
        buffer = readFile(filePath);
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        Estart = std::stod(strValue.c_str());
        index = index + 1;
        buffer[0] = buffer[0].substr(index, buffer[0].length());
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        dE = std::stod(strValue.c_str());
        index = index + 1;
        buffer[0] = buffer[0].substr(index, buffer[0].length());
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        pStart = std::stod(strValue.c_str());
        index = index + 1;
        buffer[0] = buffer[0].substr(index, buffer[0].length());
        index = buffer[0].find(",");
        strValue = buffer[0].substr(0, index);
        dp = std::stod(strValue.c_str());

        m = (int)buffer.size();
        ps[1].clear();
        for (int i = 1; i < m; i++) {
            loop = 1;
            bufferM.clear();
            while (loop) {
                if (buffer[i].find(",") != std::string::npos) {
                    index = buffer[i].find(",");
                    strValue = buffer[i].substr(0, index);
                    bufferM.push_back(std::stod(strValue.c_str()));
                    index = index + 1;
                    buffer[i] = buffer[i].substr(index, buffer[i].length());
                }
                else {
                    bufferM.push_back(std::stod(buffer[i].c_str()));
                    loop = 0;
                }
            }
            ps[1].push_back(bufferM);
        }
        m = m - 1;
        n = (int)ps[1][0].size();
        N = std::max(n, m);
        for (int i = 0; i < N; i++) {
            if (i < m) {
                energyAxis[1].push_back(Estart + dE * (double)i);
            }
            if (i < n) {
                pointingAxis[1].push_back(pStart + dp * (double)i);
            }
        }
    }

    std::vector<double> energy(int screen) { 
        return energyAxis[screen - 1];
    }

    std::vector<double> pointing(int screen) {
        return pointingAxis[screen - 1];
    }

    std::vector<std::vector<double>> parameterSpace(int screen) {
        return ps[screen - 1];
    }
};

void getOpParameters(std::vector<std::string>& settingsParameters, int& mode, double& rate, double& timeout) {
    int N = (int)settingsParameters.size();
    double buffer = 0.0;
    for (int i = 0; i < N; i++) {
        if (settingsParameters[i].find("mode") != std::string::npos) {
            size_t indexStart = settingsParameters[i].find("=") + 1;
            std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
            buffer = std::stod(strValue.c_str());
            mode = buffer;
        }
        if (settingsParameters[i].find("refreshRate") != std::string::npos) {
            size_t indexStart = settingsParameters[i].find("=") + 1;
            std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
            rate = std::stod(strValue.c_str());
        }
        if (settingsParameters[i].find("timeout") != std::string::npos) {
            size_t indexStart = settingsParameters[i].find("=") + 1;
            std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
            timeout = std::stod(strValue.c_str());
        }
    }
}

#endif