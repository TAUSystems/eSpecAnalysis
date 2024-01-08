#pragma once

#ifndef __fileIO_h__
#define __fileIO_h__

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#ifdef _WIN32
    #include <filesystem>
#else
    #include <unistd.h>
#endif
#include <string>
#include <vector>
#include <ctime>

#include <omp.h>

#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>

template<int I>
struct CvType {};
template<>
struct CvType<CV_64F> { typedef double type_t; };

class imageBW {
    double** data;
    size_t* size;
public:
    imageBW() {
        data = (double**)malloc(sizeof(double*));
        data[0] = (double*)malloc(sizeof(double));
        data[0][0] = 0.0;
        size = (size_t*)malloc(sizeof(size_t) * 2);
        size[0] = 1;
        size[1] = 1;
    }

    void resize(int Nx, int Ny) {
        size[0] = (size_t)Nx;
        size[1] = (size_t)Ny;
        data = (double**)malloc(sizeof(double*) * size[0]);

        #pragma omp parallel for
            for (int i = 0; i < (int)size[0]; i++) {
                data[i] = (double*)malloc(sizeof(double) * size[1]);
                for (int j = 0; j < (int)size[1]; j++) {
                    data[i][j] = 0.0;
                }
            }
    }

    void resize(size_t Nx, size_t Ny) {
        size[0] = Nx;
        size[1] = Ny;
        data = (double**)malloc(sizeof(double*) * size[0]);

        #pragma omp parallel for
            for (int i = 0; i < (int)size[0]; i++) {
                data[i] = (double*)malloc(sizeof(double) * size[1]);
                for (int j = 0; j < (int)size[1]; j++) {
                    data[i][j] = 0.0;
                }
            }
    }

    void resize(std::vector<size_t> winSize) {
        size[0] = winSize[0];
        size[1] = winSize[1];
        data = (double**)malloc(sizeof(double*) * size[0]);

        #pragma omp parallel for
            for (int i = 0; i < (int)size[0]; i++) {
                data[i] = (double*)malloc(sizeof(double) * size[1]);
                for (int j = 0; j < (int)size[1]; j++) {
                    data[i][j] = 0.0;
                }
            }
    }

    void resize(cv::Size imgSize) {
        size[0] = imgSize.width;
        size[1] = imgSize.height;
        data = (double**)malloc(sizeof(double*) * size[0]);

        #pragma omp parallel for
            for (int i = 0; i < (int)size[0]; i++) {
                data[i] = (double*)malloc(sizeof(double) * size[1]);
                for (int j = 0; j < (int)size[1]; j++) {
                    data[i][j] = 0.0;
                }
            }
    }

    void resize(std::vector<cv::Point2d> bounds) {
        int NsX, NfX, NsY, NfY;
        NsX = std::min(std::min((int)bounds[0].x, (int)bounds[1].x), std::min((int)bounds[2].x, (int)bounds[3].x));
        NfX = std::max(std::max((int)bounds[0].x, (int)bounds[1].x), std::max((int)bounds[2].x, (int)bounds[3].x));
        NsY = std::min(std::min((int)bounds[0].y, (int)bounds[1].y), std::min((int)bounds[2].y, (int)bounds[3].y));
        NfY = std::max(std::max((int)bounds[0].y, (int)bounds[1].y), std::max((int)bounds[2].y, (int)bounds[3].y));
        size[0] = (size_t)(NfX - NsX);
        size[1] = (size_t)(NfY - NsY);
        data = (double**)malloc(sizeof(double*) * size[0]);

        #pragma omp parallel for
            for (int i = 0; i < (int)size[0]; i++) {
                data[i] = (double*)malloc(sizeof(double) * size[1]);
                for (int j = 0; j < (int)size[1]; j++) {
                    data[i][j] = 0.0;
                }
            }
    }

    void destroy() {
        for (int i = 0; i < size[0]; i++) {
            free(data[i]);
        }
        free(data);
        free(size);

        data = (double**)malloc(sizeof(double*));
        data[0] = (double*)malloc(sizeof(double));
        data[0][0] = 0.0;
        size = (size_t*)malloc(sizeof(size_t) * 2);
        size[0] = 1;
        size[1] = 1;
    }

    void definePixel(int indexX, int indexY, double value) {
        data[indexX][indexY] = value;
    }

    void definePixel(size_t indexX, size_t indexY, double value) {
        data[indexX][indexY] = value;
    }

    void definePixel(int indexX, size_t indexY, double value) {
        data[indexX][indexY] = value;
    }

    void definePixel(size_t indexX, int indexY, double value) {
        data[indexX][indexY] = value;
    }

    double value(int indexX, int indexY) {
        return data[indexX][indexY];
    }

    size_t sizeX() {
        return size[0];
    }

    size_t sizeY() {
        return size[1];
    }

    void crop(std::vector<size_t> bounds, imageBW& output) {
        int NsX, NfX, NsY, NfY;

        if (bounds[1] > size[0]) {
            NfX = (int)size[0];
        }
        else {
            NfX = (int)bounds[1];
        }
        if (bounds[0] < 0) {
            NsX = 0;
        }
        else {
            NsX = (int)bounds[0];
        }
        if (bounds[3] > size[1]) {
            NfY = (int)size[0];
        }
        else {
            NfY = (int)bounds[3];
        }
        if (bounds[2] < 0) {
            NsY = 0;
        }
        else {
            NsY = (int)bounds[2];
        }
        int Nx = NfX - NsX;
        int Ny = NfY - NsY;
        output.resize(Nx, Ny);

        #pragma omp parallel for
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    int Ni = (int)(NsX + i);
                    int Nj = (int)(NsY + j);
                    output.definePixel(i, j, data[Ni][Nj]);
                }
            }
    }

    void crop(std::vector<int> bounds, imageBW& output) {
        int NsX, NfX, NsY, NfY;

        if (bounds[1] > size[0]) {
            NfX = (int)size[0];
        }
        else {
            NfX = (int)bounds[1];
        }
        if (bounds[0] < 0) {
            NsX = 0;
        }
        else {
            NsX = (int)bounds[0];
        }
        if (bounds[3] > size[1]) {
            NfY = (int)size[0];
        }
        else {
            NfY = (int)bounds[3];
        }
        if (bounds[2] < 0) {
            NsY = 0;
        }
        else {
            NsY = (int)bounds[2];
        }
        int Nx = NfX - NsX;
        int Ny = NfY - NsY;
        output.resize(Nx, Ny);

        
        #pragma omp parallel for
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    int Ni = (int)(NsX + i);
                    int Nj = (int)(NsY + j);
                    output.definePixel(i, j, data[Ni][Nj]);
                }
            }
    }
};

std::string getDate() {
    std::time_t timeRaw;
    std::tm timeInfo;
    char buffer[80];
    std::time(&timeRaw);
    localtime_s(&timeInfo, &timeRaw);
    std::strftime(buffer, 80, "%Y%m%d", &timeInfo);

    std::string dateString(buffer);
    return dateString;
}

std::string getPath() {
    std::string output;
    #ifdef _WIN32
        std::filesystem::path buffer = std::filesystem::current_path();
        output = buffer.string();
    #else
        char buffer[1024];
        getcwd(buffer, sizeof(buffer));
        output = buffer;
    #endif
    return output;
}

void listDir(std::string& pathDir, std::vector<std::string>& list) {
    list.clear();
    std::string file;
    bool isIMG;
    #ifdef _WIN32
        for (const auto& entry : std::filesystem::directory_iterator(pathDir)) {
            isIMG = false;
            file = entry.path().generic_string();
            if (file.find(".tiff") != std::string::npos ||
                file.find(".tif") != std::string::npos ||
                file.find(".Tiff") != std::string::npos || 
                file.find(".Tif") != std::string::npos || 
                file.find(".TIFF") != std::string::npos || 
                file.find(".TIF") != std::string::npos) {
                isIMG = true;
            }
            if (isIMG) {
                list.push_back(file);
            }
        }
    #else
        DIR* dp = opendir(pathDir);
        int dfd = dirfd(dp);
        struct dirent* dirp;
        while (readdir(dp) != NULL) {
            struct stat sb;
            fstatat(dfd, dirp->d_name, &sb, 0);
            file = pathDir + "/";
            file = file + dirp->d_name;
            isIMG = false;
            file = entry.path().generic_string();
            if (file.find(".tiff") != std::string::npos ||
                file.find(".tif") != std::string::npos ||
                file.find(".Tiff") != std::string::npos ||
                file.find(".Tif") != std::string::npos ||
                file.find(".TIFF") != std::string::npos ||
                file.find(".TIF") != std::string::npos) {
                isIMG = true;
            }
            if (isIMG) {
                list.push_back(file);
            }
        }
    #endif
}

void findRef(std::vector<std::string>& list, std::string& pathRef){
    pathRef = "-1NoRef";
    size_t N = list.size();

    for (size_t i = 0; i < N; i++) {
        if (list[i].find("refimage") != std::string::npos ||
            list[i].find("refImage") != std::string::npos || 
            list[i].find("refIMAGE") != std::string::npos || 
            list[i].find("Refimage") != std::string::npos || 
            list[i].find("RefImage") != std::string::npos || 
            list[i].find("RefIMAGE") != std::string::npos || 
            list[i].find("REFimage") != std::string::npos || 
            list[i].find("REFImage") != std::string::npos || 
            list[i].find("REFIMAGE") != std::string::npos || 
            list[i].find("refimg") != std::string::npos || 
            list[i].find("refImg") != std::string::npos || 
            list[i].find("refIMG") != std::string::npos || 
            list[i].find("Refimg") != std::string::npos || 
            list[i].find("RefImg") != std::string::npos || 
            list[i].find("RefIMG") != std::string::npos || 
            list[i].find("REFimg") != std::string::npos || 
            list[i].find("REFImg") != std::string::npos || 
            list[i].find("REFIMG") != std::string::npos ) {
            if (list[i].find(".tiff") != std::string::npos ||
                list[i].find(".tif") != std::string::npos ||
                list[i].find(".Tiff") != std::string::npos ||
                list[i].find(".Tif") != std::string::npos ||
                list[i].find(".TIFF") != std::string::npos ||
                list[i].find(".TIF") != std::string::npos) {
                pathRef = list[i];
                break;
            }
        }
    }
}

void readFile(std::string& pathFile, std::vector<std::string>& fileContent) {
    std::ifstream file;
    file.open(pathFile);
    fileContent.clear();
    std::string buffer;
    if (file.is_open()) {
        while (file) {
            std::getline(file, buffer);
            fileContent.push_back(buffer);
        }
    }
    size_t N = fileContent.size() - 1;
    fileContent.resize(N);
}

void getOpParameters(std::vector<std::string>& settingsParameters, int& mode, double& rate, double& timeout) {
    int N = (int)settingsParameters.size();
    double buffer = 0.0;
    for (int i = 0; i < N; i++) {
        if (settingsParameters[i].find("mode") != std::string::npos) {
            size_t indexStart = settingsParameters[i].find("=") + 1;
            std::string strValue = settingsParameters[i].substr(indexStart, settingsParameters[i].length() - indexStart);
            buffer = std::stod(strValue.c_str());
            mode = (int)buffer;
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

void scanNewFile(std::string path, std::vector<std::string>& refList, std::vector<int>& updateStatus) {
    std::vector<std::string> updateList;
    listDir(path, updateList);
    int nRef, nUpdate;
    nRef = (int)refList.size();
    nUpdate = (int)updateList.size();
    updateStatus.clear();
    
    if (nUpdate > nRef) {
        updateStatus.resize(nUpdate, 1);
        #pragma omp parallel for
            for (int i = 0; i < nUpdate; i++) {
                for (int j = 0; j < nRef; j++) {
                    if (updateList[i].compare(refList[j]) == 0) {
                        updateStatus[i] = 0;
                        break;
                    }
                }
            }
        refList.clear();
        refList = updateList;
    }
    else {
        updateStatus.resize(nUpdate, 0);
        refList.clear();
        refList = updateList;
    }
}

bool findFile(std::vector<std::string> list, std::string file, std::string timeStamp, int &index) {
    bool output = false;
    int N = (int)list.size();
    for (int i = 0; i < N; i++) {
        if (list[N - 1 - i].find(file) != std::string::npos) {
            if (list[N - 1 - i].find(timeStamp) != std::string::npos) {
                index = N - 1 - i;
                output = true;
                break;
            }
        }
    }
    return output;
}

void getImage(std::string file, imageBW& output) {
    cv::Mat buffer = cv::imread(file, cv::IMREAD_ANYDEPTH);
    int bitDepth;
    if (buffer.empty()) {
        std::cout << "Image Not Found: " << file << "\n";
    }

    cv::Size imgSize = buffer.size();

    if (buffer.type() == 2) {
        bitDepth = 16;
        cv::Mat buffer64f = cv::Mat::zeros(imgSize.height, imgSize.width, CV_64F);
        buffer.convertTo(buffer64f, CV_64F);

        output.resize(imgSize);
        #pragma omp parallel for
            for (int i = 0; i < imgSize.width; i++) {
                for (int j = 0; j < imgSize.height; j++) {
                    double pxValue = buffer64f.at<double>(j, i) / pow(2, bitDepth);
                    output.definePixel(i, j, pxValue);
                }
            }
    }
    else if (buffer.type() == 6) {
        bitDepth = 64;
        output.resize(imgSize);
        #pragma omp parallel for
            for (int i = 0; i < imgSize.width; i++) {
                for (int j = 0; j < imgSize.height; j++) {
                    double pxValue = buffer.at<double>(j, i);
                    output.definePixel(i, j, pxValue);
                }
            }
    }
    else{
        bitDepth = 8;
        cv::Mat buffer64f = cv::Mat::zeros(imgSize.height, imgSize.width, CV_64F);
        buffer.convertTo(buffer64f, CV_64F);

        output.resize(imgSize);
        #pragma omp parallel for
            for (int i = 0; i < imgSize.width; i++) {
                for (int j = 0; j < imgSize.height; j++) {
                    double pxValue = buffer64f.at<double>(j, i) / pow(2, bitDepth);
                    output.definePixel(i, j, pxValue);
                }
            }
    }
    
}

void getImage(cv::Mat input, imageBW& output) {
    int bitDepth;
    if (input.empty()) {
        std::cout << "Image Not Found.\n";
    }

    cv::Size imgSize = input.size();

    if (input.type() == 2) {
        bitDepth = 16;
        cv::Mat buffer64f = cv::Mat::zeros(imgSize.height, imgSize.width, CV_64F);
        input.convertTo(buffer64f, CV_64F);

        output.resize(imgSize);
        #pragma omp parallel for
            for (int i = 0; i < imgSize.width; i++) {
                for (int j = 0; j < imgSize.height; j++) {
                    double pxValue = buffer64f.at<double>(j, i) / pow(2, bitDepth);
                    output.definePixel(i, j, pxValue);
                }
            }
    }
    else if (input.type() == 6) {
        bitDepth = 64;
        output.resize(imgSize);
        #pragma omp parallel for
            for (int i = 0; i < imgSize.width; i++) {
                for (int j = 0; j < imgSize.height; j++) {
                    double pxValue = input.at<double>(j, i);
                    output.definePixel(i, j, pxValue);
                }
            }
    }
    else {
        bitDepth = 8;
        cv::Mat buffer64f = cv::Mat::zeros(imgSize.height, imgSize.width, CV_64F);
        input.convertTo(buffer64f, CV_64F);

        output.resize(imgSize);
        #pragma omp parallel for
            for (int i = 0; i < imgSize.width; i++) {
                for (int j = 0; j < imgSize.height; j++) {
                    double pxValue = buffer64f.at<double>(j, i) / pow(2, bitDepth);
                    output.definePixel(i, j, pxValue);
                }
            }
    }

}

void imgshow(imageBW image) {
    int Nx = (int)image.sizeX();
    int Ny = (int)image.sizeY();
    cv::Mat buffer(Ny, Nx, CV_64F);

    #pragma omp parallel for
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double pxBuffer = image.value(i, j);
                buffer.at<CvType<CV_64F>::type_t>(j, i) = pxBuffer;
            }
        }

    cv::namedWindow(" ", cv::WINDOW_AUTOSIZE);
    cv::imshow(" ", buffer);
    cv::waitKey(0);
}

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

        std::string rootPath;
        #pragma omp parallel for
            for (int i = 0; i < N; i++) {
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
            }
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
        std::string fileCalibration = pathCalibration + "/perspective.cal";
        std::vector<std::string> calibration;
        readFile(fileCalibration, calibration);


        int N = (int)calibration.size();

        #pragma omp parallel for
            for (int i = 0; i < N; i++) {
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
            }

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


        filePath = calPath + "/eScreenA.map";
        readFile(filePath, buffer);
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

        filePath = calPath + "/eScreenB.map";
        readFile(filePath, buffer);
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

#endif