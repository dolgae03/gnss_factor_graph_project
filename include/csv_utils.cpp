#include "csv_utils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include <limits>
#include <random>

using namespace std;


bool isValidInteger(const std::string& str) {
    std::istringstream iss(str);
    int value;
    return (iss >> value) && (iss.eof());
}

bool isValidDouble(const std::string& str) {
    if (str == "NaN") {
        return true;
    }
    std::istringstream iss(str);
    double value;
    return (iss >> value) && (iss.eof());
}

std::vector<double> generateIMUdata(size_t numEpochs, double sigma) {
    std::vector<double> data(numEpochs);

    // 난수 생성 엔진 및 분포 설정
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, sigma);  // 평균 0, 표준 편차 sigma

    // 데이터 생성
    for (int i = 0; i < numEpochs; ++i) {
        data[i] = d(gen);
    }

    return data;
}



std::vector<std::pair<int, double>> readGpsTimeCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<int, double>> data;
    std::string line;

    if (file.is_open()) {
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::stringstream lineStream(line);
            std::string cell;
            int int_value;
            double double_value;

            if (std::getline(lineStream, cell, ',')) {
                int_value = std::stoi(cell);
            }

            if (std::getline(lineStream, cell, ',')) {
                double_value = (cell == "NaN") ? std::numeric_limits<double>::quiet_NaN() : std::stod(cell);
            }

            data.emplace_back(int_value, double_value);
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return data;
}

std::vector<std::vector<std::vector<double>>> readSVPosAndVelCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<std::vector<double>>> data;
    std::string line;

    if (file.is_open()) {
        std::getline(file, line);  // 헤더 라인 무시

        while (std::getline(file, line)) {
            std::vector<std::vector<double>> epoch_data;
            std::stringstream lineStream(line);
            std::string cell;
            int satellite_index = 0;

            while (std::getline(lineStream, cell, ',')) {
                if (satellite_index % 3 == 0) {
                    // 새로운 위성 데이터를 시작합니다.
                    epoch_data.push_back(std::vector<double>(3));
                }
                epoch_data.back()[satellite_index % 3] = (cell == "NaN") ? std::numeric_limits<double>::quiet_NaN() : std::stod(cell);
                satellite_index++;
            }
            data.push_back(epoch_data);
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return data;
}

std::vector<std::vector<double>> readPseudorangeCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line;

    if (file.is_open()) {
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::vector<double> epoch_data;
            std::stringstream lineStream(line);
            std::string cell;

            while (std::getline(lineStream, cell, ',')) {
                epoch_data.push_back((cell == "NaN") ? std::numeric_limits<double>::quiet_NaN() : std::stod(cell));
            }

            data.push_back(epoch_data);
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return data;
}


std::vector<std::vector<double>> readStateCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line;

    if (file.is_open()) {
        // Skip the header line
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::vector<double> epoch_data;
            std::stringstream lineStream(line);
            std::string cell;
            int colIndex = 0;

            // Read each cell separated by commas and only capture the first 5 columns
            while (getline(lineStream, cell, ',') && colIndex < 7) {
                double value;
                if (cell == "NaN") {
                    value = std::numeric_limits<double>::quiet_NaN();
                } else {
                    value = stod(cell);
                }
                epoch_data.push_back(value);
                ++colIndex;
            }

            // Add parsed line to data vector if it has exactly 7 columns
            if (epoch_data.size() == 7) {
                data.push_back(epoch_data);
            }
        }

        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }

    return data;
}







void printSVPosAndVelCSV(const std::vector<std::vector<std::vector<double>>>& data) {
    for (const auto& epoch : data) {
        for (const auto& satellite : epoch) {
            std::cout << "x: " << satellite[0] << ", y: " << satellite[1] << ", z: " << satellite[2] << " | ";
        }
        std::cout << std::endl;
        break;
    }
}

void printPseudorangeCSV(const std::vector<std::vector<double>>& data) {
    for (const auto& epoch : data) {
        for (const auto& pr_value : epoch) {
            std::cout << pr_value << " ";
        }
        std::cout << std::endl;
        break;
    }
}

void printGpsTimeCSV(const std::vector<std::pair<int, double>>& data) {
    for (const auto& row : data) {
        std::cout << row.first << ", " << row.second << std::endl;
        break;
    }
}

void checkCSVData(const std::vector<std::vector<double>>& pr_data,
                  const std::vector<std::vector<std::vector<double>>>& sv_pos_data,
                  const std::vector<std::vector<std::vector<double>>>& sv_vel_data,
                  const std::vector<std::pair<int, double>>& time_data){ 

    std::cout << "PR Data:" << std::endl;
    printPseudorangeCSV(pr_data);

    std::cout << "SV Position Data:" << std::endl;
    printSVPosAndVelCSV(sv_pos_data);

    std::cout << "SV Velocity Data:" << std::endl;
    printSVPosAndVelCSV(sv_vel_data);

    std::cout << "Time Data:" << std::endl;
    printGpsTimeCSV(time_data);
}