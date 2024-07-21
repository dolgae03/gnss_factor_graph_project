#include "csv_utils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include <limits>

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
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::vector<std::vector<double>> epoch_data(32, std::vector<double>(3));
            std::stringstream lineStream(line);
            std::string cell;
            int satellite_index = 0;

            while (std::getline(lineStream, cell, ',')) {
                epoch_data[satellite_index / 3][satellite_index % 3] = (cell == "NaN") ? std::numeric_limits<double>::quiet_NaN() : std::stod(cell);
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