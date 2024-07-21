#ifndef CSV_UTILS_H
#define CSV_UTILS_H

#include <vector>
#include <string>

// 유효한 정수인지 확인하는 함수
bool isValidInteger(const std::string& str);

// 유효한 실수인지 확인하는 함수
bool isValidDouble(const std::string& str);

// pr.csv 파일을 읽고 데이터를 벡터로 반환하는 함수
std::vector<std::pair<int, double>> readGpsTimeCSV(const std::string& filename);

// sv_pos.csv 및 sv_vel.csv 파일을 읽고 3차원 벡터로 반환하는 함수
std::vector<std::vector<std::vector<double>>> readSVPosAndVelCSV(const std::string& filename);

// pseudorange 데이터를 읽고 벡터로 반환하는 함수
std::vector<std::vector<double>> readPseudorangeCSV(const std::string& filename);

// sv_pos.csv 및 sv_vel.csv 데이터를 출력하는 함수
void printSVPosAndVelCSV(const std::vector<std::vector<std::vector<double>>>& data);

// pseudorange 데이터를 출력하는 함수
void printPseudorangeCSV(const std::vector<std::vector<double>>& data);

// time.csv 데이터를 출력하는 함수
void printGpsTimeCSV(const std::vector<std::pair<int, double>>& data);

// CSV 데이터를 확인하는 함수
void checkCSVData(const std::vector<std::vector<double>>& pr_data,
                  const std::vector<std::vector<std::vector<double>>>& sv_pos_data,
                  const std::vector<std::vector<std::vector<double>>>& sv_vel_data,
                  const std::vector<std::pair<int, double>>& time_data);

#endif // CSV_UTILS_H