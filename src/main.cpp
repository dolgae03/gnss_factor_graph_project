
#include <gtsam/geometry/Rot2.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/dataset.h>
#include <gtsam/slam/BetweenFactor.h>
#include <fstream>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>

#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/inference/Key.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/base/OptionalJacobian.h>
#include <gtsam/base/Matrix.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/optional.hpp>

using namespace std;
using namespace gtsam;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>

// 유효한 정수인지 확인하는 함수
bool isValidInteger(const std::string& str) {
    std::istringstream iss(str);
    int value;
    return (iss >> value) && (iss.eof());
}

// 유효한 실수인지 확인하는 함수
bool isValidDouble(const std::string& str) {
    if (str == "NaN") {
        return true;
    }
    std::istringstream iss(str);
    double value;
    return (iss >> value) && (iss.eof());
}

// pr.csv 파일을 읽고 데이터를 벡터로 반환하는 함수
std::vector<std::pair<int, double>> readGpsTimeCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::pair<int, double>> data;
    std::string line;

    if (file.is_open()) {
        // 첫 줄은 헤더이므로 건너뛰기
        std::getline(file, line);

        while (std::getline(file, line)) {
            std::stringstream lineStream(line);
            std::string cell;
            int int_value;
            double double_value;

            // 첫 번째 열 (int)
            if (std::getline(lineStream, cell, ',')) {
                // std::cout << cell << std::endl;
                int_value = std::stoi(cell);
            }

            // 두 번째 열 (double)
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

// sv_pos.csv 및 sv_vel.csv 파일을 읽고 3차원 벡터로 반환하는 함수
std::vector<std::vector<std::vector<double>>> readSVPosAndVelCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<std::vector<double>>> data;
    std::string line;

    if (file.is_open()) {
        // 첫 줄은 헤더이므로 건너뛰기
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

// pseudorange 데이터를 읽고 벡터로 반환하는 함수
std::vector<std::vector<double>> readPseudorangeCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line;

    if (file.is_open()) {
        // 첫 줄은 헤더이므로 건너뛰기
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


// sv_pos.csv 및 sv_vel.csv 데이터를 출력하는 함수
void printSVPosAndVelCSV(const std::vector<std::vector<std::vector<double>>>& data) {
    for (const auto& epoch : data) {
        for (const auto& satellite : epoch) {
            std::cout << "x: " << satellite[0] << ", y: " << satellite[1] << ", z: " << satellite[2] << " | ";
        }
        std::cout << std::endl;
        break;
    }
}

// pseudorange 데이터를 출력하는 함수
void printPseudorangeCSV(const std::vector<std::vector<double>>& data) {
    for (const auto& epoch : data) {
        for (const auto& pr_value : epoch) {
            std::cout << pr_value << " ";
        }
        std::cout << std::endl;
        break;
    }
}

// time.csv 데이터를 출력하는 함수
void printGpsTimeCSV(const std::vector<std::pair<int, double>>& data) {
    for (const auto& row : data) {
        std::cout << row.first << ", " << row.second << std::endl;
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


class GpsFactor: public NoiseModelFactor1<Point3> {
    double pseudorange_;
    Point3 satellite_position_;

    public:
    GpsFactor(Key j, double pseudorange, const Point3& satellite_position, const SharedNoiseModel& model)
        : NoiseModelFactor1<Point3>(model, j), pseudorange_(pseudorange), satellite_position_(satellite_position) {}

    Vector evaluateError(const Point3& user_position, boost::optional<Matrix&> H = boost::none) const {
        // Calculate the distance between the user position and the satellite position
        double dx = user_position.x() - satellite_position_.x();
        double dy = user_position.y() - satellite_position_.y();
        double dz = user_position.z() - satellite_position_.z();
        double predicted_range = std::sqrt(dx * dx + dy * dy + dz * dz);

        // If the Jacobian matrix H is provided, calculate it
        if (H) {
        Matrix13 Hx;
        Hx << (user_position.x() - satellite_position_.x()) / predicted_range,
                (user_position.y() - satellite_position_.y()) / predicted_range,
                (user_position.z() - satellite_position_.z()) / predicted_range;
        *H = Hx;
        }

        // Return the difference between the predicted range and the pseudorange measurement
        return (Vector(1) << predicted_range - pseudorange_).finished();
    }
};


int main() {
    // 데이터 파일 경로
    std::string pr_file = "../data/pr.csv";
    std::string sv_pos_file = "../data/sv_pos.csv";
    std::string sv_vel_file = "../data/sv_vel.csv";
    std::string time_file = "../data/time.csv";

    // CSV 파일 읽기
    std::vector<std::vector<double>> pr_data = readPseudorangeCSV(pr_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data = readSVPosAndVelCSV(sv_pos_file);
    std::vector<std::vector<std::vector<double>>> sv_vel_data = readSVPosAndVelCSV(sv_vel_file);
    std::vector<std::pair<int, double>> time_data = readGpsTimeCSV(time_file);

    checkCSVData(pr_data, sv_pos_data, sv_vel_data, time_data);

    // 노이즈 모델 설정 (예제)
    auto noise = noiseModel::Isotropic::Sigma(1, 3.0); // 표준편차 1.0인 등방성 노이즈 모델

    // 팩터 그래프 생성
    NonlinearFactorGraph graph;

    // GpsFactor 추가
    for (size_t epoch = 0; epoch < pr_data.size(); ++epoch) {
        for (size_t satellite = 0; satellite < pr_data[epoch].size(); ++satellite) {
            double pr_value = pr_data[epoch][satellite];

            // NaN 값인지 확인
            if (std::isnan(pr_value)) continue;

            // 위성 위치 가져오기
            const auto& sv_position = sv_pos_data[epoch][satellite];
            Point3 sv_loc(sv_position[0], sv_position[1], sv_position[2]);

            // 팩터 그래프에 GpsFactor 추가
            graph.add(gtsam::make_shared<GpsFactor>(0, pr_value, sv_loc, noise));
        }
    }

    // 초기 사용자 위치 설정 (예제)
    Values initial;
    initial.insert(0, Point3(0, 0, 0)); // 초기 위치를 (0, 0, 0)으로 설정

    // 최적화 수행
    LevenbergMarquardtOptimizer optimizer(graph, initial);
    Values result = optimizer.optimize();

    // 결과 출력
    Point3 estimated_position = result.at<Point3>(0);
    std::cout << "Estimated user position: " << estimated_position << std::endl;

    return 0;
}
