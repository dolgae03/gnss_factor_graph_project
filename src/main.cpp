#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include "ceres/ceres.h"
#include "../include/factor.h"
#include "../include/csv_utils.h"
#include "../include/cord_utils.h"


template <typename T>
bool check_sv_data(T vector){
    for(int i=0; i<vector.size(); i++)
        if(std::isnan(vector[i]))
            return true;
    
    return false;
}


int main(int argc, char** argv) {
    std::string rover_dir = "../data/data_rover/";
    std::string station_dir = "../data/data_station/";
    std::string pr_file = "pr.csv";
    std::string dop_file = "dop.csv";
    std::string sv_pos_file = "sv_pos.csv";
    std::string sv_vel_file = "sv_vel.csv";
    std::string time_file = "time.csv";

    // CSV 파일 읽기
    std::vector<std::vector<double>> pr_data = readPseudorangeCSV(rover_dir + pr_file);
    std::vector<std::vector<double>> dop_data = readPseudorangeCSV(rover_dir + dop_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data = readSVPosAndVelCSV(rover_dir + sv_pos_file);
    std::vector<std::vector<std::vector<double>>> sv_vel_data = readSVPosAndVelCSV(rover_dir + sv_vel_file);
    std::vector<std::pair<int, double>> time_data = readGpsTimeCSV(rover_dir + time_file);

    std::vector<std::vector<double>> pr_data_station = readPseudorangeCSV(rover_dir + pr_file);
    std::vector<std::vector<double>> dop_data_station = readPseudorangeCSV(rover_dir + dop_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data_station = readSVPosAndVelCSV(rover_dir + sv_pos_file);
    std::vector<std::vector<std::vector<double>>> sv_vel_data_station = readSVPosAndVelCSV(rover_dir + sv_vel_file);
    std::vector<std::pair<int, double>> time_data_station = readGpsTimeCSV(rover_dir + time_file);

    std::vector<double * > state;
    ceres::Problem problem;

    std::vector<double> ref_location = coordinate::lla2ecef({36.372371713580250, 127.358800510185191, 91.642377777777796});

    const size_t val_num = 4;
    const size_t max_epoch = 10; //time_data.size();

    for(size_t epoch=0; epoch < max_epoch; ++epoch){
        double* current_position = new double[val_num];
        std::fill(current_position, current_position + val_num, 0.0); 

        for(size_t satellite=0; satellite < pr_data[epoch].size(); ++satellite){
            double pr_value = pr_data[epoch][satellite];
            double pr_value_station = pr_data_station[epoch][satellite];

            if (std::isnan(pr_value) || std::isnan(pr_value_station)) continue;
            if (check_sv_data(sv_pos_data[epoch][satellite])) continue;

            std::cout << "pr_value: " << pr_value << std::endl;

            factor::DiffPesudorangeFactorCostFunctor* functor = 
                new factor::DiffPesudorangeFactorCostFunctor(ref_location, sv_pos_data[epoch][satellite], pr_value-pr_value_station);

            ceres::CostFunction* cost_function =
                new ceres::AutoDiffCostFunction<factor::DiffPesudorangeFactorCostFunctor, 1, val_num>(functor);

            problem.AddResidualBlock(cost_function, nullptr, current_position);
        }

        if(epoch < 0){
            for(size_t satellite=0; satellite < pr_data[epoch].size(); ++satellite){
                double prev_dop_value = dop_data[epoch-1][satellite]; // Hz
                double next_dop_value = dop_data[epoch][satellite];
                if (std::isnan(prev_dop_value) || std::isnan(next_dop_value)) continue;
                if (check_sv_data(sv_vel_data[epoch][satellite]) || check_sv_data(sv_vel_data[epoch-1][satellite])) continue;

                double mean_dop_value = (prev_dop_value + next_dop_value) / 2;
                double time_interval = time_data[epoch].second - time_data[epoch-1].second;

                std::vector<double> sv_avg_vel_data;
                for(int i=0; i<3; i++)
                    sv_avg_vel_data.push_back((sv_vel_data[epoch][satellite][i] + sv_vel_data[epoch-1][satellite][i])/2);

                factor::DopplerFactorCostFunctor* functor = 
                    new factor::DopplerFactorCostFunctor(sv_avg_vel_data, mean_dop_value, time_interval);

                ceres::CostFunction* cost_function = 
                    new ceres::AutoDiffCostFunction<factor::DopplerFactorCostFunctor, 1, val_num, val_num>(functor);

                problem.AddResidualBlock(cost_function, nullptr, state[epoch-1], current_position);
            }
            std::cout << epoch << endl;
        }

        state.push_back(current_position);
    }
    // Solver 옵션 설정 및 실행
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.FullReport() << "\n";
    ceres::Covariance::Options cov_options;
    ceres::Covariance covariance(cov_options);

    std::vector<std::pair<const double*, const double*>> covariance_blocks;
    for (double* position : state) {
        covariance_blocks.emplace_back(position, position);
    }

    if (covariance.Compute(covariance_blocks, &problem)) {
        for (size_t epoch = 0; epoch < max_epoch; ++epoch) {
            double* position = state[epoch];
            std::cout << "Epoch(ECEF) " << epoch << ": " << position[0] << ", " << position[1] << ", " << position[2] << ", " << position[3] << "\n";
            
            // std::vector<double> ecef_position = {position[0], position[1], position[2]};

            std::vector<double> ecef_position = {position[0], position[1], position[2]};
            std::vector<double> res = coordinate::ecef2lla(ecef_position);

            std::cout << "Epoch(LLA) " << epoch << ": " << res[0] << ", " << res[1] << ", " << res[2] << "\n";

            // 공분산 출력
            double covariance_matrix[3 * 3];
            covariance.GetCovarianceBlock(position, position, covariance_matrix);
            std::cout << "Covariance Matrix:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    std::cout << covariance_matrix[3 * i + j] << " ";
                }
                std::cout << "\n";
            }

            delete[] position; // 메모리 해제
        }
    } else {
        std::cerr << "Failed to compute covariance." << std::endl;
    }


    return 0;
}