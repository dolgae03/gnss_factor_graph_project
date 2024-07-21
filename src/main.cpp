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


int main(int argc, char** argv) {
    std::string pr_file = "../data/pr.csv";
    std::string sv_pos_file = "../data/sv_pos.csv";
    std::string sv_vel_file = "../data/sv_vel.csv";
    std::string time_file = "../data/time.csv";

    // CSV 파일 읽기
    std::vector<std::vector<double>> pr_data = readPseudorangeCSV(pr_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data = readSVPosAndVelCSV(sv_pos_file);
    std::vector<std::vector<std::vector<double>>> sv_vel_data = readSVPosAndVelCSV(sv_vel_file);
    std::vector<std::pair<int, double>> time_data = readGpsTimeCSV(time_file);

    std::vector<double * > state;
    double x[3] = {100000.0, 100000.0, 100000.0};
    const double initial_x[3] = {x[0], x[1], x[2]};

    ceres::Problem problem;

    const size_t max_epoch = time_data.size();

    for(size_t epoch=0; epoch < max_epoch; ++epoch){
        double* current_position = new double[3];

        for(size_t satellite=0; satellite < pr_data[epoch].size(); ++satellite){
            double pr_value = pr_data[epoch][satellite];

            if (std::isnan(pr_value)) continue;

            factor::PesudorangeFactorCostFunctor* functor = 
                new factor::PesudorangeFactorCostFunctor(sv_pos_data[epoch][satellite], pr_value);

            ceres::CostFunction* cost_function =
                new ceres::AutoDiffCostFunction<factor::PesudorangeFactorCostFunctor, 1, 3>(functor);

            problem.AddResidualBlock(cost_function, nullptr, current_position);
        }
        state.push_back(current_position);
    }
    // Solver 옵션 설정 및 실행
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    for (size_t epoch = 0; epoch < max_epoch; ++epoch) {
        double* position = state[epoch];
        std::cout << "Epoch(ECEF) " << epoch << ": " << position[0] << ", " << position[1] << ", " << position[2] << ")\n";
        
        // std::vector<double> ecef_position = {position[0], position[1], position[2]};

        std::vector<double> ecef_position = {position[0], position[1], position[2]};
        std::vector<double> res = coordinate::ecefToLLA(ecef_position);

        std::cout << "Epoch(LLA) " << epoch << ": " << res[0] << ", " << res[1] << ", " << res[2] << ")\n";

        delete[] position; // 메모리 해제
    }


    return 0;
}