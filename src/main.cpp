#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include <map>

#include "ceres/ceres.h"
#include "../include/factors.h"
#include "../include/csv_utils.h"
#include "../include/cord_utils.h"

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;


#define DF_PR_WEIGHT (double)1.0
#define TDCP_WEIGHT (double)1.0e-4
#define CONSATANT_CLOCK_WEIGHT (double)1.0e-5


template <typename T>
bool check_sv_data(T vector){
    for(int i=0; i<vector.size(); i++)
        if(std::isnan(vector[i]))
            return true;
    
    return false;
}

bool parseCommandLineOptions(int argc, char* argv[], 
                             bool& use_df_pr, 
                             bool& use_tdcp, 
                             bool& use_clock_const, 
                             double& df_pr_weight, 
                             double& tdcp_weight, 
                             double& clock_const_weight,
                             size_t& start_epoch,
                             size_t& T,
                             std::set<int>& constellation_type,
                             std::string& constellation_name) {
    try {
        std::vector<std::string> constellations;

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("disable-df-pr", po::value<bool>(&use_df_pr)->default_value(true)->implicit_value(false), "Disable DF-PR")
            ("disable-tdcp", po::value<bool>(&use_tdcp)->default_value(true)->implicit_value(false), "Disable TDCP")
            ("disable-clock-const", po::value<bool>(&use_clock_const)->default_value(true)->implicit_value(false), "Disable Clock Const")
            ("df-pr-weight", po::value<double>(&df_pr_weight)->default_value(DF_PR_WEIGHT), "Set DF-PR weight")
            ("tdcp-weight", po::value<double>(&tdcp_weight)->default_value(TDCP_WEIGHT), "Set TDCP weight")
            ("clock-const-weight", po::value<double>(&clock_const_weight)->default_value(CONSATANT_CLOCK_WEIGHT), "Set Clock Const weight")
            ("start-epoch", po::value<size_t>(&start_epoch)->default_value(600), "Set start epoch (default 599)")
            ("T", po::value<size_t>(&T)->default_value(40), "Set T value (default 40)")
            ("constellations", po::value<std::vector<std::string>>(&constellations)->multitoken()->default_value(std::vector<std::string>{"gps"}, "gps"), 
                "Set GPS constellations (gps, bds, gal)");


        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return false;  // 도움말 출력 후 false 반환
        }

        start_epoch = start_epoch - 1;
        std::sort(constellations.begin(), constellations.end());

        std::map<std::string, int> gps_constellation_map = {
            {"gps", 0},
            {"bds", 1},
            {"gal", 2}
        };
        
        constellation_name = constellations[0];
        for (size_t i = 0; i < constellations.size(); ++i) {
            if (i > 0) {
                constellation_name += "_" + constellations[i];
            }
            
            int value = gps_constellation_map.at(constellations[i]);
            std::cout << "Processing constellation value: " << value << std::endl;

            // Check for duplicates
            if (constellation_type.find(value) != constellation_type.end()) {
                std::cerr << "Error: Constellation is duplicated" << std::endl;
                return false;
            }

            constellation_type.insert(value);
        }

    } catch (const std::out_of_range& e) {
        std::cerr << "Error: Key not found in map - " << e.what() << std::endl;
        return false;
    } catch (const po::error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }

    return true;
}


int main(int argc, char** argv) {
    std::string rover_dir = "../data/data_rover/";
    std::string station_dir = "../data/data_station/";

    std::string pr_file = "pr.csv";
    std::string ph_file = "carrier.csv";
    std::string dop_file = "dop.csv";
    std::string sv_pos_file = "sv_pos.csv";
    std::string sv_vel_file = "sv_vel.csv";
    std::string time_file = "time.csv";

    // CSV 파일 읽기
    std::vector<std::vector<double>> pr_data = readPseudorangeCSV(rover_dir + pr_file);
    std::vector<std::vector<double>> ph_data = readPseudorangeCSV(rover_dir + ph_file);
    std::vector<std::vector<double>> dop_data = readPseudorangeCSV(rover_dir + dop_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data = readSVPosAndVelCSV(rover_dir + sv_pos_file);
    std::vector<std::vector<std::vector<double>>> sv_vel_data = readSVPosAndVelCSV(rover_dir + sv_vel_file);
    std::vector<std::pair<int, double>> time_data = readGpsTimeCSV(rover_dir + time_file);

    std::vector<std::vector<double>> pr_data_station = readPseudorangeCSV(station_dir + pr_file);
    std::vector<std::vector<double>> dop_data_station = readPseudorangeCSV(station_dir + dop_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data_station = readSVPosAndVelCSV(station_dir + sv_pos_file);
    std::vector<std::vector<std::vector<double>>> sv_vel_data_station = readSVPosAndVelCSV(station_dir + sv_vel_file);
    std::vector<std::pair<int, double>> time_data_station = readGpsTimeCSV(station_dir + time_file);

    /*Define Input Factor*/

    bool use_df_pr, use_tdcp, use_clock_const;
    double df_pr_weight, tdcp_weight, clock_const_weight;

    size_t start_epoch, T;

    std::set<int> constellation_type;
    std::string constellation_name;

    if (!parseCommandLineOptions(argc, argv, 
                                use_df_pr, use_tdcp, use_clock_const, 
                                df_pr_weight, tdcp_weight, clock_const_weight, 
                                start_epoch, T,
                                constellation_type, constellation_name)) {
        return EXIT_FAILURE;  // 옵션 파싱 실패 또는 도움말 출력 시 프로그램 종료
    }

    

    // 프로그램 로직이 여기에 들어갑니다.
    // 예시 출력
    std::cout << "DF-PR enabled: " << std::boolalpha << use_df_pr << "\n";
    std::cout << "TDCP enabled: " << std::boolalpha << use_tdcp << "\n";
    std::cout << "Clock Const enabled: " << std::boolalpha << use_clock_const << "\n";
    std::cout << "DF-PR weight: " << df_pr_weight << "\n";
    std::cout << "TDCP weight: " << tdcp_weight << "\n";
    std::cout << "Clock Const weight: " << clock_const_weight << "\n";
    std::cout << "Start epoch: " << start_epoch << "\n";
    std::cout << "T: " << T << "\n";

    std::cout << "Constellation types: ";
    for (auto& type : constellation_type) {
        std::cout << type << " ";
    }
    std::cout << "\n";

    // 폴더 이름을 start_epoch, T, df_pr_weight, tdcp_weight, clock_const_weight, gps_type, gps_constellations로 구성
    std::string folder_name = "../result/" + constellation_name + "/epoch_" + std::to_string(start_epoch + 1) + "_T_" + std::to_string(T);

    if (use_df_pr)
        folder_name += "_dfpr_" + std::to_string(df_pr_weight);
    
    if (use_tdcp)
        folder_name += "_tdcp_" + std::to_string(tdcp_weight);
    
    if (use_clock_const)
        folder_name += "_clockconst_" + std::to_string(clock_const_weight);


    // 폴더 생성
    if (!fs::exists(folder_name)) {
        if (!fs::create_directories(folder_name)) {
            std::cerr << "Error: Could not create directory " << folder_name << std::endl;
            return EXIT_FAILURE;
        }
    }

    // 로그 파일 경로
    std::string log_file_ecef = folder_name + "/pos_ecef.csv";
    std::string log_file_llh = folder_name + "/pos_llh.csv";
    std::string log_file_cov = folder_name + "/error_cov.csv";

    // log results
    std::ofstream fout_ecef(log_file_ecef);
    std::ofstream fout_llh(log_file_llh);
    std::ofstream fout_cov(log_file_cov);

    fout_ecef << "Epoch, Position X(m), Position Y(m), Position Z(m), Clk Bias(s)-GPS\n";
    fout_llh << "Epoch, Latitude, Longitude, Altitude\n";
    fout_ecef  << std::fixed << std::setprecision(10);
    fout_llh  << std::fixed << std::setprecision(10);
    fout_cov  << std::fixed << std::setprecision(10);

    std::vector<double * > state;
    ceres::Problem problem;

    // std::vector<double> ref_location = coordinate::lla2ecef({36.372371713580250, 127.358800510185191, 91.642377777777796});
    std::vector<double> ref_location = {-3.119992580788137e+06, 4.086868171897103e+06, 3.761594895585738e+06};

    const size_t val_num = 6; // x, y ,z, t_gps, t_glo,
    const size_t max_epoch = start_epoch + T; 


    double* previous_position = nullptr;
    double* current_position = nullptr;

    size_t tdcp_cnt = 0, df_pr_cnt = 0, clock_const_cnt = 0;

    //Calculate Weight
    for(size_t epoch=start_epoch; epoch <= max_epoch; ++epoch){
        for(size_t satellite=0; satellite < pr_data[epoch].size(); ++satellite){
            double pr_value = pr_data[epoch][satellite];
            double pr_value_station = pr_data_station[epoch][satellite];
            double next_dop_value = dop_data[epoch][satellite];

            int satellite_type;

            if (satellite < 32) // GPS
                satellite_type = 0;
            else if(satellite < 59) //GLO
                satellite_type = 4;
            else if (satellite < 95) // GAL
                satellite_type = 1;
            else if (satellite < 158) // BDS
                satellite_type = 2;
            else
                assert(false);

            if(constellation_type.find(satellite_type) == constellation_type.end())
                continue;

            if (!std::isnan(pr_value) && !std::isnan(pr_value_station) && !std::isnan(next_dop_value) && !check_sv_data(sv_pos_data[epoch][satellite])){
                df_pr_cnt += 1;
            }

            if(epoch > start_epoch){
                double prev_dop_value = dop_data[epoch-1][satellite]; // Hz

                double prev_ph_value = ph_data[epoch-1][satellite];
                double curr_ph_value = ph_data[epoch][satellite];

                if (!std::isnan(prev_ph_value) && 
                    !std::isnan(curr_ph_value) &&
                    !check_sv_data(sv_pos_data[epoch][satellite]) && 
                    !check_sv_data(sv_pos_data[epoch-1][satellite])){

                    tdcp_cnt += 1;
                }
            }
        }
        if(epoch > start_epoch)
            clock_const_cnt += 3;
    }


    for(size_t epoch=start_epoch; epoch <= max_epoch; ++epoch){
        previous_position = current_position;
        current_position = new double[val_num];
        std::fill(current_position, current_position + val_num, 0.0); 

        for(int i=0; i<3; i++)
            current_position[i] = ref_location[i];   

        for(size_t satellite=0; satellite < pr_data[epoch].size(); ++satellite){
            int satellite_type;

            if (satellite < 32) // GPS
                satellite_type = 0;
            else if(satellite < 59) //GLO
                satellite_type = 3;
            else if (satellite < 95) // GAL
                satellite_type = 2;
            else if (satellite < 158) // BDS
                satellite_type = 1;
            else
                assert(false);

            if(constellation_type.find(satellite_type) == constellation_type.end()){
                continue;
            }
            

            double pr_value = pr_data[epoch][satellite];
            double pr_value_station = pr_data_station[epoch][satellite];
            double next_dop_value = dop_data[epoch][satellite];

            if (use_df_pr && !std::isnan(pr_value) && !std::isnan(pr_value_station) && !std::isnan(next_dop_value) && !check_sv_data(sv_pos_data[epoch][satellite])){
                factor::DiffPesudorangeFactorCostFunctor* functor = 
                    new factor::DiffPesudorangeFactorCostFunctor(ref_location, sv_pos_data[epoch][satellite], 
                                                                pr_value-pr_value_station, satellite_type, df_pr_weight / df_pr_cnt);
 
                ceres::CostFunction* cost_function =
                    new ceres::AutoDiffCostFunction<factor::DiffPesudorangeFactorCostFunctor, 1, val_num>(functor);

                problem.AddResidualBlock(cost_function, nullptr, current_position);
            }

            if(use_tdcp &&
                epoch > start_epoch){
                double prev_ph_value = ph_data[epoch-1][satellite];
                double curr_ph_value = ph_data[epoch][satellite];

                if (!std::isnan(prev_ph_value) && 
                    !std::isnan(curr_ph_value) &&
                    !check_sv_data(sv_pos_data[epoch][satellite]) && 
                    !check_sv_data(sv_pos_data[epoch-1][satellite])){

                    //Add TDCP Factor
                    factor::TDCPFactorCostFunctor* functor = 
                        new factor::TDCPFactorCostFunctor(sv_pos_data[epoch-1][satellite], sv_pos_data[epoch][satellite],
                                                          curr_ph_value - prev_ph_value,
                                                          satellite_type, tdcp_weight / tdcp_cnt);

                    ceres::CostFunction* cost_function = 
                        new ceres::AutoDiffCostFunction<factor::TDCPFactorCostFunctor, 1, val_num, val_num>(functor);

                    problem.AddResidualBlock(cost_function, nullptr, previous_position, current_position);
                }
            }
        }

        if (use_clock_const &&
            epoch > start_epoch){
            for(int i=4; i<7; i++){
            // Add Constant Clock Bias Factor
                factor::ConstantClockBiasFactorCostFunctor* functor = 
                        new factor::ConstantClockBiasFactorCostFunctor(i, clock_const_weight / clock_const_cnt);

                ceres::CostFunction* cost_function = 
                    new ceres::AutoDiffCostFunction<factor::ConstantClockBiasFactorCostFunctor, 1, val_num, val_num>(functor);
                problem.AddResidualBlock(cost_function, nullptr, previous_position, current_position);
            }
        }

        state.push_back(current_position);
    }
    // Solver 옵션 설정 및 실행
    ceres::Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    options.parameter_tolerance = 1e-12; // 더 엄격한 파라미터 수렴 조건
    options.function_tolerance = 1e-12;  // 더 엄격한 함수 값 수렴 조건
    options.gradient_tolerance = 1e-12;  // 더 엄격한 그라디언트 수렴 조건
    options.max_num_iterations = 1000;   // 최대 반복 횟수 설정

    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    ceres::Covariance::Options cov_options;
    // cov_options.algorithm_type = ceres::DENSE_SVD;
    ceres::Covariance covariance(cov_options);

    std::vector<std::pair<const double*, const double*>> covariance_blocks;
    for (double* position : state) {
        covariance_blocks.emplace_back(position, position);
    }

    bool covariance_result = covariance.Compute(covariance_blocks, &problem);

    for (size_t epoch = 0; epoch < state.size(); ++epoch) {
        double* position = state[epoch];
        std::cout << std::fixed << std::setprecision(6);
        // std::cout << "Epoch(ECEF) " << epoch << ": " << position[0] << ", " << position[1] << ", " << position[2] << ", " << position[3] << "\n";
        fout_ecef << epoch << ", " << position[0] << ", " << position[1] << ", " << position[2] << ", " << position[3] << "\n";

        std::vector<double> ecef_position = {position[0], position[1], position[2]};
        std::vector<double> res = coordinate::ecef2lla(ecef_position);

        std::cout << "Epoch(LLA) " << epoch << ": " << res[0] << ", " << res[1] << ", " << res[2] << "\n";
        fout_llh << epoch << ", " << res[0] << ", " << res[1] << ", " << res[2] << "\n";

        // if (covariance_result){
        //     double covariance_matrix[3 * 3];
        //     covariance.GetCovarianceBlock(position, position, covariance_matrix);
        //     std::cout << "Covariance Matrix:\n";
        //     for (int i = 0; i < 3; ++i) {
        //         for (int j = 0; j < 3; ++j) {
        //             std::cout << covariance_matrix[3 * i + j] << " ";
        //             fout_cov << covariance_matrix[3 * i + j] << ", ";
        //         }
        //         std::cout << "\n";
        //         fout_cov << "\n";
        //     }
        // }    

        delete[] position; // 메모리 해제
    }

    std::cout << summary.FullReport() << "\n";

    fout_ecef.close();
    fout_llh.close();
    fout_cov.close();

    return 0;
}