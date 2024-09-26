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

void calculate_factor_count(size_t& epoch, std::set<int>& constellation_type,
                            bool is_first, size_t& df_pr_cnt, size_t& tdcp_cnt, size_t& clock_const_cnt){
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

        if (!std::isnan(pr_value) && !std::isnan(pr_value_station) && 
            !std::isnan(next_dop_value) && !check_sv_data(sv_pos_data[epoch][satellite])){
            df_pr_cnt += 1;
        }

        if(!is_first){
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

    if (!is_first){
        for(int i=4; i<7; i++){
            clock_const_cnt += 1;
        }
    }
}

int test_one_factor() {

    // Define satellite positions at previous and current epochs
    std::vector<double> sv_position_prev = {-1.27378e+07, 2.32632e+07, 1.0749e+06 };        // Previous satellite position (ECEF coordinates in meters)
    std::vector<double> sv_position_curr = {-1.27381e+07, 2.32632e+07 , 1.07173e+06 };     // Current satellite position (moved slightly in Y-axis)

    double pseudorange = 2684.21;  // Pseudorange measurement in meters

    int satellite_type = 0;  // GPS L1C
    double weight = 0.909091;     // Weight of the measurement

    // Create the cost functor
    factor::TDCPFactorCostFunctor functor(sv_position_prev, sv_position_curr, pseudorange, satellite_type, weight);

    // Define previous and current user states
    // Assuming the user is at the origin at the previous state and has moved slightly in the X-axis
    double prev_state[6] = {-3119984.7525596456, 4086868.4656001395, 3761601.4722899213, -0.0000000050, 0, 0};              // [x, y, z, clock_bias_GPS, clock_bias_GALILEO, clock_bias_Beidou]
    double curr_state[6] = {-3119984.0802085600, 4086868.8650446292, 3761601.5717207277, -0.0000000055, 0, 0};          // Moved 10 meters in X, small clock bias change

    // Define residual array
    double residual[1];

    // Call the operator()
    functor(prev_state, curr_state, residual);

    // Output the residual
    std::cout << "Residual: " << residual[0] << std::endl;

    return 0;
}


int main(int argc, char** argv) {
    // test_one_factor();
    // return -1;

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

    // Solver 옵션 설정 및 실행
    ceres::Solver::Summary summary;

    ceres::Solver::Options options;
    // options.minimizer_type = ceres::LINE_SEARCH;  // MATLAB의 'quasi-newton' 알고리즘과 대응
    // options.line_search_direction_type = ceres::BFGS;  // BFGS 방법 사용
    options.minimizer_progress_to_stdout = true;  // 'Display','iter-detailed'에 대응
    options.gradient_tolerance = 1e-8;  // 'OptimalityTolerance'에 대응
    options.parameter_tolerance = 1e-10;  // 'TolX'에 대응
    options.function_tolerance = 1e-8;  // 'FunctionTolerance'에 대응
    options.gradient_check_numeric_derivative_relative_step_size = 1e-10;  // 'FiniteDifferenceStepSize'에 대응
    options.max_num_iterations = 1e+4;  // 'MaxIterations'에 대응

    std::string option_err;
    if(!options.IsValid(&option_err)){
        cout << option_err<< endl;

        exit(0);
    }

    for(size_t epoch=start_epoch; epoch < max_epoch; ++epoch){
        previous_position = current_position;
        current_position = new double[val_num];

        // init value
        std::fill(current_position, current_position + val_num, 0.0); 
        for(int i=0; i<3; i++)
            current_position[i] = ref_location[i];  

        size_t df_pr_cnt = 0, tdcp_cnt = 0, clock_const_cnt = 0;
        calculate_factor_count(epoch, constellation_type,
                               epoch == start_epoch, df_pr_cnt, tdcp_cnt, clock_const_cnt);

        // cout << epoch << " " << df_pr_cnt << " " << tdcp_cnt << " " << clock_const_cnt << endl;


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
                                                                pr_value-pr_value_station, satellite_type, df_pr_weight);
 
                ceres::CostFunction* cost_function =
                    new ceres::AutoDiffCostFunction<factor::DiffPesudorangeFactorCostFunctor, 1, val_num>(functor);


                problem.AddResidualBlock(cost_function, nullptr, current_position);
            }

                // if(satellite == 13 || satellite == 18)
                // continue;


            if(use_tdcp &&
                epoch > start_epoch){
                double prev_ph_value = ph_data[epoch-1][satellite];
                double curr_ph_value = ph_data[epoch][satellite];

                if (!std::isnan(prev_ph_value) && 
                    !std::isnan(curr_ph_value) &&
                    !check_sv_data(sv_pos_data[epoch][satellite]) && 
                    !check_sv_data(sv_pos_data[epoch-1][satellite])){
                    

                    if (epoch == 645 && satellite == 13){
                        for(int i=0; i<3; i++)
                            cout << sv_pos_data[epoch-1][satellite][i] << " ";
                        cout << endl;

                        for(int i=0; i<3; i++)
                            cout << sv_pos_data[epoch][satellite][i] << " ";
                        cout << endl;

                        cout << curr_ph_value - prev_ph_value << endl;
                        cout << tdcp_weight / tdcp_cnt << endl;
                    }

                    //Add TDCP Factor
                    factor::TDCPFactorCostFunctor* functor = 
                        new factor::TDCPFactorCostFunctor(sv_pos_data[epoch-1][satellite], sv_pos_data[epoch][satellite],
                                                          curr_ph_value - prev_ph_value,
                                                          satellite_type, tdcp_weight);

                    ceres::CostFunction* cost_function = 
                        new ceres::AutoDiffCostFunction<factor::TDCPFactorCostFunctor, 1, val_num, val_num>(functor);

                    
                    // factor::NumTDCPFactorCostFunctor* functor = 
                    //     new factor::NumTDCPFactorCostFunctor(sv_pos_data[epoch-1][satellite], sv_pos_data[epoch][satellite],
                    //                                       curr_ph_value - prev_ph_value,
                    //                                       satellite_type, tdcp_weight / tdcp_cnt);

                    // ceres::CostFunction* cost_function = 
                    //     new ceres::NumericDiffCostFunction<factor::NumTDCPFactorCostFunctor, ceres::CENTRAL, 1, val_num, val_num>(functor);

                    problem.AddResidualBlock(cost_function, nullptr, previous_position, current_position);
                }
            }
        }

        if (use_clock_const &&
            epoch > start_epoch){
            for(int i=4; i<val_num; i++){
            // Add Constant Clock Bias Factor
                factor::ConstantClockBiasFactorCostFunctor* functor = 
                        new factor::ConstantClockBiasFactorCostFunctor(i, clock_const_weight);

                ceres::CostFunction* cost_function = 
                    new ceres::AutoDiffCostFunction<factor::ConstantClockBiasFactorCostFunctor, 1, val_num, val_num>(functor);
                problem.AddResidualBlock(cost_function, nullptr, previous_position, current_position);
            }
        }

        state.push_back(current_position);
    }

    ceres::Solve(options, &problem, &summary);
    ceres::Covariance::Options cov_options;
    ceres::Covariance covariance(cov_options);

    

    std::vector<std::pair<const double*, const double*>> covariance_blocks;
    for (double* position : state) {
        covariance_blocks.emplace_back(position, position);
    }

    // bool covariance_result = covariance.Compute(covariance_blocks, &problem);

    for (size_t epoch = start_epoch; epoch < max_epoch; epoch++) {
        double* position = state[epoch - start_epoch];
        std::cout << std::fixed << std::setprecision(6);
        // std::cout << "Epoch(ECEF) " << epoch << ": " << position[0] << ", " << position[1] << ", " << position[2] << ", " << position[3] << "\n";

        fout_ecef << epoch;
        for(int i=0; i<val_num; i++)
            fout_ecef << ", " << position[i];
        fout_ecef << "\n";

        std::vector<double> ecef_position = {position[0], position[1], position[2]};
        std::vector<double> res = coordinate::ecef2lla(ecef_position);

        // std::cout << "Epoch(LLA) " << epoch << ": " << res[0] << ", " << res[1] << ", " << res[2] << "\n";
        fout_llh << epoch << ", " << res[0] << ", " << res[1] << ", " << res[2] << "\n";

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
                                                                pr_value-pr_value_station, satellite_type, df_pr_weight);
 
                double residual[1];
                functor->operator()(state[epoch-start_epoch], residual);
                
                std::cout << "pr,"<< epoch + 1 << "," << satellite + 1 << "," << residual[0] << std::endl;
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
                                                          satellite_type, tdcp_weight);

                    double residual[1];

                    functor->operator()(state[epoch - start_epoch - 1], state[epoch - start_epoch], residual);

                    // Output the residual
                    std::cout << "tdcp,"<< epoch + 1 << "," << satellite + 1 << "," << residual[0] << std::endl;
                }
            }
        }

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