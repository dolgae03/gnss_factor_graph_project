#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <random> 

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


#define DF_PR_WEIGHT (double) 1.0/sqrt(2)
#define TDCP_WEIGHT (double) 50.0/sqrt(2)
#define CONSATANT_CLOCK_WEIGHT (double) 0
#define TAU_WEIGHT (double) 100

// std::string rover_dir = "../data/rooftop4/data_rover/";
// std::string station_dir = "../data/rooftop4/data_station/";



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
                             bool& use_tau,
                             double& df_pr_weight, 
                             double& tdcp_weight, 
                             double& clock_const_weight,
                             double& tau_weight,
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
            ("disable-clock-const", po::value<bool>(&use_clock_const)->default_value(false)->implicit_value(false), "Disable Clock Const")
            ("disable-tau", po::value<bool>(&use_tau)->default_value(true)->implicit_value(false), "Disable Tau Factor")
            ("df-pr-weight", po::value<double>(&df_pr_weight)->default_value(DF_PR_WEIGHT), "Set DF-PR weight")
            ("tdcp-weight", po::value<double>(&tdcp_weight)->default_value(TDCP_WEIGHT), "Set TDCP weight")
            ("clock-const-weight", po::value<double>(&clock_const_weight)->default_value(CONSATANT_CLOCK_WEIGHT), "Set Clock Const weight")
            ("tau-weight", po::value<double>(&tau_weight)->default_value(TAU_WEIGHT), "Set Tau weight")
            ("start-epoch", po::value<size_t>(&start_epoch)->default_value(600), "Set start epoch (default 599)")
            ("T", po::value<size_t>(&T)->default_value(100), "Set T value (default 100)")
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

int runOptimization(int tau, int version, const std::string& matlab_save_dir, size_t start_epoch, size_t T, bool use_df_pr, bool use_tdcp, bool use_clock_const,  bool use_tau, double df_pr_weight, double tdcp_weight, double clock_const_weight, double tau_weight, const std::set<int>& constellation_type, const std::string& constellation_name) {

    cout << "=========================== Version " << version +1 <<" Starts ==========================="<< endl;
    std::string tau_str = "/tau_" + std::to_string(tau);
    std::string folder_name = matlab_save_dir + "/monte_carlo"+tau_str+"/epoch_" + std::to_string(start_epoch + 1) + "_T_" + std::to_string(T);
    std::string version_str = "/v" + std::to_string(version+1);  // version 1 to 100
    std::string rover_dir = "../data/monte_carlo" +tau_str + "/data_rover" + version_str;
    std::string station_dir = "../data/monte_carlo" +tau_str + "/data_base" + version_str;

    std::string pr_file = "/pr.csv";
    std::string ph_file = "/carrier.csv";
    std::string sv_pos_file = "/sv_pos.csv";
    std::string sv_vel_file = "/sv_vel.csv";
    std::string time_file = "/time.csv";

    // CSV 파일 읽기
    std::vector<std::vector<double>> pr_data = readPseudorangeCSV(rover_dir + pr_file);
    std::vector<std::vector<double>> ph_data = readPseudorangeCSV(rover_dir + ph_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data = readSVPosAndVelCSV(rover_dir + sv_pos_file);
    // std::vector<std::vector<std::vector<double>>> sv_vel_data = readSVPosAndVelCSV(rover_dir + sv_vel_file);
    std::vector<std::pair<int, double>> time_data = readGpsTimeCSV(rover_dir + time_file);

    std::vector<std::vector<double>> pr_data_station = readPseudorangeCSV(station_dir + pr_file);
    std::vector<std::vector<std::vector<double>>> sv_pos_data_station = readSVPosAndVelCSV(station_dir + sv_pos_file);
    // std::vector<std::vector<std::vector<double>>> sv_vel_data_station = readSVPosAndVelCSV(station_dir + sv_vel_file);
    std::vector<std::pair<int, double>> time_data_station = readGpsTimeCSV(station_dir + time_file);

    if (use_df_pr)
        folder_name += "_dfpr_" + std::to_string(df_pr_weight);
    
    if (use_tdcp)
        folder_name += "_tdcp_" + std::to_string(tdcp_weight);
    
    if (use_clock_const)
        folder_name += "_clockconst_" + std::to_string(clock_const_weight);

    if (use_tau)
        folder_name += "_tauWeight_" + std::to_string(tau_weight); 

    std::string folder_name_version = folder_name + version_str;
    // 폴더 생성
    if (!fs::exists(folder_name_version)) {
        if (!fs::create_directories(folder_name_version)) {
            std::cerr << "Error: Could not create directory " << folder_name_version << std::endl;
            return EXIT_FAILURE;
        }
    }

    // 로그 파일 경로
    
    std::string log_file_ecef = folder_name_version + "/pos_ecef.csv";
    std::string log_file_llh = folder_name_version + "/pos_llh.csv";
    std::string log_file_cov = folder_name_version + "/error_cov.csv";
    std::string log_file_residual_pr = folder_name_version + "/residual_pr.csv";
    std::string log_file_residual_tdcp = folder_name_version + "/residual_tdcp.csv";

    // log results
    std::ofstream fout_ecef(log_file_ecef);
    std::ofstream fout_llh(log_file_llh);
    std::ofstream fout_cov(log_file_cov);
    std::ofstream fout_residual_pr(log_file_residual_pr);
    std::ofstream fout_residual_tdcp(log_file_residual_tdcp);

    fout_ecef << "Epoch, Position X(m), Position Y(m), Position Z(m), Clk Bias(s)-GPS\n";
    fout_llh << "Epoch, Latitude, Longitude, Altitude\n";
    fout_residual_pr << "pr_residual, epoch, SV, residual" << endl;
    fout_residual_tdcp << "TDCP_residual, epoch, SV, residual" << endl;
    fout_ecef  << std::fixed << std::setprecision(10);
    fout_llh  << std::fixed << std::setprecision(10);
    fout_cov  << std::fixed << std::setprecision(10);
    fout_residual_pr  << std::fixed << std::setprecision(6);
    fout_residual_tdcp  << std::fixed << std::setprecision(6);

    std::vector<double * > state;
    std::vector<double * > position_states;
    std::vector<double * > noise_states;
    ceres::Problem problem;

    // std::vector<double> ref_location = coordinate::lla2ecef({36.372371713580250, 127.358800510185191, 91.642377777777796});
    // std::vector<double> ref_location = {-3.119992580788137e+06, 4.086868171897103e+06, 3.761594895585738e+06}; // rooftop4
    std::vector<double> ref_location = coordinate::lla2ecef({36.3727470000000, 127.357671000000, 10}); //monte-carlo

    const size_t num_var_pos = 4; // x, y ,z, t_gps, t_glo,
    const size_t num_var_meas = 7; // set arbitrarily for now
    const size_t max_epoch = start_epoch + T; 


    // Solver 옵션 설정 및 실행
    ceres::Solver::Summary summary;

    ceres::Solver::Options options;
    // options.minimizer_type = ceres::LINE_SEARCH;  // MATLAB의 'quasi-newton' 알고리즘과 대응
    // options.line_search_direction_type = ceres::BFGS;  // BFGS 방법 사용
    options.minimizer_progress_to_stdout = true;  // 'Display','iter-detailed'에 대응
    options.gradient_tolerance = 1e-12;  // 'OptimalityTolerance'를 더 엄격하게
    options.parameter_tolerance = 1e-12;  // 'TolX'를 더 엄격하게
    options.function_tolerance = 1e-12;  // 'FunctionTolerance'를 더 엄격하게
    options.gradient_check_numeric_derivative_relative_step_size = 1e-13;  // 'FiniteDifferenceStepSize'를 더 엄격하게
    options.max_num_iterations = 1e+5;  // 'MaxIterations'를 늘려 더 많은 반복 허용

    std::string option_err;
    if(!options.IsValid(&option_err)){
        cout << option_err<< endl;

        exit(0);
    }

    double* previous_position = nullptr;
    double* current_position = nullptr;
    double* previous_noise = nullptr;
    double* current_noise = nullptr;
    current_position = new double[num_var_pos];
    current_noise = new double[num_var_meas];
    std::fill(current_position, current_position + num_var_pos, 0.0); 
    std::fill(current_noise, current_noise + num_var_meas, 0.01);
    double* random_noise = nullptr;
    random_noise = new double[num_var_meas];


    for(size_t epoch=start_epoch; epoch < max_epoch; ++epoch){

        previous_position = current_position;
        previous_noise = current_noise;
        // current_position = new double[num_var_pos];
        // std::cout << endl;

        // previous_noise = current_noise;
        // current_noise = new double[num_var_meas];
        current_position = new double[num_var_pos];
        current_noise = new double[num_var_meas];
        
        for(int i=0; i<num_var_pos; i++) // initial values
            current_position[i] = previous_position[i];  
        for(int i=0; i<num_var_meas; i++)
            current_noise[i] = previous_noise[i]; 


        // init value
        // if (epoch == start_epoch) {
        //     std::fill(current_position, current_position + num_var_pos, 0.0); 
        //     std::fill(current_noise, current_noise + num_var_meas, 0.1);
        //     // previous_noise = new double[num_var_meas];
        //     // std::fill(previous_noise, previous_noise + num_var_meas, 0.1);
        //     std::fill(previous_position, previous_position + num_var_pos, 0.0); 
        //     std::fill(previous_noise, previous_noise + num_var_meas, 0.1);
        // }
        // else {
            // for(int i=0; i<num_var_pos; i++) // initial values
            //     current_position[i] = previous_position[i];  
            // for(int i=0; i<num_var_meas; i++)
            //     current_noise[i] = previous_noise[i]; 

        // }
        // std::cout << "Epoch-" << epoch << " current noise:  "<< current_noise[0] << endl;


        // if (use_tau && (epoch > start_epoch)) {        
        if (use_tau) {  

            factor::generateRandomNoise(tau, num_var_meas, random_noise);
            std::cout<< "Epoch "<< epoch << "Random Noise ";
            for (int i = 0; i < num_var_meas; i++)
                std::cout << random_noise[i] << " ";
            std::cout << std::endl;
            


            factor::TimeCorrelationFactorCostFunctor* functor = 
                new factor::TimeCorrelationFactorCostFunctor(tau, tau_weight, num_var_meas, random_noise);
            ceres::CostFunction* cost_function =
                new ceres::AutoDiffCostFunction<factor::TimeCorrelationFactorCostFunctor, 1, num_var_meas, num_var_meas>(functor);
            // problem.AddResidualBlock(cost_function, nullptr, previous_noise, current_noise);
            
            // // state.push_back(current_noise);
            noise_states.push_back(current_noise);
            
        }

        
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

            if (use_df_pr && !std::isnan(pr_value) && !std::isnan(pr_value_station) && !check_sv_data(sv_pos_data[epoch][satellite])){

                factor::DiffPesudorangeFactorCostFunctor* functor = 
                    new factor::DiffPesudorangeFactorCostFunctor(ref_location, sv_pos_data[epoch][satellite], 
                                                                pr_value-pr_value_station, satellite_type, df_pr_weight, satellite);

                ceres::CostFunction* cost_function =
                    new ceres::AutoDiffCostFunction<factor::DiffPesudorangeFactorCostFunctor, 1, num_var_pos, num_var_meas>(functor);


                problem.AddResidualBlock(cost_function, nullptr, current_position, current_noise);
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

                    ceres::CostFunction* cost_function = 
                        new ceres::AutoDiffCostFunction<factor::TDCPFactorCostFunctor, 1, num_var_pos, num_var_pos>(functor);

                    problem.AddResidualBlock(cost_function, nullptr, previous_position, current_position);
                }
            }
        }

        if (use_clock_const &&
            epoch > start_epoch){
            for(int i=4; i<num_var_pos; i++){
            // Add Constant Clock Bias Factor
                factor::ConstantClockBiasFactorCostFunctor* functor = 
                        new factor::ConstantClockBiasFactorCostFunctor(i, clock_const_weight);

                ceres::CostFunction* cost_function = 
                    new ceres::AutoDiffCostFunction<factor::ConstantClockBiasFactorCostFunctor, 1, num_var_pos, num_var_pos>(functor);
                problem.AddResidualBlock(cost_function, nullptr, previous_position, current_position);
            }
        }

        // state.push_back(current_position);
        position_states.push_back(current_position);
    }
    ceres::Solve(options, &problem, &summary);



    // Set up evaluation options
    ceres::Problem::EvaluateOptions eval_options;
    std::vector<double> residuals;
    std::vector<double> jacobian;    // Flattened Jacobian matrix
    ceres::CRSMatrix crs_jacobian;   // Compressed Row Storage (CRS) matrix for the Jacobian
    double* cost = nullptr;
    // Evaluate the problem to extract residuals and Jacobian
    problem.Evaluate(eval_options, cost, &residuals, nullptr, &crs_jacobian);


    // Print the size of the Jacobian matrix
    // std::cout << "Jacobian Matrix Size: " << std::endl;
    // std::cout << "Number of rows (residuals): " << crs_jacobian.num_rows << std::endl;
    // std::cout << "Number of columns (parameters): " << crs_jacobian.num_cols << std::endl;

    
    ceres::Covariance::Options cov_options;
    ceres::Covariance covariance(cov_options);
    
    std::vector<std::pair<const double*, const double*>> covariance_blocks;
    for (double* position : position_states) {
        covariance_blocks.emplace_back(position, position);
    }

    bool covariance_result = false;
    // bool covariance_result = covariance.Compute(covariance_blocks, &problem);
    // std::cout << " =============== Covariance_result: " << covariance_result << " ==============="<< endl;
    
    for (size_t epoch = start_epoch; epoch < max_epoch; epoch++) {
        double* noise = noise_states[epoch - start_epoch];

        std::cout << "Epoch " << epoch << ": Noise: ";
        for (int i =0; i < num_var_meas ; i++){
            std::cout << noise[i] << " ";
        }
        std::cout << endl;
    }
    for (size_t epoch = start_epoch; epoch < max_epoch; epoch++) {
        double* position = position_states[epoch - start_epoch];
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Epoch(ECEF) " << epoch << ": " << position[0] << ", " << position[1] << ", " << position[2] << ", " << position[3] << "\n";
        fout_ecef << epoch;
        for(int i=0; i<num_var_pos; i++) {
            fout_ecef << ", " << position[i];
        }
        fout_ecef << "\n";
        std::vector<double> ecef_position = {position[0], position[1], position[2]};
        std::vector<double> res = coordinate::ecef2lla(ecef_position);

        // std::cout << "Epoch(LLA) " << epoch << ": " << res[0] << ", " << res[1] << ", " << res[2] << "\n";
        // fout_llh << epoch << ", " << res[0] << ", " << res[1] << ", " << res[2] << "\n";
        double df_pr_sum = 0;
        double tdcp_sum = 0;
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
            if (use_df_pr && !std::isnan(pr_value) && !std::isnan(pr_value_station)  && !check_sv_data(sv_pos_data[epoch][satellite])){
                factor::DiffPesudorangeFactorCostFunctor functor(ref_location, sv_pos_data[epoch][satellite], 
                                                                pr_value-pr_value_station, satellite_type, df_pr_weight, satellite);
                double residual[1];
                functor(position_states[epoch-start_epoch], noise_states[epoch-start_epoch], residual);
                
                df_pr_sum += residual[0] * residual[0];
                fout_residual_pr << epoch +1 << ", " << satellite+1 << ", "<< residual[0] << endl;
                std::cout << "Epoch " << epoch << " Final Noise: " << noise_states[epoch-start_epoch][satellite]<< endl;
                // for (int i = 0; i <num_var_meas ; i ++)
                //     cout << noise_states[epoch-start_epoch][i] << " ";
                // cout << endl;
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
                    factor::TDCPFactorCostFunctor functor(sv_pos_data[epoch-1][satellite], sv_pos_data[epoch][satellite],
                                                        curr_ph_value - prev_ph_value,
                                                        satellite_type, tdcp_weight);

                    double residual[1];
                    functor(position_states[epoch - start_epoch - 1], position_states[epoch - start_epoch], residual);

                    tdcp_sum += residual[0] * residual[0];
                    // std::cout << "tdcp,"<< epoch + 1 << "," << satellite<< ","<<  residual[0] * residual[0] << std::endl;
                    // Output the residual
                    fout_residual_tdcp << epoch +1 << ", " << satellite+1 << ", "<< residual[0] << endl;
                
                    
                }
                
            }
        }
        // std::cout << "pr,"<< epoch + 1 << "," << df_pr_sum << std::endl;
        // if(use_tdcp && epoch > start_epoch)
            // std::cout << "tdcp,"<< epoch + 1 << "," << tdcp_sum << std::endl;

        
        if (covariance_result){
            double covariance_matrix[num_var_pos * num_var_pos];
            covariance.GetCovarianceBlock(position, position, covariance_matrix); 
            // std::cout << "Covariance Matrix:\n";
            for (int i = 0; i < num_var_pos; i++) {
                for (int j = 0; j < num_var_pos; j++) {
                    // std::cout << covariance_matrix[num_var_pos * i + j] << " ";
                    fout_cov << covariance_matrix[num_var_pos * i + j] << ", ";
                }
                // std::cout << "\n";
                fout_cov << "\n";
            }
        }    

    }

    std::cout << summary.FullReport() << "\n";

    fout_ecef.close();
    fout_llh.close();
    fout_cov.close();
    fout_residual_pr.close();
    fout_residual_tdcp.close();

    // Clear the state for the next iteration
    // for (double* pos : state) {
    //     delete[] pos;
    // }
    // state.clear();
    for (double* pos : position_states) {
        delete[] pos;
    }
    position_states.clear();
    for (double* noise : noise_states) {
        delete[] noise;
    }
    noise_states.clear();
    
    return 0;

}

int main(int argc, char** argv) {
    // test_one_factor();
    // return -1;

    /*Define Input Factor*/

    bool use_df_pr, use_tdcp, use_clock_const, use_tau;
    double df_pr_weight, tdcp_weight, clock_const_weight, tau_weight;
    const double CONST_C = 299792458.0;

    size_t start_epoch, T;

    std::set<int> constellation_type;
    std::string constellation_name;

    if (!parseCommandLineOptions(argc, argv, 
                                use_df_pr, use_tdcp, use_clock_const, use_tau,
                                df_pr_weight, tdcp_weight, clock_const_weight, tau_weight,
                                start_epoch, T,
                                constellation_type, constellation_name)) {
        return EXIT_FAILURE;  // 옵션 파싱 실패 또는 도움말 출력 시 프로그램 종료
    }

    

    // 프로그램 로직이 여기에 들어갑니다.
    // 예시 출력
    std::cout << "DF-PR enabled: " << std::boolalpha << use_df_pr << "\n";
    std::cout << "TDCP enabled: " << std::boolalpha << use_tdcp << "\n";
    std::cout << "Clock Const enabled: " << std::boolalpha << use_clock_const << "\n";
    std::cout << "Tau enabled: " << std::boolalpha << use_tau << "\n";
    std::cout << "DF-PR weight: " << df_pr_weight << "\n";
    std::cout << "TDCP weight: " << tdcp_weight << "\n";
    std::cout << "Clock Const weight: " << clock_const_weight << "\n";
    std::cout << "Tau weight: " << tau_weight << "\n";
    std::cout << "Start epoch: " << start_epoch << "\n";
    std::cout << "T: " << T << "\n";

    std::cout << "Constellation types: ";
    for (auto& type : constellation_type) {
        std::cout << type << " ";
    }
    std::cout << "\n";

    // 폴더 이름을 start_epoch, T, df_pr_weight, tdcp_weight, clock_const_weight, gps_type, gps_constellations로 구성
    std::string matlab_save_dir = "/mnt/c/jaeryoung/research/factor_graph/fgo_basic/error_simulation/result/result_ceres/" ;
    // std::string folder_name = matlab_save_dir + constellation_name +"/rooftop4"+ "/epoch_" + std::to_string(start_epoch + 1) + "_T_" + std::to_string(T);
    
    // double tau = 0;
    double tau = 100;
    

    int version_num = 1;

    for (int version = 0; version < version_num; version++) {
        runOptimization(tau, version, matlab_save_dir, start_epoch, T, use_df_pr, use_tdcp, use_clock_const, use_tau, df_pr_weight, tdcp_weight, clock_const_weight, tau_weight, constellation_type, constellation_name);
    }

    return 0;
}