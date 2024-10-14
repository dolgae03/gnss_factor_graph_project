#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <random> 

using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <stdexcept>
#include "ceres/ceres.h"

double GetL1Frequency(int satellite_type) {
    switch (satellite_type) {
        case 0:  // GPS L1C
            return 1575.42e6;
        case 1:  // GALILEO L1C
            return 1575.42e6;
        case 2:  // Beidou B2
            return 1561.098e6;
        default:
            throw std::invalid_argument("Unknown satellite type");
    }
}

template <typename T>
void RotateSatellitePosition(const T* user_state,
                                const std::vector<double>& satellite_position,
                                T* satellite_position_rot) {
    const double c = 299792458.0; // Speed of light in m/s
    const double omega = 7.2921151467e-5; // Earth's rotation rate in rad/s

    T diff_vector[3];
    for (int i = 0; i < 3; ++i) {
        diff_vector[i] = T(satellite_position[i]) - user_state[i];
    }

    T norm = ceres::sqrt(diff_vector[0] * diff_vector[0] +
                            diff_vector[1] * diff_vector[1] +
                            diff_vector[2] * diff_vector[2]);

    T C[3][3] = {
        {T(1),                    T(omega) * norm / T(c), T(0)},
        {T(-omega) * norm / T(c), T(1)                  , T(0)},
        {T(0)                   , T(0)                  , T(1)}
    };

    for(int i=0; i<3; i++){
        satellite_position_rot[i] = T(0);
        for(int j=0; j<3; j++){
            satellite_position_rot[i] += C[i][j] * T(satellite_position[j]);
        }
    }
}

template <typename T>
void CalculateLOS(const T* user_state,
                const T* satellite_position,
                T* los_vector) {
    T diff_vector[3];
    for (int i = 0; i < 3; ++i) {
        diff_vector[i] = satellite_position[i] - user_state[i];
    }

    T norm = ceres::sqrt(diff_vector[0] * diff_vector[0] +
                        diff_vector[1] * diff_vector[1] +
                        diff_vector[2] * diff_vector[2]);

    for (int i = 0; i < 3; ++i) 
        los_vector[i] = diff_vector[i] / norm;
}

template <typename T>
T computeRange(const T* pos1, const T* pos2) {
    T range = T(0);
    for (int i = 0; i < 3; ++i) {
        range += (pos1[i] - pos2[i]) * (pos1[i] - pos2[i]);
    }
    return ceres::sqrt(range);
}

namespace factor {
class ConstantClockBiasFactorCostFunctor {
    public:
        ConstantClockBiasFactorCostFunctor(int satellite_type, double weight)
            : satellite_type_(satellite_type), weight_(weight) {}

        template <typename T>
        bool operator()(const T* const state1, const T* const state2, T* residual) const {
            // state1과 state2의 clock bias 차이 계산
            const double c = 299792458.0;  // Speed of light in m/s
            T clock_bias1 = state1[satellite_type_];
            T clock_bias2 = state2[satellite_type_];
            
            residual[0] = (clock_bias1 - clock_bias2) *T(c) * T(weight_);
            return true;
        }

    private:
        int satellite_type_;
        double weight_;
};


class DiffPesudorangeFactorCostFunctor {
public:
    DiffPesudorangeFactorCostFunctor(const std::vector<double>& ref_position, const std::vector<double>& sv_position, double pesudorange, int satellite_type, double weight, int sv)
        : ref_position_(ref_position), sv_position_(sv_position), pesudorange_(pesudorange), satellite_type_(satellite_type), weight_(weight), sv_(sv) {}

    template <typename T>
    bool operator()(const T* state, T* residual) const {
        T rotated_sv_pos[3];
        RotateSatellitePosition(state, sv_position_, rotated_sv_pos);

        // Use computeRange for both state and reference position
        T range_to_sv = computeRange(state, rotated_sv_pos);
        T clock_bias = state[3 + satellite_type_];
        range_to_sv += clock_bias;

        // Now use ref_position_t with computeRange
        T ref_position_t[3];
        for (int i = 0; i < 3; ++i) {
            ref_position_t[i] = T(ref_position_[i]);
        }
        T range_to_ref = computeRange(ref_position_t, rotated_sv_pos);
        T estimated_pr_diff = range_to_sv - range_to_ref;
        residual[0] = (T(pesudorange_) - estimated_pr_diff) * T(sqrt(weight_));

        return true;
    }

private:
    std::vector<double> ref_position_;
    std::vector<double> sv_position_;
    double pesudorange_;
    double weight_;
    int satellite_type_;
    int sv_;
};


class DiffPesudorangeTauFactorCostFunctor {
   public:
        DiffPesudorangeTauFactorCostFunctor(const std::vector<double>& ref_position, const std::vector<double>& sv_position, double pesudorange, int satellite_type, double weight, int sv)
            : ref_position_(ref_position), sv_position_(sv_position), pesudorange_(pesudorange), satellite_type_(satellite_type), weight_(weight), sv_(sv) {}

        template <typename T>
        bool operator()(const T* state, const T* noise_state, T* residual) const {
            const double c = 299792458.0; // Speed of light in m/s
            const double omega = 7.2921151467e-5; // Earth's rotation rate in rad/s

            T rotated_sv_pos[3];
            RotateSatellitePosition(state, sv_position_, rotated_sv_pos);
            
            T range_to_sv = computeRange(state, rotated_sv_pos);
            T clock_bias = state[3 + satellite_type_];
            range_to_sv += clock_bias;

            T ref_position_T[3];
            for (int i = 0; i < 3; ++i) {
                ref_position_T[i] = T(ref_position_[i]);
            }
            T range_to_ref = computeRange(ref_position_T, rotated_sv_pos);
            T estimated_pr_diff = range_to_sv - range_to_ref + noise_state[0];
            residual[0] = (T(pesudorange_) - estimated_pr_diff) * T(sqrt(weight_));

            return true;
        }

    private:
        std::vector<double> ref_position_;
        std::vector<double> sv_position_;
        double pesudorange_;
        double weight_;
        int satellite_type_;
        int sv_;
};


// 필요한 정보듣 pr difference,     
class DopplerFactorCostFunctor {
    public:
        DopplerFactorCostFunctor(const std::vector<double>& sv_position, const std::vector<double>& sv_velocity, double doppler, double time_interval, int satellite_type, double sigma)
            : sv_position_(sv_position), sv_velocity_(sv_velocity), doppler_(doppler), time_interval_(time_interval), satellite_type_(satellite_type), sigma_(sigma) {}

        template <typename T>
        bool operator()(const T* const prev_state, const T* const next_state, T* residual) const {
            const double c = 299792458.0;  // speed of light in m/s
            T L1_frequency = T(GetL1Frequency(satellite_type_));  // L1 signal frequency in Hz
            
            T los_vector[3];
            T user_position[3];
            for (int i = 0; i < 3; ++i) {
                user_position[i] = (prev_state[i] + next_state[i]) / T(2.0);
                los_vector[i] = T(sv_position_[i]) - user_position[i];
            }

            T los_norm = ceres::sqrt(los_vector[0] * los_vector[0] +
                                    los_vector[1] * los_vector[1] +
                                    los_vector[2] * los_vector[2]);

            for (int i = 0; i < 3; ++i) 
                los_vector[i] /= los_norm;  // Normalize LOS vector

            T user_velocity[3];
            for (int i = 0; i < 3; ++i) 
                user_velocity[i] = (next_state[i] - prev_state[i]) / (T(time_interval_) + next_state[3 + satellite_type_] - prev_state[3 + satellite_type_]);
        
            T relative_velocity[3];
            for (int i = 0; i < 3; ++i) 
                relative_velocity[i] = user_velocity[i] - T(sv_velocity_[i]);

            T relative_speed = relative_velocity[0] * los_vector[0] +
                            relative_velocity[1] * los_vector[1] +
                            relative_velocity[2] * los_vector[2];

            T doppler_velocity = T(doppler_) * T(c) / L1_frequency;
    
            residual[0] = (relative_speed - doppler_velocity) / T(sigma_);

            return true;
        }

    private:
        std::vector<double> sv_position_;
        std::vector<double> sv_velocity_;
        double doppler_;
        double time_interval_;
        int satellite_type_;
        double sigma_;
};

class TDCPFactorCostFunctor {
public:
    TDCPFactorCostFunctor(const std::vector<double>& sv_position_prev,
                          const std::vector<double>& sv_position_curr,
                          double pseudorange,
                          int satellite_type,
                          double weight)
        : sv_position_prev_(sv_position_prev),
          sv_position_curr_(sv_position_curr),
          pseudorange_(pseudorange),
          satellite_type_(satellite_type),
          weight_(weight) {}

    template <typename T>
    bool operator()(const T* const prev_state, const T* const curr_state, T* residual) const {
        const double c = 299792458.0;  // Speed of light in m/s
        const double omega = 7.292115e-5;  // Earth rotation rate in rad/s

        bool debug = false;

        T L1_frequency = T(GetL1Frequency(satellite_type_));  // L1 signal frequency in Hz

        T sv_position_prev_rot[3];
        T sv_position_curr_rot[3];
  

        RotateSatellitePosition(prev_state, sv_position_prev_, sv_position_prev_rot);
        RotateSatellitePosition(curr_state, sv_position_curr_, sv_position_curr_rot);

        for (int i = 0; i < 3; ++i) {
            sv_position_prev_rot[i] = T(sv_position_prev_[i]);
            sv_position_curr_rot[i] = T(sv_position_curr_[i]);
        }     


        T los_vector_prev[3];
        T los_vector_curr[3];
        CalculateLOS(prev_state, sv_position_prev_rot, los_vector_prev);
        CalculateLOS(curr_state, sv_position_curr_rot, los_vector_curr);

        T range_prev = T(0);
        for(int i=0; i<3; i++)
            range_prev += (prev_state[i] - sv_position_prev_rot[i]) * (prev_state[i] - sv_position_prev_rot[i]);
        range_prev = ceres::sqrt(range_prev);

        T range_curr = T(0);
        for(int i=0; i<3; i++)
            range_curr += (curr_state[i] - sv_position_curr_rot[i]) * (curr_state[i] - sv_position_curr_rot[i]);
        range_curr = ceres::sqrt(range_curr);

        T clock_bias_diff = (curr_state[3 + satellite_type_] - prev_state[3 + satellite_type_]);
        T tdcp_pred = range_curr - range_prev + clock_bias_diff;
        

        // T tdcp_measure = (T(c) / L1_frequency) * T(pseudorange_);
        T tdcp_measure = T(pseudorange_);
        residual[0] = (tdcp_pred - tdcp_measure) * T(sqrt(weight_));
        return true;
    }

private:
    std::vector<double> sv_position_prev_;
    std::vector<double> sv_position_curr_;
    double pseudorange_;
    int satellite_type_;
    double weight_;
};

class NumTDCPFactorCostFunctor {
public:
    NumTDCPFactorCostFunctor(const std::vector<double>& sv_position_prev,
                          const std::vector<double>& sv_position_curr,
                          double pseudorange,
                          int satellite_type,
                          double weight)
        : sv_position_prev_(sv_position_prev),
          sv_position_curr_(sv_position_curr),
          pseudorange_(pseudorange),
          satellite_type_(satellite_type),
          weight_(weight) {}

    bool operator()(const double* const prev_state, const double* const curr_state, double* residual) const {
        const double c = 299792458.0;  // 빛의 속도 (m/s)
        const double omega = 7.292115e-5;  // 지구 자전 속도 (rad/s)

        double L1_frequency = GetL1Frequency(satellite_type_);  // L1 신호 주파수 (Hz)

        double sv_position_prev_rot[3];
        double sv_position_curr_rot[3];

        RotateSatellitePosition(prev_state, sv_position_prev_, sv_position_prev_rot);
        RotateSatellitePosition(curr_state, sv_position_curr_, sv_position_curr_rot);

        double los_vector_prev[3];
        double los_vector_curr[3];
        CalculateLOS(prev_state, sv_position_prev_rot, los_vector_prev);
        CalculateLOS(curr_state, sv_position_curr_rot, los_vector_curr);

        double D = 0.0;
        double g = 0.0;

        for (int i = 0; i < 3; ++i) {
            D += los_vector_curr[i] * sv_position_curr_rot[i] - los_vector_prev[i] * sv_position_prev_rot[i];
            g += los_vector_curr[i] * prev_state[i] - los_vector_prev[i] * prev_state[i];
        }

        double tdcp_pred = (curr_state[3 + satellite_type_] - prev_state[3 + satellite_type_]);
        for (int i = 0; i < 3; ++i) {
            tdcp_pred += -los_vector_prev[i] * (curr_state[i] - prev_state[i]);
        }

        double tdcp_measure = (c / L1_frequency) * pseudorange_ - D + g;

        residual[0] = (tdcp_pred - tdcp_measure) * sqrt(weight_);

        return true;
    }

private:
    std::vector<double> sv_position_prev_;
    std::vector<double> sv_position_curr_;
    double pseudorange_;
    int satellite_type_;
    double weight_;
};



class TimeCorrelationFactorCostFunctor {
    public:
        TimeCorrelationFactorCostFunctor(double tau, double weight)
            :tau_(tau), weight_(weight) {}

        template <typename T>
        bool operator()(const T* const previous_noise, const T* const current_noise, T* residual) const {

            T* current_noise_est = new T[1]; 
            double dt = 1.0;
            current_noise_est[0] = exp(-dt / tau_)*previous_noise[0];
            residual[0] = abs(current_noise_est[0] - current_noise[0]) * T(sqrt(weight_));

            return true;
        }

    private:
        const double tau_;
        const double weight_;

};


}

