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

namespace factor {
class ConstantClockBiasFactorCostFunctor {
    public:
        ConstantClockBiasFactorCostFunctor(int satellite_type, double weight)
            : satellite_type_(satellite_type), weight_(weight) {}

        template <typename T>
        bool operator()(const T* const state1, const T* const state2, T* residual) const {
            // state1과 state2의 clock bias 차이 계산
            T clock_bias1 = state1[satellite_type_];
            T clock_bias2 = state2[satellite_type_];
            
            residual[0] = (clock_bias1 - clock_bias2) * T(weight_);
            return true;
        }

    private:
        int satellite_type_;
        double weight_;
};


class DiffPesudorangeFactorCostFunctor {
   public:
        DiffPesudorangeFactorCostFunctor(const std::vector<double>& ref_position, const std::vector<double>& sv_position, double pesudorange, int satellite_type, double weight)
            : ref_position_(ref_position), sv_position_(sv_position), pesudorange_(pesudorange), satellite_type_(satellite_type), weight_(weight) {}

        template <typename T>
        bool operator()(const T* state, T* residual) const {
            const double c = 299792458.0; // Speed of light in m/s
            const double omega = 7.2921151467e-5; // Earth's rotation rate in rad/s

            T rotated_sv_pos[3];
            RotateSatellitePosition(state, sv_position_, rotated_sv_pos);
            
            T range_to_sv = T(0);
            for(int i=0; i<3; i++)
                range_to_sv += (state[i] - rotated_sv_pos[i]) * (state[i] - rotated_sv_pos[i]);

            T clock_bias = state[3 + satellite_type_] * T(c);
            range_to_sv = ceres::sqrt(range_to_sv) + clock_bias;


            T range_to_ref = T(0);
            for(int i=0; i<3; i++)
                range_to_ref += (T(ref_position_[i]) - rotated_sv_pos[i]) * (T(ref_position_[i])  - rotated_sv_pos[i]);

            range_to_ref = ceres::sqrt(range_to_ref);

            T estimated_pr_diff = range_to_sv - range_to_ref;

            residual[0] = (T(pesudorange_) - estimated_pr_diff) * T(weight_);

            return true;
        }

    private:
        std::vector<double> ref_position_;
        std::vector<double> sv_position_;
        double pesudorange_;
        double weight_;
        int satellite_type_;
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

        T L1_frequency = T(GetL1Frequency(satellite_type_));  // L1 signal frequency in Hz

        T sv_position_prev_rot[3];
        T sv_position_curr_rot[3];
        RotateSatellitePosition(prev_state, sv_position_prev_, sv_position_prev_rot);
        RotateSatellitePosition(curr_state, sv_position_curr_, sv_position_curr_rot);

        T los_vector_prev[3];
        T los_vector_curr[3];
        CalculateLOS(prev_state, sv_position_prev_rot, los_vector_prev);
        CalculateLOS(curr_state, sv_position_curr_rot, los_vector_curr);

        T D = T(0);
        T g = T(0);

        for (int i = 0; i < 3; ++i) {
            D += los_vector_curr[i] * T(sv_position_curr_rot[i]) - los_vector_prev[i] * T(sv_position_prev_rot[i]);
            g += los_vector_curr[i] * prev_state[i] - los_vector_prev[i] * prev_state[i];
        }

        T tdcp_pred = T(c) * (curr_state[3 + satellite_type_] - prev_state[3 + satellite_type_]);
        for (int i = 0; i < 3; ++i) //
            tdcp_pred += T(-1.0) * los_vector_curr[i] * (curr_state[i] - prev_state[i]);

        T tdcp_measure = (T(c) / L1_frequency) * T(pseudorange_) - D + g;

        residual[0] = (tdcp_pred - tdcp_measure) * T(weight_);

        return true;
    }

private:
    std::vector<double> sv_position_prev_;
    std::vector<double> sv_position_curr_;
    double pseudorange_;
    int satellite_type_;
    double weight_;
};
}