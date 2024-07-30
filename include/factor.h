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

namespace factor {

class PesudorangeFactorCostFunctor {
    public:
        PesudorangeFactorCostFunctor(const std::vector<double>& sv_position, double pesudorange, double sigma)
            : sv_position_(sv_position), pesudorange_(pesudorange), sigma_(sigma) {}

        template <typename T>
        bool operator()(const T* state, T* residual) const {
            const double c = 299792458.0; // Speed of light in m/s
            T dx = state[0] - T(sv_position_[0]);
            T dy = state[1] - T(sv_position_[1]);
            T dz = state[2] - T(sv_position_[2]);
            T clock_bias = state[3] * T(c);

            T estimated_pr = ceres::sqrt(dx * dx + dy * dy + dz * dz) + clock_bias;

            residual[0] = (estimated_pr - T(pesudorange_)) / sigma_;

            return true;
        }

    private:
        std::vector<double> sv_position_;
        double pesudorange_;
        double sigma_;
};

class DiffPesudorangeFactorCostFunctor {
    public:
        DiffPesudorangeFactorCostFunctor(const std::vector<double>& ref_position, const std::vector<double>& sv_position, double pesudorange, double sigma)
            : ref_position_(ref_position), sv_position_(sv_position), pesudorange_(pesudorange), sigma_(sigma) {}

        template <typename T>
        bool operator()(const T* state, T* residual) const {
            const double c = 299792458.0; // Speed of light in m/s

            // State parameters
            T rx = state[0];
            T ry = state[1];
            T rz = state[2];
            T clock_bias = state[3] * T(c);

            // Reference position
            T ref_x = T(ref_position_[0]);
            T ref_y = T(ref_position_[1]);
            T ref_z = T(ref_position_[2]);

            // Satellite position
            T sx = T(sv_position_[0]);
            T sy = T(sv_position_[1]);
            T sz = T(sv_position_[2]);

            // Range to receiver
            T dx = rx - sx;
            T dy = ry - sy;
            T dz = rz - sz;
            T range_to_receiver = ceres::sqrt(dx * dx + dy * dy + dz * dz) + clock_bias;

            // Range to reference position
            T dref_x = ref_x - sx;
            T dref_y = ref_y - sy;
            T dref_z = ref_z - sz;
            T range_to_ref = ceres::sqrt(dref_x * dref_x + dref_y * dref_y + dref_z * dref_z);

            // Differential pseudorange
            T estimated_pr_diff = range_to_receiver - range_to_ref;

            residual[0] = (estimated_pr_diff - T(pesudorange_)) / T(sigma_);

            return true;
        }

    private:
        std::vector<double> ref_position_;
        std::vector<double> sv_position_;
        double pesudorange_;
        double sigma_;
};

class DopplerFactorCostFunctor {
    public:
        DopplerFactorCostFunctor(const std::vector<double>& sv_position, const std::vector<double>& sv_velocity, double doppler, double time_interval, int satellite_type, double sigma)
            : sv_position_(sv_position), sv_velocity_(sv_velocity), doppler_(doppler), time_interval_(time_interval), satellite_type_(satellite_type), sigma_(sigma) {}

        template <typename T>
        bool operator()(const T* const prev_state, const T* const next_state, T* residual) const {
            const double c = 299792458.0;  // speed of light in m/s
            T L1_frequency;  // L1 signal frequency in Hz
            
            // Define frequencies for each satellite type
            switch (satellite_type_) {
                case 1:  // GPS L1C
                    L1_frequency = T(1575.42e6);
                    break;
                case 2:  // GLONASS L1C
                    // Assuming k value is provided, k = -7,…+12
                    // This example uses k = 0 for simplicity
                    L1_frequency = T(1602.0e6);
                    break;
                case 3:  // GALILEO L1C
                    L1_frequency = T(1575.42e6);
                    break;
                case 4:  // QZSS L1C
                    L1_frequency = T(1575.42e6);
                    break;
                case 5:  // Beidou B2
                    L1_frequency = T(1207.14e6);
                    break;
                default:
                    // Default to GPS L1C if an unknown satellite type is provided
                    L1_frequency = T(1575.42e6);
                    break;
            }

            // Calculate the LOS vector from the user's position to the satellite's position
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
                user_velocity[i] = (next_state[i] - prev_state[i]) / (T(time_interval_) + next_state[3] - prev_state[3]);
        
            T relative_velocity[3];
            for (int i = 0; i < 3; ++i) 
                relative_velocity[i] = user_velocity[i] - T(sv_velocity_[i]);

            // Project relative velocity onto the LOS vector to get the relative speed
            T relative_speed = relative_velocity[0] * los_vector[0] +
                            relative_velocity[1] * los_vector[1] +
                            relative_velocity[2] * los_vector[2];

            // Convert Doppler shift to velocity
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

};