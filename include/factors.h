#ifndef FACTORS_H
#define FACTORS_H

#include <vector>

namespace factor {

class ConstantClockBiasFactorCostFunctor {
public:
    ConstantClockBiasFactorCostFunctor(int satellite_type, double sigma);
    template <typename T>
    bool operator()(const T* const state1, const T* const state2, T* residual) const;

private:
    int satellite_type_;
    double sigma_;
};

class DiffPesudorangeFactorCostFunctor {
public:
    DiffPesudorangeFactorCostFunctor(const std::vector<double>& ref_position, const std::vector<double>& sv_position, double pesudorange, int satellite_type, double sigma);
    template <typename T>
    bool operator()(const T* state, T* residual) const;

private:
    std::vector<double> ref_position_;
    std::vector<double> sv_position_;
    double pesudorange_;
    double sigma_;
    int satellite_type_;
};

class DopplerFactorCostFunctor {
public:
    DopplerFactorCostFunctor(const std::vector<double>& sv_position, const std::vector<double>& sv_velocity, double doppler, double time_interval, int satellite_type, double sigma);
    template <typename T>
    bool operator()(const T* const prev_state, const T* const next_state, T* residual) const;

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
    TDCPFactorCostFunctor(const std::vector<double>& sv_position_prev, const std::vector<double>& sv_position_curr, double pseudorange, int satellite_type, double sigma);
    template <typename T>
    bool operator()(const T* const prev_state, const T* const curr_state, T* residual) const;

private:
    std::vector<double> sv_position_curr_;
    std::vector<double> sv_position_prev_;
    double pseudorange_;
    int satellite_type_;
    double sigma_;
};
}

#endif // FACTORS_H