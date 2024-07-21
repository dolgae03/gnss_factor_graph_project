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
        PesudorangeFactorCostFunctor(std::vector<double> &sv_position, double pesudorange): sv_position_(sv_position), pesudorange_(pesudorange) {}

        template <typename T>
        bool operator()(const T* state, T* residual) const {
            T dx = state[0] - T(sv_position_[0]);
            T dy = state[1] - T(sv_position_[1]);
            T dz = state[2] - T(sv_position_[2]);
            T estimated_pr = ceres::sqrt(dx * dx + dy * dy + dz * dz);

            residual[0] = T(pesudorange_) - estimated_pr;

            return true;
        }

    private:
        std::vector<double> sv_position_;
        double pesudorange_;
};

};