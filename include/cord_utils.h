#ifndef CORD_UTILS_H
#define CORD_UTILS_H

#include <vector>

namespace coordinate {
    // WGS84 타원체 상수
    const double WGS84_A = 6378137.0;         // 반경 (semi-major axis)
    const double WGS84_E = 8.1819190842622e-2; // 이심률 (eccentricity)

    // ECEF 좌표를 LLA 좌표로 변환하는 함수
    std::vector<double> ecefToLLA(const std::vector<double> &ecef_position);
}

#endif // CORD_UTILS_H