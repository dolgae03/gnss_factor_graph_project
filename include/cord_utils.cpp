#include <cmath>
#include <iostream>
#include <vector>
#include "cord_utils.h"

// WGS84 타원체 상수
const double WGS84_A = 6378137.0;         // 반경 (semi-major axis)
const double WGS84_E = 8.1819190842622e-2; // 이심률 (eccentricity)

// ECEF 좌표를 LLA 좌표로 변환하는 함수

namespace coordinate{
    std::vector<double> ecefToLLA(const std::vector<double> &ecef_position) {
        double x = ecef_position[0];
        double y = ecef_position[1];
        double z = ecef_position[2];
        double a = WGS84_A;
        double e = WGS84_E;

        double b = sqrt(a * a * (1 - e * e));
        double ep = sqrt((a * a - b * b) / (b * b));
        double p = sqrt(x * x + y * y);
        double th = atan2(a * z, b * p);

        double lat, lon, alt;

        lon = atan2(y, x);
        lat = atan2((z + ep * ep * b * pow(sin(th), 3)), (p - e * e * a * pow(cos(th), 3)));

        double N = a / sqrt(1 - e * e * sin(lat) * sin(lat));
        alt = p / cos(lat) - N;

        // 위도, 경도를 도 단위로 변환
        lat = lat * 180.0 / M_PI;
        lon = lon * 180.0 / M_PI;

        return {lat, lon, alt};
    }
};