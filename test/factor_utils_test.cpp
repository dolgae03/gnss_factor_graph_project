#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
#include <cmath>
#include "../include/factors.h"

// 테스트할 함수들 포함
double GetL1Frequency(int satellite_type);
template <typename T>
void RotateSatellitePosition(const T* user_state, const std::vector<double>& satellite_position, T* satellite_position_rot);
template <typename T>
void CalculateLOS(const T* user_state, const T* satellite_position, T* los_vector);

// GetL1Frequency에 대한 테스트
TEST(SatelliteFunctionsTest, GetL1FrequencyTest) {
    EXPECT_DOUBLE_EQ(GetL1Frequency(0), 1575.42e6);
    EXPECT_DOUBLE_EQ(GetL1Frequency(1), 1575.42e6);
    EXPECT_DOUBLE_EQ(GetL1Frequency(2), 1561.098e6);
    
    EXPECT_THROW(GetL1Frequency(3), std::invalid_argument);
}

// RotateSatellitePosition에 대한 테스트
TEST(SatelliteFunctionsTest, RotateSatellitePositionTest) {
    double user_state[3] = {1.0, 2.0, 3.0};
    std::vector<double> satellite_position = {4.0, 5.0, 6.0};
    double satellite_position_rot[3];

    RotateSatellitePosition(user_state, satellite_position, satellite_position_rot);

    // 예상 결과와 실제 결과를 비교
    EXPECT_NEAR(satellite_position_rot[0], 4.0, 1e-9);
    EXPECT_NEAR(satellite_position_rot[1], 5.000000119, 1e-9);
    EXPECT_NEAR(satellite_position_rot[2], 6.0, 1e-9);
}

// CalculateLOS에 대한 테스트
TEST(SatelliteFunctionsTest, CalculateLOSTest) {
    double user_state[3] = {1.0, 2.0, 3.0};
    double satellite_position[3] = {4.0, 5.0, 6.0};
    double los_vector[3];

    CalculateLOS(user_state, satellite_position, los_vector);

    // 예상 결과: 벡터의 방향을 확인
    EXPECT_NEAR(los_vector[0], 0.57735, 1e-5);
    EXPECT_NEAR(los_vector[1], 0.57735, 1e-5);
    EXPECT_NEAR(los_vector[2], 0.57735, 1e-5);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
