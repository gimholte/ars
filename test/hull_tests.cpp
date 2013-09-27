/*
 * all_tests.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: hoblitz
 */

#include "gtest/gtest.h"
#include <iostream>
#include "ars.h"

namespace {
class HullTest : public ::testing::Test {
protected:
    double x0, x1, h_x0, h_x1, hp_x0, hp_x1, z0;
    Hull<GammaDistribution> hull;
    void SetUp() {
        double gamma_args[] = {2.0, 2.0};
        hull.initialize(1.0, 3.0, gamma_args);
        x0 = 1.0;
        x1 = 3.0;
        h_x0 = -2.0;
        h_x1 = log(3.0) - 6.0;
        hp_x0 = 1.0 - 2.0;
        hp_x1 = 1.0 / 3.0 - 2.0;
        z0 = (h_x1 - h_x0 - x1 * hp_x1 + x0 * hp_x0) / (hp_x0 - hp_x1);
    }
};
}

TEST_F(HullTest, HullInitExtendRightPoint) {
    double gamma_args[] = {2.0, 2.0};
    hull.initialize(.2, .3, gamma_args);
    ASSERT_EQ(.6, hull.hull[1].left_x);
}

TEST_F(HullTest, HullInit) {
    // covers method initialize, initializeHullMax
    ASSERT_EQ(2, hull.num_hull_segments);
    EXPECT_DOUBLE_EQ(h_x0, hull.hull[0].h_x);
    EXPECT_DOUBLE_EQ(h_x1, hull.hull[1].h_x);
    EXPECT_DOUBLE_EQ(hp_x0, hull.hull[0].hprime_x);
    EXPECT_DOUBLE_EQ(hp_x1, hull.hull[1].hprime_x);

    EXPECT_DOUBLE_EQ(z0, hull.hull[0].z);
    EXPECT_EQ(INFINITY, hull.hull[1].z);
    EXPECT_DOUBLE_EQ(hull.upper_hull_max, -1.0);
    SUCCEED();
}

TEST_F(HullTest, IntegrateHullSegments) {
    // covers methods integrateSegment and normalizeHull
    double integral, integral1;
    double result;
    integral = (exp(h_x0 - hull.upper_hull_max - x0 * hp_x0) / hp_x0) * (exp(hp_x0 * z0) - exp(0.0));
    integral = log(integral);
    result = hull.integrateSegment(hull.hull[0], 0.0);
    EXPECT_DOUBLE_EQ(integral, result);
    EXPECT_DOUBLE_EQ(integral, hull.hull[0].raw_integral);
    EXPECT_DOUBLE_EQ(integral, hull.hull[0].raw_cumulative_integral);

    integral1 = (exp(h_x1 - hull.upper_hull_max - x1 * hp_x1) / fabs(hp_x1)) * (exp(hp_x1 * z0));
    integral1 = log(integral1);
    result = hull.integrateSegment(hull.hull[1], hull.hull[0].z);
    EXPECT_DOUBLE_EQ(integral1, result);
    EXPECT_DOUBLE_EQ(integral1, hull.hull[1].raw_integral);
    EXPECT_DOUBLE_EQ(logspaceAdd(integral, integral1), hull.hull[1].raw_cumulative_integral);

    EXPECT_EQ(0.0, hull.hull[1].cum_prob);
    EXPECT_EQ(integral - logspaceAdd(integral, integral1), hull.hull[0].cum_prob);
}

TEST_F(HullTest, BinarySearch) {
    Hull<GammaDistribution> search_test_hull;
    search_test_hull.hull[0].cum_prob = -4.0;
    search_test_hull.hull[1].cum_prob = -3.0;
    search_test_hull.hull[2].cum_prob = -2.0;
    search_test_hull.hull[3].cum_prob = -1.0;
    search_test_hull.hull[4].cum_prob = 0.0;
    search_test_hull.num_hull_segments = 5;
    int segment;

    for(int i = 0; i < 5; i++) {
        double x = -(4.0 - i + .5);
        segment = search_test_hull.argBinarySearch(x, 0, 4);
        EXPECT_EQ(i, segment);
    }

    segment = hull.argBinarySearch(log(.001), 0, 1);
    ASSERT_EQ(segment, 0);

    segment = hull.argBinarySearch(log(.999), 0, 1);
    ASSERT_EQ(segment, 1);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


