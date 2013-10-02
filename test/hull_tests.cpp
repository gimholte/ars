/*
 * all_tests.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: hoblitz
 */

#include "gtest/gtest.h"
#include "AdaptiveGibbsSampler.h"
#include "LogDensity.h"

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

class ArsTest : public ::testing::Test {
protected:
    AdaptiveGibbsSampler<GammaDistribution> gamma_sampler;
    void SetUp() {
        gamma_sampler = AdaptiveGibbsSampler<GammaDistribution>(1.0, 3.0);
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

TEST_F(HullTest, InverseCdf) {
    // checks that cdf/ inverse cdf functions are self-consistent
    // within small absolute error.
    int i, seg_idx;
    double test_percentiles[] = {.000001, .1, .2, .3, .4, .5, .6, .7, .9, .999999};
    double q, p;
    for (i = 0; i < 10; i++) {
        q = hull.inverseCdf(test_percentiles[i], seg_idx);
        EXPECT_GT(q, 0.0);
        p = hull.cdf(q);
        EXPECT_NEAR(test_percentiles[i], p, .0000000001);
    }
}

TEST_F(HullTest, InsertSegment) {
    double xnew = 1.8;
    double hxnew = hull.dist.pdf(xnew);
    hull.printHull();
    hull.insertSegment(xnew, hxnew, 1);
    hull.printHull();
    hull.insertSegment(.5, hull.dist.pdf(.5), 0);
    hull.printHull();
    hull.insertSegment(3.1, hull.dist.pdf(3.1), 3);
    hull.printHull();
}

TEST_F(ArsTest, TestInit) {
    EXPECT_EQ(1.0, gamma_sampler.x0);
    EXPECT_EQ(3.0, gamma_sampler.x1);
}

TEST_F(ArsTest, GenerateSample) {
    double gamma_args[] = {2.0, 2.0};
    RngStream rng = RngStream_CreateStream("");
    double samp;
    samp = gamma_sampler.sample(rng, gamma_args);
    EXPECT_GT(samp, 0.0);
    EXPECT_GT(gamma_sampler.x0, 0.0);
    EXPECT_GT(gamma_sampler.x1, gamma_sampler.x0);
    EXPECT_EQ(gamma_sampler.hull.num_hull_segments, 2);
    EXPECT_EQ(gamma_sampler.hull.upper_hull_max, -INFINITY);
}

TEST_F(ArsTest, GammaMomentTest) {
    double gamma_args[] = {2.9, 2.0};
    RngStream rng = RngStream_CreateStream("");
    double samp, mean = 0.0;
    int n = 1, nsamp = 100000;
    while(n <= nsamp) {
        samp = gamma_sampler.sample(rng, gamma_args);
        ASSERT_GT(samp, 0.0);
        mean = mean * (n - 1.0) / n + samp / n;
        n++;
    }
    double true_mean = gamma_args[0] / gamma_args[1];
    double true_se = sqrt(true_mean / gamma_args[1] / nsamp);
    EXPECT_GT(mean, true_mean - 2.6 * true_se);
    EXPECT_LT(mean, true_mean + 2.6 * true_se);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


