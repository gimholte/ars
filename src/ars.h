/*
 * ars.h
 *
 *  Created on: Sep 23, 2013
 *      Author: hoblitz
 *
 *  A class for generating samples from a log concave density on (0, Inf)
 *  using standard template library containers
 */

#ifndef ARS_H_
#define ARS_H_

#include <Rcpp.h>
#include <R.h>
#include <list>
#include "RngStream.h"

#define MAX_HULL_SIZE 100
#define MAX_HULL_TRIALS 1000
#define HULL_SAMPLE_ACCEPT 1
#define HULL_SAMPLE_REJECT 0

class HullSegment {
public:
    double left_x;
    double h_x;
    double hprime_x;
    double z;
    double raw_integral;
    double raw_cumulative_integral;
    double prob;
    double cum_prob;
};

template <class T>
class Hull {
    T dist;
    double upper_hull_max;
    int num_hull_segments;
    HullSegment hull[MAX_HULL_SIZE];

    void initializeHullMax();
    void normalizeHull();
    double integrateSegment(HullSegment const & segment, const double z_prev);
    void setZcoord(HullSegment & segleft, HullSegment const & segright);
    void updateZ(const int insert_idx);
    void renormalizeHull(const int start_idx);
    double inverseCDF(const double u, int & seg_idx);
    int argBinarySearch(const double log_u, int lower, int upper);
    int squeezeTest(RngStream rng, double & x_trial, double & hx_trial, int segment_idx);
public:
    Hull() {
        upper_hull_max = -INFINITY;
        num_hull_segments = 2;
    }
    int drawSample(RngStream rng, double & x_sample);
    void initialize(const double x0, const double x1, double const * const pdf_args);
    double inverseCdf(const double x);
    void insertSegment(const double x, const double h_x, const int left_idx);
};

/* base class for densities we'd like to sample*/
class LogDensity {
public:
    virtual double pdf(const double x) =0;
    virtual double pdfDeriv(const double x) =0;
    virtual ~LogDensity() =0;
    virtual void setParameters(double const * const args) =0;
};

class GammaDistribution: public LogDensity {
    double shape, rate;
public:
    void setParameters(double const * const args) {
        shape = args[0];
        rate = args[1];
    }

    double pdf(const double x) {
        return (shape - 1.0) * log(x) - x * rate;
    }

    double pdfDeriv(const double x) {
        return (shape - 1.0) / x - rate;
    }
};

template <class T>
class AdaptiveGibbsSampler {
    Hull<T> hull;
    double x0, x1;

    //int squeezeTest(RngStream rng, double & xstar, double & h_xstar);
public:
    AdaptiveGibbsSampler(double const x0_init, double const x1_init) {
        x0 = x0_init;
        x1 = x1_init;
    }
    double sample(RngStream rng, double const * const args);
};

double logspaceAdd(const double loga, const double logb);

#endif /* ARS_H_ */
