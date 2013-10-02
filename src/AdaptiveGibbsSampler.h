/*
 * AdaptiveGibbsSampelr.h
 *
 *  Created on: Sep 23, 2013
 *      Author: hoblitz
 *
 *  A class for generating samples from a log concave density on (0, Inf)
 */

#ifndef ADAPTIVEGIBBSSAMPLER_H_
#define ADAPTIVEGIBBSSAMPLER_H_

#include "Hull.h"

template <class T>
class AdaptiveGibbsSampler {
public:
    Hull<T> hull;
    double x0, x1;
    AdaptiveGibbsSampler() : x0(1.0), x1(2.0) {};
    AdaptiveGibbsSampler(double const x0_init, double const x1_init) {
        x0 = x0_init;
        x1 = x1_init;
    };

    double sample(RngStream rng, double const * const pdf_args);
};

template <class T>
double AdaptiveGibbsSampler<T>::sample(RngStream rng, double const * const pdf_args) {
    hull.initialize(x0, x1, pdf_args);
    int hull_sample_success;
    double x_star;
    hull_sample_success = hull.drawSample(rng, x_star);

    if (hull_sample_success == 0) {
        hull.printHull();
        Rcpp::stop("Maximum iterations reached in adaptive rejection sampler");
    }
    int segment_idx;
    x0 = hull.inverseCdf(.15, segment_idx);
    x1 = hull.inverseCdf(.85, segment_idx);
    hull.reset();
    return x_star;
}

#endif /* ADAPTIVEGIBBSSAMPLER_H_ */
