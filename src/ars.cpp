/*
 * ars.cpp
 *
 *  Created on: Sep 23, 2013
 *      Author: hoblitz
 */

#include "ars.h"
using namespace Rcpp;

double logspaceAdd(const double loga, const double logb) {
    if (!R_FINITE(loga))
        return logb;
    if (loga > logb)
        return logspaceAdd(logb, loga);
    return logb + log1p(exp(loga - logb));
}

template <class T>
void Hull<T>::initialize(const double x0, const double x1, double const * const pdf_args) {
    // before all else, initialize the distribution parameters
    dist.setParameters(pdf_args);

    hull[0].left_x = x0;
    hull[0].h_x = dist.pdf(x0);
    hull[0].hprime_x = dist.pdf(x0);

    hull[1].left_x = x1;
    double hpx1 = dist.pdfDeriv(x1);
    // extend right point until hull derivative is negative
    while(hpx1 >= 0.0) {
        hull[1].left_x *= 2.0;
        hpx1 = dist.pdfDeriv(hull[1].left_x);
    }
    hull[1].hprime_x = hpx1;
    hull[1].h_x = dist.pdf(hull[1].left_x);

    num_hull_segments = 2;
    // set z coordinates;
    setZcoord(hull[0], hull[1]);
    hull[1].z = INFINITY;
    initializeHullMax();
}

template <class T>
void Hull<T>::setZcoord(HullSegment & seg_left, HullSegment const & seg_right) {
    const double xj = seg_left.left_x;
    const double hxj = seg_left.h_x;
    const double hpxj = seg_left.hprime_x;

    const double xj1 = seg_right.left_x;
    const double hxj1 = seg_right.h_x;
    const double hpxj1 = seg_right.hprime_x;

    if ((hpxj - hpxj1) > 0.0) {
        seg_left.z = (hxj1 - hxj - xj1 * hpxj1 + xj * hpxj) / (hpxj - hpxj1);
    } else
        seg_left.z = (hxj + hxj1) / 2.0;
    return;
}

template <class T>
void Hull<T>::initializeHullMax() {
    const double x0 = hull[0].left_x;
    const double hx0 = hull[0].h_x;
    const double hpx0 = hull[0].hprime_x;
    const double z0 = hull[0].z;

    upper_hull_max = hpx0 > 0.0 ? hx0 + (z0 - x0) * hpx0 : hx0 - x0 * hpx0;
}

template <class T>
void Hull<T>::normalizeHull() {
    double cumulative_integral = -INFINITY;
    double segment_integral;
    double z_prev = 0.0;
    int i;
    // integrate each segment
    for(i = 0; i < num_hull_segments; i++) {
        segment_integral = integrateSegment(hull[i], z_prev);
        cumulative_integral = logspaceAdd(segment_integral, cumulative_integral);
        hull[i].raw_integral = segment_integral;
        hull[i].raw_cumulative_integral = cumulative_integral;
        z_prev = hull[i].z;
    }
    // normalize segments to sum to one
    for(i = 0; i < num_hull_segments; i++) {
        hull[i].prob = hull[i].raw_integral - cumulative_integral;
        hull[i].cum_prob = hull[i].raw_cumulative_integral - cumulative_integral;
    }
    return;
}

template <class T>
int Hull<T>::drawSample(RngStream rng, double & x_sample) {
    int num_trials = 0;
    int segment_idx;
    int test_outcome;
    double x_trial, hx_trial;
    while (num_trials < MAX_HULL_TRIALS) {
        // sample the upper hull
        x_trial = inverseCDF(RngStream_RandU01(rng), segment_idx);
        // apply squeeze test
        test_outcome = squeezeTest(rng, x_trial, hx_trial, segment_idx);
        if (test_outcome == HULL_SAMPLE_ACCEPT) {
            x_sample = x_trial;
            return 1;
        }
        // sample rejected, insert into hull
        insertSegment(x_trial, hx_trial, segment_idx);
        num_trials++;
    }
    return 0;
}

template <class T>
int Hull<T>::argBinarySearch(const double log_u, const int lower, const int upper) {
    const int mid = (lower + upper) / 2;
    if (mid == lower)
        return log_u < hull[lower].cum_prob ? lower : upper;
    if (log_u < hull[mid].cum_prob) {
        return argBinarySearch(log_u, lower, mid);
    } else {
        return argBinarySearch(log_u, mid, upper);
    }
}

template <class T>
int Hull<T>::squeezeTest(RngStream rng, double & x_trial, double & hx_trial, int segment_idx) {
    double w = log(RngStream_RandU01(rng));
    const double x = hull[segment_idx].left_x;
    const double h_x = hull[segment_idx].h_x;
    const double hp_x = hull[segment_idx].hprime_x;
    const double z = hull[segment_idx].z;
    double upr_hull_val = h_x + (x_trial - x) * hp_x;
    double lwr_hull_val;

    //compute lower hull value
    if ((x_trial < x) && (segment_idx > 0)) {
        const double x_lower = hull[segment_idx - 1].left_x;
        const double hx_lower = hull[segment_idx - 1].h_x;
        lwr_hull_val = ((x - x_trial) * h_x + (x_trial - x_lower) * hx_lower) /
                (x - x_lower);
    } else if ((x_trial > x) && (segment_idx < (num_hull_segments - 1))) {
        const double x_upper = hull[segment_idx + 1].left_x;
        const double hx_upper = hull[segment_idx + 1].h_x;
        lwr_hull_val = ((x_upper - x_trial) * hx_upper + (x_trial - x) * h_x) /
                (x_upper - x);
    } else {
        lwr_hull_val = -INFINITY;
    }

    // squeezing tests
    if (w <= lwr_hull_val - upr_hull_val) {
        return HULL_SAMPLE_ACCEPT;
    }
    hx_trial = dist.pdf(x_trial);
    if (hx_trial - upr_hull_val) {
        return HULL_SAMPLE_ACCEPT;
    }
    return HULL_SAMPLE_REJECT;
}

template <class T>
double Hull<T>::integrateSegment(HullSegment const & segment, const double z_prev) {
    const double xj = segment.left_x;
    /* subtract hull max to bound upper hull between zero and one */
    const double hxj = segment.h_x - upper_hull_max;
    const double hpxj = segment.hprime_x;
    const double zj = segment.z;

    if (hpxj == 0.0) {
        return log(zj - z_prev) + hxj;
    }

    const double pre_factor = hxj - (xj - z_prev) * hpxj - log(fabs(hpxj));
    double int_factor;
    if (hpxj > 0.0) {
        int_factor = log1p(-exp(hpxj * (z_prev - zj)));
    } else {
        int_factor = log1p(-exp(hpxj * (zj - z_prev)));
    }

    return pre_factor + int_factor;
}

template <class T>
double Hull<T>::inverseCDF(const double p, int & seg_idx) {
    seg_idx = argBinarySearch(log(p), 0, num_hull_segments - 1);
    /* seg_idx now contains the index of the hull segment such that
     * our sampled x_star is in the interval
     * (z[seg_idx - 1], z[seg_idx]]
     */
    const double x = hull[seg_idx].left_x;
    const double h_x = hull[seg_idx].h_x;
    const double hp_x = hull[seg_idx].hprime_x;
    const double z_prev = (seg_idx == 0) ? 0.0:hull[seg_idx - 1].z;
    const double cdf_prev_seg = (seg_idx == 0) ? 0.0:exp(hull[seg_idx - 1].cum_prob);
    const double p_remainder = p - cdf_prev_seg;

    const double x_star = log(p_remainder * exp(hull[num_hull_segments - 1].raw_cumulative_integral) *
            hp_x + exp((z_prev - x)* hp_x + h_x - upper_hull_max)) +
                    x*hp_x - (h_x - upper_hull_max);
    return x_star;
}

template <class T>
void Hull<T>::insertSegment(const double x_new, const double h_xnew, const int origin_idx) {
    if (num_hull_segments == MAX_HULL_SIZE)
        return;
    const double hp_xnew = dist.pdfDeriv(x_new);

    for(int k = num_hull_segments; k > origin_idx; k--) {
        hull[k] = hull[k - 1];
    }
    num_hull_segments++;
    // determine whether new segment is to left, or right, of origin segment.
    const int insert_idx = (hull[origin_idx].left_x < x_new) ? origin_idx + 1 : origin_idx;
    hull[insert_idx].left_x = x_new;
    hull[insert_idx].h_x = h_xnew;
    hull[insert_idx].hprime_x = hp_xnew;

    updateZ(insert_idx);
    renormalizeHull(insert_idx);
    return;
}

template <class T>
void Hull<T>::renormalizeHull(const int insert_idx) {
    double z_prev = 0.0;
    int idx_min = (insert_idx == 0) ? 0 : (insert_idx - 1);
    int idx_max = (insert_idx < num_hull_segments - 1) ?
            (num_hull_segments - 1) : (insert_idx + 1);
    double segment_integral, cumulative_integral = -INFINITY;

    int k;
    for(k = idx_min; k <= idx_max; k++) {
        if (k <= idx_max && k >= idx_min)
            hull[k].raw_integral = integrateSegment(hull[k], z_prev);
        segment_integral = hull[k].raw_integral;
        cumulative_integral = logspaceAdd(segment_integral, cumulative_integral);
        hull[k].raw_cumulative_integral = cumulative_integral;
        z_prev = hull[k].z;
    }

    for(k = 0; k < num_hull_segments; k++) {
        hull[k].prob = hull[k].raw_integral - cumulative_integral;
        hull[k].cum_prob = hull[k].raw_cumulative_integral - cumulative_integral;
    }

    return;
}


template <class T>
void Hull<T>::updateZ(const int insert_idx) {
    if (insert_idx == 0) {
        setZcoord(hull[insert_idx], hull[insert_idx + 1]);
    } else if (insert_idx == num_hull_segments - 1) {
        setZcoord(hull[insert_idx - 1], hull[insert_idx]);
        hull[insert_idx].z = INFINITY;
    } else {
        setZcoord(hull[insert_idx - 1], hull[insert_idx]);
        setZcoord(hull[insert_idx], hull[insert_idx + 1]);
    }
    return;
}

template <class T>
double AdaptiveGibbsSampler<T>::sample(RngStream rng, double const * const pdf_args) {
    hull.initializeHull(x0, x1, pdf_args);
    int hull_sample_success;
    double x_star;
    hull_sample_success = hull.drawSample(rng, x_star);
    if (!hull_sample_success) {
        Rcpp:stop("Maximum iterations reached in adaptive rejection sampler");
    }

    x0 = hull.inverseCDF(.15);
    x1 = hull.inverseCDF(.85);
    hull.reset();
    return x_star;
}

