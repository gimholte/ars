/*
 * LogDensityDup.h
 *
 *  Created on: Oct 2, 2013
 *      Author: hoblitz
 */

#ifndef LOGDENSITY_H_
#define LOGDENSITY_H_

#include <math.h>

class LogDensity {
public:
    virtual double pdf(const double x) =0;
    virtual double pdfDeriv(const double x) =0;
    virtual ~LogDensity() {};
    virtual void setParameters(double const * const args) =0;
};

class GammaDistribution : public LogDensity {
    double shape, rate;
public:
    void setParameters(double const * const args) {
        shape = args[0];
        rate = args[1];
    };
    double pdf(const double x) {
        return (shape - (double) 1.0) * log(x) - x * rate;
    };
    double pdfDeriv(const double x) {
        return (shape - (double) 1.0) / x - rate;
    };
};

class WeibullDistribution : public LogDensity {
    double k, ell;
public:
    void setParameters(double const * const args) {
        k = args[0];
        ell = args[1];
    }
    double pdf(const double x) {
        return (k - (double) 1.0) * log(x) - pow(x / ell, k);
    };
    double pdfDeriv(const double x) {
        return (k - (double) 1.0) / x - k * pow(x / ell, k - (double) 1.0) / ell;
    }

};


#endif /* LOGDENSITY_H_ */
