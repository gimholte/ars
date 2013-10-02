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
        return (shape - 1.0) * log(x) - x * rate;
    };
    double pdfDeriv(const double x) {
        return (shape - 1.0) / x - rate;
    };
};


#endif /* LOGDENSITY_H_ */
