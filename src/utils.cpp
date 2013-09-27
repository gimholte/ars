/*
 * utils.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: hoblitz
 */

#include "utils.h"

double logspaceAdd(const double loga, const double logb) {
    if (!R_FINITE(loga))
        return logb;
    if (loga > logb)
        return logspaceAdd(logb, loga);
    return logb + log1p(exp(loga - logb));
}


