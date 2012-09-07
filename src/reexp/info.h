/*
 * info.h
 *
 *  Created on: Dec 17, 2010
 *      Author: arauhala
 */

#ifndef INFO_H_
#define INFO_H_

#include "math.h"

namespace explib {

	inline double information(double pS) {
		return (pS == 1. || pS == 0 || isnan(pS)) ? 0 : -pS * log2(pS);
	}

	inline double information(double s1, double s2) {
		return information(s1) + information(s2);
	}

	inline double information(double s1, double s2, double s3, double s4) {
		return information(s1)
			 + information(s2)
		     + information(s3)
		     + information(s4);

	}

	inline double entropy(double p) {
		return information(p) + information(1-p);
	}

	inline double estimateInfo(double realP, double estimateP) {
		return realP == 0 || isnan(realP) ? 0 : -realP * log2(estimateP);
	}

	inline double estimateEntropy(double realP, double estimateP) {
		return estimateInfo(realP, estimateP)
             + estimateInfo(1-realP, 1-estimateP);
	}

	/*
	 * Analytical solution for p's expectation
	 * value, when assuming flat distribution
	 * a priori.
	 */
	inline double hatP(int k, int n) {
		return (k + 1) / double(n + 2);
	}

	/*
	 * Probability average can be use as simplistic
	 * biasless estimate.
	 */
	inline double p(int k, int n){
		return k / (double)n;
	}

}

#endif /* INFO_H_ */
