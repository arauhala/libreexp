/*
 * exp.h
 *
 *  Created on: Oct 21, 2010
 *      Author: arauhala
 */

#ifndef REEXP_EXP_H_
#define REEXP_EXP_H_

#include "learner.h"
#include "pred.h"


namespace reexp {

	/*
	 * Problems and have they are solved:
	 * ----------------------------------
	 *
	 * We have two demands:
	 *
	 *     * Ability to have different language versions
	 *     * Good performance requires to modification instead of copying
	 *
	 * ->
	 *
	 *     * Re-expression is done to <org, data, stats> triplets.
	 *
	 *     * Orgs can be copied. We do re-expression on a copy of original
	 *       org.
	 *
	 *     * We can get the expressions delta between original org and new org
	 *
	 *     * we can prepare another <org2, data2, stats2> triplet and apply
	 *       delta there.
	 *
	 *
	 * Problem settings:
	 *
	 *    * To have proper bit indexes, we need to combine
	 *       * expression's context
	 *       * and the expression data begin offset
	 *
	 *    * Data offset depends of data and other expressions.
	 *      In this sense, it is <data, offset> specific
	 *
	 * It seems that there is the abstract org, which is unusable as such.
	 * It can be combined with data to resolve both bit indexes and bits
	 * themselves. In this sense we have three different stages:
	 *
	 *    * org
	 *    * data
	 *    * data+org (resolved org in a sense)
	 */

/*	template <typename P>
	struct ExampleOffset {
		int sampleCopies_;
		int bitmapCopies_;
		int operator() (const cvec<P>& dim) const {
			return sampleCopies_ * dim[0] + bitmapCopies_ + b_* dim[1] * dim[2];
		}

		ExampleOffset& operator+=(const ExampleOffset<P>& offset) {
			sampleCopies_ += offset_.sampleCopies_;
			bitmapCopies_ += offset_.bitmapCopies_;
			return *this;
		}
	};

	struct ExampleProblem {
		typedef ExampleOffset<ExampleProblem> dim_to_offset;
		static const int DIM = 3;
	};*/

}

#endif /* REEXP_EXP_H_ */
