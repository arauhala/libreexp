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

	struct traits1d {
		static const int DIM = 1;
		static const int MAX_REL_VARS = 2;
	};

	struct traits2d {
		static const int DIM = 2;
		static const int MAX_REL_VARS = 2;
	};

	struct traits3d {
		static const int DIM = 3;
		static const int MAX_REL_VARS = 2;
	};

	struct traits4d {
		static const int DIM = 4;
		static const int MAX_REL_VARS = 2;
	};

	extern template class lang<traits1d>;
	extern template class lang<traits2d>;
	extern template class lang<traits3d>;
	extern template class lang<traits4d>;

	extern template class data<traits1d>;
	extern template class data<traits2d>;
	extern template class data<traits3d>;
	extern template class data<traits4d>;

	extern template struct rel_inputvars<traits1d>;
	extern template struct rel_inputvars<traits2d>;
	extern template struct rel_inputvars<traits3d>;
	extern template struct rel_inputvars<traits4d>;

	extern template struct rel_stats<traits1d>;
	extern template struct rel_stats<traits2d>;
	extern template struct rel_stats<traits3d>;
	extern template struct rel_stats<traits4d>;

	extern template class stats<traits1d>;
	extern template class stats<traits2d>;
	extern template class stats<traits3d>;
	extern template class stats<traits4d>;

	extern template class pred<traits1d>;
	extern template class pred<traits2d>;
	extern template class pred<traits3d>;
	extern template class pred<traits4d>;

	extern template void pred<traits1d>::getInputStateLogDep<std::ostream>(const rel_stats<traits1d>& rs,
											   	   	  	  	  	    	   const rel_inputvars<traits1d>& riv,
											   	   	  	  	  	    	   inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    	   std::ostream& explain) const;
	extern template void pred<traits2d>::getInputStateLogDep<std::ostream>(const rel_stats<traits2d>& rs,
											   	   	  	  	  	    	   const rel_inputvars<traits2d>& riv,
											   	   	  	  	  	           inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    	   std::ostream& explain) const;
	extern template void pred<traits3d>::getInputStateLogDep<std::ostream>(const rel_stats<traits3d>& rs,
											   	   	  	  	  	    	   const rel_inputvars<traits3d>& riv,
											   	   	  	  	  	    	   inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    	   std::ostream& explain) const;
	extern template void pred<traits4d>::getInputStateLogDep<std::ostream>(const rel_stats<traits4d>& rs,
											   	   	  	  	  	    	   const rel_inputvars<traits4d>& riv,
											   	   	  	  	  	    	   inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    	   std::ostream& explain) const;

	extern template class group_scaler<traits1d>;
	extern template class group_scaler<traits2d>;
	extern template class group_scaler<traits3d>;
	extern template class group_scaler<traits4d>;

	extern template class learner<traits1d>;
	extern template class learner<traits2d>;
	extern template class learner<traits3d>;
	extern template class learner<traits4d>;

}

#endif /* REEXP_EXP_H_ */
