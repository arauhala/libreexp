/*
 * learner.i.h
 *
 *  Created on: Nov 23, 2013
 *      Author: arau
 */

#ifndef REEXP_LEARNER_I_H_
#define REEXP_LEARNER_I_H_

#include "learner.h"
#include "stats.i.h"

namespace reexp {

	template <typename P>
	void learner<P>::scan(std::priority_queue<candidate<P> >& c) {
		const data<P>& d( stats_.data() );
		const lang<P>& l( d.lang() );
		for (int i = 0; i < l.rel_count(); i++) {
			const rel<P>& r( l.rel(i) );
			if (!r.disabled() && !excluded(r)) {
				const rel_stats<P>& r( stats_.rel(i) );
				int stateCount = r.stateCount();
				for (int j = 0; j < stateCount; j++) {
	#if 1
					double b = r.eStateBias(j);
	#else
					double b = std::abs(r.eStateBias(j));
	#endif
					if (b > threshold_) {
						candidate<P> cand(r, j, b);
						c.push(cand);
					}
				}
			}
		}
	}

	template <typename P>
	bool learner<P>::add_exp() {
		std::priority_queue<candidate<P> > cands;
		scan(cands);
		if (!cands.empty()) {
			const candidate<P>& c( cands.top() );
			const rel_stats<P>& rs = *c.rel_;
			lang_.add_exp(rs.data().rel_, c.state_);
			return true;
		}
		return false;
	}

	template <typename P>
	void learner<P>::disable_rels() {
		// optimization. updating relations' statistics is pretty much the
		// performance bottleneck for re-expression, so getting rid of
		// uninteresting relations brings major performance improvement.
		reexp::lang<P>& l = lang_;
		for (int i = 0; i < l.rel_count(); i++) {
			if (!l.rel(i).disabled()) {
				const reexp::rel_stats<P>& r( stats_.rel(i) );
				int s = r.stateCount();
				double ft = excluded(l.rel(i))?predFilterThreshold_:filterThreshold_;
				if (ft > 0) {
					while (--s >= 0) {
						if (r.eStateBias(s) >= ft) break;
					}
					if (s < 0) l.disable_rel(i);
				}
			}
		}
	}

	template <typename P>
	int learner<P>::reexpress(bool filter, int max) {
		int i = 0;
		// An optimization. We typically don't want to include
		// predicted variable in re-expression so we can skip
		// updating these relations' statistics (which is expensive)
		// disabling variable disables also relationship.
		for (int i = 0; i < exclude_.size(); i++) {
			if (var_excluded(i)) stats_.set_disabled(i, true);
		}
		if (filter) disable_rels();
		while (i != max && add_exp()) {
			++i;
			if (filter) disable_rels();
		}
		// revive disabled variables
		for (int i = 0; i < exclude_.size(); i++) {
			if (var_excluded(i)) stats_.set_disabled(i, false);
		}
		stats_.update(); // update disabled relationships
		return i;
	}


}

#endif /* REEXP_LEARNER_I_H_ */
