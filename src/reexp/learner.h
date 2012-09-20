/*
 * learner.h
 *
 *  Created on: Apr 1, 2011
 *      Author: arauhala
 */

#ifndef LEARNER_H_
#define LEARNER_H_

#include "stats.h"

#include <queue>

namespace explib {

	template <typename P>
	struct candidate {
		const rel_stats<P>* rel_;
		int state_;
		double bias_;
		candidate(const rel_stats<P>& rel, int state, double bias)
		:	rel_(&rel), state_(state), bias_(bias) {}
	};

	template <typename P>
	bool operator <(const candidate<P>& f, const candidate<P>& s) {
		return f.bias_ < s.bias_;
	}

	template <typename P>
	class learner {
		private:
			lang<P>& lang_;
			stats<P>& stats_;
			double threshold_;
			double filterThreshold_;
			double predFilterThreshold_;
			bits exclude_;

		public:

			learner(lang<P>& l, stats<P>& s, double threshold = 0, double filterThreshold = -1, double predFilterTreshold = -1)
			: lang_(l), stats_(s), threshold_(threshold) {
				filterThreshold_ = filterThreshold < 0 ? threshold_ : filterThreshold;
				predFilterThreshold_ = predFilterTreshold;
			}

			void exclude(int var) {
				if (var >= exclude_.size()) {
					exclude_.resize(var+1);
				}
				exclude_[var] = true;
			}
			void exclude(const explib::bits& b) {
				if (b.size() >= exclude_.size()) {
					exclude_.resize(b.size());
				}
				exclude_ |= b;
			}
			void unexclude(int var) {
				if (var < exclude_.size()) {
					exclude_[var] = false;
				}
			}
			void unexclude(const explib::bits& b) {
				if (b.size() >= exclude_.size()) {
					exclude_.resize(b.size());
				}
				exclude_.andNeg(b);
			}

			bool var_excluded(int vid) const {
				return vid < exclude_.size() && exclude_[vid];
			}

			bool excluded(const rel<P>& r) const {
				int varCount = r.varCount();
				for (int j = 0; j < varCount; j++) {
					int vid = r.entries()[j].var_->id();
					if (var_excluded(vid)) {
						return true;
					}
				}
				return false;
			}

			void scan(std::priority_queue<candidate<P> >& c) {
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

			bool add_exp() {
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

			void disable_rels() {
				// optimization. updating relations' statistics is pretty much the
				// performance bottleneck for re-expression, so getting rid of
				// uninteresting relations brings major performance improvement.
				explib::lang<P>& l = lang_;
				for (int i = 0; i < l.rel_count(); i++) {
					if (!l.rel(i).disabled()) {
						const explib::rel_stats<P>& r( stats_.rel(i) );
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

			int reexpress(bool filter = true, int max = -1) {
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

	};

}


#endif /* LEARNER_H_ */
