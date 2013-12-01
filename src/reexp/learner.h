/*
 * learner.h
 *
 *  Created on: Apr 1, 2011
 *      Author: arauhala
 */

#ifndef REEXP_LEARNER_H_
#define REEXP_LEARNER_H_

#include "stats.h"

#include <queue>

namespace reexp {

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
			void exclude(const reexp::bits& b) {
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
			void unexclude(const reexp::bits& b) {
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
					if (var_excluded(vid)) return true;
				}
				return false;
			}

			void scan(std::priority_queue<candidate<P> >& c);

			bool add_exp();

			void disable_rels();

			int reexpress(bool filter = true, int max = -1);

	};

}


#endif /* REEXP_LEARNER_H_ */
