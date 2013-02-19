/*
 * knn.h
 *
 *  Created on: Jul 31, 2012
 *      Author: arau
 */

#ifndef REEXP_KNN_H_
#define REEXP_KNN_H_

#include "stats.h"

namespace explib {

	template <typename P>
	class knn {
		private:
			const stats<P>& stats_;
			explib::bits pred_; // predicted vars
			cvec<P> selection_; // selected dims
			std::vector<std::pair<double, int> > influence_;
			int k_;

		public:

			knn(const stats<P>& stats, const explib::bits& pred, const cvec<P>& selection, int k)
			: stats_(stats), pred_(pred), selection_(selection), influence_(), k_(k) {

				selection_ = 1;
				for (int i = 0; (select_>>i); ++i) {
					if ((select_>>i)&1) {
						selection_ *= stats.data().dim()[i];
					}
				}

				int vc = stats_.data().lang().var_count();
				for (int i = 0; i < vc; ++i) {
					if (i >= pred.size() || !pred[i]) {
						influence_.push_back(
							std::pair<double, int>(calc_influence(i, pred_), i));
					}
				}
				std::sort(influence_.begin(), influence_.end());
			}

			std::vector<std::pair<double, int>> find(const data<P>& d, const cvec<P>& at) {
				std::vector<double> distances;
				distances.resize(selection_.volume());
				// start from the end
				for (int i = influence_.size()-1; i >= 0; --i) {
				}
			}

			std::vector<double> p(const data<P>& d, const cvec<P>& at) {

			}

			double calc_influence(int from, const explib::bits& to) {
				const explib::stats<P>& s = stats_;
				const explib::lang<P>& l = s.data().lang();
				double rv = 0;

				for (int i = 0; i < l.rel_count(); ++i) {
					const explib::rel<P>& r = l.rel(i);
					bool inflFound = 0;
					bool efound = 0;
					for (size_t j = 0; j < r.entries().size(); ++j) {
						int id = r.entries()[j].var_->id();
						if (id == from) inflFound = true;
						if (id < to.size() && to[id]) efound = true;
					}
					if (inflFound && efound) {
						const explib::rel_stats<P>& rs = s.rel(i);
						rv += rs.n() * (rs.eNaiveEntropy() - rs.eEntropy());
					}
				}
				return rv;
			}

	};

}



#endif /* REEXP_KNN_H_ */
