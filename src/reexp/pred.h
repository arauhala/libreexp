/*
 * pred.h
 *
 *  Created on: Aug 16, 2011
 *      Author: arau
 */

#ifndef REEXP_PRED_H_
#define REEXP_PRED_H_


#include "stats.h"

namespace reexp {

	inline int pow3i(int v) {
		int rv = 1;
		while (v-- > 0) rv *= 3;
		return rv;
	}
	template <typename P>
	struct rel_inputvars {
		const reexp::rel<P>& r_;
		int predicted_;
		rel_inputvars(const reexp::rel<P>& r,
					  int predicted)
		: 	r_(r), predicted_(predicted) {}
		inline int stateCount() const {
			return pow3i(r_.varCount()-1);
		}
		enum DefState {
			FALSE = 0,
			TRUE = 1,
			UNDEFINED = 2
		};
		/**
		 * Transforms relation state into rel_inputvars state.
		 * dmask contains the definition bits, while smask contains
		 * the state bits.
		 */
		int maskState(int dmask, int smask) const;

		inline int varDefState(int var, int state) const {
			if (var > predicted_) var--;
			return (state / pow3i(var)) % 3;
		}
		bool includesRelState(int istate, int rstate) const;
	};

	struct inputstate_logdep {
		inline inputstate_logdep(int max_rels)
		: d_(2*size_t(pow3i(max_rels-1))) {}
		inline double& depStateV(int state) {
			return d_[state*2];
		}
		inline double& depStateNotV(int state) {
			return d_[(state*2)+1];
		}
		std::vector<double> d_;
	};

	struct fake_output {
		template <typename T>
		inline fake_output& operator<<(T) {
			return *this;
		}
	};

	template <typename P>
	class pred {
		private:
			const reexp::stats<P>& stats_;
			double prioriWeight_;
			bool useLogDepB_;

		public:
			const reexp::stats<P>& stats() const { return stats_; }

			pred(const reexp::stats<P>& s, double prioriWeight = 2., bool useLogDepB = false)
			: stats_(s), prioriWeight_(prioriWeight), useLogDepB_(useLogDepB) {}


#if 0
			/**
			 * Probabilistic re-expression function. Input data
			 * must be in original language format.
			 */
			std::vector<double> reexpress(const reexp::data<P>& d);
#endif


			template <typename Out>
			void getInputStateLogDepA(const rel_stats<P>& rs,
									  const rel_inputvars<P>& riv,
									  inputstate_logdep& inputStateLogDep,
									  Out& explain) const;

			template <typename Out>
			void getInputStateLogDepB(const rel_stats<P>& rs,
									  const rel_inputvars<P>& riv,
									  inputstate_logdep& inputStateLogDep,
									  Out& explain) const;


			template <typename Out>
			void getInputStateLogDepB(int n,
									  const std::vector<int>& varFreqs,
									  const std::vector<int>& stateFreqs,
									  const rel_inputvars<P>& riv,
									  inputstate_logdep& inputStateLogDep,
									  Out& explain) const;
			template <typename Out>
			void getInputStateLogDep(const reexp::rel_stats<P>& rs,
									 const rel_inputvars<P>& riv,
									 inputstate_logdep& inputStateLogDep,
									 Out& explain) const;


			template <typename Out>
			std::vector<double> bitP(const data<P>& d,
									 int varid,
									 Out& explain) const;

			/**
			 * Bit by bit predicting
			 */
			inline std::vector<double> bitP(const data<P>& d, int varid) const {
				fake_output out;
				return bitP(d, varid, out);
			}


			/**
			 * Row predicting. This method of predicting is especially optimized
			 * for situation, where relation has many-to-single bits mapping.
			 */
			template <typename Out = fake_output>
			std::vector<double> rowP(const data<P>& d, int varid, Out explain = Out()) const;

			std::vector<double> p(const data<P>& d, int varid) const;

			void info(const std::vector<double>& ps,
					  const data<P>& d,
					  int varid,
					  double& totalinfo,
					  double& entryinfo) const;

			void info(const data<P>&d,
					  int varid,
					  double& totalinfo,
					  double& entryinfo) const {
				info(p(d, varid), d, varid, totalinfo, entryinfo);
			}

	};

	template <typename P>
	class group_scaler {
		private:
			std::vector<double> scales_;
			std::vector<double> avers_;

		public:
			group_scaler(const pred<P>& pr, int var, int target_group_sz);

			void scale(std::vector<double>& ps);
	};

}

#endif /* REEXP_PRED_H_ */
