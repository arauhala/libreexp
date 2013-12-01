/*
 * stats.h
 *
 *  Created on: Dec 15, 2010
 *      Author: arauhala
 */

#ifndef REEXP_STATS_H_
#define REEXP_STATS_H_

#include "data.h"
#include "info.h"

#include <stdexcept>

namespace reexp {

	template <typename P>
	class stats;

	template <typename P>
	class var_stats {
		public:
			var_stats(const data_var<P>& var)
			:	var_(var), disabled_(false), version_(-1) {
				const reexp::ctx<P>& c = var.var_.ctx();
				for (rowdim_= 0; c.v_[rowdim_] < 0; ++rowdim_) {}
			}
			int freq() const {
				return freq_;
			}
			int n() const {
				return n_;
			}
			bool disabled() const {
				return disabled_;
			}
			void set_disabled(bool disabled) {
				disabled_ = disabled;
			}
			void row_update() {
				n_ = var_.defined().popcount();
				bits states;;
				states = var_.states();
				states &= var_.defined();
				freq_ = states.popcount();
			}
			void states_undefined(const reexp::bits& undefs) {
				row_update(); // no optimization needed
			}
			void update() {
				if (var_.version_ > version_) {
					row_update();
					//bit_update();
					version_ = var_.version_;
				}
			}
			void bit_update() {
				freq_ = 0;
				n_ = 0;
				// could be heavily optimized by treating row by row
				cvec<P> dim( var_.dim() );
				typedef const_cond_bit_ref cbit;
				for (dim_iterator<P> i(dim); i; ++i) {
					const cvec<P>& at(*i);
					cbit bit = var_[at];
					if (bit) {
						if (*bit) freq_++;
						n_++;
					}
				}
			}
			double p() const {
				return reexp::p(freq_, n_);
			}
			double eP() const {
				return reexp::hatP(freq_, n_);
			}
			double eEntropy() const {
				return reexp::entropy(eP());
			}
			double entropy() const {
				return reexp::entropy(p());
			}
			double information() const {
				return entropy() * n_;
			}
			const data_var<P>& data() const {
				return var_;
			}
		private:
			// statistics
			const data_var<P>& var_;
			int freq_;
			int n_;

			// optimizations
			bool disabled_;
			int version_;
			int rowdim_;
	};


	/**
	 * Relation statistics for some specific data set
	 */
	template <typename P>
	struct rel_stats {
		private:
			data_rel<P> rel_;
			std::vector<int> varFreqs_;
			std::vector<int> stateFreqs_;
			int n_;
			// optimization related
			int rowdim_;
			std::bitset<P::MAX_REL_VARS> varrows_;
			int version_;

		public:
			rel_stats(const data_rel<P>& rel)
			:	rel_(rel),
			 	varFreqs_(),
			 	stateFreqs_(),
			 	n_(0),
			 	rowdim_(0),
			 	varrows_(),
			 	version_(-1) {
				int vars = rel.entries().size();
				varFreqs_.resize(vars);
				stateFreqs_.resize(1<<vars);
				const std::vector<rel_entry<P> >& e = rel.entries();
				const reexp::ctx<P>& c = rel.rel().ctx();
				for (; c.v_[rowdim_] < 0; ++rowdim_) {}
				varrows_.reset();
				for (size_t i = 0; i < e.size(); ++i) {
					if (e[i].var_->ctx().v_[rowdim_] >= 0) {
						varrows_[i] = true;
					}
				}
				if (!varrows_.count()) {
					//rowdim_ = -1;
					throw std::runtime_error("no variables have bit rows for this relation?");
				}
			}
			void reset() {
				for (size_t i = 0; i < varFreqs_.size(); i++) {
					varFreqs_[i] = 0;
				}
				for (size_t i = 0; i < stateFreqs_.size(); i++) {
					stateFreqs_[i] = 0;
				}
				n_ = 0;
			}

			const std::vector<int>& varFreqs() const {
				return varFreqs_;
			}
			const std::vector<int>& stateFreqs() const {
				return stateFreqs_;
			}
			int n() const {
				return n_;
			}
			const data_rel<P>& data() const {
				return rel_;
			}
			int stateCount() const {
				return rel_.rel_.stateCount();
			}
			int varCount() const {
				return rel_.rel_.varCount();
			}
			inline double varPriori(int var) const {
				return 0.5;
			}
			inline double prioriWeight() const {
				return 1;
			}
			double varP(int var) const {
				return varFreqs_[var] / (double)n_;
			}
			double eVarP(int var) const {
				return (varFreqs_[var]+varPriori(var)) / (double)(n_ + prioriWeight());
			}
			double stateP(int state) const {
				return stateFreqs_[state] / (double)n_;
			}
			double eStateP(int state) const {
				return (stateFreqs_[state]+rel_.rel().statePrioriP(state)*prioriWeight()) / (double)(n_+prioriWeight());
			}
			double stateNaiveP(int state) const {
				double rv = 1;
				for (size_t i = 0; i < varFreqs_.size(); i++) {
					bool s = state & (1 << i);
					rv *= s ? varP(i) : 1. - varP(i);
				}
				return rv;
			}
			double eStateNaiveP(int state) const {
				double rv = 1;
				for (size_t i = 0; i < varFreqs_.size(); i++) {
					bool s = state & (1 << i);
					rv *= s ? eVarP(i) : 1. - eVarP(i);
				}
				return rv;
			}
			double stateInfo(int state) const {
				double rv = stateP(state);
				if (rv > 0) rv = log(rv) / log(2);
				return -rv;
			}
			double eStateInfo(int state) const {
				double rv = eStateP(state);
				if (rv > 0) rv = log(rv) / log(2);
				return -rv;
			}
			double stateEntropy(int state) const {
				return n_?stateP(state) * stateInfo(state):0;
			}
			double eStateEntropy(int state) const {
				return eStateP(state) * eStateInfo(state);
			}
			double entropy() const {
				double rv = 0;
				for (int i = 0; i < stateCount(); ++i) rv += stateEntropy(i);
				return rv;
			}
			double eEntropy() const {
				double rv = 0;
				for (int i = 0; i < stateCount(); ++i) rv += eStateEntropy(i);
				return rv;
			}
			double stateNaiveInfo(int state) const {
				double rv = stateNaiveP(state);
				if (rv > 0) rv = log(rv) / log(2);
				return -rv;
			}
			double eStateNaiveInfo(int state) const {
				double rv = eStateNaiveP(state);
				if (rv > 0) rv = log(rv) / log(2);
				return -rv;
			}
			double stateNaiveEntropy(int state) const {
				return stateNaiveP(state) * stateNaiveInfo(state);
			}
			double eStateNaiveEntropy(int state) const {
				return eStateNaiveP(state) * eStateNaiveInfo(state);
			}
			double naiveEntropy() const {
				double rv = 0;
				for (int i = 0; i < stateCount(); ++i) rv += stateNaiveEntropy(i);
				return rv;
			}
			double eNaiveEntropy() const {
				double rv = 0;
				for (int i = 0; i < stateCount(); ++i) rv += eStateNaiveEntropy(i);
				return rv;
			}
			double stateBias(int state) const {
				return stateFreqs_[state] * (stateNaiveInfo(state) - stateInfo(state));
			}
			double eStateBias(int state) const {
				return stateFreqs_[state] * (eStateNaiveInfo(state) - eStateInfo(state));
			}

			void row_update();


			void states_sparse_undefined(const reexp::data_var<P>& dv,
								   	     const std::vector<int>& dirtychunks,
								   	     const reexp::bits& undefs);

			/**
			 * This doesn't work, if the undefined variable is present
			 * multiple times in the relation.
			 *
			 * There are basically few ways to avoid this.
			 *    One way would be in preparing separate delta row
			 *    in cases, where there are 2 variables.
			 */
			void states_undefined(const reexp::data_var<P>& dv,
								  const reexp::bits& undefs,
								  bool version);

			int targetVersion(reexp::stats<P>& stats) const;

			void update(reexp::stats<P>& stats);

			void bit_update();
	};

	/**
	 * exp_rel_stats provides statistics for the version of a relation
	 * that the expression was based on.
	 *
	 * FIXME: Even if someone changes the data (by hand) this entry
	 *        will not be updated.
	 */
	template <typename P>
	class exp_rel_stats {
		private:
			const exp<P>* exp_;
			std::vector<int> varFreqs_;
			std::vector<int> stateFreqs_;
			int n_;

		public:
			inline const std::vector<int>& varFreqs() const {
				return varFreqs_;
			}
			inline const std::vector<int>& stateFreqs() const {
				return stateFreqs_;
			}
			inline int n() const {
				return n_;
			}
		public:
			exp_rel_stats(const exp<P>* e);

			void update(const data<P>& d);
	};

	template <typename P>
	class stats : public data_obs<P> {
		private:
			const reexp::data<P>& data_;
			util::dense_vector<var_stats<P> > vars_;
			util::dense_vector<rel_stats<P> > rels_;
			util::dense_vector<exp_rel_stats<P> > exp_rels_;
		public:
			stats(reexp::data<P>& data)
			: data_(data), vars_(), rels_(), exp_rels_() {
				data.set_obs(*this);
				const lang<P>& lang = data_.lang();
				for (int i = 0; i < lang.var_count(); i++) {
					var_added(data_.var(i));
				}
				for (int i = 0; i < lang.rel_count(); i++) {
					rel_added(data.rel(i));
				}
			}
			void set_disabled(int varid, bool state) {
				vars_[varid].set_disabled(state);
			}
			void rel_added(const data_rel<P>& rel) {
				rels_.push_back(rel_stats<P>(rel));
				rels_.back().update(*this);
			}
			void var_added(const data_var<P>& var) {
				vars_.push_back(var_stats<P>(var));
				const exp<P>* e = dynamic_cast<const exp<P>*>(&var.var_);
				exp_rels_.push_back(exp_rel_stats<P>(e));
				exp_rels_.back().update(data_);
				update();
			}
			void states_undefined(const data_var<P>& var,
								  const reexp::bits& undefs,
								  int version) {
				const reexp::var<P>& v = var.var();

/*				std::vector<int> dirtychunks;
				for (size_t i = 0; i < undefs.chunks().size(); ++i) {
					if (undefs.chunks()[i]) dirtychunks.push_back(i);
				}*/

				for (const reexp::rel<P>* r : v.rels()) {
//					rels_[r->id()].states_sparse_undefined(var, dirtychunks, undefs);
					rels_[r->id()].states_undefined(var, undefs, version);
				}
			}
			void update() {
				for (int i = 0; i < vars_.size(); i++) {
					vars_[i].update();
				}
				for (int i = 0; i < rels_.size(); i++) {
					rels_[i].update(*this);
				}
			}
			double naiveInfo() const {
				double rv = 0;
				for (int i = 0; i < vars_.size(); i++) {
					rv += vars_[i].information();
				}
				return rv;
			}
			const reexp::data<P>& data() const {
				return data_;
			}
			const var_stats<P>& var(int i) const {
				return vars_[i];
			}
			const rel_stats<P>& rel(int i) const {
				return rels_[i];
			}
			const exp_rel_stats<P>& exp_rel(int i) const {
				return exp_rels_[i];
			}
	};


/*
	template <typename P>
	const_cond_bit_ref var<P>::at(const data<P>& data, const cvec<P>& at) const {
		return data.bits()[index(data, at)];
	}

	template <typename P>
	cond_bit_ref var<P>::at(data<P>& data, const cvec<P>& at) const {
		return data.bits()[index(data, at)];
	}*/


}

#endif /* REEXP_STATS_H_ */
