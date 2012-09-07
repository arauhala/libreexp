/*
 * stats.h
 *
 *  Created on: Dec 15, 2010
 *      Author: arauhala
 */

#ifndef STATS_H_
#define STATS_H_

#include "data.h"
#include "info.h"

#include <stdexcept>

#define USE_ROW_UPDATE

namespace explib {

	template <typename P>
	class stats;

	template <typename P>
	class var_stats {
		public:
			var_stats(const data_var<P>& var)
			:	var_(var), disabled_(false), version_(-1) {
				const explib::ctx<P>& c = var.var_.ctx();
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
			void states_undefined(const explib::bits& undefs) {
				row_update(); // no optimization needed
			}
			void update() {
				if (var_.version_ > version_) {
#ifdef USE_ROW_UPDATE
					row_update();
#else
					bit_update();
#endif
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
				return explib::p(freq_, n_);
			}
			double eP() const {
				return explib::hatP(freq_, n_);
			}
			double eEntropy() const {
				return explib::entropy(eP());
			}
			double entropy() const {
				return explib::entropy(p());
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
	 * Relation statistics for some specif{ic data set
	 */
	template <typename P>
	struct rel_stats {
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
				const explib::ctx<P>& c = rel.rel().ctx();
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
/*			inline double statePriori(int state) const {
				for (size_t i = 0; i < varFreqs_.size(); ++i) {

				}
				return 0.25;
			}*/
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

			void row_update2() {
				reset();
				const explib::rel<P>& r = rel_.rel();
				if (r.disabled()) {
					return; // skip
				}

				explib::cvec<P> dim = rel_.dim();
				explib::bits defs(dim.volume());
				explib::bits var_states(defs.size());

				// 1. first solve the defined vector
				//
				defs.fill(true);
				for (size_t v = 0; v < rel_.entries().size(); ++v) {
					rel_.setand_var_def_bits(v, defs);
				}
				n_ = defs.popcount();

				// 2. init the state bit vectors with the defined vector
				//
				explib::bits rel_states[1<<P::MAX_REL_VARS];
				for (int s = 0; s < stateCount(); ++s) {
					rel_states[s] = defs;
				}

				// 3. read state vectors for each variable and resolve state vectors
				//
				for (int v = 0; v < varCount(); ++v) {
					var_states.copy(defs); // init
					rel_.setand_var_state_bits(v, var_states);
					varFreqs_[v] = var_states.popcount();
					for (int s = 0; s < stateCount(); ++s) {
						if (r.varState(v, s)) {
							rel_states[s] &= var_states;
						} else {
							rel_states[s].andNeg(var_states);
						}
					}
				}

				// 4. read statistics from state vectors
				//
				for (int s = 0; s < stateCount(); ++s) {
					stateFreqs_[s] = rel_states[s].popcount();
				}
			}

			void row_update() {

				// few thoughts:
				//  * row is only exist to the lowest
				//    context variable that is defined for variable
				//  * sometimes the 'row' dimension is not
				//    for the first context variable
				//  * sometimes the row dimension for all variables
				//    is not the same.
				//  * Is it even possible that there is no row?

				reset();
				const explib::rel<P>& r = rel_.rel();
				if (r.disabled()) {
					return; // skip
				}
				typedef cond_bit_ref cbit;
				cvec<P> d( rel_.dim() );

				dim_row_iterator<P> i( d, rowdim_ );
				int vars = rel_.entries().size();
				int states = stateCount();
				int rowl = i.length();
				bits defined(rowl);
				bits statef[1<<P::MAX_REL_VARS];
				bits rows[P::MAX_REL_VARS]; // support at max 8 vars in rel
				for (int j = 0; j < states; ++j) {
					statef[j].resize(rowl);
				}
				for (int j = 0; j < P::MAX_REL_VARS; ++j) {
					rows[j].resize(rowl);
				}
				for (; i; ++i) {
					defined.fill(true);
					// make copy of the row and figure out
					// defined bits
					cvec<P> begin( i.begin() );
					for (int j = 0; j < vars; ++j) {
						if (varrows_[j]) {
							cond_bits_ref row = rel_.var(j).bitrow(begin, rowl);
							rows[j] = row;
							defined &= rows[j];
							rows[j] = *row;
						} else {
							cond_bit_ref bit = rel_.var(j)[begin];
							if (!bit) defined.fill(false);
							rows[j].fill(*bit);
						}
					}

					for (int j = 0; j < states; ++j) {
						statef[j] = defined;
					}
					for (int j = 0; j < vars; ++j) {
						bits& row = rows[j];
						row &= defined;
						for (int k = 0; k < states; ++k) {
							if (r.varState(j, k)) {
								statef[k] &= row;
							} else {
								statef[k].andNeg(row); // &= ~
							}
						}
					}
					n_ += defined.popcount();
					for (int k = 0; k < vars; ++k) {
						varFreqs_[k] += rows[k].popcount();
					}
					for (int k = 0; k < states; ++k) {
						//std::cout<<"state "<<k<<" "<<vector_todensestring(statef[k])<<"\n";
						stateFreqs_[k] += statef[k].popcount();
					}
				}
			}

			void states_sparse_undefined(const explib::data_var<P>& dv,
								   	     const std::vector<int>& dirtychunks,
								   	     const explib::bits& undefs) {

				// NOTE: DOES NOT WORK
				// NOTE: SPEED BOTTLENECK IS ELSEWHERE (after delta optimization)
				const explib::rel<P>& r = rel_.rel();
				if (r.disabled()) return; // skip
				const explib::var<P>& v = dv.var();
				const explib::cvec<P>& vdim = dv.dim();

				typedef cond_bit_ref cbit;
				cvec<P> d( rel_.dim() );

				int vshift = 0;
				for (const rel_entry<P>& e : rel_.entries()) {
					if (e.var_ == &v) {
						vshift = vdim.offset(e.shift_);
					}
				}

				const explib::bits* statebits[P::MAX_REL_VARS];
				const explib::bits* defbits[P::MAX_REL_VARS];
				explib::cvec<P> dims[P::MAX_REL_VARS];
				int shifts[P::MAX_REL_VARS];

				for (size_t i = 0; i < rel_.entries().size(); ++i) {
					const explib::data_var<P>& dv = rel_.data_.var(rel_.entries()[i].var_->id());
					statebits[i] = &dv.states();
					defbits[i] = &dv.defined();
					dims[i] = dv.dim();
					shifts[i] = dims[i].offset(rel_.entries()[i].shift_);
				}

				for (int cidx : dirtychunks) {
					bchunk undefc = undefs[cidx];
					cvec<P> begin( d.at(cidx * bchunk_bsize - vshift) );

					// TODO: handle rows and row boundaries. begin may be even
					// outside the row. If row is longer than the chunk, chunk
					// needs to be split.
					//
					bchunk varchunks[P::MAX_REL_VARS];
					for (size_t v = 0; v < rel_.entries().size(); ++v) {
						if (rel_.entries()[v].var_ != &dv.var()) {
							int chunkbegin = dims[v].offset(begin) + shifts[v];
							undefc &= defbits[v]->chunk_from(chunkbegin);
						}
					}
					for (size_t v = 0; v < rel_.entries().size(); ++v) {
						int chunkbegin = dims[v].offset(begin) + shifts[v];
						varchunks[v] = statebits[v]->chunk_from(chunkbegin) & undefc;
						varFreqs_[v] -= __builtin_popcountl(varchunks[v]);
					}
					for (size_t s = 0; s < size_t(stateCount()); ++s) {
						bchunk statechunk = undefc;
						for (size_t v = 0; v < rel_.entries().size(); ++v) {
							if (r.varState(v, s)) {
								statechunk &= varchunks[v];
							} else {
								statechunk &= ~varchunks[v];
							}
						}
						stateFreqs_[s] -= __builtin_popcountl(statechunk);
					}
					n_ -= __builtin_popcountl(undefc);
				}
			}

			/**
			 * This doesn't work, if the undefined variable is present
			 * multiple times in the relation.
			 *
			 * There are basically few ways to avoid this.
			 *    One way would be in preparing separate delta row
			 *    in cases, where there are 2 variables.
			 */
			void states_undefined(const explib::data_var<P>& dv,
								  const explib::bits& undefs,
								  bool version) {
				const explib::rel<P>& r = rel_.rel();
				if (r.disabled()) return; // skip
				const explib::var<P>& v = dv.var();
				const explib::cvec<P>& vdim = dv.dim();

				typedef cond_bit_ref cbit;
				cvec<P> d( rel_.dim() );

				int vshift = 0;
				int changedvars = 0;
				for (const rel_entry<P>& e : rel_.entries()) {
					if (e.var_ == &v) {
						vshift = vdim.offset(e.shift_);
						changedvars++;
					}
				}
				if (changedvars > 1) {
					return;
				}
				const explib::bits* statebits[P::MAX_REL_VARS];
				const explib::bits* defbits[P::MAX_REL_VARS];
				explib::cvec<P> dims[P::MAX_REL_VARS];
				int shifts[P::MAX_REL_VARS];

				for (size_t i = 0; i < rel_.entries().size(); ++i) {
					const explib::data_var<P>& dv =
						rel_.data_.var(rel_.entries()[i].var_->id());
					statebits[i] = &dv.states();
					defbits[i] = &dv.defined();
					dims[i] = dv.dim();
					shifts[i] = dims[i].offset(rel_.entries()[i].shift_);
				}

				// let's update the dirty rows, and the dirty rows only
				//
				for (dim_row_iterator<P> i( d, rowdim_ ); i; ) {
					// first of all, how to find dirty rows?
					//

					// make copy of the row and figure out defined bits
					//
					cvec<P> begin( i.begin() );
					int rowbegin = vdim.offset(begin) + vshift;

					explib::bchunk_istream2<false_tail_fill> in(undefs.chunks().data(), rowbegin, i.length());
					int at = 0;
					while (in) {
						bchunk undefc;
						in>>undefc;

						if (undefc) {
							bchunk varchunks[P::MAX_REL_VARS];
							for (size_t v = 0; v < rel_.entries().size(); ++v) {
								if (rel_.entries()[v].var_ != &dv.var()) {
									int chunkbegin = dims[v].offset(begin) + shifts[v] + at;
									undefc &= defbits[v]->chunk_from(chunkbegin);
								}
							}
							for (size_t v = 0; v < rel_.entries().size(); ++v) {
								int chunkbegin = dims[v].offset(begin) + shifts[v] + at;
								varchunks[v] = statebits[v]->chunk_from(chunkbegin) & undefc;
								varFreqs_[v] -= __builtin_popcountl(varchunks[v]);
							}
							for (size_t s = 0; s < size_t(stateCount()); ++s) {
								bchunk statechunk = undefc;
								for (size_t v = 0; v < rel_.entries().size(); ++v) {
									if (r.varState(v, s)) {
										statechunk &= varchunks[v];
									} else {
										statechunk &= ~varchunks[v];
									}
								}
								stateFreqs_[s] -= __builtin_popcountl(statechunk);
							}
							n_ -= __builtin_popcountl(undefc);
						}
						at += bchunk_bsize;
					}
					++i;
				}
				version_ = version;
			}

			int targetVersion(explib::stats<P>& stats) const;

			void update(explib::stats<P>& stats) {
				int v = targetVersion(stats);
				if (v > version_ && rel_.dim().volume()) {
#ifdef USE_ROW_UPDATE
					row_update2();
					//row_update();
#else
					bit_update();
#endif
					version_ = v;
				}
			}

			void bit_update() {
				reset();
				if (rel_.rel().disabled()) {
					return; // skip
				}
				typedef cond_bit_ref cbit;
				cvec<P> d( rel_.dim() );

				dim_iterator<P> i( d );

				int vars = rel_.entries().size();

				for (; i; ++i) {
					cvec<P> at( *i );
					bool defined = true;
					for (int j = 0; j < vars; j++) {
						cbit bit = rel_.var(j)[at];
						defined &= bit;
					}
					if (defined) {
						int state = 0;
						for (int j = 0; j < vars; j++) {
							cbit bit = rel_.var(j)[at];
							if (*bit) {
								varFreqs_[j]++;
								state |= 1<<j;
							}
						}
						stateFreqs_[state]++;
						n_++;
					}
				}

			}
		private:
			data_rel<P> rel_;
			std::vector<int> varFreqs_;
			std::vector<int> stateFreqs_;
			int n_;
			// optimization related
			int rowdim_;
			std::bitset<P::MAX_REL_VARS> varrows_;
			int version_;
	};

	template <typename P>
	class stats : public data_obs<P> {
		public:
			stats(explib::data<P>& data)
			: data_(data) {
				data.set_obs(*this);
				const lang<P>& lang = data_.lang();
				for (int i = 0; i < lang.var_count(); i++) {
					vars_.push_back(var_stats<P>(data_.var(i)));
				}
				for (int i = 0; i < lang.rel_count(); i++) {
					rels_.push_back(rel_stats<P>(data.rel(i)));
				}
				update();
			}
			void set_disabled(int varid, bool state) {
				vars_[varid].set_disabled(state);
			}
			void rel_added(const data_rel<P>& rel) {
				rels_.push_back(rel_stats<P>(rel));
				rels_.back().update(*this); // NOTE: the real bottleneck
			}
			void var_added(const data_var<P>& var) {
				vars_.push_back(var_stats<P>(var));
				update();
			}
			void states_undefined(const data_var<P>& var,
								  const explib::bits& undefs,
								  int version) {
				const explib::var<P>& v = var.var();

/*				std::vector<int> dirtychunks;
				for (size_t i = 0; i < undefs.chunks().size(); ++i) {
					if (undefs.chunks()[i]) dirtychunks.push_back(i);
				}*/

				for (const explib::rel<P>* r : v.rels()) {
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
			const explib::data<P>& data() const {
				return data_;
			}
			const var_stats<P>& var(int i) const {
				return vars_[i];
			}
			const rel_stats<P>& rel(int i) const {
				return rels_[i];
			}
		private:
			const explib::data<P>& data_;
			util::dense_vector<var_stats<P> > vars_;
			util::dense_vector<rel_stats<P> > rels_;
	};

	template <typename P>
	int rel_stats<P>::targetVersion(explib::stats<P>& stats) const {
		int v = -1;
		typedef typename std::vector<rel_entry<P> >::const_iterator iter;
		for (iter i = rel_.entries().begin(); i != rel_.entries().end(); ++i) {
			const explib::var_stats<P>& vs = stats.var(i->var_->id());
			if (vs.disabled()) return -1;
			v = std::max(v, vs.data().version_);
		}
		return v;
	}

	template <typename P>
	void filter_rels(explib::lang<P>& l, const explib::stats<P>& s, double threshold) {
		for (int i = 0; i < l.rel_count(); i++) {
			if (!l.rel(i).disabled()) {
				const explib::rel_stats<P>& r( s.rel(i) );
				int s = r.stateCount();
				while (--s >= 0) {
					if (r.eStateBias(s) >= threshold) break;
				}
				if (s < 0) l.disable_rel(i);
			}
		}
	}

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

#endif /* STATS_H_ */
