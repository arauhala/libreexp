/*
 * stats.i.h
 *
 *  Created on: Nov 22, 2013
 *      Author: arau
 */

#ifndef REEXP_STATS_I_H_
#define REEXP_STATS_I_H_

#include "reexp/stats.h"
#include "reexp/data.i.h"

namespace reexp {

	template <typename P>
	void rel_stats<P>::row_update() {
		reset();
		const reexp::rel<P>& r = rel_.rel();
		if (r.disabled()) {
			return; // skip
		}

		reexp::cvec<P> dim = rel_.dim();
		reexp::bits defs(dim.volume());
		reexp::bits var_states(defs.size());

		// 1. first solve the defined vector
		//
		defs.fill(true);
		for (size_t v = 0; v < rel_.entries().size(); ++v) {
			rel_.setand_var_def_bits(v, defs);
		}
		n_ = defs.popcount();

		// 2. init the state bit vectors with the defined vector
		//
		reexp::bits rel_states[1<<P::MAX_REL_VARS];
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

	template <typename P>
	void rel_stats<P>::states_sparse_undefined(const reexp::data_var<P>& dv,
								 	 	   const std::vector<int>& dirtychunks,
								 	 	   const reexp::bits& undefs) {

		// NOTE: DOES NOT WORK
		// NOTE: SPEED BOTTLENECK IS ELSEWHERE (after delta optimization)
		const reexp::rel<P>& r = rel_.rel();
		if (r.disabled()) return; // skip
		const reexp::var<P>& v = dv.var();
		const reexp::cvec<P>& vdim = dv.dim();

		cvec<P> d( rel_.dim() );

		int vshift = 0;
		for (const rel_entry<P>& e : rel_.entries()) {
			if (e.var_ == &v) {
				vshift = vdim.offset(e.shift_);
			}
		}

		const reexp::bits* statebits[P::MAX_REL_VARS];
		const reexp::bits* defbits[P::MAX_REL_VARS];
		reexp::cvec<P> dims[P::MAX_REL_VARS];
		int shifts[P::MAX_REL_VARS];

		for (size_t i = 0; i < rel_.entries().size(); ++i) {
			const reexp::data_var<P>& dv = rel_.data_.var(rel_.entries()[i].var_->id());
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
	template <typename P>
	void rel_stats<P>::states_undefined(const reexp::data_var<P>& dv,
						  	  	  	const reexp::bits& undefs,
						  	  	  	bool version) {
		const reexp::rel<P>& r = rel_.rel();
		if (r.disabled()) return; // skip
		const reexp::var<P>& v = dv.var();
		const reexp::cvec<P>& vdim = dv.dim();

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
		const reexp::bits* statebits[P::MAX_REL_VARS];
		const reexp::bits* defbits[P::MAX_REL_VARS];
		reexp::cvec<P> dims[P::MAX_REL_VARS];
		int shifts[P::MAX_REL_VARS];

		for (size_t i = 0; i < rel_.entries().size(); ++i) {
			const reexp::data_var<P>& dv =
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

			reexp::bchunk_istream2<false_tail_fill> in(undefs.chunks().data(), rowbegin, i.length());
			int at = 0;
			while (in) {
				bchunk undefc;
				in>>undefc;

				if (undefc) {
					bchunk varchunks[P::MAX_REL_VARS];
					for (int i = 0; i < P::MAX_REL_VARS; ++i) varchunks[i] = 0; // make error go away

					size_t rsz = rel_.entries().size();
					for (size_t v = 0; v < rsz; ++v) {
						if (rel_.entries()[v].var_ != &dv.var()) {
							int chunkbegin = dims[v].offset(begin) + shifts[v] + at;
							undefc &= defbits[v]->chunk_from(chunkbegin);
						}
					}
					for (size_t v = 0; v < rsz; ++v) {
						int chunkbegin = dims[v].offset(begin) + shifts[v] + at;
						varchunks[v] = statebits[v]->chunk_from(chunkbegin) & undefc;
						varFreqs_[v] -= __builtin_popcountl(varchunks[v]);
					}
					for (size_t s = 0; s < size_t(stateCount()); ++s) {
						bchunk statechunk = undefc;
						for (size_t v = 0; v < rsz; ++v) {
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

	template <typename P>
	int rel_stats<P>::targetVersion(reexp::stats<P>& stats) const {
		int v = -1;
		typedef typename std::vector<rel_entry<P> >::const_iterator iter;
		for (iter i = rel_.entries().begin(); i != rel_.entries().end(); ++i) {
			const reexp::var_stats<P>& vs = stats.var(i->var_->id());
			if (vs.disabled()) return -1;
			v = std::max(v, vs.data().version_);
		}
		return v;
	}


	template <typename P>
	void rel_stats<P>::update(reexp::stats<P>& stats) {
		int v = targetVersion(stats);
		if (v > version_ && rel_.dim().volume()) {
			row_update();
			//bit_update();
			version_ = v;
		}
	}

	template <typename P>
	void rel_stats<P>::bit_update() {
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

	template <typename P>
	exp_rel_stats<P>::exp_rel_stats(const exp<P>* e)
	: exp_(e),
	  varFreqs_(),
	  stateFreqs_(),
	  n_() {
		if (exp_)
		{
			const rel<P>& r = exp_->rel();
			int vars = r.varCount();
			varFreqs_.resize(vars);
			stateFreqs_.resize(1<<vars);
		}
	}


	template <typename P>
	void exp_rel_stats<P>::update(const data<P>& d) {
		if (exp_)
		{
			const rel<P>& r = exp_->rel();
			data_rel<P> dr = d.rel(r.id());

			// FIXME: redundancy with rel::row_update()

			for (size_t i = 0; i < varFreqs_.size(); i++) {
				varFreqs_[i] = 0;
			}
			for (size_t i = 0; i < stateFreqs_.size(); i++) {
				stateFreqs_[i] = 0;
			}

			reexp::cvec<P> dim = dr.dim();
			int vol = dim.volume();
			reexp::bits defs(vol);
			reexp::bits var_states(defs.size());

			// 1. first copy the defined vector from cache
			//
			defs = d.exp_rel_defs(exp_->id());
			n_ = defs.popcount();

			// 2. init the state bit vectors with the defined vector
			//
			reexp::bits rel_states[1<<P::MAX_REL_VARS];
			for (int s = 0; s < r.stateCount(); ++s) {
				rel_states[s] = defs;
			}

			// 3. read state vectors for each variable and resolve state vectors
			//
			for (int v = 0; v < r.varCount(); ++v) {
				var_states.copy(defs); // init
				dr.setand_var_state_bits(v, var_states);
				varFreqs_[v] = var_states.popcount();
				for (int s = 0; s < r.stateCount(); ++s) {
					if (r.varState(v, s)) {
						rel_states[s] &= var_states;
					} else {
						rel_states[s].andNeg(var_states);
					}
				}
			}

			// 4. read statistics from state vectors
			//
			for (int s = 0; s < r.stateCount(); ++s) {
				stateFreqs_[s] = rel_states[s].popcount();
			}
		}
	}


	template <typename P>
	void filter_rels(reexp::lang<P>& l, const reexp::stats<P>& s, double threshold) {
		for (int i = 0; i < l.rel_count(); i++) {
			if (!l.rel(i).disabled()) {
				const reexp::rel_stats<P>& r( s.rel(i) );
				int s = r.stateCount();
				while (--s >= 0) {
					if (r.eStateBias(s) >= threshold) break;
				}
				if (s < 0) l.disable_rel(i);
			}
		}
	}

}


#endif /* STATS_I_H_ */
