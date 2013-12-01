#ifndef REEXP_DATA_I_H_
#define REEXP_DATA_I_H_

#include "reexp/data.h"

namespace reexp {


	template <typename P>
	void data_var<P>::set_implicit(const cvec<P>& at, int cause) {
		int i = index(at);

		if (i < cause) {
			cond_bits& bits = bits_;
			cbit b = bits[i];
			b = false;
			if (*b) {
				const reexp::exp<P>* e = dynamic_cast<const reexp::exp<P>*>(&var_);
				if (e) {
					// in case state was determined to be true, also the
					// expressed variables' states becomes implictly known
					const std::vector<rel_entry<P> >& re = e->rel().entries();
					for (size_t j = 0; j < re.size(); j++) {
						data_.var(re[j].var_->id()).set_implicit(at + re[j].shift_, cause);
					}
				}
			}
			for (auto i = var_.deps().begin();
				 i != var_.deps().end();
				 ++i) {
				if (i->var_->id() < data_.var_count()) {
					// when re-expressing new data, all data vars may not be added yet
					reexp::data_var<P>& v(data_.var(i->var_->id()));
					reexp::cvec<P> limits = v.var().ctx().dim(data_.dim());
					bool skip = false;
					for (int j = 0; j < P::DIM; ++j) {
						int l = at[j] + i->shift_[j];
						if (l < 0 || l >= limits[j]) {
							skip = true;
							break;
						}
					}
					if (!skip) {
						cvec<P> a(at + i->shift_);
						int j = v.index(a);
						if (j < cause) {
							bits[j] = false;
						}
					}
				}
			}
		}
	}

	template <typename P>
	void data_var<P>::touch_deps(int version) {
		if (version_ != version) { // prevents recursion
			touch(version);
			const reexp::exp<P>* e = dynamic_cast<const reexp::exp<P>*>(&var_);
			if (e) {
				// in case state was determined to be true, also the
				// expressed variables' states becomes implictly known
				const std::vector<rel_entry<P> >& re = e->rel().entries();
				for (size_t j = 0; j < re.size(); j++) {
					size_t id = re[j].var_->id();
					if (id < data_.var_count()) {
						data_.var(id).touch_deps(version);
					}
				}
			}
			for (auto i = var_.deps().begin();
				 i != var_.deps().end();
				 ++i) {
				size_t id = i->var_->id();
				if (id < data_.var_count()) {
					data_.var(id).touch(version);
				}
			}
		}
	}

	template <typename P>
	void data_rel<P>::setand_var_def_bits(size_t var, reexp::bits& defs) const {
		const reexp::cvec<P> d = dim();
		int rowdim = d.rowdim();
		assert(defs.size() >= d.volume());

		const rel_entry<P>& e = rel_.entries()[var];
		const data_var<P>& dv = data_.var(e.var_->id());

		int rowlength = d[rowdim];
		int destidx = 0;

		reexp::dim_row_iterator<P> i(d);
		const reexp::var<P>& v = dv.var_;

		if (v.ctx().v_[rowdim] < 0) {
			while (i) {
				cvec<P> at = i.begin() + e.shift_;
				if (!dv[at]) {
					defs.from(destidx, rowlength).fill(false);
				}
				destidx += rowlength;
				++i;
			}
		} else {
			const reexp::cvec<P> srcdim(dv.dim());
			const bits& srcdefs = dv.defined();
			int srcshift = srcdim.offset(e.shift_);
			while (i) {
				int srcidx = srcdim.offset(i.begin())+srcshift;
				defs.from(destidx, rowlength) &= srcdefs.from(srcidx, rowlength);
				destidx += rowlength;
				++i;
			}
		}
	}

	template <typename P>
	void data_rel<P>::setand_var_state_bits(int var, reexp::bits& states) const {
		const reexp::cvec<P> d = dim();
		int rowdim = d.rowdim();
		assert(states.size() >= d.volume());

		const rel_entry<P>& e = rel_.entries()[var];
		const data_var<P>& dv = data_.var(e.var_->id());

		int rowlength = d[rowdim];
		int destidx = 0;

		reexp::dim_row_iterator<P> i(d);
		const reexp::var<P>& v = dv.var_;

		if (v.ctx().v_[rowdim] < 0) {
			while (i) {
				cvec<P> at = i.begin() + e.shift_;
				if (!*dv[at]) {
					states.from(destidx, rowlength).fill(false);
				}
				destidx += rowlength;
				++i;
			}
		} else {
			const reexp::cvec<P> srcdim(dv.dim());
			const bits& srcstates = dv.states();
			int srcshift = srcdim.offset(e.shift_);
			while (i) {
				int srcidx = srcdim.offset(i.begin())+srcshift;
				states.from(destidx, rowlength) &= srcstates.from(srcidx, rowlength);
				destidx += rowlength;
				++i;
			}
		}

	}

	template <typename P>
	data<P>::data(reexp::lang<P>& lang,
					   const reexp::cvec<P>& dim)
	: lang_(lang), dim_(dim), obs_() {
		lang.set_obs(*this);
		for (int i = 0; i < lang.orig_count(); i++) {
			var_added(lang.orig(i));
		}
	}

	template <typename P>
	void data<P>::apply_exps() {
		for (int i = vars_.size(); i < lang_.var_count(); i++) {
			const reexp::exp<P>* e(
				dynamic_cast<const reexp::exp<P>*>(&lang_.var(i)));
			if (e) var_added(*e);
		}
	}

	// needed by serialization
	template <typename P>
	void data<P>::apply_state(const reexp::exp<P>& exp, cvec<P> at, bool state) {
		data_var<P>& dvar(var(exp.id()));
		*dvar[at] = state;
		if (state) {
			for (const reexp::implmask<P>& mask : exp.expmasks()) {
				data_var<P>& mvar = var(mask.var_->id());
				apply_mask(at,
						   mask,
						   mvar.dim(),
						   mvar.defined(),
						   false);
			}
		}
	}


	template <typename P>
	void data<P>::apply(const reexp::exp<P>& exp) {
		data_var<P>& dvar(var(exp.id()));
		data_rel<P> rel( *this, exp.rel() );
		bits& rel_defs( exp_rel_defs_[exp.id()]);

		cvec<P> dim( dvar.dim() );
		int esz = rel.entries().size();
		int state = exp.state();

		// store the exp rel defs bits vector. this is needed for populating
		// statistics needed for re-expression with probabilistic data
		rel_defs.resize(dim.volume());
		rel_defs.fill(true);
		for (int i = 0; i < esz; ++i) {
			rel.setand_var_def_bits(i, rel_defs);
		}

		reexp::bits& states = dvar.states();

		const std::vector<reexp::implmask<P>>& applied_masks( exp.expmasks() );
		// mark affected variables as dirty:
		int v = util::next_version();
		for (const reexp::implmask<P>& m : applied_masks) {
			vars_[m.var_->id()].version_ = v;
		}

		states.fill(true); // by default

		bits firstDelta; // def delta for the given variable
		firstDelta.resize(var(rel.entries()[0].var_->id()).defined().size());

		// first figure out what states could be true, if we don't take into
		// account the possibility that states may exclude each other
		//
		int rowdim = dim.rowdim();
		int rowlength = dim[rowdim];
		for (int j = esz-1; j >= 0; --j) {
			bool varstate = rel.rel_.varState(j, state);
			const rel_entry<P>& e = rel.entries()[j];

			const data_var<P>& dv = var(e.var_->id());
			const reexp::bits& srcdef = dv.defined();
			const reexp::bits& srcstates = dv.states();

			int destidx = 0;

			int srcshift = dv.dim().offset(e.shift_);
			reexp::dim_row_iterator<P> i(dim);
			while (i) {
				int srcidx = dv.dim().offset(i.begin())+srcshift;
				states.from(destidx, rowlength) &= srcdef.from(srcidx, rowlength);
				if (j == 0) {
					firstDelta.from(srcidx, rowlength) |= states.from(destidx, rowlength);
				}
				if (varstate) {
					states.from(destidx, rowlength) &= srcstates.from(srcidx, rowlength);
				} else {
					states.from(destidx, rowlength).andNeg(srcstates.from(srcidx, rowlength));
				}
				destidx += rowlength;
				++i;
			}
		}

		/*
		 * this is the mask that is applied on the applied expression
		 * it is a bit trickier, because we are determining the expression
		 * states at the same time we are applying the mask.
		 */
		const reexp::implmask<P>* expmask = 0;
		for (const reexp::implmask<P>& mask : applied_masks) {
			if (mask.var_->id() == exp.id()) {
				expmask = &mask;
				break;
			}
		}

		// second: exclude states and mark implicit values using the mask
		//         for the applied expression.
		//
		reexp::bits& def = dvar.defined();
		def.fill(true);
		int chunkcount = states.chunks().size();
		for (int i = 0; i < chunkcount; ++i) {
			bchunk ichunk = ~def.chunks()[i];
			bchunk schunk = states.chunks()[i];
			if (schunk || ichunk) {
				bchunk mask = 1;
				for (size_t j = 0; j < bchunk_bsize; ++j) {
					if (ichunk & mask) {
						size_t w = i*bchunk_bsize + j;
						if (w >= size_t(def.size())) break;
						def[w] = true; // turn it back true;
						states[w] = false; // set state false instead;

						ichunk = ~def.chunks()[i]; // restore chunks
						schunk = states.chunks()[i];
					} else if (schunk & mask) {
						data_var<P>& dv = vars_[expmask->var_->id()];
						int w = i*bchunk_bsize + j;
						cvec<P> at = dim.at(w);
						apply_mask(at, *expmask, dv.dim(), dv.defined(), false);

						def[w] = true; // turn it back true;
						// revive chunks
						ichunk = ~def.chunks()[i];
						schunk = states.chunks()[i];
					}
					mask <<= 1;
				}
			}
		}

		// third, apply masks on other variables. this operation basically
		//       mark all the states implicit determined by the expression states.
		//
		bits defdelta; // def delta for the given variable
		for (const reexp::implmask<P>& mask : applied_masks) {
			if (&mask != expmask) {
				data_var<P>& dv = vars_[mask.var_->id()];
				defdelta.resize(dv.defined().size());
				defdelta.fill(false);
				// following can likely be optimized even further:
				for (int i = 0; i < chunkcount; ++i) {
					bchunk schunk = states.chunks()[i];
					if (schunk) {
						bchunk m = 1;
						for (size_t j = 0; j < bchunk_bsize; ++j) {
							if (schunk & m) {
								cvec<P> at = dim.at(i*bchunk_bsize + j);
								apply_mask(at, mask, dv.dim(), defdelta, true);
							}
							m<<=1;
						}
					}
				}
				defdelta &= dv.defined();
				if (mask.var_ == rel.entries()[0].var_) {
					defdelta |= firstDelta;
				}
				dv.defined().andNeg(defdelta);
				if (obs_) obs_->states_undefined(dv, defdelta, v);
			}
		}

	}
}


#endif //REEXP_DATA_I_H_
