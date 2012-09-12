/*
 * data.h
 *
 *  Created on: Dec 15, 2010
 *      Author: arauhala
 */

#ifndef DATA_H_
#define DATA_H_

#include "lang.h"

namespace explib {

	/**
	 * Variable at some specific information sample
	 */
	template <typename P>
	struct data_var {
		typedef cond_bit_ref cbit;
		typedef const_cond_bit_ref const_cbit;
		data_var(data<P>& data, const explib::var<P>& var)
		:	data_(data),
			var_(var),
			bits_(),
			version_() {
			bits_.resize(dim().volume());
			bits_.defined().fill(true);
		}
		inline cvec<P> dim() const {
			return var_.ctx().dim(data_.dim());
		}
		inline int index(const cvec<P>& at) const {
			return var_.ctx().offset(data_.dim(), at);
		}
		inline cond_bits_ref bitrow(const cvec<P>& at, int rowl) const {
			return bits_.from(index(at), rowl);
		}
		inline cond_bits& bits() {
			return bits_;
		}
		inline const cond_bits& bits() const {
			return bits_;
		}
		inline explib::bits& defined() {
			return bits_.defined();
		}
		inline const explib::bits& defined() const {
			return bits_.defined();
		}
		inline explib::bits& states() {
			return bits_.states();
		}
		inline const explib::bits& states() const {
			return bits_.states();
		}
		inline const_cbit operator[](const cvec<P>& at) const {
			return bits_[index(at)];
		}
		inline cbit operator[](const cvec<P>& at) {
			return bits_[index(at)];
		}

		/**
		 * basically, if e -> x, and x becomes implicitly known,
		 * then e becomes implicitly known, except that...
		 */
		void set_implicit(const cvec<P>& at, int cause) {
			int i = index(at);

			if (i < cause) {
				cond_bits& bits = bits_;
				cbit b = bits[i];
				b = false;
				if (*b) {
					const explib::exp<P>* e = dynamic_cast<const explib::exp<P>*>(&var_);
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
						explib::data_var<P>& v(data_.var(i->var_->id()));
						explib::cvec<P> limits = v.var().ctx().dim(data_.dim());
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

		void touch(int version) {
			version_ = version;
		}

		void touch_deps(int version) {
			if (version_ != version) { // prevents recursion
				touch(version);
				const explib::exp<P>* e = dynamic_cast<const explib::exp<P>*>(&var_);
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

		const explib::var<P>& var() const {
			return var_;
		}

		data<P>& data_;
		const explib::var<P>& var_;
		cond_bits bits_;
		int version_;
	};

	template <typename P>
	struct shifted_data_var {
		shifted_data_var(data_var<P>& var, const cvec<P>& shift)
		:	var_(var), shift_(shift) {}
		inline cond_bits_ref bitrow(const cvec<P>& at, int rowl) const {
			return var_.bitrow(at + shift_, rowl);
		}
		inline int index(const cvec<P>& at) const {
			return var_.index(at + shift_);
		}
		inline const_cond_bit_ref operator[](const cvec<P>& at) const {
			return var_[at + shift_];
		}
		inline cond_bit_ref operator[](const cvec<P>& at) {
			return var_[at + shift_];
		}
		inline void set_implicit(const cvec<P>& at, int cause) {
			var_.set_implicit(at + shift_, cause);
		}
		data_var<P>& var_;
		const cvec<P>& shift_;
	};

	template <typename P>
	class data;

	template <typename P>
	struct data_rel {
	private:
		data_rel() {}
	public:
		typedef shifted_data_var<P> rvar;
		data_rel(data<P>& d, const explib::rel<P>& r)
		:	data_(d), rel_(r) {}
		inline explib::cvec<P> dim() const {
			return rel_.ctx().dim(data_.dim());
		}
		inline rvar var(int i) const {
			return shifted_data_var<P>(data_.var(rel_.entries()[i].var_->id()),
									   rel_.entries()[i].shift_);
		}
		inline const std::vector<rel_entry<P> >& entries() const {
			return rel_.entries();
		}
		const explib::rel<P>& rel() const {
			return rel_;
		}

		void setand_var_def_bits(size_t var, explib::bits& defs) const {
			const explib::cvec<P> d = dim();
			int rowdim = d.rowdim();
			assert(defs.size() >= d.volume());

			const rel_entry<P>& e = rel_.entries()[var];
			const data_var<P>& dv = data_.var(e.var_->id());

			int rowlength = d[rowdim];
			int destidx = 0;

			explib::dim_row_iterator<P> i(d);
			const explib::var<P>& v = dv.var_;

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
				const explib::cvec<P> srcdim(dv.dim());
				const bits& srcdefs = dv.defined();
				int srcshift = srcdim.offset(e.shift_);
				while (i) {
					int srcidx = srcdim.offset(i.begin())+srcshift;
					defs.from(destidx, rowlength) &= srcdefs.from(srcidx, rowlength);
					destidx += rowlength;
					++i;
				}
			}
/*			int srcrowlength = dv.dim()[rowdim];
			int srcidx = srcdim.offset(e.shift_);
			while (destidx < defs.size()) {
				defs.from(destidx, rowlength) &= srcdefs.from(srcidx, rowlength);
				destidx += rowlength;
				srcidx += srcrowlength;
			}*/
		}

		void setand_var_state_bits(int var, explib::bits& states) const {
			const explib::cvec<P> d = dim();
			int rowdim = d.rowdim();
			assert(states.size() >= d.volume());

			const rel_entry<P>& e = rel_.entries()[var];
			const data_var<P>& dv = data_.var(e.var_->id());

			int rowlength = d[rowdim];
			int destidx = 0;

			explib::dim_row_iterator<P> i(d);
			const explib::var<P>& v = dv.var_;

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
				const explib::cvec<P> srcdim(dv.dim());
				const bits& srcstates = dv.states();
				int srcshift = srcdim.offset(e.shift_);
				while (i) {
					int srcidx = srcdim.offset(i.begin())+srcshift;
					states.from(destidx, rowlength) &= srcstates.from(srcidx, rowlength);
					destidx += rowlength;
					++i;
				}
			}



/*			int srcshift = srcdim.offset(e.shift_);
			explib::dim_row_iterator<P> i(d);
			while (i) {
				int srcidx = srcdim.offset(i.begin())+srcshift;
				states.from(destidx, rowlength) &= srcstates.from(srcidx, rowlength);
				destidx += rowlength;
				++i;
			}*/

/*			int srcrowlength = dv.dim()[rowdim];
			int srcidx = srcdim.offset(e.shift_);
			while (destidx < states.size()) {
				states.from(destidx, rowlength) &= srcstates.from(srcidx, rowlength);
				destidx += rowlength;
				srcidx += srcrowlength;
			}*/
		}
		explib::data<P>& data_;
		const explib::rel<P>& rel_;
	};


	template <typename P>
	class data_obs {
		public:
			virtual void rel_added(const data_rel<P>& rel) = 0;
			virtual void var_added(const data_var<P>& var) = 0;
			virtual void states_undefined(const data_var<P>& var,
										  const explib::bits& undefs,
										  int version) = 0;
	};

	template <typename P>
	class data : public lang_obs<P> {
		public:
			data(explib::lang<P>& lang, const explib::cvec<P>& dim)
			: lang_(lang), dim_(dim), obs_() {
				lang.set_obs(*this);
				for (int i = 0; i < lang.orig_count(); i++) {
					var_added(lang.orig(i));
				}
			}
 			// This is required, if the language already contains
			// expressions. This will add space for the expression variables
			// and do re-expression for data. DO NOT call twice.
			// It is important that each expression will be applied
			// only once.
			//
			void apply_exps() {
				for (int i = vars_.size(); i < lang_.var_count(); i++) {
					const explib::exp<P>* e(
						dynamic_cast<const explib::exp<P>*>(&lang_.var(i)));
					if (e) var_added(*e);
				}
			}
			void set_obs(explib::data_obs<P>& obs) {
				obs_ = &obs;
			}
			void rel_added(const rel<P>& rel) {
				data_rel<P> dv(*this, rel);
				if (obs_) obs_->rel_added(dv);
			}
			void var_added(const var<P>& var) {
				vars_.push_back(data_var<P>(*this, var));
				//bits_.resize(bits_.size() + vars_.back().bsize_);
				const exp<P>* e = dynamic_cast<const exp<P>*>(&var);
				if (e) apply(*e);
				if (obs_) obs_->var_added(vars_.back());
			}
			const explib::lang<P>& lang() const {
				return lang_;
			}
			explib::lang<P>& lang() {
				return lang_;
			}
			const cvec<P>& dim() const {
				return dim_;
			}
			void set_dim(const cvec<P>& d) {
				dim_ = d;
			}
			const data_var<P>& var(int i) const {
				assert(i < vars_.size());
				return vars_[i];
			}
			data_var<P>& var(int i) {
				assert(i < vars_.size());
				return vars_[i];
			}
			data_rel<P> rel(int i) {
				return data_rel<P>(*this, lang_.rel(i));
			}
			data_rel<P> rel(int i) const {
				return data_rel<P>(const_cast<data&>(*this), lang_.rel(i));
			}
			size_t var_count() const {
				return vars_.size();
			}
			// size in bits
			size_t size() const {
				size_t rv = 0;
				for (int i = 0; i < vars_.size(); ++i) {
					rv += vars_[i].bits().size();
				}
				return rv;
			}
		private:
			void apply(const explib::exp<P>& exp);
//			void row_apply(const explib::exp<P>& exp);

			explib::lang<P>& lang_;
			cvec<P> dim_;
			//cond_bits bits_;
			util::arrays_list<explib::data_var<P> > vars_;
			data_obs<P>* obs_;

	};



	/**
	 * The redundancy elimination step is actually very
	 * challenging. It was done on Java side successfully,
	 * but unfortunately the mechanism was clumsy & expensive
	 * computation wise.
	 *
	 * Basically, whether X is true or false, some amount of
	 * expressed variable states becomes implicitly known.
	 * Because the states becomes implicitly known, some
	 * states of expressions (somewhere) can be implicitly
	 * determined and will not need explicit declaration.
	 *
	 * With explicits & implicits, the critical thing to
	 * remember is the definition order.
	 *
	 * Few considerations:
	 *
	 *    * Also relationship statistics should be updated
	 *
	 *    * All expressions related to the named variable
	 *      should be updated.
	 *
	 *    * By redundancy check, expressions can only be
	 *      determined to be false
	 *
	 *    * Is one-time scan or iterative point modification
	 *      cheaper?
	 *
	 *    * Is it possible to have some kind of pre-constructed
	 *    	filter, that can be applied. E.g. it could be
	 *    	just a set of index offsets (or cvec's), which
	 *    	points out dependent locations in the data.
	 *
	 * Redundancy elimination rules:
	 *
	 *    * expression
	 *    * known last
	 *    * implicit -> false expressions become impossible
	 *
	 * Filter mechanism
	 *
	 *    * if state v of V becomes implicit
	 *    * E determines v
	 *    * v implicit, E -> v, E -> implicit
	 *
	 * Filter usage:
	 *
	 * int eidx = ...;
	 * int vidx = ...;
	 * filter& filter = ...;
	 * for (int i = 0; i < filter.size(); i++) {
	 *     int e2idx= filter[vidx];
	 *     if (e2idx < vidx) {
	 *         bits_[e2idx] = false; // make implicit
	 *     	   // recursion
	 *     }
	 * }
	 *
	 * Design?
	 *
	 * filter {
	 *    exp& dep_;
	 *    cvec<P> offset_;
	 *    operator[] {}
	 * }
	 *
	 */

	template <typename P>
	void apply_mask(const cvec<P>& at, const explib::implmask<P>& mask, const cvec<P>& targetdim, explib::bits& targetdef, bool deltaMode) {
		const bitmatrix<P>& m = mask.mask_;
		cvec<P> targetbegin = at;
		targetbegin += m.ndim_.shift_;

		int rowdim = m.ndim_.dim_.rowdim();

		dim_row_iterator<P> i(m.ndim_.dim_);
		for (int j = rowdim; j < P::DIM; ++j) {
			if (targetbegin[j] < 0) {
				i.settable_begin()[j] = -targetbegin[j];
			}
			if (i.settable_dim()[j] > targetdim[j] - targetbegin[j]) {
				i.settable_dim()[j] = targetdim[j] - targetbegin[j];
			}
		}
		int rowlength = i.settable_dim()[rowdim] - i.settable_begin()[rowdim];

		for (; i; ++i) {
			const cvec<P>& sourcebegin = i.begin();

			int sourceIdx = m.ndim_.dim_.offset(sourcebegin);
			int targetIdx = targetdim.offset(sourcebegin + targetbegin);

			if (deltaMode) { // mark only delta
				targetdef.between(targetIdx, rowlength) |=
					m.bits_.from(sourceIdx, rowlength);
			} else {	// apply directly
				targetdef.between(targetIdx, rowlength).andNeg(
					m.bits_.from(sourceIdx, rowlength));
			}
		}

	}


	template <typename P>
	void data<P>::apply(const explib::exp<P>& exp) {
		data_var<P>& dvar(var(exp.id()));
		data_rel<P> rel( *this, exp.rel() );
		cvec<P> dim( dvar.dim() );
		int esz = rel.entries().size();
		int state = exp.state();


		typedef typename data_var<P>::cbit cbit;

		explib::bits& states = dvar.states();

		const std::vector<explib::implmask<P>>& applied_masks( exp.expmasks() );
#define DELTAS
		// mark affected variables as dirty:
		int v = util::next_version();
		for (const explib::implmask<P>& m : applied_masks) {
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
			const explib::bits& srcdef = dv.defined();
			const explib::bits& srcstates = dv.states();

			int destidx = 0;

/*			int srcshift = srcdim.offset(e.shift_);
			int srcrowlength = dv.dim()[rowdim];
			int shift_offset = dv.dim().offset(e.shift_);
			explib::dim_row_iterator<P> i(d);
			while (i) {
				int srcidx = srcdim.offset(i.begin())+srcshift;
				defs.from(destidx, rowlength) &= srcdefs.from(srcidx, rowlength);
				destidx += rowlength;
				++i;
			}*/
			int srcshift = dv.dim().offset(e.shift_);
			explib::dim_row_iterator<P> i(dim);
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
#if 0
			int srcidx = shift_offset;
			while (destidx < states.size()) {
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
				srcidx += srcrowlength; // note this likely doesn't work for DIM > 2, delta between rows varies
			}
#endif
		}

		const explib::implmask<P>* expmask = 0;
		for (const explib::implmask<P>& mask : applied_masks) {
			if (mask.var_->id() == exp.id()) {
				expmask = &mask;
				break;
			}
		}

		// second: exclude states and mark implicit values
		//
		explib::bits& def = dvar.defined();
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
		//printf("exp %d\n", int(expmask->var_->id()));

		bits defdelta; // def delta for the given variable
		for (const explib::implmask<P>& mask : applied_masks) {
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
#ifdef DELTAS
				if (obs_) obs_->states_undefined(dv, defdelta, v);
#endif
			}
		}
#if 0

		// now we (should) have the expression bits known
		// calculate deltas for applying effect

		for (int i = 0; i < chunkcount; ++i) {
			bchunk ichunk = ~def.chunks()[i];
			bchunk schunk = states.chunks()[i];
			if (schunk || ichunk) {
				bchunk mask = 1;
				for (int j = 0; j < bchunk_bsize; ++j) {
					if (ichunk & mask) {
						int w = i*bchunk_bsize + j;
						if (w >= def.size()) break;
						def[w] = true; // turn it back true;
						states[w] = false; // set state false instead;
						ichunk = ~def.chunks()[i]; // restore chunks
						schunk = states.chunks()[i];
					} else if (schunk & mask) {
						int w = i*bchunk_bsize + j;
						cvec<P> c = dim.at(w);

						for (const explib::implmask<P>& mask : applied_masks) {
							const bitmatrix<P>& m = mask.mask_;
							data_var<P>& dv = vars_[mask.var_->id()];
							cvec<P> targetdim = dv.dim();
							cvec<P> targetbegin = c;
							targetbegin += m.ndim_.shift_;

							int rowdim = m.ndim_.dim_.rowdim();

							dim_row_iterator<P> i(m.ndim_.dim_);
							for (int j = rowdim; j < P::DIM; ++j) {
								if (targetbegin[j] < 0) {
									i.settable_begin()[j] = -targetbegin[j];
								}
								if (i.settable_dim()[j] > targetdim[j] - targetbegin[j]) {
									i.settable_dim()[j] = targetdim[j] - targetbegin[j];
								}
							}
							int rowlength = i.settable_dim()[rowdim] - i.settable_begin()[rowdim];

							explib::bits& targetdef = dv.defined();
							for (; i; ++i) {
								const cvec<P>& sourcebegin = i.begin();

								int sourceIdx = m.ndim_.dim_.offset(sourcebegin);
								int targetIdx = targetdim.offset(sourcebegin + targetbegin);

								targetdef.between(targetIdx, rowlength).andNeg(
									m.bits_.from(sourceIdx, rowlength));
							}
						}
						def[w] = true; // turn it back true;
						// revive chunks
						ichunk = ~def.chunks()[i];
						schunk = states.chunks()[i];
					}
					mask <<= 1;
				}
			}
		}
#endif

	}

}

#endif /* DATA_H_ */
