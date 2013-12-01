/*
 * data.h
 *
 *  Created on: Dec 15, 2010
 *      Author: arauhala
 */

#ifndef REEXP_DATA_H_
#define REEXP_DATA_H_

#include "lang.h"

namespace reexp {

	/**
	 * Variable at some specific information sample
	 */
	template <typename P>
	struct data_var {
		typedef cond_bit_ref cbit;
		typedef const_cond_bit_ref const_cbit;

		data<P>& data_;
		const reexp::var<P>& var_;
		cond_bits bits_;
		int version_;

		inline data_var(data<P>& data, const reexp::var<P>& var)
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
		inline reexp::bits& defined() {
			return bits_.defined();
		}
		inline const reexp::bits& defined() const {
			return bits_.defined();
		}
		inline reexp::bits& states() {
			return bits_.states();
		}
		inline const reexp::bits& states() const {
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
		void set_implicit(const cvec<P>& at, int cause);

		inline void touch(int version) {
			version_ = version;
		}

		void touch_deps(int version);

		inline const reexp::var<P>& var() const {
			return var_;
		}
	};

	template <typename P>
	struct shifted_data_var {
		inline shifted_data_var(data_var<P>& var, const cvec<P>& shift)
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
	public:
		typedef shifted_data_var<P> rvar;
		inline data_rel(data<P>& d, const reexp::rel<P>& r)
		:	data_(d), rel_(r) {}
		inline reexp::cvec<P> dim() const {
			return rel_.ctx().dim(data_.dim());
		}
		inline rvar var(int i) const {
			return shifted_data_var<P>(data_.var(rel_.entries()[i].var_->id()),
									   rel_.entries()[i].shift_);
		}
		inline const std::vector<rel_entry<P> >& entries() const {
			return rel_.entries();
		}
		inline const reexp::rel<P>& rel() const {
			return rel_;
		}

		void setand_var_def_bits(size_t var, reexp::bits& defs) const;

		void setand_var_state_bits(int var, reexp::bits& states) const;

		reexp::data<P>& data_;
		const reexp::rel<P>& rel_;
	};


	template <typename P>
	class data_obs {
		public:
			virtual ~data_obs() {}
			virtual void rel_added(const data_rel<P>& rel) = 0;
			virtual void var_added(const data_var<P>& var) = 0;
			virtual void states_undefined(const data_var<P>& var,
										  const reexp::bits& undefs,
										  int version) = 0;
	};

	template <typename P>
	class data : public lang_obs<P> {
		public:
			data(reexp::lang<P>& lang, const reexp::cvec<P>& dim);
 			// This is required, if the language already contains
			// expressions. This will add space for the expression variables
			// and do re-expression for data. DO NOT call twice.
			// It is important that each expression will be applied
			// only once.
			//
			void apply_exps();
			void set_obs(reexp::data_obs<P>& obs) {
				obs_ = &obs;
			}
			void rel_added(const rel<P>& rel) {
				data_rel<P> dv(*this, rel);
				if (obs_) obs_->rel_added(dv);
			}
			void var_added(const var<P>& var) {
				vars_.push_back(data_var<P>(*this, var));
				exp_rel_defs_.push_back(bits());
				//bits_.resize(bits_.size() + vars_.back().bsize_);
				const exp<P>* e = dynamic_cast<const exp<P>*>(&var);
				if (e) apply(*e);
				if (obs_) obs_->var_added(vars_.back());
			}
			const reexp::lang<P>& lang() const {
				return lang_;
			}
			reexp::lang<P>& lang() {
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
			const bits& exp_rel_defs(int i) const {
				return exp_rel_defs_[i];
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
			void apply_state(const reexp::exp<P>& exp, cvec<P> at, bool state);
		private:
			void apply(const reexp::exp<P>& exp);

			reexp::lang<P>& lang_;
			cvec<P> dim_;
			util::arrays_list<reexp::data_var<P> > vars_;
			std::vector<bits> exp_rel_defs_; //
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
	void apply_mask(const cvec<P>& at,
					const reexp::implmask<P>& mask,
					const cvec<P>& targetdim,
					reexp::bits& targetdef,
					bool deltaMode) {
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

}

#endif /* REEXP_DATA_H_ */
