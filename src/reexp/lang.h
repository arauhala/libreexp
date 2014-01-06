/*
 * lang.h
 *
 *  Created on: Dec 15, 2010
 *      Author: arauhala
 */

#ifndef REEXP_LANG_H_
#define REEXP_LANG_H_

#include "ctx.h"
#include "bits.h"
#include "util.h"
#include <stdint.h>
#include <values.h>

#include <unordered_set>


namespace reexp {

	template <typename P>
	class lang;

	template <typename P>
	class data;

	template <typename P>
	class rel;

	template <typename P>
	class var;

	template <typename P>
	class orig;

	template <typename P>
	struct dep_var {
		public:
			reexp::cvec<P> shift_;
			reexp::var<P>* var_;
			dep_var(const cvec<P>& s, var<P>& v)
			: shift_(s), var_(&v) {}
	};

	/**
	 * specifies shifted dimensions. dim_ part describes the width, height etc., while
	 * shift_ informs the shift.
	 */
	template <typename P>
	struct ndim {
		reexp::cvec<P> shift_;
		reexp::cvec<P> dim_;
		ndim() : shift_(), dim_() {}
		ndim(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) : shift_(shift), dim_(dim) {}
		bool fit(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim);
	};

	inline void bitref_assign_and(bits_ref dest, const_bits_ref src);

	inline void bitref_assign_and_neg(bits_ref dest, const_bits_ref src);

	template <typename P, typename _rowop>
	void generic_blit(const cvec<P>& destdim, const cvec<P>& destshift, reexp::bits& dest,
			     	  const cvec<P>& srcdim, const cvec<P>& srcshift, const reexp::bits& src,
			     	  const cvec<P>& clip, _rowop op);

	template <typename P>
	struct bitmatrix {
		reexp::ndim<P> ndim_;
		reexp::bits bits_;
		bitmatrix() : ndim_(), bits_() {}
		bitmatrix(const ndim<P>& ndim)
		:   ndim_(ndim),
		    bits_(ndim_.dim_.volume()) {}
		bitmatrix(const ndim<P>& nd, const bits& b)
		:   ndim_(nd),
		    bits_(b) {}
		bitmatrix(const cvec<P>& d)
		:   ndim_(cvec<P>(), d),
		    bits_(ndim_.dim_.volume()) {}
		bitmatrix(const cvec<P>& d, const bits& b)
		:   ndim_(cvec<P>(), d),
		    bits_(b) {}
		bool operator[](const reexp::cvec<P>& at) const {
			return bits_[ndim_.dim_.offset(at)];
		}
		bit_ref operator[](const reexp::cvec<P>& at) {
			return bits_[ndim_.dim_.offset(at - ndim_.shift_)];
		}
		void unite(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim);
		bool fit(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim);
		void resizing_blit(const reexp::cvec<P>& shift, const bitmatrix& v);
		void resizing_set(const cvec<P>& shift, bool v);
		void blit(const cvec<P>& shift, const bitmatrix& v);
		void blitAnd(const cvec<P>& shift, const bitmatrix& v);
		void blitAndNeg(const cvec<P>& shift, const bitmatrix& v);
		void blitAndNegBefore(const cvec<P>& shift, const bitmatrix& v, const cvec<P>& before);
		bool overlap(const cvec<P>& shift, const bitmatrix& v) const;
	};

	// mask marking that which variables are implicitly known
	template <typename P>
	struct implmask {
		const reexp::var<P>* var_;
		bitmatrix<P> mask_;
		inline implmask(const reexp::var<P>* var)
		: var_(var), mask_() {
		}
		inline implmask() : var_(), mask_() {}
		inline implmask(const implmask& o) : var_(o.var_), mask_(o.mask_) {}
		inline void unite(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
			mask_.unite(shift, dim);
		}
		inline void blit(const cvec<P>& shift, const implmask& v) {
			mask_.blit(shift, v.mask_);
		}
		inline void blitAndNeg(const cvec<P>& shift, const implmask& v) {
			mask_.blitAndNeg(shift, v.mask_);
		}
		inline void blitAndNegBefore(const cvec<P>& shift, const implmask& v, const cvec<P>& before) {
			mask_.blitAndNegBefore(shift, v.mask_, before);
		}
		inline void resizing_blit(const cvec<P>& shift, const implmask& v) {
			mask_.resizing_blit(shift, v.mask_);
		}
		inline void resizing_set(const cvec<P>& shift, bool v) {
			mask_.resizing_set(shift, v);
		}
	};

	template <typename P>
	class var {
		protected:
			size_t id_;
			reexp::ctx<P> ctx_;
			std::vector<rel<P>*> rels_;
			std::vector<dep_var<P>> deps_;
			std::vector<implmask<P> > rootmasks_; // used for checking the overlap
			std::vector<implmask<P> > expmasks_; // used for marking excluded exps
			bool disabled_;
			double prioriP_;

			// utility function
			static implmask<P>& have_implmask(std::vector<implmask<P> >& vm,
								  	  	  	  const reexp::var<P>& v);
		public:

			var(const reexp::ctx<P>& ctx, double prioriP);
			virtual ~var();
			int hash() const;
			double prioriP() const;
			void setPrioriP(double prioriP);
			bool disabled() const;
			void setDisabled(bool disabled);
			const std::vector<implmask<P> >& rootmasks() const;
			const std::vector<implmask<P> >& expmasks() const;
			const reexp::ctx<P>& ctx() const;
			const std::vector<rel<P>*>& rels() const;
			void add_rel(rel<P>* rel);
			void remove_rel(rel<P>* rel);
			void set_id(int id);
			size_t id() const;

			/**
			 * In principle, the dependencies should be filled.
			 */
			void add_dep(const reexp::cvec<P>& shift,
						 reexp::var<P>& v);
			const std::vector<dep_var<P>>& deps() const;
		};

	/**
	 * Variable. Is used to organize bits.
	 */
	template <typename P>
	class orig : public var<P> {
		public:
			orig(const ctx<P>& c) : var<P>(c, 0.5f) {}
			void init();
	};


	template <typename P>
	struct rel_entry {
		cvec<P> shift_;
		var<P>* var_;
		rel_entry(const cvec<P>& shift, var<P>& var)
		: shift_(shift),
		  var_(&var) {}
		bool get(const data<P>& data, const cvec<P>& at) const {
			return var_->get(data, at + shift_);
		}
	};
	template <typename P>
	bool operator < (const reexp::rel_entry<P>& e1,
					 const reexp::rel_entry<P>& e2);

	template <typename P>
	class rel {
		private:
			int id_;
			bool disabled_;
			util::hash_code hash_;
			reexp::ctx<P> ctx_;
			std::vector<rel_entry<P> > e_;
		public:
			rel(int id, const reexp::ctx<P>& c);
			double statePrioriP(int state) const;
			void sort_entries();
			void finalize();
			void disable();
			bool disabled() const;
			int id() const;
			void add_var(const cvec<P>& shift, var<P>& var);
			const std::vector<rel_entry<P> >& entries() const;
			const reexp::ctx<P>& ctx() const;
			int stateCount() const;
			int varCount() const;
			static bool varState(int var, int state);
			int hash() const;
			bool operator==(reexp::rel<P>& r) const;
	};

	template <typename P>
	struct rel_ptr_hash {
		typedef int(*func_type)(rel<P>*);
		static int hash(rel<P>* p) {
			return p->hash();
		}
	};

	/**
	 * Expression
	 *
	 */
	template <typename P>
	class exp : public var<P> {
		private:
			const reexp::rel<P>& rel_;
			int state_;
			std::vector<reexp::rel<P>*> rels_;
		public:
			exp(const reexp::rel<P>& rel, int state);
			int min_var_idx() const;
			// rootmasks need to be renamed root masks
			// state of rel needs to be checked, if true, use root masks
			// if false, blit directly.
			void fill_rootmasks();
			void fill_expmasks();
			void init(int id);
			void add_deps(const cvec<P>& shift, const reexp::exp<P>& at);
			int hash() const;
			static bool is_overlap(const reexp::var<P>& v1,
								   	   	  const cvec<P>& shift,
								   	   	  const reexp::var<P>& v2);
			static bool is_overlap(const std::vector<reexp::rel_entry<P> >& rvars);
/*			static bool is_overlap(bitmatrix<P>& matrix,
								   std::vector<int>& vars,
								   const std::vector<reexp::rel_entry<P> >& rvars);*/
			void gen_rels(lang<P>& lang);
			cvec<P> ctx_dim(const data<P>& data) const;
			/** offset of pixel as compared to exp data offset*/
			int offset(const cvec<P>& dim, const cvec<P>& at) const;
			int state() const;
			void add_rel(reexp::rel<P>* rel);
			const reexp::rel<P>& rel() const;
			const std::vector<reexp::rel<P>*>& rels() const;
	};

	template <typename P>
	class lang_obs {
		public:
			virtual void var_added(const var<P>& exp) = 0;
			virtual void rel_added(const rel<P>& exp) = 0;
	};


	/**
	 * Organization of data. There should be one organization per
	 * problem domain. This object also easily reveals semantic
	 * meaning of data. Organization of data allows and helps learning.
	 */
	template <typename P>
	class lang {
		private:
		public: // needed for measuring performance
			/** Organizes bits in variables */
			util::arrays_list<reexp::orig<P> > origs_;
			/** Expressions */
			util::arrays_list<reexp::exp<P> > exps_;
			/** Organizes variables in relations */
			util::arrays_list<reexp::rel<P> > rels_;
			/**/
			lang_obs<P>* obs_;

		public:
			lang();
			lang(const lang& l);
			const reexp::orig<P>& orig(int i) const;
			const reexp::exp<P>& exp(int i) const;
			const reexp::rel<P>& rel(int i) const;
			void disable_rel(int i);
			int enabled_rel_count() const;
			int rel_count() const;
			reexp::var<P>& var(int i);
			const reexp::var<P>& var(int i) const;
			int orig_count() const;
			int exp_count() const;
			int var_count() const;
			const reexp::exp<P>& exp_back() const;
			void add_orig(const reexp::orig<P>& o);
			/** after allocating relation it still needs to be either finalized or wasted */
			reexp::rel<P>& alloc_rel(const reexp::ctx<P>& ctx);
			/** note this can be used for allocated relation before it is wasted / finalized
			 * still, the relation entries need to be sorted before calling this method.
			 * finalization will sort the relation anyway */
			bool is_rel_back_unique();
			void waste_rel();
			void rel_done();
			void add_exp(const reexp::rel<P>& rel, int state);
			void set_obs(lang_obs<P>& o);
			void unset_obs();
	};



}

#endif /* REEXP_LANG_H_ */
