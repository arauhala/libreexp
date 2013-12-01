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
	 * Dimensions
	 */
	template <typename P>
	struct ndim {
		reexp::cvec<P> shift_;
		reexp::cvec<P> dim_;
		ndim() : shift_(), dim_() {}
		bool fit(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
			bool changed = false;
			for (int i = 0; i < P::DIM; ++i) {
				if (dim[i]) {
					int delta = shift[i] -shift_[i];
					if (delta < 0) {
						shift_[i] += delta;
						dim_[i] -= delta;
						changed = true;
					}
					int up = shift[i] + dim[i] - shift_[i];
					if (up > dim_[i]) {
						dim_[i] = up;
						changed = true;
					}
				}
			}
			return changed;
		}
	};

	template <typename P>
	struct bitmatrix {
		reexp::ndim<P> ndim_;
		reexp::bits bits_;
		bitmatrix() : ndim_(), bits_() {}
		bitmatrix(const ndim<P>& ndim)
		:   ndim_(ndim),
		    bits_(ndim_.dim_.volume()) {}
		bool operator[](const reexp::cvec<P>& at) const {
			return bits_[ndim_.dim_.offset(at)];
		}
		bit_ref operator[](const reexp::cvec<P>& at) {
			return bits_[ndim_.dim_.offset(at - ndim_.shift_)];
		}
		void unite(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
			ndim_.fit(shift, dim);
			if (ndim_.dim_.volume() > bits_.size()) {
				bits_.resize(ndim_.dim_.volume());
			}
		}
		bool fit(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
			reexp::ndim<P> ndim(ndim_);
			bool resize = ndim.fit(shift, dim);
			if (resize) {
				bitmatrix<P> resized(ndim);
				resized.blit(cvec<P>(), *this);
				*this = resized;
			}
			return resize;
		}
		void resizing_blit(const reexp::cvec<P>& shift, const bitmatrix& v) {
			fit(shift + v.ndim_.shift_, v.ndim_.dim_);
			blit(shift, v);
		}
		void resizing_set(const cvec<P>& shift, bool v) {
			fit(shift, cvec<P>::unit_vector());
			int offset = ndim_.dim_.offset(shift - ndim_.shift_);
			bits_[offset] = v;
		}
		void blit(const cvec<P>& shift, const bitmatrix& v) {
			cvec<P> s = shift + v.ndim_.shift_ - ndim_.shift_ ;
			for (dim_row_iterator<P> i(v.ndim_.dim_); i; ++i) {
				const cvec<P>& b = i.begin();
				int bidx = v.ndim_.dim_.offset(b);
				int len = i.length();
				int targetIdx = ndim_.dim_.offset(b + s);

				bits_.from(targetIdx, len) |= v.bits_.from(bidx, len);
			}
		}
		void blitAndNeg(const cvec<P>& shift, const bitmatrix& v) {
			cvec<P> s = shift + v.ndim_.shift_ - ndim_.shift_;
			for (dim_row_iterator<P> i(v.ndim_.dim_); i; ++i) {
				const cvec<P>& b = i.begin();
				int bidx = v.ndim_.dim_.offset(b);
				int targetIdx = ndim_.dim_.offset(b + s);
				int len = i.length();

				bits_.from(targetIdx, len).andNeg(v.bits_.from(bidx, len));
			}
		}

		void blitAndNegBefore(const cvec<P>& shift, const bitmatrix& v, const cvec<P>& before) {
			int beforeidx = ndim_.dim_.offset(before - ndim_.shift_);
			cvec<P> s = shift + v.ndim_.shift_ - ndim_.shift_ ;
			for (dim_row_iterator<P> i(v.ndim_.dim_); i; ++i) {
				const cvec<P>& b = i.begin();
				int bidx = v.ndim_.dim_.offset(b);
				int targetIdx = ndim_.dim_.offset(b + s);
				int end = std::min(beforeidx, targetIdx + i.length());
				int len = end - targetIdx;
				if (len <= 0) break;

				bits_.from(targetIdx, len).andNeg(v.bits_.from(bidx, len));
			}
		}

		inline bool overlap(const cvec<P>& shift, const bitmatrix& v) const {

			// duplication, very similar algorithm is found in data<P>::apply_exp
			cvec<P> targetbegin = shift + v.ndim_.shift_ - ndim_.shift_ ;

			int rowdim = v.ndim_.dim_.rowdim();

			// here things will not really work, when we add dimensions
			dim_row_iterator<P> i(v.ndim_.dim_);
			for (int j = rowdim; j < P::DIM; ++j) {
				if (targetbegin[j] < 0) {
					i.settable_begin()[j] = -targetbegin[j];
				}
				if (i.settable_dim()[j] > ndim_.dim_[j] - targetbegin[j]) {
					i.settable_dim()[j] = ndim_.dim_[j] - targetbegin[j];
				}
			}
			int rowlength = i.settable_dim()[rowdim] - i.settable_begin()[rowdim];

			for (; i; ++i) {
				const cvec<P>& sourcebegin = i.begin();

				int sourceIdx = v.ndim_.dim_.offset(sourcebegin);
				int targetIdx = ndim_.dim_.offset(sourcebegin + targetbegin);

				if (bits_.from(targetIdx, rowlength).true_overlap(v.bits_.from(sourceIdx, rowlength))) {
					return true;
				}
			}
			return false;
		}
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
								  	  	  	  const reexp::var<P>& v) {
				size_t found = 0;
				for (; found < vm.size(); ++found) {
					if (vm[found].var_ == &v) break;
				}
				if (found == vm.size()) vm.push_back(implmask<P>(&v));
				return vm[found];
			}

		public:

			var(const reexp::ctx<P>& ctx, double prioriP)
			: id_(),
			  ctx_(ctx),
			  rels_(),
			  deps_(),
			  rootmasks_(),
			  expmasks_(),
			  disabled_(false),
			  prioriP_(prioriP) {}
			virtual ~var() {}
			int hash() const{
				return id;
			}
			double prioriP() const {
				return prioriP_;
			}
			void setPrioriP(double prioriP) {
				prioriP_ = prioriP;
			}
			bool disabled() const {
				return disabled_;
			}
			void setDisabled(bool disabled) {
				disabled_ = disabled;
			}
			inline const std::vector<implmask<P> >& rootmasks() const {
				return rootmasks_;
			}
			inline const std::vector<implmask<P> >& expmasks() const {
				return expmasks_;
			}
			const reexp::ctx<P>& ctx() const {
				return ctx_;
			}
			const std::vector<rel<P>*>& rels() const {
				return rels_;
			}
			void add_rel(rel<P>* rel) {
				if (std::find(rels_.begin(), rels_.end(), rel) == rels_.end()) {
					rels_.push_back(rel);
				}
			}
			void remove_rel(rel<P>* rel) {
				if (std::find(rels_.begin(), rels_.end(), rel) == rels_.end()) {
					rels_.push_back(rel);
				}
			}
			void set_id(int id) {
				id_ = id;
			}
			size_t id() const {
				return id_;
			}
			/**
			 * In principle, the dependencies should be filled.
			 */
			void add_dep(const reexp::cvec<P>& shift,
						 reexp::var<P>& v) {

				deps_.push_back(reexp::dep_var<P>(shift, v));
				// here goes the magic
			}
			const std::vector<dep_var<P>>& deps() const {
				return deps_;
			}
		};

	/**
	 * Variable. Is used to organize bits.
	 */
	template <typename P>
	class orig : public var<P> {
		public:
			orig(const ctx<P>& c) : var<P>(c, 0.5f) {}
			void init() {
				var<P>::rootmasks_.push_back(implmask<P>(this));
				var<P>::rootmasks_.back().resizing_set(reexp::cvec<P>(), true);
			}
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
	class rel {
		private:
			int id_;
			bool disabled_;
			util::hash_code hash_;
			reexp::ctx<P> ctx_;
			std::vector<rel_entry<P> > e_;
		public:
			rel(int id, const reexp::ctx<P>& c)
			: id_(id), disabled_(false), hash_(), ctx_(c), e_() {
				for (int i = 0; i < P::DIM; ++i) {
					hash_ << (ctx_.v_[i] < 0);
				}
			}
			double statePrioriP(int state) const {
				double p = 1;
				for (size_t i = 0; i < e_.size(); ++i) {
					double vp = e_[i].var_->prioriP();
					p *= varState(i, state)?vp:1.-vp;
				}
				return p;
			}
			void finalize() {
				for (auto i = e_.begin(); i != e_.end(); ++i) {
					i->var_->add_rel(this);
				}
			}
			void disable() {
				disabled_ = true;
				for (auto i = e_.begin(); i != e_.end(); ++i) {
					i->var_->remove_rel(this);
				}
			}
			bool disabled() const {
				return disabled_;
			}
			int id() const {
				return id_;
			}
			void add_var(const cvec<P>& shift, var<P>& var) {
				hash_ << shift.hash();
				hash_ << var.id();
				e_.push_back(rel_entry<P>(shift, var));
				ctx_.introduce(var.ctx(), shift);
			}
			inline const std::vector<rel_entry<P> >& entries() const {
				return e_;
			}
			const reexp::ctx<P>& ctx() const {
				return ctx_;
			}
			int stateCount() const {
				return 1 << e_.size();
			}
			int varCount() const {
				return e_.size();
			}
			static bool varState(int var, int state) {
				return ((state>>var)&0x1);
			}
			int hash() const {
				return hash_();
			}
			bool operator==(reexp::rel<P>& r) const {
				if (hash_() != r.hash_()) return false;
				if (ctx_.v_ != r.ctx_.v_) return false;
				if (e_.size() != r.e_.size()) return false;
				for (size_t i = 0; i < e_.size(); i++) {
					if (e_[i].var_ != r.e_[i].var_) {
						return false;
					}
				}
				return true;
			}
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
			exp(const reexp::rel<P>& rel, int state)
			:	var<P>(rel.ctx(), rel.statePrioriP(state)), rel_(rel), state_(state), rels_() {
			}
			int min_var_idx() const {
				size_t minvar = SIZE_MAX;
				for (const rel_entry<P>& e : rel_.entries()) {
					minvar = std::min(minvar, e.var_->id());
				}
				return minvar;
			}
			// rootmasks need to be renamed root masks
			// state of rel needs to be checked, if true, use root masks
			// if false, blit directly.
			void fill_rootmasks() {
				for (size_t i = 0; i < rel_.entries().size(); ++i) {
					var<P>* v = rel_.entries()[i].var_;
					const cvec<P>& shift = rel_.entries()[i].shift_;
					if (rel_.varState(i, state_)) {
						for (const implmask<P>& o : v->rootmasks()) {
							implmask<P>& imask = var<P>::have_implmask(var<P>::rootmasks_, *o.var_);
							imask.resizing_blit(shift, o);
						}
					} else {
						implmask<P>& imask = var<P>::have_implmask(var<P>::rootmasks_, *v);
						imask.resizing_set(shift, true);
					}
				}
			}
			void fill_expmasks() {
				const std::vector<implmask<P>>& rmasks = var<P>::rootmasks();
				// few points: 1) the bits before the first variable don't need
				// to be made implicit. 2) each variable after first variable
				// may reduce its fill_expmask of this variable for bits before
				// that variable.

				size_t minvar = min_var_idx();

				for (const implmask<P>& m : rmasks) {
					const var<P>& v = *m.var_;
					const std::vector<dep_var<P>>& deps = v.deps();
					const bits& b = m.mask_.bits_;
					const cvec<P>& shift = m.mask_.ndim_.shift_;
					const cvec<P>& d = m.mask_.ndim_.dim_;
					if (v.id() >= minvar) {
						implmask<P>& imask = var<P>::have_implmask(var<P>::expmasks_, v);
						imask.resizing_blit(reexp::cvec<P>(), m);
					}
					dim_iterator<P> i(d);
					int offset = 0;
					while (i) {
						if (b[offset]) {
							cvec<P> at = *i;
							at += shift;
							for (const dep_var<P>& dep : deps) {
								const var<P>& dvar = *dep.var_;
								if (dvar.id() >= minvar) {
									implmask<P>& imask = var<P>::have_implmask(var<P>::expmasks_, dvar);
									imask.resizing_set(at + dep.shift_, true);
								}
							}
						}
						++i; offset++;
					}
				}
				// this is mostly optimization
				for (size_t i = 0; i < rel_.entries().size(); ++i) {
					const rel_entry<P>& e = rel_.entries()[i];
					const var<P>& v = *e.var_;
					const cvec<P>& s = e.shift_;
					if (rel_.varState(i, state_)) {
						for (const implmask<P>& rvarmask : v.expmasks()) {
							if (rvarmask.var_->id() >= minvar &&
								rvarmask.var_->id() <= v.id()) {
								implmask<P>& mymask = var<P>::have_implmask(var<P>::expmasks_, *rvarmask.var_);
								if (rvarmask.var_->id() == v.id()) {
									mymask.blitAndNegBefore(s, rvarmask, s); // TODO: blit only pixels before the target exp
								} else {
									mymask.blitAndNeg(s, rvarmask);
								}
							}
						}
					}
//					implmask<P>& m = have_implmask(var<P>::expmasks_, v);
//					m.resizing_set(s, true);
				}
			}
			void init(int id) {
				const std::vector<rel_entry<P> >& re = rel_.entries();
				// note this may not be the correct way to adjust priories
				for (size_t i = 0; i < re.size(); ++i) {
					const rel_entry<P>& e = re[i];
					var<P>& v = *e.var_;
					bool state = rel_.varState(i, state_);
					double stateP = state?v.prioriP():1-v.prioriP();
					stateP = (stateP - var<P>::prioriP_) / (1 - var<P>::prioriP_);
					v.setPrioriP(state?stateP:1-stateP);
				}

				var<P>::set_id(id);
				add_deps(cvec<P>(), *this);
				fill_rootmasks();
				fill_expmasks();
			}
			void add_deps(const cvec<P>& shift, const reexp::exp<P>& at) {
				const std::vector<rel_entry<P> >& re = at.rel_.entries();
				for (size_t i = 0; i < re.size(); ++i) {
					const rel_entry<P>& e = re[i];
					cvec<P> sh = shift - e.shift_;
					bool s = at.rel_.varState(i, at.state_);
					e.var_->add_dep(sh, *this);
					// add also the states
					reexp::exp<P>* ex = dynamic_cast<reexp::exp<P>*>(e.var_);
					if (s && ex) add_deps(sh, *ex); // go recursive
				}
			}
			virtual int hash() const {
				return rel_.hash() + state_;
			}

			static inline bool is_overlap(const reexp::var<P>& v1,
								   	   	  const cvec<P>& shift,
								   	   	  const reexp::var<P>& v2) {
				for (const implmask<P>& im : v1.rootmasks()) {
					for (const implmask<P>& jm : v2.rootmasks()) {
						if (im.var_ == jm.var_) {
							if (im.mask_.overlap(shift, jm.mask_)) {
								return true;
							}
							break;
						}
					}
				}
				return false;
			}

			static bool is_overlap(const std::vector<reexp::rel_entry<P> >& rvars) {
				for (size_t i = 1; i < rvars.size(); ++i) {
					const reexp::rel_entry<P>& ie( rvars[i] );
					for (size_t j = 0; j < i; j++) {
						const reexp::rel_entry<P>& je( rvars[j] );
						reexp::cvec<P> delta(je.shift_);
						delta -= ie.shift_;
						if (is_overlap(*ie.var_, delta, *je.var_)) {
							return true;
						}
					}
				}
				return false;
			}

			static bool is_overlap(bitmatrix<P>& matrix,
								   std::vector<int>& vars,
								   const std::vector<reexp::rel_entry<P> >& rvars) {
				std::fill(vars.begin(), vars.end(), 0);
				for (size_t i = 0; i < rvars.size(); i++) {
					const reexp::rel_entry<P>& ie( rvars[i] );
					for (const implmask<P>& o : ie.var_->rootmasks()) {
						matrix.unite(ie.shift_, o.mask_.ndim_.dim_);
						vars[o.var_->id()]++;
					}
				}
				for (size_t v = 0; v < vars.size(); ++v) {
					matrix.bits_.fill(false);
					for (size_t i = 0; vars[v] && i < rvars.size(); i++) {
						const reexp::rel_entry<P>& ie( rvars[i] );
						for (const implmask<P>& o : ie.var_->rootmasks()) {
							if (o.var_->id() == v) {
								if (matrix.blit(ie.shift_, o.mask_)) {
									return true;
								}
								vars[v]--;
								break;
							}
						}
					}
				}
				return false;
			}

			void gen_rels(lang<P>& lang) {
				const std::vector<reexp::rel_entry<P> >& vars( rel_.entries() );

				typedef std::unordered_set<reexp::rel<P>*,
										   typename rel_ptr_hash<P>::func_type> rel_set;
				rel_set	rs(10, &rel_ptr_hash<P>::hash);
				for (auto i = vars.begin(); i != vars.end(); i++) {
					const std::vector<reexp::rel<P>*>& v( i->var_->rels() );
					for (auto r = v.begin(); r != v.end(); r++) {
						if (!(*r)->disabled()) rs.insert(*r);
					}
				}
				/*
				 * So what we are doing is that we are generating combinations
				 * from the expression and relations related to expression
				 * components.
				 */

				/*
				 * Big loop goes through the relations.
				 */
				for (auto i = rs.begin(); i != rs.end(); i++) {
					reexp::rel<P>& r = **i;
					std::vector<reexp::rel_entry<P> > rvars( r.entries() );
					std::vector<int> up; // we use this to mark current combination
					up.resize(rvars.size());

					/*
					 * This loop goes through combinations that can
					 * be generated from relation r.
					 */
					while (true) {
						/*
						 * So we seek to find from associated relation's variables an entry, which
						 * is included within this expression. E.g. if e(x) = <a(x), b(x+s)> and there is
						 * relation r(a(x), Z), we can try to generate r(e(x), Z).
						 *
						 * So go forth and try to find a match!
						 */

						repeat:
						for (size_t j = 0; j < rvars.size(); j++) {
							reexp::rel_entry<P>& re( rvars[j] );
							reexp::var<P>* replaced = r.entries()[j].var_;
							for (size_t k = up[j]; k < vars.size(); k++) {
								const reexp::rel_entry<P>& te( vars[k] );

								if (replaced == te.var_) {
									// relation variable matches expression variable
									// let's replace it with the expression
									re.var_ = this;
									re.shift_ = r.entries()[j].shift_ - te.shift_;

									// reset the previously replaced relation variables
									// (these may be reset with following iteration
									for (size_t l = 0; l < j; l++) {
										up[l] = 0;
										rvars[l] = r.entries()[l];
									}
									up[j] = k+1;

									if (!is_overlap(rvars)) {
										// normalization
										reexp::cvec<P> norm;
										for (int l = 0; l < P::DIM; l++) {
											norm[l] = 1000; // replace with MAX_INT
										}

										for (size_t k = 0; k < rvars.size(); k++) {
											reexp::cvec<P>& v = rvars[k].shift_;
											const reexp::ctx<P>& ctx = rvars[k].var_->ctx();
											for (int l = 0; l < P::DIM; l++) {
												if (ctx.v_[l] >= 0) {
													norm[l] = std::min(norm[l], v[l]);
												}
											}
										}
										for (size_t k = 0; k < rvars.size(); k++) {
											reexp::cvec<P>& v = rvars[k].shift_;
											const reexp::ctx<P>& ctx = rvars[k].var_->ctx();
											for (int l = 0; l < P::DIM; l++) {
												if (ctx.v_[l] >= 0) {
													v[l] -= norm[l];
												} else {
													v[l] = 0;
												}
											}
										}

										// add redundancy calculations later
										reexp::rel<P>& nr(lang.alloc_rel(r.ctx()));
										for (size_t k = 0; k < rvars.size(); k++) {
											nr.add_var(rvars[k].shift_, *rvars[k].var_);
										}
										if (lang.is_rel_back_unique()) {
											lang.rel_done();
										} else {
											lang.waste_rel();
										}
									}
									goto repeat;
								}
							}
						}
						break;
					}
				}
			}
			cvec<P> ctx_dim(const data<P>& data) const {
				return rel_.ctx().dim(data.dim());
			}
			/** offset of pixel as compared to exp data offset*/
			int offset(const cvec<P>& dim, const cvec<P>& at) const {
				return rel_.ctx().dim(dim).offset(at);
			}
			int state() const {
				return state_;
			}
			void add_rel(reexp::rel<P>* rel) {
				if (std::find(rels_.begin(), rels_.end(), rel) == rels_.end()) {
					rels_.push_back(rel);
				}
			}
			const reexp::rel<P>& rel() const {
				return rel_;
			}
			const std::vector<reexp::rel<P>*>& rels() const {
				return rels_;
			}

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
			lang() : origs_(), exps_(), rels_(), obs_() {}
			lang(const lang& l) : origs_(), exps_(), rels_(), obs_() {
				if (l.exp_count()) {
					throw std::runtime_error("constructing lang with new exps not yet supported");
				}
				for (int i = 0; i < l.orig_count(); ++i) {
					add_orig(reexp::orig<P>(l.orig(i).ctx()));
					origs_.back().setPrioriP(l.orig(i).prioriP());
				}
				for (int i = 0; i < l.rel_count(); ++i) {
					const reexp::rel<P>& orig = l.rel(i);
					reexp::rel<P>& rl(alloc_rel(orig.ctx()));
					for (const rel_entry<P>& e : orig.entries()) {
						rl.add_var(e.shift_, var(e.var_->id()));
					}
					rel_done();
				}
			}
			const reexp::orig<P>& orig(int i) const {
				assert(i >= 0 && i < origs_.size());
				return origs_[i];
			}
			const reexp::exp<P>& exp(int i) const {
				assert(i >= 0 && i < exps_.size());
				return exps_[i];
			}
			const reexp::rel<P>& rel(int i) const {
				assert(i >= 0 && i < rels_.size());
				return rels_[i];
			}
			void disable_rel(int i) {
				rels_[i].disable();
			}
			int enabled_rel_count() const {
				int rv = 0;
				for (int i = 0; i < rels_.size(); ++i) {
					rv += !rels_[i].disabled();
				}
				return rv;
			}
			int rel_count() const {
				return rels_.size();
			}
			reexp::var<P>& var(int i) {
				assert(i >= 0 && i < var_count());
				if (i < origs_.size()) {
					return origs_[i];
				}
				i -= origs_.size();
				return exps_[i];
			}
			const reexp::var<P>& var(int i) const {
				assert(i >= 0 && i < var_count());
				if (i < origs_.size()) {
					return origs_[i];
				}
				i -= origs_.size();
				return exps_[i];
			}
			int orig_count() const {
				return origs_.size();
			}
			int exp_count() const {
				return exps_.size();
			}
			int var_count() const {
				return origs_.size() + exps_.size();
			}
			const reexp::exp<P>& exp_back() const {
				return exps_.back();
			}
			void add_orig(const reexp::orig<P>& o) {
				int id = origs_.size();
				origs_.push_back(o);
				origs_.back().set_id(id);
				origs_.back().init();
				if (obs_) obs_->var_added(origs_.back());
			}
			reexp::rel<P>& alloc_rel(const reexp::ctx<P>& ctx) {
				rels_.push_back(reexp::rel<P>(rels_.size(), ctx));
				return rels_.back();
			}
			bool is_rel_back_unique() {
				reexp::rel<P>& rel = rels_.back(); // is this duplicate or unique?
				for (int i = 0; i < rels_.size()-1; ++i) {
					if (rel == rels_[i]) return false;
				}
				return true;
			}
			void waste_rel() {
//				rels_.back().disable();
				rels_.pop_back();
			}
			void rel_done() {
				rels_.back().finalize();
				if (obs_) obs_->rel_added(rels_.back());
			}
			void add_exp(const reexp::rel<P>& rel, int state) {
				int id = origs_.size() + exps_.size();
				exps_.push_back(reexp::exp<P>(rel, state));
				exps_.back().init(id);
				if (obs_) obs_->var_added(exps_.back());
				exps_.back().gen_rels(*this);
			}
			void set_obs(lang_obs<P>& o) {
				obs_ = &o;
			}
			void unset_obs() {
				obs_ = 0;
			}
	};



}

#endif /* REEXP_LANG_H_ */
