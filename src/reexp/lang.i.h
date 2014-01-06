/*
 * lang.i.h
 *
 *  Created on: Jan 6, 2014
 *      Author: arau
 */

#ifndef LANG_I_H_
#define LANG_I_H_

#include "lang.h"

namespace reexp {

	template <typename P>
	bool ndim<P>::fit(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
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

	void bitref_assign_and(bits_ref dest, const_bits_ref src) {
		dest &= src;
	}

	void bitref_assign_and_neg(bits_ref dest, const_bits_ref src) {
		dest.andNeg(src);
	}

	template <typename P, typename _rowop>
	void generic_blit(const cvec<P>& destdim, const cvec<P>& destshift, reexp::bits& dest,
					  const cvec<P>& srcdim, const cvec<P>& srcshift, const reexp::bits& src,
					  const cvec<P>& clip, _rowop op) {
		// I'm not convinced this is the fastest possible implementation, especially with short rows
		// most of the time (always?) srcidx and destdim are incremented by constant deltas per iteration
		for (dim_row_iterator<P> i(clip); i; ++i) {
			const cvec<P>& at = i.begin();
			int srcidx = srcdim.offset(at + srcshift);
			int destidx = destdim.offset(at + destshift);
			int len = i.length();

			op(dest.from(destidx, len), src.from(srcidx, len));
		}
	}


	template <typename P>
	void bitmatrix<P>::unite(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
		ndim_.fit(shift, dim);
		if (ndim_.dim_.volume() > bits_.size()) {
			bits_.resize(ndim_.dim_.volume());
		}
	}
	template <typename P>
	bool bitmatrix<P>::fit(const reexp::cvec<P>& shift, const reexp::cvec<P>& dim) {
		reexp::ndim<P> ndim(ndim_);
		bool resize = ndim.fit(shift, dim);
		if (resize) {
			bitmatrix<P> resized(ndim);
			resized.blit(cvec<P>(), *this);
			*this = resized;
		}
		return resize;
	}
	template <typename P>
	void bitmatrix<P>::resizing_blit(const reexp::cvec<P>& shift, const bitmatrix& v) {
		fit(shift + v.ndim_.shift_, v.ndim_.dim_);
		blit(shift, v);
	}
	template <typename P>
	void bitmatrix<P>::resizing_set(const cvec<P>& shift, bool v) {
		fit(shift, cvec<P>::unit_vector());
		int offset = ndim_.dim_.offset(shift - ndim_.shift_);
		bits_[offset] = v;
	}
	template <typename P>
	void bitmatrix<P>::blit(const cvec<P>& shift, const bitmatrix& v) {
		cvec<P> s = shift + v.ndim_.shift_ - ndim_.shift_ ;
		for (dim_row_iterator<P> i(v.ndim_.dim_); i; ++i) {
			const cvec<P>& b = i.begin();
			int bidx = v.ndim_.dim_.offset(b);
			int len = i.length();
			int targetIdx = ndim_.dim_.offset(b + s);

			bits_.from(targetIdx, len) |= v.bits_.from(bidx, len);
		}
	}
	template <typename P>
	void bitmatrix<P>::blitAnd(const cvec<P>& shift, const bitmatrix& v) {
		cvec<P> s = shift + v.ndim_.shift_ - ndim_.shift_;

		for (dim_row_iterator<P> i(v.ndim_.dim_); i; ++i) {
			const cvec<P>& b = i.begin();
			int bidx = v.ndim_.dim_.offset(b);
			int targetIdx = ndim_.dim_.offset(b + s);
			int len = i.length();

			bits_.from(targetIdx, len) &= (v.bits_.from(bidx, len));
		}
	}

	template <typename P>
	void bitmatrix<P>::blitAndNeg(const cvec<P>& shift, const bitmatrix& v) {
		cvec<P> s = shift + v.ndim_.shift_ - ndim_.shift_;
		for (dim_row_iterator<P> i(v.ndim_.dim_); i; ++i) {
			const cvec<P>& b = i.begin();
			int bidx = v.ndim_.dim_.offset(b);
			int targetIdx = ndim_.dim_.offset(b + s);
			int len = i.length();

			bits_.from(targetIdx, len).andNeg(v.bits_.from(bidx, len));
		}
	}

	template <typename P>
	void bitmatrix<P>::blitAndNegBefore(const cvec<P>& shift, const bitmatrix& v, const cvec<P>& before) {
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

	template <typename P>
	bool bitmatrix<P>::overlap(const cvec<P>& shift, const bitmatrix& v) const {

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


	// utility function

	template <typename P>
	implmask<P>& var<P>::have_implmask(std::vector<implmask<P> >& vm,
									  const reexp::var<P>& v) {
		size_t found = 0;
		for (; found < vm.size(); ++found) {
			if (vm[found].var_ == &v) break;
		}
		if (found == vm.size()) vm.push_back(implmask<P>(&v));
		return vm[found];
	}

	template <typename P>
	var<P>::var(const reexp::ctx<P>& ctx, double prioriP)
			: id_(),
			  ctx_(ctx),
			  rels_(),
			  deps_(),
			  rootmasks_(),
			  expmasks_(),
			  disabled_(false),
			  prioriP_(prioriP) {}

	template <typename P>
	var<P>::~var() {}

	template <typename P>
	int var<P>::hash() const{
		return id_;
	}
	template <typename P>
	double var<P>::prioriP() const {
		return prioriP_;
	}
	template <typename P>
	void var<P>::setPrioriP(double prioriP) {
		prioriP_ = prioriP;
	}
	template <typename P>
	bool var<P>::disabled() const {
		return disabled_;
	}
	template <typename P>
	void var<P>::setDisabled(bool disabled) {
		disabled_ = disabled;
	}
	template <typename P>
	const std::vector<implmask<P> >& var<P>::rootmasks() const {
		return rootmasks_;
	}
	template <typename P>
	const std::vector<implmask<P> >& var<P>::expmasks() const {
		return expmasks_;
	}
	template <typename P>
	const reexp::ctx<P>& var<P>::ctx() const {
		return ctx_;
	}
	template <typename P>
	const std::vector<rel<P>*>& var<P>::rels() const {
		return rels_;
	}
	template <typename P>
	void var<P>::add_rel(rel<P>* rel) {
		if (std::find(rels_.begin(), rels_.end(), rel) == rels_.end()) {
			rels_.push_back(rel);
		}
	}
	template <typename P>
	void var<P>::remove_rel(rel<P>* rel) {
		if (std::find(rels_.begin(), rels_.end(), rel) == rels_.end()) {
			rels_.push_back(rel);
		}
	}
	template <typename P>
	void var<P>::set_id(int id) {
		id_ = id;
	}
	template <typename P>
	size_t var<P>::id() const {
		return id_;
	}
	/**
	 * In principle, the dependencies should be filled.
	 */
	template <typename P>
	void var<P>::add_dep(const reexp::cvec<P>& shift,
				 reexp::var<P>& v) {

		deps_.push_back(reexp::dep_var<P>(shift, v));
		// here goes the magic
	}
	template <typename P>
	const std::vector<dep_var<P>>& var<P>::deps() const {
		return deps_;
	}

	template <typename P>
	void orig<P>::init() {
		var<P>::rootmasks_.push_back(implmask<P>(this));
		var<P>::rootmasks_.back().resizing_set(reexp::cvec<P>(), true);
	}


	template <typename P>
	bool operator < (const reexp::rel_entry<P>& e1,
					 const reexp::rel_entry<P>& e2) {
		if (e1.var_->id() < e2.var_->id()) return true; // e1 is less than e2
		if (e1.var_->id() > e2.var_->id()) return false; // e1 is not less than e2
		for (int i = P::DIM-1; i >= 0; --i) {
			if (e1.shift_[i] < e2.shift_[i]) return true; // e1 bit is before e2 bit
			if (e1.shift_[i] > e2.shift_[i]) return false;// e2 bit is before e1 bit
		}
		return false; // equal
	}

	template <typename P>
	rel<P>::rel(int id, const reexp::ctx<P>& c)
	: id_(id), disabled_(false), hash_(), ctx_(c), e_() {
		// reorder variables so that they are in order that reflect the variables' bit order when serializing
		for (int i = 0; i < P::DIM; ++i) {
			hash_ << (ctx_.v_[i] < 0);
		}
	}
	template <typename P>
	double rel<P>::statePrioriP(int state) const {
		double p = 1;
		for (size_t i = 0; i < e_.size(); ++i) {
			double vp = e_[i].var_->prioriP();
			p *= varState(i, state)?vp:1.-vp;
		}
		return p;
	}
	template <typename P>
	void rel<P>::sort_entries() {
		std::sort(e_.begin(), e_.end());
	}
	template <typename P>
	void rel<P>::finalize() {
		sort_entries();
		for (auto i = e_.begin(); i != e_.end(); ++i) {
			i->var_->add_rel(this);
		}
	}
	template <typename P>
	void rel<P>::disable() {
		disabled_ = true;
		for (auto i = e_.begin(); i != e_.end(); ++i) {
			i->var_->remove_rel(this);
		}
	}
	template <typename P>
	bool rel<P>::disabled() const {
		return disabled_;
	}
	template <typename P>
	int rel<P>::id() const {
		return id_;
	}
	template <typename P>
	void rel<P>::add_var(const cvec<P>& shift, var<P>& var) {
		hash_ << shift.hash();
		hash_ << var.id();
		e_.push_back(rel_entry<P>(shift, var));
		ctx_.introduce(var.ctx(), shift);
	}
	template <typename P>
	const std::vector<rel_entry<P> >& rel<P>::entries() const {
		return e_;
	}
	template <typename P>
	const reexp::ctx<P>& rel<P>::ctx() const {
		return ctx_;
	}
	template <typename P>
	int rel<P>::stateCount() const {
		return 1 << e_.size();
	}
	template <typename P>
	int rel<P>::varCount() const {
		return e_.size();
	}
	template <typename P>
	bool rel<P>::varState(int var, int state) {
		return ((state>>var)&0x1);
	}
	template <typename P>
	int rel<P>::hash() const {
		return hash_();
	}
	template <typename P>
	bool rel<P>::operator==(reexp::rel<P>& r) const {
		if (hash_() != r.hash_()) return false;
		if (ctx_.v_ != r.ctx_.v_) return false;
		if (e_.size() != r.e_.size()) return false;
		for (size_t i = 0; i < e_.size(); i++) {
			if (e_[i].var_ != r.e_[i].var_
			 || e_[i].shift_ != r.e_[i].shift_) {
				return false;
			}
		}
		return true;
	}

	template <typename P>
	exp<P>::exp(const reexp::rel<P>& rel, int state)
	:	var<P>(rel.ctx(), rel.statePrioriP(state)), rel_(rel), state_(state), rels_() {
	}
	template <typename P>
	int exp<P>::min_var_idx() const {
		size_t minvar = SIZE_MAX;
		for (const rel_entry<P>& e : rel_.entries()) {
			minvar = std::min(minvar, e.var_->id());
		}
		return minvar;
	}
	// rootmasks need to be renamed root masks
	// state of rel needs to be checked, if true, use root masks
	// if false, blit directly.
	template <typename P>
	void exp<P>::fill_rootmasks() {
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
	template <typename P>
	void exp<P>::fill_expmasks() {
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
	template <typename P>
	void exp<P>::init(int id) {
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
	template <typename P>
	void exp<P>::add_deps(const cvec<P>& shift, const reexp::exp<P>& at) {
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
	template <typename P>
	int exp<P>::hash() const {
		return rel_.hash() + state_;
	}

	template <typename P>
	bool exp<P>::is_overlap(const reexp::var<P>& v1,
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

	template <typename P>
	bool exp<P>::is_overlap(const std::vector<reexp::rel_entry<P> >& rvars) {
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

	/*
	template <typename P>
	bool exp<P>::is_overlap(bitmatrix<P>& matrix,
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
	}*/

	template <typename P>
	void exp<P>::gen_rels(lang<P>& lang) {
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
								// 1. normalization of shifts
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

								// 2. create new normalized relation
								reexp::rel<P>& nr(lang.alloc_rel(r.ctx()));
								for (size_t k = 0; k < rvars.size(); k++) {
									nr.add_var(rvars[k].shift_, *rvars[k].var_);
								}
								nr.sort_entries();

								// 3. if the relation is new, keep it, otherwise scrap
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
	template <typename P>
	cvec<P> exp<P>::ctx_dim(const data<P>& data) const {
		return rel_.ctx().dim(data.dim());
	}
	/** offset of pixel as compared to exp data offset*/
	template <typename P>
	int exp<P>::offset(const cvec<P>& dim, const cvec<P>& at) const {
		return rel_.ctx().dim(dim).offset(at);
	}
	template <typename P>
	int exp<P>::state() const {
		return state_;
	}
	template <typename P>
	void exp<P>::add_rel(reexp::rel<P>* rel) {
		if (std::find(rels_.begin(), rels_.end(), rel) == rels_.end()) {
			rels_.push_back(rel);
		}
	}
	template <typename P>
	const reexp::rel<P>& exp<P>::rel() const {
		return rel_;
	}
	template <typename P>
	const std::vector<reexp::rel<P>*>& exp<P>::rels() const {
		return rels_;
	}

	template <typename P>
	lang<P>::lang() : origs_(), exps_(), rels_(), obs_() {}
	template <typename P>
	lang<P>::lang(const lang& l) : origs_(), exps_(), rels_(), obs_() {
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
	template <typename P>
	const reexp::orig<P>& lang<P>::orig(int i) const {
		assert(i >= 0 && i < origs_.size());
		return origs_[i];
	}
	template <typename P>
	const reexp::exp<P>& lang<P>::exp(int i) const {
		assert(i >= 0 && i < exps_.size());
		return exps_[i];
	}
	template <typename P>
	const reexp::rel<P>& lang<P>::rel(int i) const {
		assert(i >= 0 && i < rels_.size());
		return rels_[i];
	}
	template <typename P>
	void lang<P>::disable_rel(int i) {
		rels_[i].disable();
	}
	template <typename P>
	int lang<P>::enabled_rel_count() const {
		int rv = 0;
		for (int i = 0; i < rels_.size(); ++i) {
			rv += !rels_[i].disabled();
		}
		return rv;
	}
	template <typename P>
	int lang<P>::rel_count() const {
		return rels_.size();
	}
	template <typename P>
	reexp::var<P>& lang<P>::var(int i) {
		assert(i >= 0 && i < var_count());
		if (i < origs_.size()) {
			return origs_[i];
		}
		i -= origs_.size();
		return exps_[i];
	}
	template <typename P>
	const reexp::var<P>& lang<P>::var(int i) const {
		assert(i >= 0 && i < var_count());
		if (i < origs_.size()) {
			return origs_[i];
		}
		i -= origs_.size();
		return exps_[i];
	}
	template <typename P>
	int lang<P>::orig_count() const {
		return origs_.size();
	}
	template <typename P>
	int lang<P>::exp_count() const {
		return exps_.size();
	}
	template <typename P>
	int lang<P>::var_count() const {
		return origs_.size() + exps_.size();
	}
	template <typename P>
	const reexp::exp<P>& lang<P>::exp_back() const {
		return exps_.back();
	}
	template <typename P>
	void lang<P>::add_orig(const reexp::orig<P>& o) {
		int id = origs_.size();
		origs_.push_back(o);
		origs_.back().set_id(id);
		origs_.back().init();
		if (obs_) obs_->var_added(origs_.back());
	}
	template <typename P>
	reexp::rel<P>& lang<P>::alloc_rel(const reexp::ctx<P>& ctx) {
		rels_.push_back(reexp::rel<P>(rels_.size(), ctx));
		return rels_.back();
	}
	template <typename P>
	bool lang<P>::is_rel_back_unique() {
		reexp::rel<P>& rel = rels_.back(); // is this duplicate or unique?
		for (int i = 0; i < rels_.size()-1; ++i) {
			if (rel == rels_[i]) return false;
		}
		return true;
	}
	template <typename P>
	void lang<P>::waste_rel() {
//				rels_.back().disable();
		rels_.pop_back();
	}
	template <typename P>
	void lang<P>::rel_done() {
		rels_.back().finalize();
		if (obs_) obs_->rel_added(rels_.back());
	}
	template <typename P>
	void lang<P>::add_exp(const reexp::rel<P>& rel, int state) {
		int id = origs_.size() + exps_.size();
		exps_.push_back(reexp::exp<P>(rel, state));
		exps_.back().init(id);
		if (obs_) obs_->var_added(exps_.back());
		exps_.back().gen_rels(*this);
	}
	template <typename P>
	void lang<P>::set_obs(lang_obs<P>& o) {
		obs_ = &o;
	}
	template <typename P>
	void lang<P>::unset_obs() {
		obs_ = 0;
	}

}


#endif /* LANG_I_H_ */
