/*
 * io.h
 *
 *  Created on: Nov 20, 2013
 *      Author: arau
 */

#ifndef IO_H_
#define IO_H_

#include "reexp/stats.h"
#include "reexp/arithmetic.h"
#include <functional>

namespace reexp {

	template <typename P>
	class index_over_var_bits {
	private:
		int size_;
		std::vector<int> offsets_;

	public:
		index_over_var_bits(const data<P>& d) : size_(), offsets_(){
			offsets_.resize(d.var_count()+1);
			for (size_t i = 0; i < d.var_count(); i++) {
				offsets_[i] = size_;
				size_ += d.var(i).bits().size();
			}
			offsets_[d.var_count()] = size_;
		}
		int from_var_bit(int varid, int bitindex) const {
			return offset(varid) + bitindex;
		}
		bool to_var_bit(int index, int& varid, int& bitindex) const {
			if (index >= size_) return false;
			varid = 0;
			while (index >= offsets_[varid+1]) varid++;
			bitindex = index - offsets_[varid];
			return true;
		}
		int offset(int varid) const {
			return offsets_[varid];
		}
		int size() const {
			return size_;
		}
	};

	template <typename P>
	class io {
	private:
		const reexp::stats<P>& s_;

	public:
		io(const reexp::stats<P>& s) : s_(s) {}
		void write_states(arithmetic_bit_ostream& o, const reexp::data<P>& d) {
			for (int i = d.var_count()-1; i >= 0; --i) {
				const reexp::data_var<P>& v = d.var(i);
				int int_p = as_int_p(s_.var(i).eP());
				for (int j = v.bits().size()-1; j >= 0; --j) {
					if (v.bits()[j]) { // if defined, write
						o._write(int_p, *v.bits()[j]);
					}
				}
			}
		}
/*		void mark_expressed(reexp::data<P>& d,
				            const reexp::exp<P>& e,
							cvec<P> at,
							const std::vector<int>& knownOffsets,
							std::vector<int>& knownAt,
							int kat) {
			const reexp::rel<P>& r = e.rel();
			int reidx = 0;
			for (const reexp::rel_entry<P>& re : r.entries()) {
				reexp::data_var<P>& rdv = d.var(re.var_->id());
				cvec<P> rat = at + re.shift_;
				rdv[rat] = false;
				*rdv[rat] = r.varState(reidx++, e.state());
				int idx = rdv.index(rat);
				int rvknownoffset = knownOffsets[re.var_->id()];
				int rkat = rvknownoffset + idx;
				knownAt[rkat] = kat;
				const reexdp::exp<P>* rexp = dynamic_cast<const reexp::exp<P>*>(re.var_);
				if (rexp) {
					mark_expressed(d, *rexp, rat, knownOffsets, knownAt, kat);
				}
			}
		};*/


		// knownAt should mark the moment, when certain bit was marked undefined/determined
		// for the first time (smallest index) when re-expressing. note that reading
		// is done in reverse order compared to re-expressing and writing.
		//

		void read_states(arithmetic_bit_istream& in, reexp::data<P>& d) {
			index_over_var_bits<P> indexspace(d);
			std::vector<int> knownAt;
			read_states(in, d, indexspace, knownAt);
		}

		void read_states(arithmetic_bit_istream& in, reexp::data<P>& d, const index_over_var_bits<P>& indexspace, std::vector<int>& knownAt) {
			knownAt.resize(indexspace.size());
			for (int& k : knownAt) k = INT_MAX;
			for (size_t i = 0; i < d.var_count(); i++) {
				d.var(i).defined().fill(true);
				d.var(i).states().fill(false);
			}
			for (int i = d.var_count()-1; i >= 0; --i) {
				reexp::data_var<P>& dv = d.var(i);
				int vknownoffset = indexspace.offset(i);
				int int_p = as_int_p(s_.var(i).eP());
				const reexp::var<P>& v = d.lang().var(i);
				const reexp::exp<P>* vexp = dynamic_cast<const reexp::exp<P>*>(&v);
				cvec<P> vdim = dv.dim();

				std::vector<const exp<P>*> sexcl;
				for (int j = i+1; j < d.var_count(); ++j) {
					const exp<P>* e = dynamic_cast<const reexp::exp<P>*>(&d.lang().var(j));
					if (e && e->rel().entries()[0].var_ == &v) {
						sexcl.push_back(e);
					}
				}

				for (int j = dv.bits().size()-1; j >= 0; --j) {
					cvec<P> at = vdim.at(j);
					int kat = vknownoffset + j;
					for (const exp<P>* e : sexcl) {
						const rel<P>& erel = e->rel();
						cvec<P> eat = at - erel.entries()[0].shift_;
						const data_var<P>& edv = d.var(e->id());
						cvec<P> edim = edv.dim();
						if (edim.contains(eat) && !*edv[eat]) {
							bool exclude = true;
							int ekat = indexspace.from_var_bit(e->id(), edv.index(eat));
							for (int k = 1; k < erel.entries().size(); ++k) {
								const rel_entry<P>& re = erel.entries()[1];
								cvec<P> evat = eat + re.shift_;
								const data_var<P>& evdv = d.var(re.var_->id());
								exclude &= (*evdv[evat] == erel.varState(k, e->state()));
								int evkat = indexspace.from_var_bit(re.var_->id(), evdv.index(evat));
								exclude &= knownAt[evkat] > ekat;
							}
							if (exclude) {
								auto bit = dv.bits()[j];
								if (bit) {
									bit = false; // implicitly known
									*bit = !erel.varState(0, e->state());
								}
								// this bit was marked determined at the moment of ekat
								knownAt[kat] = std::min(ekat, knownAt[kat]);
							}
						}
					}
					auto b = dv.bits()[j];
					if (b) *b = in._read(int_p); // if defined read the state
					if (*b) {
						for (const implmask<P>& m : v.expmasks()) {
							int mid = m.var_->id();
							reexp::data_var<P>& mdv = d.var(mid);
							int mvknownoffset = indexspace.offset(mid);
							const bitmatrix<P>& mask = m.mask_;
							for (dim_iterator<P> di(mask.ndim_.dim_); di; ++di) {
								cvec<P> mat = *di + mask.ndim_.shift_ + at;
								if (mdv.dim().contains(mat) && m.mask_[*di]) {
									int idx = mdv.index(mat);
									int mkat = mvknownoffset + idx;
									if (mkat < kat) {
										auto mbit = mdv.bits()[idx];
										if (mbit) {
											mbit = false; // mark implicitly known
											*mbit = false; // assume it is excluded
										}
										knownAt[mkat] = kat;
									}
								}
							}
						}
						if (vexp) { // mark expressed states
//							mark_expressed(d, *vexp, at, knownOffsets, knownAt, kat);
							const reexp::rel<P>& r = vexp->rel();
							int reidx = 0;
							for (const reexp::rel_entry<P>& re : r.entries()) {
								reexp::data_var<P>& rdv = d.var(re.var_->id());
								cvec<P> rat = at + re.shift_;
								*rdv[rat] = r.varState(reidx++, vexp->state());
							}
						}
					}
				}
			}
		}

#if 0
		/**
		 * NOTE: the data need to have dimensions set before hand
		 *       the serialization does not serialize the data dimensions,
		 *       but only the binary states.
		 *
		 * NOTE: this is probably pretty slow algorithm (still depending of
		 *       the data), unfortunately
		 */
		void read_states2(arithmetic_bit_istream& in, reexp::data<P>& d) {

			// needed for reverting state exclusion rule
			//
			int totalbitcount = 0;
			std::vector<int> knownOffsets(d.var_count());

			// prepare the data
			for (size_t i = 0; i < d.var_count(); i++) {
				knownOffsets.push_back(totalbitcount);
				totalbitcount += d.var(i).bits().size();
				d.var(i).defined().fill(true);
				d.var(i).states().fill(false);
			}
			std::vector<int> knownAt(totalbitcount);
			std::vector<size_t> later;
			for (int i = d.var_count()-1; i >= 0; --i) {
				reexp::data_var<P>& v = d.var(i);
				cvec<P> vdim = v.dim();
				int vknownoffset = knownOffsets[i];
				int int_p = as_int_p(s_.var(i).eP());
				later.reset();

				// this is all & good, but it doesn't really work for expressions
				// where both variables are the same.
				for (size_t e = i+1 ; e < d.var_count(); ++e) {
					reexp::data_var<P>& ed = d.var(e);
					const reexp::exp<P>* ev =
						dynamic_cast<const reexp::exp<P>*>(&ed.var());
					if (ev) { // is an expression
						const reexp::rel<P>& er = ev->rel();
						if (&v.var_ == er.entries()[0].var_) {
							later.push_back(e);
#if 0
							bool process_later = false;
							// sanity check, we cannot rely on expression states that haven't been read yet
							for (size_t i = 1; i < er.entries().size(); ++i) {
								if (&ed.var_ == er.entries()[i].var_) {
									// there goes batch processing & bit optimizations :(
									process_later = true;
									break;
								}
							}
							if (process_later) {
								later.push_back(e);
							} else {
								// fast algorithm, which produces valid results only
								// under some strict condition
								cvec<P> nvshift = -er.entries()[0].shift_;
								reexp::bitmatrix<P> undef(ndim<P>(nvshift, ed.dim()));
								undef.bits_.fill(true); // mark all true
	/*							reexp::bitmatrix<P> expstates(ed.dim(), ed.states());
								undef.blitAndNeg(cvec<P>(), expstates); // if expstate is false, mark true*/
								generic_blit(vdim, cvec<P>(), undef.bits_,
											 ed.dim(), -nvshift, ed.states(),
											 ed.dim(), bitref_assign_and);

								for (size_t i = 1; i < er.entries().size(); ++i) {
									reexp::data_var<P>& rv = d.var(er.entries()[i].var_->id());
									const reexp::cvec<P>& rs = er.entries()[i].shift_;
									if (er.varState(i, ev->state())) {
										generic_blit(vdim, cvec<P>(), undef.bits_,
													 rv.dim(), rs - nvshift, rv.states(),
													 ed.dim(), bitref_assign_and);

	//									undef.blitAnd(rs, varstates); // if expstate is false, mark true
									} else {
										generic_blit(vdim, cvec<P>(), undef.bits_,
													 rv.dim(), rs - nvshift, rv.states(),
													 ed.dim(), bitref_assign_and_neg);
	//									undef.blitAndNeg(rs, varstates); // if expstate is false, mark true
									}
								}
								v.defined().andNeg(undef.bits_);
							}
#endif
						}
					}
				}
				// there has to be a better way
				// (maintaining a lookup for 'active' state exclusions)
				// bitmap for each expression (?)
				// consider the rule:
				//   a & b <-> e 		(a the time of e)
				//   -> e! && b -> a
				//

				for (int j = v.states().size()-1; j >= 0; --j) {
					if (later.size()) { // slow
						cvec<P> vat = vdim.at(j);

						for (size_t e : later) {
							reexp::data_var<P>& ed = d.var(e);
							const reexp::exp<P>& ev =
								dynamic_cast<const reexp::exp<P>&>(ed.var());
							const reexp::rel<P>& er = ev->rel();
							reexp::cvec<P> vshift = er.entries()[0].shift_;
							reexp::cvec<P> eat = vat - vshift;
							reexp::cvec<P> edim = ed.dim();
							int eidx = edim.offset(eat);
							int eimploff = knownOffsets[e] + eidx;
							if (edim.contains(eat) && !ed[eat] && knownAt[vknownoffset]) {
								bool exclude = true;
								for (size_t i = 1; i < er.entries().size(); ++i) {
									size_t rid = er.entries()[i].var_->id();
									reexp::data_var<P> rd = d.var(rid);
									reexp::cvec<P> rat = eat + er.entries()[i].shift_;
									int ridx = rd.index(rat);
									auto bit = rd[rat];
									exclude &= (bit && knownAt[knownOffsets[rid] + irdx] && *bit == er.varState(i, ev.state()));
								}
								if (exclude) {
									knownAt[vknownoffset + j] = eimploff;
									v.defined()[j] = false;
								}
							}
						}
					}
					auto b = v.bits()[j];
					if (b) { // if defined, read the state
						*b = in._read(int_p);
						knownAt[vknownoffset+j] = vknownoffsetoffset+j;
						// apply the state mask
						const reexp::exp<P>* e = dynamic_cast<const reexp::exp<P>*>(&v.var_);
						if (e) {
							// todo: add implicit offsets
							d.apply_state(*e, v.dim().at(j), *b);
						}
					}
				}
			}

		}
#endif
	};

}

#endif /* IO_H_ */
