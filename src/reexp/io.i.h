/*
 * io.i.h
 *
 *  Created on: Jan 16, 2014
 *      Author: arau
 */

#ifndef IO_I_H_
#define IO_I_H_

#include "reexp/io.h"
#include "reexp/stats.i.h"

namespace reexp {

	template <typename P>
	index_over_var_bits<P>::index_over_var_bits(const data<P>& d) : size_(), offsets_(){
		offsets_.resize(d.var_count()+1);
		for (size_t i = 0; i < d.var_count(); i++) {
			offsets_[i] = size_;
			size_ += d.var(i).bits().size();
		}
		offsets_[d.var_count()] = size_;
	}
	template <typename P>
	int index_over_var_bits<P>::from_var_bit(int varid, int bitindex) const {
		return offset(varid) + bitindex;
	}
	template <typename P>
	bool index_over_var_bits<P>::to_var_bit(int index, int& varid, int& bitindex) const {
		if (index >= size_) return false;
		varid = 0;
		while (index >= offsets_[varid+1]) varid++;
		bitindex = index - offsets_[varid];
		return true;
	}
	template <typename P>
	int index_over_var_bits<P>::offset(int varid) const {
		return offsets_[varid];
	}
	template <typename P>
	int index_over_var_bits<P>::size() const {
		return size_;
	}

	template <typename P>
	io<P>::io(const reexp::stats<P>& s) : s_(s) {}

	template <typename P>
	void io<P>::write_states(arithmetic_bit_ostream& o, const reexp::data<P>& d) {
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

	template <typename P>
	void io<P>::read_states(arithmetic_bit_istream& in, reexp::data<P>& d) {
		index_over_var_bits<P> indexspace(d);
		std::vector<int> knownAt;
		read_states(in, d, indexspace, knownAt);
	}

	template <typename P>
	void io<P>::read_states(arithmetic_bit_istream& in, reexp::data<P>& d, const index_over_var_bits<P>& indexspace, std::vector<int>& knownAt) {
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

}



#endif /* IO_I_H_ */
