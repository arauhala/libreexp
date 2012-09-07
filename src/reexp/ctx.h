/*
 * ctx.h
 *
 *  Created on: Dec 15, 2010
 *      Author: arauhala
 */

#ifndef CTX_H_
#define CTX_H_

#include <stdarg.h>
#include <algorithm>
#include <bitset>
#include <iostream>
#include <assert.h>

#include "util.h"

namespace explib {

	/**
	 * So, we basically have following kind of thing
	 *
	 *  * point (cvec)
	 *  * space (cvec called dim)
	 *     * subspace (cvec called dim)
	 *   * selection, which picks some subspace out of space (called ctx)
	 *     * Selects dimensions, e.g. two out of four
	 *     * Limits the length of dimension, e.g. by 1 or 2
	 */

	/*
	 * context vector
	 */
	template <typename P>
	struct cvec {
		int v_[P::DIM];
		void clear() {
			for (int i = 0; i < P::DIM; i++) {
				(*this)[i] = 0;
			}
		}
		static cvec<P> unit_vector() {
			cvec<P> rv;
			for (int i = 0; i < P::DIM; i++) {
				rv[i] = 1;
			}
			return rv;
		}
		inline int size() const {
			return P::DIM;
		}
		inline bool operator!=(cvec<P>& v) const {
			for (int i = 0; i < P::DIM; i++) {
				if (v_[i] != v.v_[i]) return true;
			}
			return false;
		}
		inline int& operator[](int i) {
			assert(i >= 0 && i < P::DIM);
			return v_[i];
		}
		inline int operator[](int i) const {
			assert(i >= 0 && i < P::DIM);
			return v_[i];
		}
		inline bool operator!() const {
			for (int i = 0; i < P::DIM; i++) {
				if (v_[i]) return false;
			}
			return true;
		}

		inline int volume() const {
			int rv = 1;
			for (int i = 0; i < P::DIM; i++) {
				if (v_[i] != -1) {
					rv *= v_[i];
				}
			}
			return rv;
		}

		cvec(std::initializer_list<int> list) {
			int i = 0;
			for (int v : list) {
				(*this)[i++] = v;
			}
			while (i < P::DIM) (*this)[i++] = 0;
		}
		cvec(int first, ...) {
			(*this)[0] = first;
			va_list vl;
			va_start(vl, first);
			for (int i = 1; i < P::DIM; i++) {
				(*this)[i] = va_arg(vl, int);
			}
			va_end(vl);
		}
		cvec() {
			for (int i = 0; i < P::DIM; i++) {
				(*this)[i] = 0;
			}
		}
		inline cvec<P> operator+(const cvec<P>& v) const {
			cvec<P> rv;
			for (int i = 0; i < P::DIM; i++) {
				rv.v_[i] = v_[i] + v[i];
			}
			return rv;
		}
		inline cvec<P> operator-(const cvec<P>& v) const {
			cvec<P> rv;
			for (int i = 0; i < P::DIM; i++) {
				rv.v_[i] = v_[i] - v[i];
			}
			return rv;
		}
		cvec<P>& operator+=(const cvec<P>& v) {
			for (int i = 0; i < P::DIM; i++) {
				(*this)[i] += v[i];
			}
			return *this;
		}
		cvec<P>& operator-=(const cvec<P>& v) {
			for (int i = 0; i < P::DIM; i++) {
				(*this)[i] -= v[i];
			}
			return *this;
		}
		cvec<P> operator-() const {
			cvec<P> rv;
			for (int i = 0; i < P::DIM; i++) {
				rv[i] = -(*this)[i];
			}
			return rv;
		}
		int offset(const cvec<P>& v) const {
			int rv = 0;
			int exp = 1;
			for (int i = 0; i < P::DIM; i++) {
				if ((*this)[i] >= 0) {
					rv += v[i] * exp;
					exp *= (*this)[i];
				}
			}
			return rv;
		}
		cvec<P> at(int at) const {
			cvec<P> rv;
			for (int i = 0; i < P::DIM && at; i++) {
				int v = (*this)[i];
				if (v != -1) {
					rv[i] = at % v;
					at /= v;
				}
			}
			return rv;
		}
		int rowdim() const {
			for (int i = 0; i < P::DIM; i++) {
				if ((*this)[i] > 0) return i;
			}
			return -1;
		}
		int hash() const {
			util::hash_code h;
			for (int i = 0; i < P::DIM; i++) {
				h << (*this)[i];
			}
			return h();
		}
	};

	template <typename P>
	std::ostream& operator<<(std::ostream& out,
							 const explib::cvec<P>& cvec) {
		out<<"[";
		for (int i = 0; i < P::DIM; i++) {
			if (i) out<<", ";
			out<<cvec[i];
		}
		out<<"]";
		return out;
	}


	template <typename P>
	struct dim_iterator {
		dim_iterator(const cvec<P>& dim)
		:	dim_(dim), i_(), eoi_(false) {
			for (int i = 0; i < dim.size(); i++) {
				if (dim[i] == 0) eoi_ = true;
			}
		}
		const cvec<P>& operator* () const {
			return i_;
		}
		operator bool() const {
			return !eoi_;
		}
		dim_iterator& operator++() {
			for (int i = 0; i < P::DIM; i++) {
				if (i_[i] + 1 < dim_[i]) {
					i_[i]++;
					return *this;
				} else {
					i_[i] = 0;
				}
			}
			eoi_ = true;
			return *this;
		}
		cvec<P> dim_;
		cvec<P> i_;
		bool eoi_;
	};

/*	template <typename P>
	struct window_row_idx_iterator {
		private:
			int row_begin_;
			const cvec<P>& window_;
			const cvec<P>& dim_;
			int rowdim_;
		public:
			window_row_idx_iterator(const cvec<P>& window, const cvec<P>& dim)
			: row_begin_(), window_(window), dim_(dim), rowdim_(dim.rowdim()) {}
			int row_begin() const {
				return row_begin_;
			}
			int row_length() const {
				return window_[rowdim_];
			}
			window_row_idx_iterator& operator++() {
				int div = 1;
				for (int i = rowdim_+1; i < P::DIM; ++i) {
					div *= dim_;
					int v = row_begin_ % div;
					if (v < window_[i]) {

					} else {
						row_begin -= window_[i];
					}
				}
				return *this;
			}
	};*/

	template <typename P>
	struct dim_row_iterator {
		dim_row_iterator(const cvec<P>& dim)
		:	dim_(dim), i_(), eoi_(false), rowdim_(dim.rowdim()) {
			for (int i = 0; i < dim.size(); i++) {
				if (dim[i] == 0) eoi_ = true;
			}
		}
		dim_row_iterator(const cvec<P>& dim, int rowdim)
		:	dim_(dim), i_(), eoi_(false), rowdim_(rowdim) {
			for (int i = 0; i < dim.size(); i++) {
				if (dim[i] == 0) eoi_ = true;
			}
		}

		cvec<P>& settable_dim() {
			return dim_;
		}
		cvec<P>& settable_begin() {
			return i_;
		}
		const cvec<P>& begin() const {
			return i_;
		}
		// returns row length
		int length() const {
			return dim_[rowdim_];
		}
		operator bool() const {
			return !eoi_;
		}
		dim_row_iterator& operator++() {
			if (P::DIM > 1) {
				for (int i = rowdim_+1; i < P::DIM; i++) {
					if (i_[i] + 1 < dim_[i]) {
						i_[i]++;
						return *this;
					} else {
						i_[i] = 0;
					}
				}
			}
			eoi_ = true;
			return *this;
		}
		cvec<P> dim_;
		cvec<P> i_;
		bool eoi_;
		int rowdim_;
	};

	template <typename P>
	cvec<P> max(const cvec<P>& v1, const cvec<P>& v2) {
		cvec<P> rv;
		for (int i = 0; i < P::DIM; i++) {
			rv[i] = std::max(v1[i], v2[i]);
		}
		return rv;
	}

	template <typename P>
	struct ctx {
		cvec<P> v_;
		ctx(const cvec<P>& v) : v_(v) {}
		int hash() const {
			return v_.hash();
		}
		void introduce(const ctx<P>& c, const cvec<P>& v) {
			for (int i = 0; i < P::DIM; ++i) {
				if (v_[i] >= 0 && c.v_[i] >= 0) {
					v_[i] = std::max(v_[i], v[i] + c.v_[i]);
				}
			}
		}

		/**
		 * Should be equivalent with
		 * dim(limits).offset(at)
		 */
		int offset(const cvec<P>& limits, const cvec<P>& at) const {
			int rv = 0;
			int exp = 1;
			for (int i = 0; i < P::DIM; i++) {
				if (v_[i] >= 0) {
					rv += at[i] * exp;
					exp *= limits[i]-v_[i];
				}
			}
			// assert(rv == dim(limits).offset(at));
			return rv;
		}
		cvec<P> dim(const cvec<P>& limits) const {
			cvec<P> rv(limits);
			for (int i = 0; i < P::DIM; ++i) {
				if (v_[i] >= 0) {
					if (rv[i] > v_[i]) {
						rv[i] -= v_[i];
					} else {
						rv[i] = 0;
					}
				} else {
					rv[i] = -1;
				}
			}
			return rv;
		}
	};

}

#endif /* CTX_H_ */
