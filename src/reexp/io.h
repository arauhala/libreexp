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

	/**
	 * provides an indexing scheme, where it is possible to identify
	 * an arbitrary bit in an arbitrary variable in a single index.
	 * the indexing scheme is bound language/data dimensions specific.
	 * variable index / bit index pair can be freely transformed into
	 * 'unified' index and back
	 */
	template <typename P>
	class index_over_var_bits {
	private:
		int size_;
		std::vector<int> offsets_;

	public:
		index_over_var_bits(const data<P>& d);
		int from_var_bit(int varid, int bitindex) const;
		bool to_var_bit(int index, int& varid, int& bitindex) const;
		int offset(int varid) const;
		int size() const;
	};

	template <typename P>
	class io {
	private:
		const reexp::stats<P>& s_;

	public:
		io(const reexp::stats<P>& s);
		void write_states(arithmetic_bit_ostream& o, const reexp::data<P>& d);
		void read_states(arithmetic_bit_istream& in, reexp::data<P>& d);
		/*
		 * knownAt marks the moment, when certain bit was marked determined/undefined.
		 * this field is exposed, because it was very useful in testing.
		 */
		void read_states(arithmetic_bit_istream& in,
						 reexp::data<P>& d,
						 const index_over_var_bits<P>& indexspace,
						 std::vector<int>& knownAt);
	};

}

#endif /* IO_H_ */
