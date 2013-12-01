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

namespace reexp {


	template <typename P>
	class io {
	private:
		reexp::stats<P> s_;

	public:

		void write_states(arithmetic_bit_ostream& o, const reexp::data<P>& d) {
			for (int i = d.var_count()-1; i >= 0; --i) {
				const reexp::data_var<P>& v = d.var(i);
				int int_p = as_int_probability(s_.var(i).eP());
				for (size_t j = 0; j < v.bits().size(); ++j) {
					if (v.bits()[j]) { // if defined, write
						o._write(int_p, *v.bits()[j]);
					}
				}
			}
		}

		/**
		 * NOTE: the data need to have dimensions set before hand
		 *       the serialization does not serialize the data dimensions,
		 *       but only the binary states.
		 */
		void read_states(arithmetic_bit_istream& o, reexp::data<P>& d) {
			// prepare the data
			for (size_t i = 0; i < d.var_count(); i++) {
				d.var(i).defined().fill(true);
				d.var(i).states().fill(false);
			}
			for (int i = d.var_count()-1; i >= 0; --i) {
				reexp::data_var<P>& v = d.var(i);
				int int_p = as_int_probability(s_.var(i).eP());
				for (size_t j = 0; j < v.bits().size(); ++j) {
					auto b = v.bits()[j];
					if (b) { // if defined, read the state
						*b = o._read(int_p);
						// apply the state mask
					}
				}
			}

		}



	};

}

#endif /* IO_H_ */
