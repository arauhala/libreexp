/*
 * arithmetic.h
 *
 *  Created on: Nov 18, 2013
 *      Author: arau
 */

#ifndef ARITHMETIC_H_
#define ARITHMETIC_H_

#include "reexp/bits.h"

// my very own arithmetic encoding implementation, may not be very optimal
// uses only 30 bits, because Java signed bit and need for overflow bit

namespace reexp {

	static const int ONE_BITS = 30; // leave sign bit unused
	static const int ONE = 1<<ONE_BITS; // leave sign bit unused

	// converts float to integer probability representation
	static int as_int_p(double number) {
		return (int)(ONE * (1-number));
	}

	class arithmetic_bit_istream {
	private:
		bit_istream& in_;
		int offset = 0;
		int size = ONE;
		int select = 0;
		int ssize = ONE;
	public:
		inline arithmetic_bit_istream(bit_istream& in) : in_(in), offset(0), size(ONE), select(0), ssize(ONE) {}
		bool _read(int prob);
		inline bool read(double prob) {
			return _read(as_int_p(prob));
		}

	};

	class arithmetic_bit_ostream {
		private:
			bit_ostream& out_;
			long size = ONE;
			int offset = 0;

		public:
			inline arithmetic_bit_ostream(bit_ostream& o) : out_(o), size(ONE), offset(0) {}
			void _write(int prob, bool bit);
			inline void write(double prob, bool bit) {
				_write(as_int_p(prob), bit);
			}
			void finish();
	};
}


#endif /* ARITHMETIC_H_ */
