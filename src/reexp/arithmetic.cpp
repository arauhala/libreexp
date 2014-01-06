/*
 * arithmetic.cpp
 *
 *  Created on: Nov 19, 2013
 *      Author: arau
 */


#include "reexp/arithmetic.h"

namespace reexp {

	static const int HALF = 1<<(ONE_BITS-1); // leave sign bit unused
	static const int MASK = ~((1<<31) | ONE); // Contains only 31 bits

	bool arithmetic_bit_istream::_read(int prob) {
		// Loose select
		while (true) {
			if (offset >= 0 && size + offset <= HALF) { // right border is under cut
			} else if (offset >= HALF && size + offset <= ONE) { // left border is above cut
				offset -= HALF;
			} else {
				break;
			}
			offset<<=1;
			size<<=1;
			select=(select<<1)&MASK;
			ssize<<=1; // this should not be able to raise above 32 bit
		}

		int cut = (int)(offset + (((long)prob*(long)size) >> (ONE_BITS)));

		while (true) {
			while (true) {
				// Side select
				if (cut <= select) {
					size -= (cut - offset);
					offset = (int)cut;
					return true;
				} else if (cut >= select + ssize) {
					size = cut - offset;
					return false;
				} else {
					break;
				}
			}

			// Strict select
			bool b; in_>>b;
			ssize>>=1;
			if (b) {
				select |= ssize;
			}
		}
	}

	void arithmetic_bit_ostream::_write(int prob, bool bit) {
		// Cut
		int cut = (int)((size * (long)prob)>>ONE_BITS);
		if (cut == size || !cut) {
			throw std::runtime_error("this shouldn't happen");
		}

		// Side select
		if (bit) {
			size -= cut;
			offset += cut;
		} else {
			size = cut;
		}
		if (!size) {
			throw std::runtime_error("this shouldn't happen");
		}

		// Loose select
		while (true) {
			if (!size) {
				throw std::runtime_error("this shouldn't happen");
			}

			if (offset >= 0 && size + offset <= HALF) { // right border is under cut
				out_<<false;
			} else if (offset >= HALF && size + offset <= ONE) { // left border is above cut
				offset -= HALF;
				out_<<true;
			} else {
				break;
			}
			size <<= 1;
			offset <<= 1;
		}
	}

	void arithmetic_bit_ostream::finish() {
		// strict selecting
		while (true) {
			if (offset + size < ONE) {
				out_<<false;
			} else if (offset > 0) {
				offset -= HALF;
				out_<<true;
			} else if (offset <= 0 && offset + size >= ONE) {
				break;
			} else {
				throw new std::runtime_error("should not happen");
			}
			size <<= 1;
			offset <<= 1;
		}
	}
}
