/*
 * util.cpp
 *  Created on: Oct 11, 2011
 *
 *      Author: arau
 */

#include "util.h"

namespace reexp {
	namespace util {

		int next_version_ = 0;

		int next_version() {
			return ++next_version_;
		}
	}
}
