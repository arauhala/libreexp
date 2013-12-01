/*
 * templates.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: arau
 */

#include "reexp/reexp.h"
#include "reexp/pred.i.h"
#include "tester.h"

namespace reexp {

	template std::vector<double> pred<traits3d>::rowP<test_tool&>(const data<traits3d>& d, int varid, test_tool& explain) const;
}
