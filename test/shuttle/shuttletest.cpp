/*
 * shuttletest.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: arau
 */

#if 0
#include "reexp/all.h"

#include "tester.h"
#include <stdexcept>
#include <iostream>
#include "evaluation/evaluation.h"

namespace {

}

void addshuttletest(TestRunner& runner) {
	runner.add("shuttle/setup", {"func"}, &setup_test);
}

#endif
