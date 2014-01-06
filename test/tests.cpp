/*
 * test.cpp
 *
 *  Created on: Jul 3, 2011
 *      Author: arau
 */

#include <iostream>
#include "tester.h"
#include <string.h>

void addgermancredittest(TestRunner& runner);
void addtexttest(TestRunner& runner);
void addvarstest(TestRunner& runner);
void addroomstest(TestRunner& runner);
void addutilstest(TestRunner& runner);
void addoptdigitstest(TestRunner& runner);
void addoptdigits2test(TestRunner& runner);
void addsimpletest(TestRunner& runner);
void addimagetest(TestRunner& runner);
void addroomstest(TestRunner& runner);
void addstatlogbitstest(TestRunner& runner);
void adddnasplicetest(TestRunner& runner);

using namespace std;

int main(int argc, char** argv) {
	TestRunner runner("libreexp test");

	addgermancredittest(runner);
	addvarstest(runner);
	addtexttest(runner);
	addimagetest(runner);
	addroomstest(runner);
	addutilstest(runner);
	addoptdigitstest(runner);
	addoptdigits2test(runner);
	addsimpletest(runner);
	addstatlogbitstest(runner);
	adddnasplicetest(runner);
/*	adddigitstest(runner);*/

	return runner.exec(argc, argv);
}
