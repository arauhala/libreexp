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
void addsimpletest(TestRunner& runner);
void adddigitstest(TestRunner& runner);
void addimagetest(TestRunner& runner);
void addroomstest(TestRunner& runner);
void addstatlogbitstest(TestRunner& runner);

using namespace std;

int main(int argc, char** argv) {
	TestRunner runner;

	addgermancredittest(runner);
	addvarstest(runner);
	addtexttest(runner);
	addimagetest(runner);
	addroomstest(runner);
	addutilstest(runner);
	addoptdigitstest(runner);
	addsimpletest(runner);
	addstatlogbitstest(runner);
/*	adddigitstest(runner);*/

	bool list = false, print = false, verbose = false, path = false;
	set<string> keys;
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-l") == 0) {
			list = true;
		} else if (strcmp(argv[i], "-p") == 0) {
			print = true;
		} else if (strcmp(argv[i], "--path") == 0) {
			path = true;
		} else if (strcmp(argv[i], "-v") == 0) {
		    verbose = true;
		} else {
			keys.insert(argv[i]);
		}
	}
	if (list) {
		runner.map(keys, [](const TestEntry& t) {
			cout.width(30);
			cout.setf(_S_left, _S_left);
			cout<<t.name_<<"{";
			bool first = true;
			for (const string& k : t.tags_) {
				if (!first) cout<<", ";
				cout<<k; first = false;
			}
			cout<<"}\n";
		});
	} else if (print) {
		runner.map(keys, [&runner](const TestEntry& t) {
			cout<<t.name_<<": ";
			std::ifstream in(runner.expfilepath(t.name_));
			if (in) {
				std::cout<<"\n";
				int c;
				while (in && (c = in.get()) >= 0) {
					std::cout<<char(c);;
				}
			} else {
				cout<<"(not defined)\n";
			}
		});
	} else if (path) {
		runner.map(keys, [&runner](const TestEntry& t) {
			cout<<runner.expfilepath(t.name_)<<"\n";
		});
	} else {
		cout<<"explib test.."<<endl<<endl;

		TimeSentry ms;
		bool ok = runner.run(keys, verbose);
		cout<<endl;
		cout<<ms.ms()<<"ms.. ";

		if (ok) {
			cout<<"done.\n";
		} else {
			cout<<"failed.\n";
		}
	}
}
