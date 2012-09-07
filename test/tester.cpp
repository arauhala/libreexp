/*
 * tester.cpp
 *
 *  Created on: Jul 6, 2011
 *      Author: arau
 */

#include "tester.h"

//checkl_t checkl;
using namespace std;

ReportOutput::ReportOutput(TestTool& tool, linemod_t mod)
: tool_(tool), mod_(mod) {
	if (tool_.line_.str().size()) {
		tool_<<"\n";
	}
}
ReportOutput::~ReportOutput() {
	if (tool_.line_.str().size()) {
		(*this)<<"\n";
	}
}

TestTool::TestTool(const string& test, bool& ok, bool verbose)
:	expfile_(test + "_exp.txt"),
	outfile_(test + "_out.txt"),
	recexpfile_(test + "_rec_exp.txt"),
	recoutfile_(test + "_rec_out.txt"),
	exp_(expfile_),
	out_(outfile_),
	rec_(recoutfile_),
	line_(),
	linenumber_(),
	verbose_(verbose),
	cleanline_(),
	fail_(),
	ok_(ok),
	time_(),
	records_() {
	cout<<test<<".. ";
	frozen_ = exp_;
	if (!frozen_) {
		verbose_ = true;
		cout<<"test not yet frozen";
	}
	if (verbose_) { cout<<endl<<endl; cleanline_ = true; }
	std::ifstream rec(recexpfile_);
	if (rec) {
		while (rec) {
			int keys = 0;
			rec>>keys;
			record_entry e;
			for (int i = 0; i < keys; ++i) {
				std::string key;
				rec>>key;
				e.first.insert(key);
			}
			e.first.insert("run:exp");
			rec>>e.second;
			records_.push_back(e);
		}
	}
}

TestTool::~TestTool() {
	int ms = time_.ms();
	cout<<ms<<"ms.. ";
	if (verbose_) cout<<endl;

	rec_.close();
	out_.close();

	if (records_.size() == 0) remove(recoutfile_.c_str());
	if (frozen_ && !fail_) {
		cout<<"ok. "<<endl;
		remove(outfile_.c_str());
	} else {
		if (frozen_) cout<<"failed. ";
		while (true) {
			if (frozen_) cout<<"re[f]reeze";
			else cout<<" [f]reeze";
			cout<<" or view [d]iff? ";
			string answer;
			cin>>answer;
			if (answer == "f") {
				cout<<endl;
				int err = rename(outfile_.c_str(), expfile_.c_str());
				if (!err) {
					if (records_.size() > 0) {
						err = rename(recoutfile_.c_str(), recexpfile_.c_str());
					} else {
						remove(recexpfile_.c_str());
					}
					if (!err) {
						cout<<"exp file created."<<endl;
					}
				}
				if (err) {
					cout<<"failed creating exp file."<<endl;
				}
				break;
			} else if (answer == "d") {
				std::ostringstream cmd;
				cmd<<"meld ";
				cmd<<expfile_<<" ";
				cmd<<outfile_;
				int rv = system(cmd.str().c_str());
				if (rv) cout<<"no meld?\n";
				cout<<"\n";
			} else {
				break;
			}
		}
		cout<<endl;
	}
	ok_ &= !fail_;
}

TestTool& TestTool::operator<<(linemod_t mod) {
	string l = line_.str();
	line_.str("");
	linenumber_++;
	out_ << l << endl;
	string e;
	bool fail = false;
	if (mod == faill ||
	    (frozen_ && (!getline(exp_, e)
	    		 || (mod == checkl && l != e)))) {
		fail = fail_ = true;
	}

	if (verbose_ || fail || mod == reportl) {
		if (!cleanline_) { cout<<endl; cleanline_ = true; };
		cout<<"[";
		cout<<(linenumber_/1000)%10;
		cout<<(linenumber_/100)%10;
		cout<<(linenumber_/10)%10;
		cout<<(linenumber_%10);
		cout<<"] ";
		switch (mod) {
			case reportl:  cout<<"R "; break;
			case faill: cout<<"F "; break;
			default: cout<<"  ";
		}
		cout<<l<<endl;
		if (mod == checkl && fail) {
			cout<<" [!=]    "<<e<<endl;
		}
	}
	return *this;
}

void TestTool::record(const std::set<std::string>& tags, double value) {
	rec_<<tags.size()<<" ";
	for (const std::string& tag : tags) {
		rec_<<tag<<" ";
	}
	rec_<<value<<"\n";
	records_.push_back(record_entry(tags, value));
	records_.back().first.insert("run:out");
}

TestRunner::TestRunner() : tests_() {}

void TestRunner::add(const char* name, const std::set<std::string>& tags, testfunc func) {
	tests_.push_back(TestEntry(name, tags, func));
	std::string ntag( name );
	while (true) {
		tests_.back().tags_.insert(ntag);
		size_t sep = ntag.find_last_of("/");
		if (sep == std::string::npos) break;
		ntag = ntag.substr(0, sep);
	}
}

bool TestRunner::run(const std::set<std::string>& keys, bool verbose) {
	bool ok = true;
	map(keys, [&](const TestEntry& e) {
		std::string test;
		test += "test/";
		test += e.name_;
		TestTool t(test, ok, verbose);
		(*e.func_)(t);
	});
	return ok;
}

std::string TestRunner::expfilepath(const std::string& testcase) {
	return sup()<<"test/"<<testcase<<"_exp.txt";
}
