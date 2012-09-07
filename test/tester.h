/*
 * testsuite.h
 *
 *  Created on: Jul 5, 2011
 *      Author: arau
 */

#ifndef TESTER_H
#define TESTER_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <time.h>
#include <sys/time.h>

/*struct checkl_t {}; // end of line
extern checkl_t checkl;
struct ignorel_t {}; // end of line
extern ignorel_t ignorel;
struct warnl_t {}; // end of line
extern warnl_t warnl;
struct errorl_t {}; // end of line
extern errorl_t errorl;*/

enum linemod_t {
	checkl = 0,
	ignorel,
	reportl,
	faill,
};

/*class TimeSentry {
private:
	clock_t begin_;

public:
	TimeSentry()
	: begin_(clock()) {
	}

	size_t ms() {
		size_t end = clock();
		return ((end-begin_)*1000) / CLOCKS_PER_SEC;
	}
	size_t us() {
		size_t end = clock();
		return ((end-begin_)*1000000) / CLOCKS_PER_SEC;
	}
};*/

class TimeSentry {
private:
	timeval begin_;

public:
	TimeSentry() : begin_() {
		reset();
	}
	void reset() {
		gettimeofday(&begin_, 0);
	}

	size_t us() {
		timeval end;
		gettimeofday(&end, 0);
		long secdiff = end.tv_sec - begin_.tv_sec;
		long usecdiff = end.tv_usec - begin_.tv_usec;
		usecdiff += 1000000 * secdiff;
		return usecdiff;
	}
	size_t ms() {
		return us() / 1000;
	}
};

typedef std::pair<std::set<std::string>, double> record_entry;
/*
class IReporter<T> {
	public:
		virtual void feed(record_entry& e) = 0;
		virtual T report() = 0;
};*/

struct Average {
	typedef double return_type;
	double sum_;
	int n_;
	Average() : sum_(), n_() {}
	inline void feed(const record_entry& e) {
		sum_ += e.second;
		n_++;
	}
	inline double report() {
		return sum_ / n_;
	}
};

template <typename T, int N>
struct Multiplied {
	typedef double return_type;
	T t_;
	Multiplied(const T& t = T()) : t_(t) {}
	inline void feed(const record_entry& e) {
		t_.feed(e);
	}
	inline double report() {
		return t_.report() * N;
	}
};

template <typename T, typename RV = typename T::type>
struct Filtered {
	typedef RV return_type;
	const std::string& str_;
	T v_;
	Filtered(const std::string& str, const T& v)
	: str_(str), v_(v) {}
	inline void feed(const record_entry& e) {
		if (e.first.count(str_)) {
			v_.feed(e);
		}
	}
	inline RV report() {
		return v_.report();
	}
};

class TestTool;

class ReportOutput {
	private:
		TestTool& tool_;
		linemod_t mod_;
	public:
		ReportOutput(TestTool& tool, linemod_t mod);
		~ReportOutput();

		template <typename T>
		ReportOutput& operator<<(const T& c);

};

class TestTool {

	private:

		friend class ReportOutput;

		std::string expfile_;
		std::string outfile_;
		std::string recexpfile_;
		std::string recoutfile_;
		std::ifstream exp_;
		std::ofstream out_;
		std::ofstream rec_; // record file
		std::ostringstream line_;

		int linenumber_;
		bool verbose_;
		bool cleanline_;
		bool frozen_;
		bool fail_;
		bool& ok_;
		TimeSentry time_;

		std::vector<record_entry> records_;

	public:

		TestTool(const std::string& test, bool& ok, bool verbose = false);
		~TestTool();

		TestTool& operator<<(linemod_t mod);

		void record(const std::set<std::string>& tags, double value);
		template <typename R, typename T = typename R::return_type>
		T report(R r) const {
			for (const record_entry& e : records_) {
				r.feed(e);
			}
			return r.report();
		}

		ReportOutput ignored() {
			return ReportOutput(*this, ignorel);
		}

		ReportOutput reported() {
			return ReportOutput(*this, reportl);
		}

		ReportOutput failed() {
			return ReportOutput(*this, faill);
		}

		template <typename T>
		TestTool& operator<<(const T& c) {
			std::ostringstream buf;
			buf<<c;
			std::string s( buf.str() );
			for (size_t i = 0; i < s.size(); ++i) {
				if (s[i] == '\n') {
					(*this)<<checkl;
				} else {
					line_<<s[i];
				}
			}
			return *this;
		}

};

template <typename T>
ReportOutput& ReportOutput::operator<<(const T& c) {
	std::ostringstream buf;
	buf<<c;
	std::string s( buf.str() );
	for (size_t i = 0; i < s.size(); ++i) {
		if (s[i] == '\n') {
			tool_<<mod_;
		} else {
			tool_.line_<<s[i];
		}
	}
	return *this;
}

typedef void (*testfunc)(TestTool& t);

struct TestEntry {
	std::string name_;
	std::set<std::string> tags_;
	testfunc func_;
	TestEntry(const char* name,
			  const std::set<std::string>& tags,
			  testfunc func)
	: name_(name),
	  tags_(tags),
	  func_(func) {}
};

class TestRunner {
	public:
		TestRunner();
		void add(const char* name, const std::set<std::string>& tags, testfunc func);
		template <typename F>
		void map(const std::set<std::string>& keys, F f) {
			for (const TestEntry& e : tests_) {
				bool run = true;
				for (std::string key : keys) {
					if (e.tags_.count(key) == 0) {
						run = false;
						break;
					}
				}
				if (run) f(e);
			}
		}
		bool run(const std::set<std::string>& keys, bool verbose = false);
		std::string expfilepath(const std::string& testcase);
	private:
		std::vector<TestEntry> tests_;
};

// mostly syntatctic sugar, fast way for constructing string
class sup {
private:
	std::ostringstream str_;
public:
	inline sup() : str_() {}
	template <typename T>
	inline sup& operator<<(T t) {
		str_<<t;
		return *this;
	}
	inline operator std::string () {
		return str_.str();
	}
};

// helps to create tag sets
inline std::set<std::string> operator+(const std::set<std::string>& set, const std::string& t) {
	std::set<std::string> rv;
	for (auto m : set) rv.insert(m);
	rv.insert(t);
	return rv;
}

class Table;

struct TableFormat {
	const Table& t_;
	int precision_;
	TableFormat(const Table& t, int precision) : t_(t), precision_(precision) {}
};

class Table {
	private:
		std::vector<std::string> xlabels_;
		std::vector<std::string> ylabels_;
		std::string yprefix_;
		std::vector<double> values_;
	public:
		Table(const std::vector<std::string>& xlabels,
			  const std::vector<std::string>& ylabels,
			  const std::string& yprefix)
		: 	xlabels_(xlabels),
		  	ylabels_(ylabels),
		  	yprefix_(yprefix),
		  	values_() {
			values_.resize(xlabels_.size() * ylabels_.size());
		}
		const std::vector<std::string>& xlabels() const {
			return xlabels_;
		}
		const std::vector<std::string>& ylabels() const {
			return ylabels_;
		}
		Table operator-(const Table& table) const {
			if (xlabels_ != table.xlabels_
	  		 || ylabels_ != table.ylabels_
	  		 || yprefix_ != table.yprefix_) {
				throw std::runtime_error("incompatible tables for operator-");
			}
			Table rv(xlabels_, ylabels_, yprefix_);
			rv.values_ = values_;
			for (size_t i = 0; i < values_.size(); ++i) {
				rv.values_[i] -= table.values_[i];
			}
			return rv;
		}
		double& at(int x, int y) {
			return values_[y * xlabels_.size() + x];
		}
		const double& at(int x, int y) const {
			return values_[y * xlabels_.size() + x];
		}
		double& at(const std::string& x, const std::string& y) {
			size_t xi = 0, yi = 0;
			for( ;xi < xlabels_.size(); xi++) {
				if (xlabels_[xi] == x) break;
			}
			for( ;yi < ylabels_.size(); yi++) {
				if (ylabels_[yi] == y) break;
			}
			return at(xi, yi);
		}
		TableFormat formatted(int precision) const {
			return TableFormat(*this, precision);
		}
		const Table& format(std::ostream& o, int precision = 3) const {
			std::ostringstream buf;
			buf.setf(std::ios::fixed,std::ios::floatfield);
			buf.precision(precision);
			buf.setf(std::ios::left, std::ios::adjustfield);

			buf.width(12);
			buf<<yprefix_;
			for (size_t x = 0; x < xlabels_.size(); ++x) {
				buf.width(12);
				buf<<xlabels_[x];
			}
			buf<<std::endl;
			for (size_t y = 0; y < ylabels_.size(); ++y) {
				buf.width(12);
				buf<<ylabels_[y];
				for (size_t x = 0; x < xlabels_.size(); ++x) {
					buf.width(12);
					buf<<at(x, y);
				}
				buf<<std::endl;
			}
			o<<buf.str();
			return *this;
		}
		const Table& operator>>(std::ostream& o) const {
			return format(o);
		}
		std::string toplot(int xscale, int height, double minbase = 0) const {
			std::ostringstream buf;
			buf.setf(std::ios::fixed,std::ios::floatfield);
			buf.setf(std::ios::left, std::ios::adjustfield);

			double max = 0, min = minbase;
			for (size_t y = 0; y < ylabels_.size(); ++y) {
				for (size_t x = 0; x < xlabels_.size(); ++x) {
					max = std::max(max, at(x, y));
					min = std::min(min, at(x, y));
				}
			}
			double vheight = (max - min)*1.1;

			buf.precision(1);
			size_t lsize = 4, l = 0;
/*			for (size_t i = 0; i < ylabels_.size(); i++) {
				lsize = std::max(lsize, ylabels_[i].size());
			}
			lsize = std::min(size_t(5), lsize); lsize = std::max(size_t(xscale), lsize);*/
			for (size_t i = 0; i < ylabels_.size(); i++) {
				size_t begin = xscale*i;
				while (l < begin) {l++; buf<<' '; }
				if (l == begin) {
					int sz = std::min(size_t(lsize), ylabels_[i].size());
					if (l + sz > xscale*ylabels_.size()) break;
					buf<<ylabels_[i];
					buf<<' ';
					l += sz+1;
				}
			}
			buf<<'\n';

			buf.precision(3);
			for (int y = height -1; y >= 0; y--) {
				for (size_t i = 0; i < ylabels_.size(); i++) {
					char c = '.';
					for (size_t j = 0; j < xlabels_.size(); j++) {
						double value = at(j, i);
						double relat = (value - min) / vheight;
						int discat = int(relat * height);
						if (discat == y) {
							c = char(j + 'A');
						}
					}
					for (int j = 0; j < xscale; ++j) buf<<c;
				}
				buf<<"  ";
				buf<<(((double(y))*vheight)/height+min);
				buf<<'\n';
			}
			buf<<"\n";
			for (size_t j = 0; j < xlabels_.size(); j++) {
				buf<<char(j+'A')<<"   "<<xlabels_[j]<<"\n";
			}


			return buf.str();
		}
};

inline std::ostream& operator<<(std::ostream& o, const Table& t) {
	t>>o;
	return o;
}

inline std::ostream& operator<<(std::ostream& o, const TableFormat& t) {
	t.t_.format(o, t.precision_);
	return o;
}

inline bool subset(const std::set<std::string>& set,
				   const std::set<std::string>& subset) {
	for (auto s : subset) {
		if (set.count(s) == 0) return false;
	}
	return true;
}

inline bool cuttail(const std::string& path,
					const std::string& head,
					std::string& tail) {
	if (path.substr(0, head.size()) == head) {
		tail = path.substr(head.size());
		return true;
	}
	return false;
}

inline bool cuttail(const std::set<std::string>& paths,
					const std::string& head,
					std::string& tail) {
	for (auto p : paths) {
		if (cuttail(p, head, tail)) return true;
	}
	return false;
}

template <typename T>
class ToTable {
		typedef Table return_type;
	private:
		std::set<std::string> filter_;
		std::string xprefix_;
		std::string yprefix_;
		std::map<std::set<std::string>, T> entries_;
		std::vector<std::string> xlabels_;
		std::vector<std::string> ylabels_;
	public:
		ToTable(const std::set<std::string>& filter,
				const std::string& xprefix,
				const std::string& yprefix)
		: filter_(filter), xprefix_(xprefix), yprefix_(yprefix), entries_(), xlabels_(), ylabels_() {}
		void feed(const record_entry& e) {
			if (subset(e.first, filter_)) {
				std::string xlabel;
				std::string ylabel;
				if (cuttail(e.first, xprefix_, xlabel)) {
					if (cuttail(e.first, yprefix_, ylabel)) {
						std::set<std::string> key = {xprefix_+xlabel, yprefix_+ylabel};
						entries_[key].feed(e);
						if (std::find(xlabels_.begin(), xlabels_.end(), xlabel)
						    == xlabels_.end()) {
							xlabels_.push_back(xlabel);
						}
						if (std::find(ylabels_.begin(), ylabels_.end(), ylabel)
						    == ylabels_.end()) {
							ylabels_.push_back(ylabel);
						}
					}
				}
			}
		}
		Table report() {
			Table t(xlabels_, ylabels_, yprefix_);
			for (auto x : t.xlabels()) {
				std::string fx = xprefix_;
				fx += x;
				for (auto y : t.ylabels()) {
					std::string fy = yprefix_;
					fy += y;
					for (auto r : entries_) {
						if (r.first.count(fx) && r.first.count(fy)) {
							t.at(x, y) = r.second.report();
							break;
						}
					}
				}
			}
			return t;
		}


};


class Plot {
	public:
		Plot();
		~Plot();
	private:

};


#endif
