/*
 * RecognitionTest.cpp
 *
 *  Created on: Oct 22, 2010
 *      Author: arauhala
 */

#include "reexp/all.h"

#include "tester.h"

#include "stdio.h"

#include <iostream>
#include <fstream>

namespace {

	using namespace reexp;

	typedef reexp::traits2d dnasplice_problem;

	enum {
		varid_low,
		varid_high,
		varid_count
	};

	enum {
		cvarid_sample,
		cvarid_index
	};

	void setup_lang(lang<dnasplice_problem>& l) {
		typedef dnasplice_problem p;
		for (int i = 0; i < varid_count ; ++i) {
			l.add_orig(orig<p>(cvec<p>(0,0)));
		}

		for (int i = 0; i < varid_count ; ++i) {
			for (int j = 0; j < varid_count ; ++j) {
				reexp::rel<p>& rl(l.alloc_rel(reexp::cvec<p>(0, 0)));
				rl.add_var(reexp::cvec<p>(0, 0),
						   l.var(i));
				rl.add_var(reexp::cvec<p>(0, 1),
						   l.var(j));
				l.rel_done();

				reexp::rel<p>& rl2(l.alloc_rel(reexp::cvec<p>(0, 0)));
				rl2.add_var(reexp::cvec<p>(0, 0),
						   l.var(i));
				rl2.add_var(reexp::cvec<p>(0, 2),
						   l.var(j));
				l.rel_done();
			}
		}
	}

	void setup_data(lang<dnasplice_problem>& l, data<dnasplice_problem>& d) {
		typedef dnasplice_problem p;
		std::ifstream in("test/dnasplice/splice.data");
		std::string ln;
		int s = 0; // sample
		while (std::getline(in, ln)) {
			if (ln.size()) {
				auto at = ln.begin();
				auto comma = std::find(at, ln.end(), ',');

				std::string clazz;
				for (; at != comma; ++at) {
					if (*at != ' ') clazz += *at;
				}
				at = std::find(comma+1, ln.end(), ',');

				std::string dna;
				for (++at; at != ln.end(); ++at) {
					if (*at != ' ') dna += *at;
				}
				if (dna.size() != 60) {
					throw std::runtime_error("reading dna failed");
				}
				for (int i = 0; i < 60; ++i) {
					int let = 0;
					switch (dna[i]) {
						case 'A': let=0; break;
						case 'C': let=1; break;
						case 'T': let=2; break;
						case 'G': let=3; break;
					}
					*d.var(varid_low)[cvec<p>(s, i)] = (let&1) > 0;
					*d.var(varid_high)[cvec<p>(s, i)] = (let>>1) > 0;
					d.var(varid_low)[cvec<p>(s, i)] = true;
					d.var(varid_high)[cvec<p>(s, i)] = true;
				}
				s++;
			}
		}
	}

	pinfo get_info() {
		pinfo i;
		i.cvnames_.push_back("s"); // sample
		i.cvnames_.push_back("i"); // index
		i.vnames_.push_back("low");
		i.vnames_.push_back("high");
		return i;
	}

	void setup_test(test_tool& t) {
		typedef dnasplice_problem p;
		lang<p> l;
		setup_lang(l);
		data<p> d(l, cvec<p>(3190, 60));
		setup_data(l, d);

		stats<p> s(d);

		pinfo i = get_info();

		stats_info<p> si(i, s);
		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"naive info: "<<s.ndl()<<" bits\n";
		t<<"            "<<(s.ndl()/d.dim()[0])<<" bits/sample\n";
		t<<"            "<<(s.ndl()/d.dim().volume())<<" bits/dna\n";
	}

	void reexp_test(test_tool& t) {
		typedef dnasplice_problem p;
		lang<p> l;
		setup_lang(l);
		data<p> d(l, cvec<p>(3190, 60));
		setup_data(l, d);

		stats<p> s(d);

		t<<"naive info: "<<s.ndl()<<" bits\n";
		t<<"            "<<(s.ndl()/d.dim()[0])<<" bits/sample\n";
		t<<"            "<<(s.ndl()/d.dim().volume())<<" bits/dna\n\n";

		t<<"reexpressing:\n";

		learner<p> ln(l, s, 2000, 3);

		int add = 5;
		int exps = 0;
		while (true) {
			double preinfo = (s.ndl()/d.dim().volume());
			int added = ln.reexpress(true, add);
			exps += added;
			t<<"delta: "<<preinfo<<"->"<<(s.ndl()/d.dim().volume())<<" bits/dna\n\n";
			if (added < add) break;
		}

		t<<exps<<" expressions added\n\n";

		pinfo i = get_info();
		stats_info<p> si(i, s);

		t<<si.vars_tostring();

/*		t<<"vars:\n\n";
		t<<si.vars_tostring();*/

		s.update();

		t<<"naive info: "<<s.ndl()<<" bits\n";
		t<<"            "<<(s.ndl()/d.dim()[0])<<" bits/sample\n";
		t<<"            "<<(s.ndl()/d.dim().volume())<<" bits/dna\n";
	}

}

void adddnasplicetest(TestRunner& runner) {
	runner.add("dnasplice/setup", {"func"}, &setup_test);
	runner.add("dnasplice/reexp", {"func"}, &reexp_test);
}
