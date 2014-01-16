/*
 * simple.cpp
 *
 *  Created on: Aug 23, 2011
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"
#include "exptesttools.h"

namespace {

	typedef reexp::traits2d simple_problem;


	const char* image[] = {
		"XX...",
		".XX..",
		".....",
		"..XX.",
		"...XX"
	};

	static const int Width = 5;

	static const int Height = 5;

	static const int VarCount = 1;

	// context variables
	namespace cvarid {
		static const int x = 0;
		static const int y = 1;
	}
	// bit variables (used to group bits)
	namespace varid {
		static const int pixel = 0;
	}
	// relations (used to organize bits)
	namespace relid {
		static const int left_right = 0;
		static const int up_down = 1;
	}
	const char* varnames[] =  {
		"px"
	};
	const char* cvarnames[] =  {
		"x",
		"y"
	};

	template <typename P>
	void populate(reexp::data<P>& data) {
		reexp::data_var<P>& p = data.var(varid::pixel);
		reexp::cvec<P> at(0, 0);

		for (int j = 0; j < Height; j++) {
			at[cvarid::y] = j;
			for (int k = 0; k < Width; k++) {
				at[cvarid::x] = k;
				p[at] = true;
				*p[at] = image[j][k] == 'X';
			}
		}
	}

	reexp::cvec<simple_problem> simple_dim() {
		reexp::cvec<simple_problem> v;
		v[cvarid::x] = Width;
		v[cvarid::y] = Height;
		return v;
	}

	template <typename P>
	void setup_lang(reexp::lang<P>& lang) {
		reexp::ctx<P> pixel_ctx(reexp::cvec<P>(0, 0, 0));

		lang.add_orig(reexp::orig<P>(pixel_ctx));

		reexp::rel<P>& rl(lang.alloc_rel(pixel_ctx)); // right left
		rl.add_var(reexp::cvec<P>(0, 0),
				   lang.var(varid::pixel));
		rl.add_var(reexp::cvec<P>(1, 0),
				  lang.var(varid::pixel));
		lang.rel_done();

		reexp::rel<P>& ud(lang.alloc_rel(pixel_ctx)); // up down
		ud.add_var(reexp::cvec<P>(0, 0),
				   lang.var(varid::pixel));
		ud.add_var(reexp::cvec<P>(0, 1),
				   lang.var(varid::pixel));
		lang.rel_done();
	}

	template <typename P>
	void setup_reg(reexp::lang<P>& lang, reexp::data<P>& data) {
		setup_lang(lang);
		populate(data);
	}

	template <typename P>
	void setup_names(reexp::pinfo& info) {
		for (int i = 0; i < VarCount; ++i) {
			info.vnames_.push_back(varnames[i]);
		}
		for (int i = 0; i < P::DIM; ++i) {
			info.cvnames_.push_back(cvarnames[i]);
		}
	}

	template <typename P>
	void print_overview(test_tool& t, reexp::stats_info<P>& si, reexp::data<P>& data) {
		t<<"vars:\n\n";
		t<<si.vars_tostring();
		t<<"rels:\n\n";
		t<<si.rels_tostring();

		for (int i = 0; i < data.lang().var_count(); ++i) {
			t<<"px"<<i<<"\n\n";

			exptest::print_pixels<P>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
			t<<"\n\n";
		}
	}

	void setuptest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::lang_info<p> li(i, lang);

		t<<"vars:\n\n";
		t<<li.vars_tostring()<<"\n\n";
		t<<"rels:\n\n";
		t<<li.rels_tostring()<<"\n\n";
		t<<"pixels:\n\n";
		exptest::print_pixels<p>(t,
								 data.var(varid::pixel),
								 cvarid::x,
								 cvarid::y);
	}

	void statstest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		exptest::print_pixels<p>(t,
								 data.var(varid::pixel),
								 cvarid::x,
								 cvarid::y);

		t<<"\n";

		reexp::stats_info<p> si(i, stats);

		print_overview(t, si, data);
	}

	template <typename P>
	void test_overlap(test_tool& t,
					  reexp::lang_info<P>& li,
					  const reexp::var<P>& v1,
					  const reexp::cvec<P>& shift,
					  const reexp::var<P>& v2) {
		t<<li.var_tostring(v1)<<"\n"<<li.var_tostring(v1, shift)<<"\n";
		if (reexp::exp<P>::is_overlap(v1, shift, v2)) {
			t<<"overlap";
		} else {
			t<<"independent";
		}
		t<<"\n\n";
	}

	void overlaptest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());
		setup_reg<p>(lang, data);

		lang.add_exp(lang.rel(relid::left_right), 3);
		const reexp::exp<p> le = lang.exp_back();

		reexp::pinfo pi;
		setup_names<p>(pi);
		reexp::lang_info<p> li(pi, lang);

		test_overlap(t, li, le, reexp::cvec<p>(-2, 0), le);
		test_overlap(t, li, le, reexp::cvec<p>(-1, 0), le);
		test_overlap(t, li, le, reexp::cvec<p>(0, 0), le);
		test_overlap(t, li, le, reexp::cvec<p>(1, 0), le);
		test_overlap(t, li, le, reexp::cvec<p>(2, 0), le);

		test_overlap(t, li, le, reexp::cvec<p>(-2, 1), le);
		test_overlap(t, li, le, reexp::cvec<p>(-1, 1), le);
		test_overlap(t, li, le, reexp::cvec<p>(0, 1), le);
		test_overlap(t, li, le, reexp::cvec<p>(1, 1), le);
		test_overlap(t, li, le, reexp::cvec<p>(2, 1), le);
	}


	void scantest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);
		t<<si.scan_tostring();
	}


	void learningtest(test_tool& t, int n) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);

		reexp::learner<p> learner(lang, stats);
		double e = stats.ndl();
		t<<"scan:\n\n"<<si.scan_tostring(3)<<"\n";
		for (int i = 0; i < n; ++i) {
			if (!learner.add_exp()) break;
			t<<"exp added\n\n";
			t<<"information "<<e<<" -> "<<stats.ndl()<<"\n\n";
			e = stats.ndl();
			t<<"scan:\n\n"<<si.scan_tostring(3)<<"\n";
		}
		print_overview<p>(t, si, data);
	}

	void learning1test(test_tool& t) {
		learningtest(t, 1);
	}
	void learning2test(test_tool& t) {
		learningtest(t, 2);
	}
	void learning3test(test_tool& t) {
		learningtest(t, 3);
	}

	void tryallexpstest(test_tool& t) {
		for (int r = 0; r < 2; ++r) {
			for (int s = 0; s < 4; ++s) {
				typedef simple_problem p;

				reexp::lang<p> lang;
				reexp::data<p> data(lang, simple_dim());

				setup_reg<p>(lang, data);

				reexp::stats<p> stats(data);

				reexp::pinfo i;
				setup_names<p>(i);
				reexp::stats_info<p> si(i, stats);

				double e = stats.ndl();
				lang.add_exp(lang.rel(r), s);
				t<<"exp added\n\n";

				t<<"vars:\n\n";
				t<<si.vars_tostring();

				t<<"information "<<e<<" -> "<<stats.ndl()<<"\n\n";

				if (stats.ndl() > e) {
					t<<"entropy increased!\n\n";
				}
			}
		}
	}

	void entropytest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);

		exptest::print_pixels<p>(t,
								 data.var(varid::pixel),
								 cvarid::x,
								 cvarid::y);

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		double e = stats.ndl();
		lang.add_exp(lang.rel(relid::left_right), 0);
		t<<"exp added\n\n";

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"information "<<e<<" -> "<<stats.ndl()<<"\n\n";

		if (stats.ndl() > e) {
			t<<"entropy increased!\n\n";
		}

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);


		}
	}

	void fourbittest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		setup_lang<p>(lang);

		reexp::data<p> data(lang, reexp::cvec<p>(4, 1));
		const char bitmap[] = "...X";

		{
			reexp::data_var<p>& px = data.var(varid::pixel);
			for (int i = 0; i < 4; i++) {
				reexp::cvec<p> at(i, 0);
				px[at] = true;
				*px[at] = bitmap[i] == 'X';
			}
		}

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
		}

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		double e = stats.ndl();
		lang.add_exp(lang.rel(relid::left_right), 0);
		t<<"exp added\n\n";

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"information "<<e<<" -> "<<stats.ndl()<<"\n\n";

		if (stats.ndl() > e) {
			t<<"entropy increased!\n\n";
		}

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
		}
	}

	void fourbitallexptest(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		setup_lang<p>(lang);

		reexp::data<p> data(lang, reexp::cvec<p>(4, 1));
		const char bitmap[] = "...X";

		{
			reexp::data_var<p>& px = data.var(varid::pixel);
			for (int i = 0; i < 4; i++) {
				reexp::cvec<p> at(i, 0);
				px[at] = true;
				*px[at] = bitmap[i] == 'X';
			}
		}

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
		}

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		double e = stats.ndl();
		lang.add_exp(lang.rel(relid::left_right), 2);
		lang.add_exp(lang.rel(relid::left_right), 0);
		lang.add_exp(lang.rel(relid::left_right), 1);
		lang.add_exp(lang.rel(relid::left_right), 3);
		t<<"4 exps added\n\n";

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"information "<<e<<" -> "<<stats.ndl()<<"\n\n";

		if (stats.ndl() > e) {
			t<<"entropy increased!\n\n";
		}

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
		}
	}

	void print_exp_rel_stats(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);
		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::learner<p> learner(lang, stats);
		reexp::stats_info<p> si(i, stats);

		learner.reexpress(0, 3);

		for (int i = 0; i < lang.exp_count(); ++i) {
			const reexp::exp<p>& e = lang.exp(i);
			const reexp::exp_rel_stats<p>& ers =
				stats.exp_rel(e.id());

			const reexp::bits& rel_defs = data.exp_rel_defs(e.id());

			t<<"exp: "<<si.var_tostring(e.id())<<":\n";
			t<<"exp rel defs: "<<reexp::vector_todensestring(rel_defs)<<"\n";
			t<<"    f: "<<rel_defs.popcount()<<" / "<<rel_defs.size()<<"\n";

			t<<"n: "<<ers.n();
			int expstate = e.state();
			int expstatef = ers.stateFreqs()[expstate];
			t<<"\nexpf: "<<expstatef<<"\n";
			t<<"\nvars:  ";
			for (size_t i = 0; i < ers.varFreqs().size(); ++i) {
				if (i) t<<", ";
				t<<ers.varFreqs()[i];
			}
			t<<"\nstats:  ";
			for (size_t i = 0; i < ers.stateFreqs().size(); ++i) {
				if (i) t<<", ";
				t<<ers.stateFreqs()[i];
			}

			reexp::rel<p> r = e.rel();

			for (size_t i  = 0; i < ers.varFreqs().size(); ++i) {
				bool vstate = r.varState(i, expstate);
				int vstatef = ers.varFreqs()[i];
				if (!vstate) vstatef = ers.n() - vstatef;

				// frequentist approach
				double pExpOnV = (expstatef) / (double)vstatef;

				// this is one plausible estimate
				double hatPExpOnV = (expstatef + 1) / (double)(vstatef + 2);

				t<<"\np(E|"<<si.lang_info().var_tostring(*r.entries()[i].var_)<<")="<<pExpOnV;
				t<<"\nhat p(E|"<<si.lang_info().var_tostring(*r.entries()[i].var_)<<")="<<hatPExpOnV;
			}

			t<<"\n\n";
		}
	}
#if 0

	void probreexp(test_tool& t) {
		typedef simple_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);
		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::learner<p> learner(lang, stats);
//		learner.reexpress(0, 1);

		reexp::stats_info<p> si(i, stats);

		reexp::pred<p> pred(stats);

		reexp::data<p> data2(lang, reexp::cvec<p>(5, 1));
		data2.var(0).defined().fill(true);
		data2.var(0).defined()[2] = false;
		data2.var(0).defined()[4] = false;
		data2.var(0).states()[1] = true;

		std::vector<double> ps = pred.reexpress(data2);
		for (size_t i = 0; i < ps.size(); ++i) {
			if (i) t<<", ";
			t<<ps[i];
		}
		t<<"\n";
	}
#endif

}

void addsimpletest(TestRunner& runner) {
#if 0
	runner.add("simple/probreexp", 		{"func"},   &probreexp);
#endif
	runner.add("simple/exprelstats", 	{"func"},   &print_exp_rel_stats);
	runner.add("simple/fourbitallexp", 	{"func"},   &fourbitallexptest);
	runner.add("simple/fourbit", 		{"func"},   &fourbittest);
	runner.add("simple/entropy", 		{"func"},   &entropytest);
	runner.add("simple/tryallexps", 	{"func"},   &tryallexpstest);
	runner.add("simple/learning_1", 	{"func"},   &learning1test);
	runner.add("simple/learning_2", 	{"func"},   &learning2test);
	runner.add("simple/learning_3", 	{"func"},   &learning3test);
	runner.add("simple/scan", 			{"func"},   &scantest);
	runner.add("simple/stats", 			{"func"},   &statstest);
	runner.add("simple/overlap", 		{"func"},   &overlaptest);
	runner.add("simple/setup", 			{"func"},   &setuptest);
}
