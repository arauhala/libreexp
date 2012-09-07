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

	struct simple_problem {
		static const int DIM = 2; // two context variables
		static const int MAX_REL_VARS = 2; // two max relation variables
	};

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
	void populate(explib::data<P>& data) {
		explib::data_var<P>& p = data.var(varid::pixel);
		explib::cvec<P> at(0, 0);

		for (int j = 0; j < Height; j++) {
			at[cvarid::y] = j;
			for (int k = 0; k < Width; k++) {
				at[cvarid::x] = k;
				p[at] = true;
				*p[at] = image[j][k] == 'X';
			}
		}
	}

	explib::cvec<simple_problem> simple_dim() {
		explib::cvec<simple_problem> v;
		v[cvarid::x] = Width;
		v[cvarid::y] = Height;
		return v;
	}

	template <typename P>
	void setup_lang(explib::lang<P>& lang) {
		explib::ctx<P> pixel_ctx(explib::cvec<P>(0, 0, 0));

		lang.add_orig(explib::orig<P>(pixel_ctx));

		explib::rel<P>& rl(lang.alloc_rel(pixel_ctx)); // right left
		rl.add_var(explib::cvec<P>(0, 0),
				   lang.var(varid::pixel));
		rl.add_var(explib::cvec<P>(1, 0),
				  lang.var(varid::pixel));
		lang.rel_done();

		explib::rel<P>& ud(lang.alloc_rel(pixel_ctx)); // up down
		ud.add_var(explib::cvec<P>(0, 0),
				   lang.var(varid::pixel));
		ud.add_var(explib::cvec<P>(0, 1),
				   lang.var(varid::pixel));
		lang.rel_done();
	}

	template <typename P>
	void setup_reg(explib::lang<P>& lang, explib::data<P>& data) {
		setup_lang(lang);
		populate(data);
	}

	template <typename P>
	void setup_names(explib::pinfo& info) {
		for (int i = 0; i < VarCount; ++i) {
			info.vnames_.push_back(varnames[i]);
		}
		for (int i = 0; i < P::DIM; ++i) {
			info.cvnames_.push_back(cvarnames[i]);
		}
	}

	template <typename P>
	void print_overview(TestTool& t, explib::stats_info<P>& si, explib::data<P>& data) {
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

	void setuptest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		explib::pinfo i;
		setup_names<p>(i);

		explib::lang_info<p> li(i, lang);

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

	void statstest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);

		exptest::print_pixels<p>(t,
								 data.var(varid::pixel),
								 cvarid::x,
								 cvarid::y);

		t<<"\n";

		explib::stats_info<p> si(i, stats);

		print_overview(t, si, data);
	}

	template <typename P>
	void test_overlap(TestTool& t,
					  explib::lang_info<P>& li,
					  const explib::var<P>& v1,
					  const explib::cvec<P>& shift,
					  const explib::var<P>& v2) {
		t<<li.var_tostring(v1)<<"\n"<<li.var_tostring(v1, shift)<<"\n";
		if (explib::exp<P>::is_overlap(v1, shift, v2)) {
			t<<"overlap";
		} else {
			t<<"independent";
		}
		t<<"\n\n";
	}

	void overlaptest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, simple_dim());
		setup_reg<p>(lang, data);

		lang.add_exp(lang.rel(relid::left_right), 3);
		const explib::exp<p> le = lang.exp_back();

		explib::pinfo pi;
		setup_names<p>(pi);
		explib::lang_info<p> li(pi, lang);

		test_overlap(t, li, le, explib::cvec<p>(-2, 0), le);
		test_overlap(t, li, le, explib::cvec<p>(-1, 0), le);
		test_overlap(t, li, le, explib::cvec<p>(0, 0), le);
		test_overlap(t, li, le, explib::cvec<p>(1, 0), le);
		test_overlap(t, li, le, explib::cvec<p>(2, 0), le);

		test_overlap(t, li, le, explib::cvec<p>(-2, 1), le);
		test_overlap(t, li, le, explib::cvec<p>(-1, 1), le);
		test_overlap(t, li, le, explib::cvec<p>(0, 1), le);
		test_overlap(t, li, le, explib::cvec<p>(1, 1), le);
		test_overlap(t, li, le, explib::cvec<p>(2, 1), le);
	}


	void scantest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);
		explib::stats_info<p> si(i, stats);
		t<<si.scan_tostring();
	}


	void learningtest(TestTool& t, int n) {
		typedef simple_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);
		explib::stats_info<p> si(i, stats);

		explib::learner<p> learner(lang, stats);
		double e = stats.naiveInfo();
		t<<"scan:\n\n"<<si.scan_tostring(3)<<"\n";
		for (int i = 0; i < n; ++i) {
			if (!learner.add_exp()) break;
			t<<"exp added\n\n";
			t<<"information "<<e<<" -> "<<stats.naiveInfo()<<"\n\n";
			e = stats.naiveInfo();
			t<<"scan:\n\n"<<si.scan_tostring(3)<<"\n";
		}
		print_overview<p>(t, si, data);
	}

	void learning1test(TestTool& t) {
		learningtest(t, 1);
	}
	void learning2test(TestTool& t) {
		learningtest(t, 2);
	}
	void learning3test(TestTool& t) {
		learningtest(t, 3);
	}

	void tryallexpstest(TestTool& t) {
		for (int r = 0; r < 2; ++r) {
			for (int s = 0; s < 4; ++s) {
				typedef simple_problem p;

				explib::lang<p> lang;
				explib::data<p> data(lang, simple_dim());

				setup_reg<p>(lang, data);

				explib::stats<p> stats(data);

				explib::pinfo i;
				setup_names<p>(i);
				explib::stats_info<p> si(i, stats);

				double e = stats.naiveInfo();
				lang.add_exp(lang.rel(r), s);
				t<<"exp added\n\n";

				t<<"vars:\n\n";
				t<<si.vars_tostring();

				t<<"information "<<e<<" -> "<<stats.naiveInfo()<<"\n\n";

				if (stats.naiveInfo() > e) {
					t<<"entropy increased!\n\n";
				}
			}
		}
	}

	void entropytest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, simple_dim());

		setup_reg<p>(lang, data);

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);
		explib::stats_info<p> si(i, stats);

		exptest::print_pixels<p>(t,
								 data.var(varid::pixel),
								 cvarid::x,
								 cvarid::y);

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		double e = stats.naiveInfo();
		lang.add_exp(lang.rel(relid::left_right), 0);
		t<<"exp added\n\n";

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"information "<<e<<" -> "<<stats.naiveInfo()<<"\n\n";

		if (stats.naiveInfo() > e) {
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

	void fourbittest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		setup_lang<p>(lang);

		explib::data<p> data(lang, explib::cvec<p>(4, 1));
		const char bitmap[] = "...X";

		{
			explib::data_var<p>& px = data.var(varid::pixel);
			for (int i = 0; i < 4; i++) {
				explib::cvec<p> at(i, 0);
				px[at] = true;
				*px[at] = bitmap[i] == 'X';
			}
		}

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);
		explib::stats_info<p> si(i, stats);

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
		}

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		double e = stats.naiveInfo();
		lang.add_exp(lang.rel(relid::left_right), 0);
		t<<"exp added\n\n";

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"information "<<e<<" -> "<<stats.naiveInfo()<<"\n\n";

		if (stats.naiveInfo() > e) {
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

	void fourbitallexptest(TestTool& t) {
		typedef simple_problem p;

		explib::lang<p> lang;
		setup_lang<p>(lang);

		explib::data<p> data(lang, explib::cvec<p>(4, 1));
		const char bitmap[] = "...X";

		{
			explib::data_var<p>& px = data.var(varid::pixel);
			for (int i = 0; i < 4; i++) {
				explib::cvec<p> at(i, 0);
				px[at] = true;
				*px[at] = bitmap[i] == 'X';
			}
		}

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);
		explib::stats_info<p> si(i, stats);

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<"\n"<<si.lang_info().var_tostring(i)<<"\n\n";
			exptest::print_pixels<p>(t,
									 data.var(i),
									 cvarid::x,
									 cvarid::y);
		}

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		double e = stats.naiveInfo();
		lang.add_exp(lang.rel(relid::left_right), 2);
		lang.add_exp(lang.rel(relid::left_right), 0);
		lang.add_exp(lang.rel(relid::left_right), 1);
		lang.add_exp(lang.rel(relid::left_right), 3);
		t<<"4 exps added\n\n";

		t<<"vars:\n\n";
		t<<si.vars_tostring();

		t<<"information "<<e<<" -> "<<stats.naiveInfo()<<"\n\n";

		if (stats.naiveInfo() > e) {
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

}

void addsimpletest(TestRunner& runner) {
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
