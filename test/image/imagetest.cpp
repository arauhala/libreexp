/*
 * imagetest.cpp
 *
 *  Created on: Jan 31, 2012
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"
#include "exptesttools.h"

#include <fstream>

namespace {

	struct image_problem {
		static const int DIM = 2; // two context variables
		static const int MAX_REL_VARS = 2; // two max relation variables
	};

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
	void populate(explib::data<P>& data, std::istream& in, int w, int h) {
		explib::data_var<P>& p = data.var(varid::pixel);
		explib::cvec<P> at(0, 0);

		std::string line;

		for (int j = 0; j < h; j++) {
			at[cvarid::y] = j;
			in>>line;
			for (int k = 0; k < w; k++) {
				at[cvarid::x] = k;
				p[at] = true;
				*p[at] = (line[k] == '#');
			}
		}
	}

	template <typename P>
	void setup_lang(explib::lang<P>& lang) {
		explib::ctx<P> pixel_ctx(explib::cvec<P>(0, 0, 0));

		lang.add_orig(explib::orig<P>(pixel_ctx));

		explib::rel<P>& rl(lang.alloc_rel(pixel_ctx)); // left right
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

		explib::rel<P>& lrd(lang.alloc_rel(pixel_ctx)); // left right down
		lrd.add_var(explib::cvec<P>(0, 0),
				   lang.var(varid::pixel));
		lrd.add_var(explib::cvec<P>(1, 1),
				  lang.var(varid::pixel));
		lang.rel_done();

		explib::rel<P>& lru(lang.alloc_rel(pixel_ctx)); // left right up
		lru.add_var(explib::cvec<P>(0, 1),
				    lang.var(varid::pixel));
		lru.add_var(explib::cvec<P>(1, 0),
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

	void setup_test(TestTool& t) {
		std::ifstream in("test/image/letters.txt");

		int w, h;
		in>>w>>h;

		t<<"image w:"<<w<<" h:"<<h<<"\n\n";

		typedef image_problem p;

		explib::cvec<p> dim(w, h);
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);

		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		explib::stats<p> stats(data);

		explib::pinfo names;
		setup_names<p>(names);

		explib::stats_info<p> si(names, stats);
		t<<"vars:\n\n";
		t<<si.vars_tostring();
		t<<"rels:\n\n";
		t<<si.rels_tostring();
		t<<"pixels:\n\n";
		exptest::print_pixels<p>(t,
								 data.var(varid::pixel),
								 cvarid::x,
								 cvarid::y);
	}

	void test_image(TestTool& t, const char* file, double threshold) {
		std::ifstream in(file);

		int w, h;
		in>>w>>h;

		t<<"image w:"<<w<<" h:"<<h<<" threshold:"<<threshold<<"\n\n";

		typedef image_problem p;

		explib::cvec<p> dim(w, h);
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);

		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		explib::stats<p> stats(data);
		std::vector<double> varorigPs;
		for (int i = 0; i < lang.var_count(); ++i) {
			varorigPs.push_back(stats.var(i).p());
		}
		explib::learner<p> learner(lang, stats, threshold, 0.35*threshold);

		double infoBefore = stats.naiveInfo();
		TimeSentry ms;
		int exps = learner.reexpress(true);
		t.ignored()<<"rexpression took "<<ms.ms()<<" ms.\n\n";
/*		while (true) {
			int e = learner.reexpress(true, 3);
			exps += e;
			t<<e<<" exps added. ";
			int rels = 0;
			for (int j = 0; j < lang.rel_count(); ++j) {
				if (!lang.rel(j).disabled()) rels++;
			}
			t<<rels<<" rels. \n";
			if (!e) break;

			if (exps > 100) break;
		}*/
		t<<exps<<" exps formed\n\n";
		double infoAfter = stats.naiveInfo();

		explib::pinfo names;
		setup_names<p>(names);

		explib::stats_info<p> si(names, stats);

		t<<"vars:\n\n";
		t<<si.drawn_vars_tostring_byexpbias(varorigPs, cvarid::x, cvarid::y);

		t<<"entropy "<<infoBefore<<" -> "<<infoAfter<<"\n";
	}

	void letters_test(TestTool& t) {
		test_image(t, "test/image/letters.txt", 9);
	}

	void invaders_test(TestTool& t) {
		test_image(t, "test/image/invaders.txt", 25);
	}

	typedef std::vector<std::pair<explib::var<image_problem>*,
							      explib::bits> >& effect_mask;

	void overlap_test(TestTool& t) {
		std::ifstream in("test/image/letters.txt");

		int w, h;
		in>>w>>h;

		t<<"image w:"<<w<<" h:"<<h<<"\n\n";

		typedef image_problem p;

		explib::cvec<p> dim(w, h);
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);

		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		explib::stats<p> stats(data);
		explib::learner<p> learner(lang, stats, 20, 20);

		learner.reexpress(true);

		explib::pinfo names;
		setup_names<p>(names);

		explib::stats_info<p> si(names, stats);

		for (int i = 0; i < lang.exp_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.exp(i))<<"->\n";
			t<<si.lang_info().drawn_var_rootmasks_to_string(lang.exp(i))<<"\n\n";
		}
	}


	void deps_test(TestTool& t) {
		std::ifstream in("test/image/invaders.txt");

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		explib::cvec<p> dim(w, h);
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		explib::stats<p> stats(data);
		explib::learner<p> learner(lang, stats, 10, 10);
		for (int i = 0; i < 40; ++i) {
			std::ostringstream buf;
			buf<<"exp:"<<i;
			{
				TimeSentry time;
				std::priority_queue<explib::candidate<p> > cands;
				learner.scan(cands);
				if (cands.empty()) break;
				const explib::candidate<p>& c( cands.top() );
				// lang.add_exp(c.rel_->data().rel_, c.state_);
				const explib::rel<p>& rel = c.rel_->data().rel_;
				int state = c.state_;
				int id = lang.origs_.size() + lang.exps_.size();
				lang.exps_.push_back(explib::exp<p>(rel, state));
				lang.exps_.back().init(id);
				{	// re-expression
					TimeSentry time2;
					data.var_added(lang.exps_.back());
					t.record({"property:data ms", buf.str()}, time2.ms());
				}
				{
					int relsBefore = lang.rel_count();
					TimeSentry time2;
					lang.exps_.back().gen_rels(lang);
					t.record({"property:gen rels ms", buf.str()}, time2.ms());
					int relsAfter = lang.rel_count();
					t.record({"property:rels added", buf.str()}, relsAfter - relsBefore);
					t.record({"property:postgen rel", buf.str()}, lang.enabled_rel_count());
				}
				{
					TimeSentry time2;
					learner.disable_rels();
					t.record({"property:filter ms", buf.str()}, time2.ms());
					t.record({"property:filtered rl", buf.str()}, lang.enabled_rel_count());
				}
				t.record({"property:rexp ms", buf.str()}, time.ms());
			}

			const std::vector<explib::implmask<p>>& masks( lang.exp_back().expmasks() );
			// mark affected variables as dirty:
			int v = explib::util::next_version();
			for (const explib::implmask<p>& m : masks) {
				data.var(m.var_->id()).version_ = v;
			}
			int influence = 0;
			for (size_t j = 0; j < data.var_count(); ++j) {
				if (data.var(j).version_ == v) influence++;
			}
			t.record({"property:influence", buf.str()}, masks.size());

			{
				TimeSentry time;
				stats.update();
				long ms = time.ms();
				t.record({"property:update ms", 		buf.str()},	      ms);
			}
		}

		Table table(
			t.report(ToTable<Average>({}, "property:", "exp:")));
		t.ignored()<<table;
	}

	void implmask_test(TestTool& t) {
		std::ifstream in("test/image/letters.txt");

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		explib::cvec<p> dim(w, h);
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		explib::stats<p> stats(data);
		explib::learner<p> learner(lang, stats, 20, 20);
		learner.reexpress(true, 5);

		explib::pinfo names;
		setup_names<p>(names);

		explib::stats_info<p> si(names, stats);

		t<<"vars:\n\n";
		t<<si.lang_info().drawn_vars_tostring(cvarid::x, cvarid::y);

		t<<"masks:\n\n";

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.var(i))<<"->\n\n";
			t<<si.lang_info().drawn_var_implmasks_to_string(lang.var(i))<<"\n\n";
		}
	}

	void applybits_test(TestTool& t, const char* file, int exps) {
		std::ifstream in(file);

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		explib::cvec<p> dim(w, h);
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		explib::stats<p> stats(data);
		explib::learner<p> learner(lang, stats, 20, 20);

		explib::pinfo names;
		setup_names<p>(names);

		explib::stats_info<p> si(names, stats);
		const explib::lang_info<p>& li(si.lang_info());

		t<<si.rels_tostring()<<"\n";

		double infoBefore = stats.naiveInfo();

		learner.reexpress(true, exps);

		for (int i = 0; i < lang.var_count(); ++i) {
			const explib::var<p>& v = lang.var(i);

			t<<si.var_tostring(i);
			t<<"\n\n";
			t<<li.drawn_var_tostring(cvarid::x, cvarid::y, v);
			t<<"\n";
			t<<si.drawn_data_var_tostring(cvarid::x, cvarid::y, data.var(v.id()));
			t<<"\n\n";
		}

		t<<si.rels_tostring()<<"\n";

		double infoAfter = stats.naiveInfo();
		t<<"entropy "<<infoBefore<<" -> "<<infoAfter<<"\n";
	}

	void applybits_1exp_test(TestTool& t) {
		applybits_test(t, "test/image/letters.txt", 1);
	}

	void applybits_2exp_test(TestTool& t) {
		applybits_test(t, "test/image/letters.txt", 2);
	}

	void applybits_3exp_test(TestTool& t) {
		applybits_test(t, "test/image/letters.txt", 3);
	}

	void inv_applybits_1exp_test(TestTool& t) {
		applybits_test(t, "test/image/invaders.txt", 1);
	}

	void inv_applybits_2exp_test(TestTool& t) {
		applybits_test(t, "test/image/invaders.txt", 2);
	}
	void inv_applybits_3exp_test(TestTool& t) {
		applybits_test(t, "test/image/invaders.txt", 3);
	}
	void inv_applybits_4exp_test(TestTool& t) {
		applybits_test(t, "test/image/invaders.txt", 4);
	}
	void inv_applybits_11exp_test(TestTool& t) {
		applybits_test(t, "test/image/invaders.txt", 11);
	}

}

void addimagetest(TestRunner& runner) {
	runner.add("image/setup", 	{"func"}, 	&setup_test);
	runner.add("image/invaders",{"func"},   &invaders_test);
	runner.add("image/letters", {"func"}, 	&letters_test);
	runner.add("image/deps", 	{"func"}, 	&deps_test);
	runner.add("image/overlap",	{"func"}, 	&overlap_test);
	runner.add("image/implmask",{"func"}, 	&implmask_test);
	runner.add("image/applybits_1exp",{"func"}, &applybits_1exp_test);
	runner.add("image/applybits_2exp",{"func"}, &applybits_2exp_test);
	runner.add("image/applybits_3exp",{"func"}, &applybits_3exp_test);
	runner.add("image/inv_applybits_1exp",{"func"}, &inv_applybits_1exp_test);
	runner.add("image/inv_applybits_2exp",{"func"}, &inv_applybits_2exp_test);
	runner.add("image/inv_applybits_3exp",{"func"}, &inv_applybits_3exp_test);
	runner.add("image/inv_applybits_4exp",{"func"}, &inv_applybits_4exp_test);
	runner.add("image/inv_applybits_11exp",{"func"}, &inv_applybits_11exp_test);
}
