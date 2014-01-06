/*
 * imagetest.cpp
 *
 *  Created on: Jan 31, 2012
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"
#include "exptesttools.h"

#include "reexp/io.h"

#include <fstream>

namespace {

	typedef reexp::traits2d image_problem;

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
	void populate(reexp::data<P>& data, std::istream& in, int w, int h) {
		reexp::data_var<P>& p = data.var(varid::pixel);
		reexp::cvec<P> at(0, 0);

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
	void setup_lang(reexp::lang<P>& lang) {
		reexp::ctx<P> pixel_ctx(reexp::cvec<P>(0, 0, 0));

		lang.add_orig(reexp::orig<P>(pixel_ctx));

		reexp::rel<P>& rl(lang.alloc_rel(pixel_ctx)); // left right
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

		reexp::rel<P>& lrd(lang.alloc_rel(pixel_ctx)); // left right down
		lrd.add_var(reexp::cvec<P>(0, 0),
				   lang.var(varid::pixel));
		lrd.add_var(reexp::cvec<P>(1, 1),
				  lang.var(varid::pixel));
		lang.rel_done();

		reexp::rel<P>& lru(lang.alloc_rel(pixel_ctx)); // left right up
		lru.add_var(reexp::cvec<P>(1, 0),
				    lang.var(varid::pixel));
		lru.add_var(reexp::cvec<P>(0, 1),
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

	void setup_test(test_tool& t) {
		std::ifstream in("test/image/letters.txt");

		int w, h;
		in>>w>>h;

		t<<"image w:"<<w<<" h:"<<h<<"\n\n";

		typedef image_problem p;

		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);

		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);

		reexp::pinfo names;
		setup_names<p>(names);

		reexp::stats_info<p> si(names, stats);
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

	void test_image(test_tool& t, const char* file, double threshold) {
		std::ifstream in(file);

		int w, h;
		in>>w>>h;

		t<<"image w:"<<w<<" h:"<<h<<" threshold:"<<threshold<<"\n\n";

		typedef image_problem p;

		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);

		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		std::vector<double> varorigPs;
		for (int i = 0; i < lang.var_count(); ++i) {
			varorigPs.push_back(stats.var(i).p());
		}
		reexp::learner<p> learner(lang, stats, threshold, 0.35*threshold);

		double infoBefore = stats.naiveInfo();
		time_sentry ms;
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

		reexp::pinfo names;
		setup_names<p>(names);

		reexp::stats_info<p> si(names, stats);

		t<<"vars:\n\n";
		t<<si.drawn_vars_tostring_byexpbias(varorigPs, cvarid::x, cvarid::y);

		t<<"entropy "<<infoBefore<<" -> "<<infoAfter<<"\n";
	}

	void letters_test(test_tool& t) {
		test_image(t, "test/image/letters.txt", 9);
	}

	void invaders_test(test_tool& t) {
		test_image(t, "test/image/invaders.txt", 25);
	}

	typedef std::vector<std::pair<reexp::var<image_problem>*,
							      reexp::bits> >& effect_mask;

	void overlap_test(test_tool& t) {
		std::ifstream in("test/image/letters.txt");

		int w, h;
		in>>w>>h;

		t<<"image w:"<<w<<" h:"<<h<<"\n\n";

		typedef image_problem p;

		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);

		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, 20, 20);

		learner.reexpress(true);

		reexp::pinfo names;
		setup_names<p>(names);

		reexp::stats_info<p> si(names, stats);

		for (int i = 0; i < lang.exp_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.exp(i))<<"->\n";
			t<<si.lang_info().drawn_var_rootmasks_to_string(lang.exp(i))<<"\n\n";
		}
	}


	void deps_test(test_tool& t) {
		std::ifstream in("test/image/invaders.txt");

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, 10, 10);
		for (int i = 0; i < 40; ++i) {
			std::ostringstream buf;
			buf<<"exp:"<<i;
			{
				time_sentry time;
				std::priority_queue<reexp::candidate<p> > cands;
				learner.scan(cands);
				if (cands.empty()) break;
				const reexp::candidate<p>& c( cands.top() );
				// lang.add_exp(c.rel_->data().rel_, c.state_);
				const reexp::rel<p>& rel = c.rel_->data().rel_;
				int state = c.state_;
				int id = lang.origs_.size() + lang.exps_.size();
				lang.exps_.push_back(reexp::exp<p>(rel, state));
				lang.exps_.back().init(id);
				{	// re-expression
					time_sentry time2;
					data.var_added(lang.exps_.back());
					t.record({"property:data ms", buf.str()}, time2.ms());
				}
				{
					int relsBefore = lang.rel_count();
					time_sentry time2;
					lang.exps_.back().gen_rels(lang);
					t.record({"property:gen rels ms", buf.str()}, time2.ms());
					int relsAfter = lang.rel_count();
					t.record({"property:rels added", buf.str()}, relsAfter - relsBefore);
					t.record({"property:postgen rel", buf.str()}, lang.enabled_rel_count());
				}
				{
					time_sentry time2;
					learner.disable_rels();
					t.record({"property:filter ms", buf.str()}, time2.ms());
					t.record({"property:filtered rl", buf.str()}, lang.enabled_rel_count());
				}
				t.record({"property:rexp ms", buf.str()}, time.ms());
			}

			const std::vector<reexp::implmask<p>>& masks( lang.exp_back().expmasks() );
			// mark affected variables as dirty:
			int v = reexp::util::next_version();
			for (const reexp::implmask<p>& m : masks) {
				data.var(m.var_->id()).version_ = v;
			}
			int influence = 0;
			for (size_t j = 0; j < data.var_count(); ++j) {
				if (data.var(j).version_ == v) influence++;
			}
			t.record({"property:influence", buf.str()}, masks.size());

			{
				time_sentry time;
				stats.update();
				long ms = time.ms();
				t.record({"property:update ms", 		buf.str()},	      ms);
			}
		}

		table table(
			t.report(to_table<average>({}, "property:", "exp:")));
		t.ignored()<<table;
	}

	void implmask_test(test_tool& t) {
		std::ifstream in("test/image/letters.txt");

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, 20, 20);
		learner.reexpress(true, 5);

		reexp::pinfo names;
		setup_names<p>(names);

		reexp::stats_info<p> si(names, stats);

		t<<"vars:\n\n";
		t<<si.lang_info().drawn_vars_tostring(cvarid::x, cvarid::y);

		t<<"masks:\n\n";

		for (int i = 0; i < lang.var_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.var(i))<<"->\n\n";
			t<<si.lang_info().drawn_var_implmasks_to_string(lang.var(i))<<"\n\n";
		}
	}

	void applybits_test(test_tool& t, const char* file, int exps) {
		std::ifstream in(file);

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, 20, 20);

		reexp::pinfo names;
		setup_names<p>(names);

		reexp::stats_info<p> si(names, stats);
		const reexp::lang_info<p>& li(si.lang_info());

		t<<si.rels_tostring()<<"\n";

		double infoBefore = stats.naiveInfo();

		learner.reexpress(true, exps);

		for (int i = 0; i < lang.var_count(); ++i) {
			const reexp::var<p>& v = lang.var(i);

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

	void applybits_1exp_test(test_tool& t) {
		applybits_test(t, "test/image/letters.txt", 1);
	}

	void applybits_2exp_test(test_tool& t) {
		applybits_test(t, "test/image/letters.txt", 2);
	}

	void applybits_3exp_test(test_tool& t) {
		applybits_test(t, "test/image/letters.txt", 3);
	}

	void inv_applybits_1exp_test(test_tool& t) {
		applybits_test(t, "test/image/invaders.txt", 1);
	}

	void inv_applybits_2exp_test(test_tool& t) {
		applybits_test(t, "test/image/invaders.txt", 2);
	}
	void inv_applybits_3exp_test(test_tool& t) {
		applybits_test(t, "test/image/invaders.txt", 3);
	}
	void inv_applybits_4exp_test(test_tool& t) {
		applybits_test(t, "test/image/invaders.txt", 4);
	}
	void inv_applybits_11exp_test(test_tool& t) {
		applybits_test(t, "test/image/invaders.txt", 11);
	}

	void io_test(test_tool& t, const char* file, bool verbose, double threshold, int exps = -1) {
		std::ifstream in(file);

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, threshold, 0.35*threshold);
		learner.reexpress(true, exps);


		reexp::pinfo names;
		setup_names<p>(names);
		reexp::stats_info<p> si(names, stats);

		t<<"variablesl\n";
		t<<si.vars_tostring();

		{
			reexp::bits buffer(64*1024*1024); // big enough buffer
			reexp::bit_ostream bout = buffer.ostream(0);
			reexp::arithmetic_bit_ostream out(bout);

			reexp::io<p> io(stats);
			io.write_states(out, data);
			out.finish();

			t<<"re-expressed image was serialized in "<<bout.pos()<<" bits, \n";
			t<<"while naive description length was "<<stats.naiveInfo()<<".\n";


			reexp::bit_istream bin = buffer.istream(0);
			reexp::arithmetic_bit_istream in(bin);
			reexp::data<p> data2(lang, dim);
			data2.prepare_exp_vars();
			io.read_states(in, data2);

			t<<bin.pos()<<" bits read from encoded blob\n";

			t<<"\n";

			if (verbose) t<<"restored image, variable by variable:\n";

			int ddiff = 0;
			int sdiff = 0;
			for (int i = 0; i < lang.var_count(); ++i) {
//				t<<si.lang_info().drawn_var_expmasks_to_string(lang.var(i));
				reexp::bits b;
				b = data.var(i).defined();
				b ^= data2.var(i).defined();
				int dd = b.popcount();
				b = data.var(i).states();
				b ^= data2.var(i).states();
				int sd = b.popcount();
				ddiff += dd;
				sdiff += sd;

				b = data.var(i).states();
				b ^= data2.var(i).states();
				sdiff += b.popcount();

				if (verbose || sd || dd) {
					t<<"variable "<<si.var_tostring(i);

					if (sd || dd) {
						t<<" failed."<<faill;
					} else {
						t<<" ok.\n";
					}
					t<<dd<<" errors in defined.\n";

					t<<sd<<" errors in states.\n";

					t<<"\noriginal:\n\n";
					t<<si.drawn_data_var_tostring(cvarid::x, cvarid::y, data.var(i));

					t<<"\nrestored:\n\n";
					t<<si.drawn_data_var_tostring(cvarid::x, cvarid::y, data2.var(i));

					t<<"\n";
					t<<"\n\n";
				}
			}

			if (sdiff || ddiff) {
				t<<"failed"<<faill;
			} else {
				t<<"ok.\n";
			}
			t<<ddiff<<" errors in defined.\n";
			t<<sdiff<<" errors in states.\n";
		}
	}

	void generate_test(test_tool& t, const char* file, bool verbose, double threshold, int exps = -1) {
		std::ifstream in(file);

		int w, h;
		in>>w>>h;
		typedef image_problem p;
		reexp::cvec<p> dim(w, h);
		reexp::lang<p> lang;
		reexp::data<p> data(lang, dim);
		setup_lang<p>(lang);
		populate<p>(data, in, w, h);

		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, threshold, 0.35*threshold);
		learner.reexpress(true, exps);


		reexp::pinfo names;
		setup_names<p>(names);
		reexp::stats_info<p> si(names, stats);

		t<<"variablesl\n";
		t<<si.vars_tostring();

		{
			reexp::bits buffer(1024*1024); // big enough buffer

			// fill the buffer with garbage
			for (int i = 0; i < buffer.size(); ++i) {
				buffer[i] = rand() % 2;
			}
			reexp::bit_istream bin = buffer.istream(0);
			reexp::arithmetic_bit_istream in(bin);
			reexp::data<p> data2(lang, dim);
			data2.prepare_exp_vars();
			reexp::io<p> io(stats);
			io.read_states(in, data2);

			t<<bin.pos()<<" bits read from encoded blob\n";

			t<<"\n";

			t<<"\n generated image:\n\n";
			t<<si.drawn_data_tostring(reexp::cvec<p>(), cvarid::x, cvarid::y, data2);
		}
	}

	void simple_io_test(test_tool& t) {
		io_test(t, "test/image/simple.txt", true, 2, 2);
	}

	void verbose_invaders_io_test(test_tool& t) {
		io_test(t, "test/image/invaders.txt", true, 25, -1);
	}

	void invaders_io_test(test_tool& t) {
		io_test(t, "test/image/invaders.txt", false, 25, -1);
	}

	void generate_invaders_test(test_tool& t) {
		generate_test(t, "test/image/invaders.txt", false, 25, -1);
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
	runner.add("image/io",{"func"}, 	&simple_io_test);
	runner.add("image/verbose_invaders_io",{"func"}, 	&verbose_invaders_io_test);
	runner.add("image/invaders_io",{"func"}, 	&invaders_io_test);
	runner.add("image/generate_invaders",{"func"}, 	&generate_invaders_test);
}
