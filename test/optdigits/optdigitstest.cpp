

#include "reexp/all.h"

#include "tester.h"
#include "exptesttools.h"

#include "stdio.h"

#include <iostream>
#include <algorithm>

#include "optdigits.h"

// bit variables (used to group bits)

namespace optdigits {
	const char* TRA_DATA_FILE = "test/optdigits/optdigits-orig.tra";
	int TRA_DATA_FILE_SAMPLES = 1934;
	const char* CV_DATA_FILE = "test/optdigits/optdigits-orig.cv";
	int CV_DATA_FILE_SAMPLES = 946;
	const char* WDEP_DATA_FILE = "test/optdigits/optdigits-orig.wdep";
	int WDEP_DATA_FILE_SAMPLES = 943;
	// use this one testing
	const char* WINDEP_DATA_FILE = "test/optdigits/optdigits-orig.windep";
	int WINDEP_DATA_FILE_SAMPLES = 1797;

	namespace varid {
		static const int pixel = 10;
	}
}

namespace {
	using namespace optdigits;

	const int DefaultThreshold = 2000;
	const int FilterThreshold = 1000;
	//
	// Each data entry is 7*3 + 9 = 30 bits
	//

	static const int Pack = 2;

	static const int Width = OriginalWidth / Pack;

	static const int Height = OriginalHeight / Pack;

	static const int DefaultSampleCount = 500;

	static const int DigitCount = 10;

	static const int VarCount = 11;

	// context variables
	namespace cvarid {
		static const int x = 0;
		static const int y = 1;
		static const int sample = 2;
	}
	// relations (used to organize bits)
	namespace relid {
		static const int left_right = 0;
		static const int up_down = 1;
	}
	const char* varnames[] =  {
		"d0",
		"d1",
		"d2",
		"d3",
		"d4",
		"d5",
		"d6",
		"d7",
		"d8",
		"d9",
		"px"
	};
	const char* cvarnames[] =  {
		"x",
		"y",
		"n" // sample
	};

	struct optdigits_problem {
		static const int DIM = 3;
		static const int MAX_REL_VARS = 2;
	};

	reexp::cvec<optdigits_problem> optdigits_dim(int samples) {
		return reexp::cvec<optdigits_problem>(Width, Height, samples);
	}

	template <typename P>
	void populate(reexp::data<P>& data, int samples, std::ifstream& in, int offset) {
		reexp::data_var<P>& p = data.var(varid::pixel);
		reexp::cvec<P> at(0, 0, 0);

		std::string line;

		for (int i = 0; i < DATA_FILE_HEADER_LINES; ++i) {
			std::getline(in, line);
		}

		reexp::bitmap pic;

		for (int i = 0; i < samples; i++) {
			// read
			for (int j = 0; j < OriginalHeight; j++) {
				std::getline(in, line);
				for (int k = 0; k < OriginalWidth; k++) {
					pic.set(k, j, (line[k] == '1'));
				}
			}

			std::getline(in, line);
			int number = line[1]-'0';

			// prepare an entry
			at[cvarid::sample] = i + offset;
			for (int x = 0; x < Width; x++) {
				at[cvarid::x] = x;
				for (int y = 0; y < Height; y++) {
					at[cvarid::y] = y;
					int popcount = 0;
					for (int x2 = 0; x2 < Pack; x2++) {
						for (int y2 = 0; y2 < Pack; y2++) {
							if (pic.get(x*Pack + x2, y*Pack + y2)) popcount++;
						}
					}
					p[at] = true;
					*p[at] = (popcount*2 >= Pack*Pack);
				}
			}

			for (int j = 0; j < DigitCount; j++) {
				reexp::data_var<P>& digit = data.var(varid::digit0 + j);
				digit[at] = true;
				*digit[at] = (number == j);
			}
		}
	}

	template <typename P>
	int populate_from_file(reexp::data<P>& data, int samples, const char* name, int offset = 0) {
		std::ifstream in(name);
		populate<P>(data, samples, in, offset);
		return offset + samples;
	}


	template <typename P>
	void setup_lang(reexp::lang<P>& lang) {
		reexp::ctx<P> pixel_ctx(reexp::cvec<P>(0, 0, 0));
		reexp::ctx<P> shape_ctx(reexp::cvec<P>(-1, -1, 0));
		reexp::ctx<P> full_ctx(reexp::cvec<P>(0, 0, 0));

		for (int i = 0; i < DigitCount; ++i) {
			lang.add_orig(reexp::orig<P>(shape_ctx));
		}
		lang.add_orig(reexp::orig<P>(pixel_ctx));

		reexp::rel<P>& rl(lang.alloc_rel(pixel_ctx)); // left right
		rl.add_var(reexp::cvec<P>(0, 0, 0),
				   lang.var(varid::pixel));
		rl.add_var(reexp::cvec<P>(1, 0, 0),
				  lang.var(varid::pixel));
		lang.rel_done();

		reexp::rel<P>& ud(lang.alloc_rel(pixel_ctx)); // up down
		ud.add_var(reexp::cvec<P>(0, 0, 0),
				   lang.var(varid::pixel));
		ud.add_var(reexp::cvec<P>(0, 1, 0),
				   lang.var(varid::pixel));
		lang.rel_done();

/*		explib::rel<P>& ldr(lang.alloc_rel(pixel_ctx)); // left down right
		ldr.add_var(explib::cvec<P>(0, 0, 0),
				    lang.var(varid::pixel));
		ldr.add_var(explib::cvec<P>(1, 1, 0),
				    lang.var(varid::pixel));
		lang.rel_done();

		explib::rel<P>& dlr(lang.alloc_rel(pixel_ctx)); // down left right
		dlr.add_var(explib::cvec<P>(0, 1, 0),
				    lang.var(varid::pixel));
		dlr.add_var(explib::cvec<P>(1, 0, 0),
				    lang.var(varid::pixel));
		lang.rel_done();*/

		for (int i = 0; i < DigitCount; ++i) {
			reexp::rel<P>& ps(lang.alloc_rel(full_ctx)); // pixel-shape
			ps.add_var(reexp::cvec<P>(0, 0, 0), // pixel
					   lang.var(varid::pixel));
			ps.add_var(reexp::cvec<P>(0, 0, 0), // shape
					   lang.var(varid::digit0 + i));
			lang.rel_done();
		}
	}

	template <typename P>
	void setup_problem(reexp::lang<P>& lang, reexp::data<P>& data) {
		setup_lang(lang);
		populate_from_file(data, DefaultSampleCount, TRA_DATA_FILE, 0);
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
	void setup_learner(reexp::learner<P>& learner) {
		for (int i = varid::digit0; i <= varid::digit9; ++i) {
			learner.exclude(i);
		}
	}

	void setup_test(test_tool& t) {
		typedef optdigits_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(DefaultSampleCount));

		setup_problem<p>(lang, data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::lang_info<p> li(i, lang);

		t<<"vars:\n\n";
		t<<li.vars_tostring()<<"\n\n";
		t<<"rels:\n\n";
		t<<li.rels_tostring()<<"\n\n";
		t<<"pixels:\n\n";

		for (int i = 0; i < 5; ++i) {
			reexp::cvec<p> at(0, 0, i);
			for (int j = 0; j < DigitCount; ++j) {
				if (*data.var(varid::digit0+j)[at]) {
					t<<j<<":\n\n";
					break;
				}
			}

			exptest::print_pixels<p>(t,
									 data.var(varid::pixel),
									 cvarid::x,
									 cvarid::y,
									 at);
			t<<"\n";
		}
	}

	void learning_test(test_tool& t) {
		typedef optdigits_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(DefaultSampleCount));

		setup_problem<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);


		reexp::stats_info<p> si(i, stats);
		t<<"scan:\n\n"<<si.scan_tostring(3, 2)<<"\n";

		reexp::learner<p> learner(lang, stats, DefaultThreshold, FilterThreshold);
		setup_learner<p>(learner);

		{
			time_sentry time;
			while (true) {
				float before = stats.naiveInfo();
				if (!learner.add_exp()) break;
				t<<"expression added."<<checkl;
				float after = stats.naiveInfo();
				t<<"info "<<before<<" -> "<<after<<checkl<<checkl;

				t<<"scan:\n\n"<<si.scan_tostring(3, 2)<<"\n";
			}
///			printf("\nreexpression took %dms.\n", int(time.ms()));
		}

		t<<checkl<<"learning done."<<checkl<<checkl;

		t<<"scan:\n\n"<<si.scan_tostring(3, 2)<<"\n";

		t<<"vars:\n\n"<<si.vars_tostring()<<"\n";
		t<<"rels:\n\n"<<si.rels_tostring()<<"\n";
	}

#if 0

	template <typename P>
	void evaluate(test_tool& t,
			      const reexp::data<P>& data,
			      const reexp::stats<P>& stats,
			      const std::set<std::string>& tags) {
		reexp::pred<P> pred(stats);
		std::vector<double> ps[10];

		time_sentry timer;
		for (int i = 0; i < 10; i++) {
			ps[i] = std::move(pred.p(data, i+varid::digit0));
		}
		long us = timer.us();
		int samples = ps[0].size();
		t.record(tags+"perf:pred us", us/samples);

		std::vector<int> decisions;
		std::vector<int> correct;
		decisions.resize(samples);
		correct.resize(samples);
		const reexp::data_var<P>* vars[10];

		for (int i = 0; i < 10; i++) {
			vars[i] = &data.var(i + varid::digit0);
		}

		timer.reset();
		double infos[10];
		int oks[10];
		int ns[10];
		for (int i = 0; i < 10; ++i) {
			oks[i] = 0; infos[i] = 0; ns[i] = 0;
		}
		reexp::cvec<P> at;
		int ok = 0;
		double info = 0;
		for (int s = 0; s < samples; s++) {
			int digit = -1;
			// find out the right answer
			at[cvarid::sample] = s;
			for (int i = 0; i < 10; i++) {
				if (*((*vars[i])[at])) {
					correct[s] = digit = i;
					break;
				}
			}
			// normalize
			double totalP = 0;
			for (int d = 0; d < 10; d++) {
				totalP += ps[d][s];
			}
			for (int d = 0; d < 10; d++) {
				ps[d][s] /= totalP;
			}
			// find out the peak probability & calculate information
			int top = -1;
			double topP = -1;
			double i = 0;
			for (int d = 0; d < 10; d++) {
				if (ps[d][s] > topP) {
					top = d;
					topP = ps[d][s];
				}
			}
			i = -log(ps[digit][s]) / log(2);
			info += i;
			infos[digit] += i;
			ns[digit]++;
			decisions[s] = top;
			if (digit == top) {
				ok++;
				oks[digit]++;
			}
		}
		us = timer.us();
		t.record(tags+"perf:norm us", us/samples);

		info /= samples;

		t.record(tags+"prop:entropy", info);
		t.record(tags+"prop:accuracy", (ok/(double)samples));
		t.record(tags+"prop:error", ((samples-ok)/(double)samples));
	}
#endif
	template <typename P>
	std::string predictions_tostring(const reexp::data<P>& data, const reexp::stats<P>& stats, bool details = false) {
		std::ostringstream buf;
		reexp::pred<P> pred(stats);
		std::vector<double> ps[10];
		for (int i = 0; i < 10; i++) {
			ps[i] = std::move(pred.p(data, i+varid::digit0));
		}
		int samples = ps[0].size();
		std::vector<int> decisions;
		std::vector<int> correct;
		decisions.resize(samples);
		correct.resize(samples);
		const reexp::data_var<P>* vars[10];

		for (int i = 0; i < 10; i++) {
			vars[i] = &data.var(i + varid::digit0);
		}
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);

		double infos[10];
		int oks[10];
		int ns[10];
		for (int i = 0; i < 10; ++i) {
			oks[i] = 0; infos[i] = 0; ns[i] = 0;
		}
		reexp::cvec<P> at;
		int ok = 0;
		double info = 0;
		for (int s = 0; s < samples; s++) {
			int digit = -1;
			// find out the right answer
			at[cvarid::sample] = s;
			for (int i = 0; i < 10; i++) {
				if (*((*vars[i])[at])) {
					correct[s] = digit = i;
					break;
				}
			}

			// normalize
			double totalP = 0;
			for (int d = 0; d < 10; d++) {
				totalP += ps[d][s];
			}
			for (int d = 0; d < 10; d++) {
				ps[d][s] /= totalP;
			}
			// find out the peak probability & calculate information
			int top = -1;
			double topP = -1;
			double i = 0;
			for (int d = 0; d < 10; d++) {
				if (ps[d][s] > topP) {
					top = d;
					topP = ps[d][s];
				}
			}
			i = -log(ps[digit][s]) / log(2);
			info += i;
			infos[digit] += i;
			ns[digit]++;
			decisions[s] = top;
			bool o = false;
			if (digit == top) {
				ok++;
				oks[digit]++;
				o = true;
			}
			if (details) {
				buf<<"s"<<((s/100)%10)<<((s/10)%10)<<(s%10)<<" d"<<correct[s]<<" ";
				if (o) {
					buf<<"==";
				} else {
					buf<<"!=";
				}
				buf<<decisions[s]<<" p="<<ps[digit][s]<<" i="<<i<<"\n";
			}
		}
		info /= samples;
		double baseInfo = -log(0.1) / log(2);

		if (details) {
			for (int i = 0; i < 10; i++) {
				infos[i] /= samples;
				buf<<"d"<<i<<" ok:"<<oks[i]<<"/"<<ns[i]<<" a:"<<(oks[i]/(double)ns[i])<<" i:"<<infos[i]<<"\n";
			}
		}

		buf<<"entropy	  "<<info<<" / "<<baseInfo<<"\n";
		buf<<"accuracy    "<<(ok/(double)samples)<<"\n";
		buf<<"error       "<<((samples-ok)/(double)samples)<<"\n";


		return buf.str();
	}

	void predict_test(test_tool& t) {
		typedef optdigits_problem p;

		int samples = 100;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(samples));

		setup_lang<p>(lang);
		populate_from_file<p>(data, samples, TRA_DATA_FILE);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::stats_info<p> si(i, stats);

		reexp::learner<p> learner(lang, stats, 50);
		setup_learner<p>(learner);

		double before = stats.naiveInfo();
		int exps = learner.reexpress();
		t<<exps<<" expressions added.\n\n";

		t<<"information "<<before<<" -> "<<stats.naiveInfo()<<"\n\n";

		t<<"scan:\n\n"<<si.scan_tostring(3)<<"\n";
		t<<"vars:\n\n"<<si.vars_tostring()<<"\n";
//		t<<"rels:\n\n"<<si.rels_tostring()<<"\n";

		t<<"predictions:\n\n"<<si.preds_tostring(varid::digit0, varid::digit9, true);

		t<<predictions_tostring(data, stats);
	}


	void detailedpredict_test(test_tool& t) {
		typedef optdigits_problem p;

		int samples = 100;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(samples));

		setup_lang<p>(lang);
		populate_from_file<p>(data, samples, TRA_DATA_FILE);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::stats_info<p> si(i, stats);

		reexp::learner<p> learner(lang, stats, 50);
		setup_learner<p>(learner);

		double before = stats.naiveInfo();
		int exps = learner.reexpress(true, 1);
		t<<exps<<" expressions added.\n\n";

		t<<"information "<<before<<" -> "<<stats.naiveInfo()<<"\n\n";

		t<<"scan:\n\n"<<si.scan_tostring(3)<<"\n";
		t<<"vars:\n\n"<<si.vars_tostring()<<"\n";
		t<<"rels:\n\n"<<si.rels_tostring()<<"\n";

//		t<<"predictions:\n\n"<<si.preds_tostring(varid::digit0, varid::digit9, true);

		reexp::pred<p> pred(stats);
		for (int i = 0; i < 10; i++) {
			typedef test_tool& ref_type;
			pred.rowP<ref_type>(data, i+varid::digit0, t);
		}

		t<<predictions_tostring(data, stats, false);
	}

	void visuals_test(test_tool& t) {
		typedef optdigits_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(DefaultSampleCount));

		setup_problem<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::stats_info<p> si(i, stats);
		const reexp::lang_info<p>& li = si.lang_info();

		reexp::learner<p> learner(lang, stats, DefaultThreshold, FilterThreshold);
		setup_learner<p>(learner);

		{
			time_sentry time;
			int exps = learner.reexpress(true);
			t<<exps<<" expressions added.\n\n";
//			printf("reexpression took %dms.\n", int(time.ms()));
		}

		t<<"exps:\n\n";
		for (int i = varid::digit9+1; i < lang.var_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.var(i))<<"\n";
		}

		t<<"predicting:\n\n"; // <<si.lang_info().drawn_vars_tostring(cvarid::x, cvarid::y);

		for (int i = 0; i < 10; ++i) {
			reexp::var<p>& var(lang.var(varid::digit0 + i));
			t<<si.lang_info().var_tostring(var)<<":\n\n";

			std::priority_queue<reexp::candidate<p>> cands;
			si.var_scan(var, cands);

			int top = 10;
			while (cands.size() && top) {
				const reexp::candidate<p>& c = cands.top();
				if (c.state_ == 3) {
					const reexp::rel<p>& rel = c.rel_->data().rel();
					const reexp::var<p>& cvar = *rel.entries()[0].var_;

					t<<"#"<<rel.id()<<" "<<c.bias_<<" "<<li.rel_tostring(rel, c.state_)<<"\n\n";
					t<<li.drawn_var_tostring(cvarid::x, cvarid::y, cvar)<<"\n";
					top--;
				}
				cands.pop();
			}
		}
		t<<predictions_tostring(data, stats);
//		t<<"rels:\n\n"<<si.lang_info().drawn_rels_tostring(cvarid::x, cvarid::y);
	}
#if 0
	template <typename P>
	void do_evaluate(test_tool& t,
					 reexp::lang<P>& lang,
					 reexp::data<P>& data,
					 reexp::stats<P>& stats,
					 reexp::learner<P>& learner) {
		int tsamples = 200;
		std::set<std::string> tags = {sup()<<"exps: "<<lang.exp_count()};
		reexp::data<P> tdata(lang, optdigits_dim(tsamples));
		populate<P>(tdata, tsamples, true);
		tdata.apply_exps();
		t.record(tags+"gen:info", stats.naiveInfo());
		evaluate(t, data, stats, tags+"data:train");
		evaluate(t, tdata, stats, tags+"data:test");
		lang.set_obs(data); // return observer
	}
#endif

	void byexps_test(test_tool& t) {
		typedef optdigits_problem p;

		int samples = TRA_DATA_FILE_SAMPLES;
				    /*+ CV_DATA_FILE_SAMPLES
				    + WDEP_DATA_FILE_SAMPLES;*/

//		int tsamples = WINDEP_DATA_FILE_SAMPLES;
		int tsamples = CV_DATA_FILE_SAMPLES;
//		int samples = 500;

		reexp::lang<p> lang;
		setup_lang<p>(lang);

		reexp::data<p> tdata(lang, optdigits_dim(tsamples));
//		populate_from_file(tdata, WINDEP_DATA_FILE_SAMPLES, WINDEP_DATA_FILE, 0);
		populate_from_file(tdata, CV_DATA_FILE_SAMPLES, CV_DATA_FILE, 0);

		reexp::data<p> data(lang, optdigits_dim(samples));
		int at = 0;
		at = populate_from_file(data, TRA_DATA_FILE_SAMPLES, TRA_DATA_FILE, at);
/*		at = populate_from_file(data, CV_DATA_FILE_SAMPLES, CV_DATA_FILE, at);
		at = populate_from_file(data, WDEP_DATA_FILE_SAMPLES, WDEP_DATA_FILE, at);*/

		reexp::stats<p> stats(data);

		double th = 350;
		reexp::learner<p> learner(lang, stats, th, 0.4*th, 0);
		setup_learner<p>(learner);

		int expsPerStep = 5;
		int exps = 0;
		do_evaluate(t, lang, data, tdata, stats, learner, cvarid::sample);

		while (exps < 150) {
			time_sentry timer;
			int added = learner.reexpress(true, expsPerStep);
			if (added) {
				long us = timer.us();
				std::set<std::string> tags = {sup()<<"exps: "<<lang.exp_count()};
				t.record(tags+"perf:reexp us", us/added);
				exps += added;
				tdata.apply_exps();
				do_evaluate(t, lang, data, tdata, stats, learner, cvarid::sample);
			}
			// do measuring here
			if (added < expsPerStep) break;
		}
		t<<exps<<" expressions added.\n\n";

		std::set<std::string> tags;
		tags.insert("run:out");

		table exptable(
			t.report(to_table<average>(tags+"data:test", "prop:", "exps:")));
		t<<"prediction quality: (test)\n\n"<<exptable<<"\n";

		table exptable3(
			t.report(to_table<average>(tags+"data:train", "prop:", "exps:")));
		t<<"prediction quality: (train)\n\n"<<exptable3<<"\n";

		table exptable4(
			t.report(to_table<average>(tags, "gen:", "exps:")));
		t<<"train sample compression: \n\n"<<exptable4<<"\n";

		table exptable2(
			t.report(to_table<average>(tags, "perf:", "exps:")));
		t.ignored()<<"performance:\n\n"<<exptable2<<"\n";

		t.ignored()<<"reexp us/exp:\n"
				   <<t.report(to_table<average>(tags, "perf:reexp us", "exps:"))
					  .to_plot(3, 20)<<"\n";

		t.ignored()<<"pred us:\n"
				   <<t.report(to_table<average>(tags+"perf:pred us", "data:", "exps:"))
					  .to_plot(3, 20)<<"\n";
		t<<"\nentropy:\n\n"<<t.report(to_table<average>(tags+"prop:entropy", "data:", "exps:"))
			.to_plot(3, 20)<<"\n";

		t<<"\naccuracy:\n\n"<<t.report(to_table<average>(tags+"prop:accuracy", "data:", "exps:"))
			.to_plot(3, 20)<<"\n";

	}

	void big_test(test_tool& t) {
		typedef optdigits_problem p;

		int samples = 945;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(samples));

		setup_lang<p>(lang);
		populate_from_file<p>(data, samples, TRA_DATA_FILE);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::learner<p> learner(lang, stats, 100);
		setup_learner<p>(learner);

		double before = stats.naiveInfo();
		int exps = learner.reexpress(true);
		t<<exps<<" expressions added.\n\n";
		t<<"information "<<before<<" -> "<<stats.naiveInfo()<<"\n\n";

		reexp::stats_info<p> si(i, stats);
		const reexp::lang_info<p>& li = si.lang_info();

		t<<"exps:\n\n";
		for (int i = varid::digit9+1; i < lang.var_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.var(i))<<"\n";
		}

		t<<"predicting:\n\n"; // <<si.lang_info().drawn_vars_tostring(cvarid::x, cvarid::y);

		for (int i = 0; i < 10; ++i) {
			reexp::var<p>& var(lang.var(varid::digit0 + i));
			t<<si.lang_info().var_tostring(var)<<":\n\n";

			std::priority_queue<reexp::candidate<p>> cands;
			si.var_scan(var, cands);

			int top = 10;
			while (cands.size() && top) {
				const reexp::candidate<p>& c = cands.top();
				if (c.state_ == 3) {
					const reexp::rel<p>& rel = c.rel_->data().rel();
					const reexp::var<p>& cvar = *rel.entries()[0].var_;

					t<<"#"<<rel.id()<<" "<<c.bias_<<" "<<li.rel_tostring(rel, c.state_)<<"\n\n";
					t<<li.drawn_var_tostring(cvarid::x, cvarid::y, cvar)<<"\n";
					top--;
				}
				cands.pop();
			}
		}

		t<<si.drawn_data_tostring(reexp::cvec<p>(), cvarid::x, cvarid::y, cvarid::sample, 40)<<"\n";

		t<<predictions_tostring(data, stats, true);

		int tsamples = 200;
		t<<"prepared test data... ";
		reexp::data<p> tdata(lang, optdigits_dim(tsamples));
		t<<"ok.\npopulating it... ";
		populate_from_file<p>(data, tsamples, CV_DATA_FILE);
		t<<"ok.\nre-expressing it... ";
		tdata.apply_exps(); // this is behind the problem (!!!)
		t<<"ok.\nmaking predictions... \n";
		t<<predictions_tostring(tdata, stats, false);
		t<<"done.\n";
	}

	void measure_test(test_tool& t) {
		int samples = DefaultSampleCount;

		typedef optdigits_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, optdigits_dim(samples));

		setup_lang<p>(lang);
		populate_from_file<p>(data, samples, TRA_DATA_FILE);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);

		reexp::stats_info<p> si(i, stats);
		reexp::learner<p> learner(lang, stats, 65);
		setup_learner<p>(learner);

		t<<predictions_tostring(data, stats, false);

		double before = stats.naiveInfo();
		for (int i = 0; i < 5; i++) {
			int exps = learner.reexpress(true, 5);
			t<<exps<<" expressions added.\n\n";
			t<<predictions_tostring(data, stats, false);
		}

		t<<"information "<<before<<" -> "<<stats.naiveInfo()<<"\n\n";

		t<<"exps:\n\n";
		for (int i = varid::digit9+1; i < lang.var_count(); ++i) {
			t<<si.lang_info().drawn_var_tostring(cvarid::x, cvarid::y, lang.var(i))<<"\n";
		}

		t<<predictions_tostring(data, stats, true);

	}
}


void addoptdigitstest(TestRunner& runner) {
	runner.add("optdigits/byexps", 			{"real"},   &byexps_test);
	runner.add("optdigits/big", 			{"real"},   &big_test);
	runner.add("optdigits/measure", 		{"func"},   &measure_test);
	runner.add("optdigits/visuals",			{"func"},   &visuals_test);
	runner.add("optdigits/predict", 		{"func"},   &predict_test);
	runner.add("optdigits/detailedpredict", {"func"},   &detailedpredict_test);
	runner.add("optdigits/learning",		{"func"},   &learning_test);
	runner.add("optdigits/setup", 			{"func"},   &setup_test);
//	big_test(ok);
/*	measure_test(ok);
	visuals_test(ok);
	predict_test(ok);
	learning_test(ok);
	setup_test(ok);*/
}
