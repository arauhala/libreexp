/*
 * optdigits2test.cpp

 *
 *  Created on: Sep 11, 2012
 *      Author: arau
 */


#include "reexp/all.h"

#include "tester.h"
#include "exptesttools.h"

#include <stdio.h>

#include <iostream>
#include <algorithm>

#include "optdigits.h"

// bit variables (used to group bits)
namespace optdigits {
	namespace varid {
		static const int first_pixel = 10;
	}
}

namespace {

	using namespace optdigits;

	//
	// Each data entry is 7*3 + 9 = 30 bits
	//
	static const int Pack = 2;

	static const int Width = OriginalWidth / Pack;

	static const int Height = OriginalHeight / Pack;

	static const int PixelBits = 1;

	static const int PixelBitBound[] = {
		2,
	};


//	static int PixelCount = Width * Height;

	static const int DigitCount = 10;

	namespace cvarid {
		static const int sample = 0;
	}


	const char* dig_varnames[] =  {
		"d0",
		"d1",
		"d2",
		"d3",
		"d4",
		"d5",
		"d6",
		"d7",
		"d8",
		"d9"
	};
	const char* cvarnames[] =  {
		"s" // sample
	};


	struct optdigits_problem {
		static const int DIM = 1;
		static const int MAX_REL_VARS = 2;
	};

	int pixel_varid(int x, int y, int bit) {
		return varid::first_pixel + PixelBits * (x + y*Width) + bit;
	}

	explib::cvec<optdigits_problem> optdigits_dim(int samples) {
		return explib::cvec<optdigits_problem>(samples);
	}

	template <typename P>
	int populate(explib::data<P>& data, int samples, std::ifstream& in, int offset = 0) {
		std::string line;

		for (int i = 0; i < DATA_FILE_HEADER_LINES; ++i) {
			std::getline(in, line);
		}

		explib::bitmap pic;

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

			explib::cvec<P> at(offset+i);
			for (int x = 0; x < Width; x++) {
				for (int y = 0; y < Height; y++) {
					int popcount = 0;
					for (int x2 = 0; x2 < Pack; x2++) {
						for (int y2 = 0; y2 < Pack; y2++) {
							if (pic.get(x*Pack + x2, y*Pack + y2)) popcount++;
						}
					}
					for (int j = 0; j < PixelBits; ++j) {
						explib::data_var<P>& p = data.var(pixel_varid(x, y, j));
						p[at] = true;
						//*p[at] = (popcount*2 >= Pack*Pack);
						*p[at] = (popcount >= PixelBitBound[j]);
					}
				}
			}

			for (int j = 0; j < DigitCount; j++) {
				explib::data_var<P>& digit = data.var(varid::digit0 + j);
				digit[at] = true;
				*digit[at] = (number == j);
			}
		}
		return offset + samples;
	}

	template <typename P>
	int populate_from_file(explib::data<P>& data, int samples, const char* name, int offset) {
		std::ifstream in(name);
		populate<P>(data, samples, in, offset);
		return offset + samples;
	}

	template <typename P>
	void setup_lang(explib::lang<P>& lang) {
		explib::ctx<P> ctx(explib::cvec<P>(0));

		for (int i = 0; i < DigitCount; ++i) {
			lang.add_orig(explib::orig<P>(ctx));
		}
		for (int x = 0; x < Width; ++x) {
			for (int y = 0; y < Height; ++y) {
				for (int b = 0; b < PixelBits; ++b) {
					lang.add_orig(explib::orig<P>(ctx));
				}
			}
		}

		for (int x = 0; x < Width; ++x) {
			for (int y = 0; y < Height; ++y) {
				for (int b = 0; b < PixelBits; ++b) {
					explib::var<P>& pixel = lang.var(pixel_varid(x, y, b));
					for (int b2 = 0; b2 <= b; ++b2) {
						if (b2 != b) {
							explib::var<P>& pixel2 = lang.var(pixel_varid(x, y, b2));
							explib::rel<P>& rl(lang.alloc_rel(ctx)); // left right
							rl.add_var(explib::cvec<P>(0), pixel);
							rl.add_var(explib::cvec<P>(0), pixel2);
							lang.rel_done();
						}
						if (x+1 <Width) {
							explib::var<P>& rightpixel = lang.var(pixel_varid(x+1, y, b2));
							explib::rel<P>& rl(lang.alloc_rel(ctx)); // left right
							rl.add_var(explib::cvec<P>(0), pixel);
							rl.add_var(explib::cvec<P>(0), rightpixel);
							lang.rel_done();
							if (b2 != b) {
								explib::var<P>& leftpixel2 = lang.var(pixel_varid(x, y, b2));
								explib::var<P>& rightpixel2 = lang.var(pixel_varid(x+1, y, b));
								explib::rel<P>& rl2(lang.alloc_rel(ctx)); // left right
								rl2.add_var(explib::cvec<P>(0), leftpixel2);
								rl2.add_var(explib::cvec<P>(0), rightpixel2);
								lang.rel_done();
							}
						}
						if (y+1 < Height) {
							explib::var<P>& downpixel = lang.var(pixel_varid(x, y+1, b2));
							explib::rel<P>& ud(lang.alloc_rel(ctx)); // up down
							ud.add_var(explib::cvec<P>(0), pixel);
							ud.add_var(explib::cvec<P>(0), downpixel);
							lang.rel_done();
							if (b2 != b) {
								explib::var<P>& downpixel2 = lang.var(pixel_varid(x, y, b2));
								explib::var<P>& uppixel2 = lang.var(pixel_varid(x, y+1, b));
								explib::rel<P>& ud2(lang.alloc_rel(ctx)); // left right
								ud2.add_var(explib::cvec<P>(0), downpixel2);
								ud2.add_var(explib::cvec<P>(0), uppixel2);
								lang.rel_done();
							}
						}
						if (x+1 < Width && y+1 < Height) {
							explib::var<P>& rightup = lang.var(pixel_varid(x+1, y, b));
							explib::var<P>& leftdown = lang.var(pixel_varid(x, y+1, b2));
							explib::rel<P>& r1(lang.alloc_rel(ctx)); // left right
							r1.add_var(explib::cvec<P>(0), rightup);
							r1.add_var(explib::cvec<P>(0), leftdown);
							lang.rel_done();

							explib::var<P>& rightdown = lang.var(pixel_varid(x+1, y+1, b));
							explib::var<P>& leftup = lang.var(pixel_varid(x, y, b2));
							explib::rel<P>& r2(lang.alloc_rel(ctx)); // left right
							r2.add_var(explib::cvec<P>(0), rightdown);
							r2.add_var(explib::cvec<P>(0), leftup);
							lang.rel_done();
						}
					}
					for (int i = 0; i < DigitCount; ++i) {
						explib::rel<P>& ps(lang.alloc_rel(ctx)); // pixel-shape
						ps.add_var(explib::cvec<P>(0), // pixel
								   pixel);
						ps.add_var(explib::cvec<P>(0), // shape
								   lang.var(varid::digit0 + i));
						lang.rel_done();
					}
				}
			}
		}
	}

	template <typename P>
	void setup_names(explib::pinfo& info) {
		for (int i = 0; i < varid::first_pixel; ++i) {
			info.vnames_.push_back(dig_varnames[i]);
		}
		for (int x = 0; x < Width; ++x) {
			for (int y = 0; y < Height; ++y) {
				info.vnames_.push_back(sup()<<"pixel"<<x<<","<<y);
			}
		}
		for (int i = 0; i < P::DIM; ++i) {
			info.cvnames_.push_back(cvarnames[i]);
		}
	}

	template <typename P>
	void setup_learner(explib::learner<P>& learner) {
		for (int i = varid::digit0; i <= varid::digit9; ++i) {
			learner.exclude(i);
		}
	}

	void setup_test(TestTool& t) {
		typedef optdigits_problem p;
		int samples = 1000;

		explib::lang<p> lang;
		explib::data<p> data(lang, optdigits_dim(samples));

		setup_lang<p>(lang);
		std::ifstream in(TRA_DATA_FILE);
		populate(data, samples, in);

		explib::pinfo i;
		setup_names<p>(i);

		explib::lang_info<p> li(i, lang);

		t<<"vars:\n\n";
		t<<li.vars_tostring()<<"\n\n";
		t<<"rels:\n\n";
		t<<li.rels_tostring()<<"\n\n";
	}


	void learning_test(TestTool& t) {
		typedef optdigits_problem p;
		int samples = 1934;

		explib::lang<p> lang;
		explib::data<p> data(lang, optdigits_dim(samples));

		setup_lang<p>(lang);
		std::ifstream in(TRA_DATA_FILE);
		populate(data, samples, in);

		explib::stats<p> stats(data);

		explib::pinfo i;
		setup_names<p>(i);

		explib::stats_info<p> si(i, stats);

		t<<"scan:\n\n"<<si.scan_tostring(3, 2)<<"\n";

		double threshold = 200;
		explib::learner<p> learner(lang, stats, threshold, 0.4*threshold, 5);
		setup_learner<p>(learner);

		int exps = learner.reexpress();
		/*
		{
			TimeSentry time;
			while (true) {
				float before = stats.naiveInfo();
				if (!learner.add_exp()) break;
		//		t<<"expression added.\n";
				float after = stats.naiveInfo();
				t<<"info "<<before<<" -> "<<after<<"\n\n";

				t<<"scan:\n\n"<<si.scan_tostring(3, 2)<<"\n";
			}
///			printf("\nreexpression took %dms.\n", int(time.ms()));
		}*/

		t<<"\nlearning done.\n\n";

		t<<exps<<" exps:\n\n";

		for (int i = lang.orig_count(); i < lang.var_count(); ++i) {
			t<<si.var_tostring(i);
		}

/*		t<<"scan:\n\n"<<si.scan_tostring(3, 2)<<"\n";

		t<<"vars:\n\n"<<si.vars_tostring()<<"\n";
		t<<"rels:\n\n"<<si.rels_tostring()<<"\n";*/
	}


	void byexps_test(TestTool& t) {
		typedef optdigits_problem p;

		int samples = TRA_DATA_FILE_SAMPLES+CV_DATA_FILE_SAMPLES+WDEP_DATA_FILE_SAMPLES;
		int tsamples = WINDEP_DATA_FILE_SAMPLES;
		static const double rel_filter = 0;

		explib::lang<p> lang;
		setup_lang<p>(lang);

		explib::data<p> tdata(lang, optdigits_dim(tsamples));
		populate_from_file<p>(tdata, WINDEP_DATA_FILE_SAMPLES, WINDEP_DATA_FILE, 0);

		explib::data<p> data(lang, optdigits_dim(samples));
		int at = 0;
		at = populate_from_file<p>(data, TRA_DATA_FILE_SAMPLES, TRA_DATA_FILE, at);
		at = populate_from_file<p>(data, CV_DATA_FILE_SAMPLES, CV_DATA_FILE, at);
		at = populate_from_file<p>(data, WDEP_DATA_FILE_SAMPLES, WDEP_DATA_FILE, at);

		explib::stats<p> stats(data);

		double th = 30;
		explib::learner<p> learner(lang, stats, th, 0.4*th, rel_filter);
		setup_learner<p>(learner);

		int expsPerStep = 50;
		int exps = 0;
		do_evaluate(t, lang, data, tdata, stats, learner, cvarid::sample);

		while (exps < 1000) {
			TimeSentry timer;
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

		Table exptable(
			t.report(ToTable<Average>(tags+"data:test", "prop:", "exps:")));
		t<<"prediction quality: (test)\n\n"<<exptable<<"\n";

		Table exptable3(
			t.report(ToTable<Average>(tags+"data:train", "prop:", "exps:")));
		t<<"prediction quality: (train)\n\n"<<exptable3<<"\n";

		Table exptable4(
			t.report(ToTable<Average>(tags, "gen:", "exps:")));
		t<<"train sample compression: \n\n"<<exptable4<<"\n";

		Table exptable2(
			t.report(ToTable<Average>(tags, "perf:", "exps:")));
		t.ignored()<<"performance:\n\n"<<exptable2<<"\n";

		t.ignored()<<"reexp us/exp:\n"
				   <<t.report(ToTable<Average>(tags, "perf:reexp us", "exps:"))
					  .toplot(3, 20)<<"\n";

		t.ignored()<<"pred us:\n"
				   <<t.report(ToTable<Average>(tags+"perf:pred us", "data:", "exps:"))
					  .toplot(3, 20)<<"\n";
		t<<"\nentropy:\n\n"<<t.report(ToTable<Average>(tags+"prop:entropy", "data:", "exps:"))
			.toplot(3, 20)<<"\n";

/*		t<<"\naccuracy:\n\n"<<t.report(ToTable<Average>(tags+"prop:accuracy", "data:", "exps:"))
			.toplot(3, 20)<<"\n";*/

		t<<"\nerror:\n\n"<<t.report(ToTable<Average>(tags+"prop:error", "data:", "exps:"))
			.toplot(3, 20)<<"\n";

	}


}

// Unlike in optdigits test suite, optdigits2 uses separate variable for
// each pixel
void addoptdigits2test(TestRunner& runner) {
	runner.add("optdigits/o2_setup", 			{"func"},   &setup_test);
	runner.add("optdigits/o2_learning", 		{"func"},   &learning_test);
	runner.add("optdigits/o2_exps", 			{"func"},  	&byexps_test);
}