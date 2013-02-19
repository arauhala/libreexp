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

namespace {

	using namespace std;

	//
	// Each data entry is 7*3 + 9 = 30 bits
	//
	static int Width = 3;

	static int Height = 5;

	static int SampleCount = 10;

	static int DigitCount = 10;

	static int VarCount = 11;

	// context variables
	namespace cvarid {
		static const int x = 0;
		static const int y = 1;
		static const int sample = 2;
	}
	// bit variables (used to group bits)
	namespace varid {
		static const int pixel = 0;
		static const int digit0 = 1;
		static const int digit1 = 2;
		static const int digit2 = 3;
		static const int digit3 = 4;
		static const int digit4 = 5;
		static const int digit5 = 6;
		static const int digit6 = 7;
		static const int digit7 = 8;
		static const int digit8 = 9;
		static const int digit9 = 10;
	}
	// relations (used to organize bits)
	namespace relid {
		static const int left_right = 0;
		static const int up_down = 1;
	}
	const char* varnames[] =  {
		"px",
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
		"x",
		"y",
		"n" // sample
	};

	struct recognition_problem {
		static const int DIM = 3;
		static const int MAX_REL_VARS = 2;
	};

	reexp::cvec<recognition_problem> sample_dim() {
		reexp::cvec<recognition_problem> dim;
		dim[cvarid::x] = Width;
		dim[cvarid::y] = Height;
		dim[cvarid::sample] = SampleCount;
		return dim;
	}


	//
	// 3x7
	//

	const char* sample0[] = {
		"XXX",
		"X.X",
		"X.X",
		"X.X",
		"XXX"
	};

	const char* sample1[] = {
		".X.",
		".X.",
		".X.",
		".X.",
		".X."
	};
	const char* sample2[] = {
		"XXX",
		"..X",
		"XXX",
		"X..",
		"XXX"
	};
	const char* sample3[] = {
		"XXX",
		"..X",
		"XXX",
		"..X",
		"XXX"
	};
	const char* sample4[] = {
		"X.X",
		"X.X",
		"XXX",
		"..X",
		"..X"
	};

	const char* sample5[] = {
		"XXX",
		"X..",
		"XXX",
		"..X",
		"XXX"
	};

	const char* sample6[] = {
		"XXX",
		"X..",
		"XXX",
		"X.X",
		"XXX"
	};

	const char* sample7[] = {
		"XXX",
		"..X",
		"..X",
		"..X",
		"..X"
	};

	const char* sample8[] = {
		"XXX",
		"X.X",
		"XXX",
		"X.X",
		"XXX"
	};

	const char* sample9[] = {
		"XXX",
		"X.X",
		"XXX",
		"..X",
		"..X"
	};

	struct sample {
		int number_;
		const char** bitmap_;
	};

	sample samples[] = {
		{0, sample0},
		{1, sample1},
		{2, sample2},
		{3, sample3},
		{4, sample4},
		{5, sample5},
		{6, sample6},
		{7, sample7},
		{8, sample8},
		{9, sample9}
	};

	int xyoffset(int x, int y) {
		return y*3 + x;
	}

	template <typename P>
	void populate(reexp::data<P>& data) {
		reexp::data_var<P>& p = data.var(varid::pixel);
		reexp::cvec<P> at(0, 0, 0);

		for (int i = 0; i < SampleCount; i++) {
			sample& sa = samples[i];
			at[cvarid::sample] = i;
			for (int j = 0; j < DigitCount; j++) {
				reexp::data_var<P>& digit = data.var(varid::digit0 + j);
				digit[at] = true;
				*digit[at] = (sa.number_ == j);
			}
			for (int j = 0; j < Height; j++) {
				at[cvarid::y] = j;
				for (int k = 0; k < Width; k++) {
					at[cvarid::x] = k;
					p[at] = true;
					*p[at] = sa.bitmap_[j][k] == 'X';
				}
			}
		}
	}

	template <typename P>
	void setup_lang(reexp::lang<P>& lang) {
		reexp::ctx<P> pixel_ctx(reexp::cvec<P>(0, 0, 0));
		reexp::ctx<P> shape_ctx(reexp::cvec<P>(-1, -1, 0));
		reexp::ctx<P> full_ctx(reexp::cvec<P>(0, 0, 0));

		lang.add_orig(reexp::orig<P>(pixel_ctx));
		for (int i = 0; i < DigitCount; ++i) {
			lang.add_orig(reexp::orig<P>(shape_ctx));
		}

		reexp::rel<P>& rl(lang.alloc_rel(pixel_ctx)); // right left
		rl.add_var(reexp::cvec<P>(0, 0, 0, 0),
				   lang.var(varid::pixel));
		rl.add_var(reexp::cvec<P>(1, 0, 0, 0),
				  lang.var(varid::pixel));
		lang.rel_done();

		reexp::rel<P>& ud(lang.alloc_rel(pixel_ctx)); // up down
		ud.add_var(reexp::cvec<P>(0, 0, 0, 0),
				   lang.var(varid::pixel));
		ud.add_var(reexp::cvec<P>(0, 1, 0, 0),
				   lang.var(varid::pixel));
		lang.rel_done();

		for (int i = 0; i < DigitCount; ++i) {
			reexp::rel<P>& ps(lang.alloc_rel(full_ctx)); // pixel-shape
			ps.add_var(reexp::cvec<P>(0, 0, 0, 0), // pixel
					   lang.var(varid::pixel));
			ps.add_var(reexp::cvec<P>(0, 0, 0, 0), // shape
					   lang.var(varid::digit0 + i));
			lang.rel_done();
		}
	}

	template <typename P>
	void setup_reg(reexp::lang<P>& lang, reexp::data<P>& data) {
		setup_lang(lang);
		populate(data);
	}

	void setup_pinfo(reexp::pinfo& info) {
		for (int i = 0; i < VarCount; ++i) {
			info.vnames_.push_back(varnames[i]);

		}
		for (int i = 0; i < recognition_problem::DIM; ++i) {
			info.cvnames_.push_back(cvarnames[i]);

		}
	}

	template <typename P>
	void setup_learner(reexp::learner<P>& learner) {
		for (int i = varid::digit0; i <= varid::digit9; ++i) {
			learner.exclude(i);
		}
	}

	template <typename P>
	void print_pixels(test_tool& t, const reexp::data_var<P>& pixels) {
		reexp::cvec<P> d = pixels.dim();
		for (int k = 0; k < d[cvarid::sample]; k++) {
			for (int i = 0; i < d[cvarid::y]; i++) {
				for (int j = 0; j < d[cvarid::x]; j++) {
					reexp::cvec<P> at(j, i, k, 0);
					auto b = pixels[at];
					if (b) {
						if (*b) {
							t<<"X";
						} else {
							t<<".";
						}
					} else {
						t<<"?";
					}
				}
				t<<checkl;
			}
			t<<checkl;
		}
	}

	template <typename P>
	void print_offsets(test_tool& t, reexp::data<P>& data, int var, int x, int y, int sample) {
		reexp::cvec<P> at(x, y, sample);
		reexp::ctx<P> ctx( data.lang().var(var).ctx() );
		reexp::cvec<P> dim = data.dim();
		const reexp::data_var<P>& dvar = data.var(var);
		t<<"var id "<<var<<" dim "<<ctx.dim(dim)<<checkl;
		t<<"offset a for "<<at<<" = "<<ctx.offset(dim, at)<<checkl;
		t<<"offset b for "<<at<<" = "<<ctx.dim(dim).offset(at)<<checkl;
		t<<"index for    "<<at<<" = "<<dvar.index(at)<<checkl;
		t<<"def for      "<<at<<" = "<<bool(dvar[at])<<checkl;
		t<<"value for    "<<at<<" = "<<bool(*dvar[at])<<checkl;
		t<<checkl;
	}

	void run_ctx_test(bool& ok) {
		test_tool t("test/digits/ctx", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;

		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::var<p>& pixel( lang.var(varid::pixel) );
		reexp::var<p>& shape( lang.var(varid::digit0) );

		t<<"shape has "<<shape.ctx().dim(data.dim()).volume()<<" bits."<<checkl<<checkl;

		print_offsets<p>(t, data, varid::digit0, 0, 0, 0);
		print_offsets<p>(t, data, varid::digit0, 1, 0, 0);
		print_offsets<p>(t, data, varid::digit0, 0, 1, 0);
		print_offsets<p>(t, data, varid::digit0, 0, 0, 1);
		print_offsets<p>(t, data, varid::digit1, 0, 0, 0);
		print_offsets<p>(t, data, varid::digit0, 0, 0, 9);
		print_offsets<p>(t, data, varid::digit9, 0, 0, 0);
		print_offsets<p>(t, data, varid::digit9, 0, 0, 9);
		print_offsets<p>(t, data, varid::digit9, 2, 4, 9);

		t<<"pixel has "<<pixel.ctx().dim(data.dim()).volume()<<" bits."<<checkl<<checkl;

		print_offsets<p>(t, data, varid::pixel, 0, 0, 0);
		print_offsets<p>(t, data, varid::pixel, 0, 0, 0);
		print_offsets<p>(t, data, varid::pixel, 0, 0, 1);
		print_offsets<p>(t, data, varid::pixel, 0, 1, 0);
		print_offsets<p>(t, data, varid::pixel, 1, 0, 0);

		print_offsets<p>(t, data, varid::pixel, 2, 0, 0);
		print_offsets<p>(t, data, varid::pixel, 0, 4, 0);
		print_offsets<p>(t, data, varid::pixel, 0, 0, 9);
		print_offsets<p>(t, data, varid::pixel, 0, 0, 0);
		print_offsets<p>(t, data, varid::pixel, 2, 4, 9);
	}

	void print_bits(test_tool& t, reexp::cond_bits& bits, int offset, int n) {
		for (int i = 0; i < n; i++) {
			if (!(i % 10)&&i) {
				t<<checkl;
			}
			auto b = bits[offset+i];
			if (b) {
				if (*b) {
					t<<"1";
				} else {
					t<<"0";
				}
			} else {
				t<<"?";
			}
		}
		t<<checkl;
	}

	void run_setup_test(bool& ok) {
		test_tool t("test/digits/setup", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::data_var<p>& pixels(data.var(varid::pixel));
		int psz = pixels.dim().volume();
		t<<"pixels, "<<psz<<" bits"<<checkl;
		print_bits(t, pixels.bits(), 0, pixels.bits().size());
		t<<checkl;

		for (int i = 0; i < DigitCount; ++i) {
			reexp::data_var<p>& digit(data.var(varid::digit0+i));
			int dsz = digit.dim().volume();
			t<<"digit"<<i<<", "<<dsz<<" bits"<<checkl;
			print_bits(t, digit.bits(), 0, digit.bits().size());
			t<<checkl;
		}
	}

	template <typename P>
	void print_var_name(test_tool& t, const reexp::var<P>& var) {
		switch (var.id()) {
			case 0: t<<"px"; break;
			case varid::digit0:
			case varid::digit1:
			case varid::digit2:
			case varid::digit3:
			case varid::digit4:
			case varid::digit5:
			case varid::digit6:
			case varid::digit7:
			case varid::digit8:
			case varid::digit9:
				t<<"d"<<(var.id()-varid::digit0); break;
			default: t<<"e"<<(var.id()-varid::digit9-1); break;
		}
	}

	template <typename P>
	void print_var(test_tool& t, const reexp::var<P>& var) {
		switch (var.id()) {
			case 0: t<<"px"<<checkl; break;
			case varid::digit0:
			case varid::digit1:
			case varid::digit2:
			case varid::digit3:
			case varid::digit4:
			case varid::digit5:
			case varid::digit6:
			case varid::digit7:
			case varid::digit8:
			case varid::digit9:
				t<<"d"<<(var.id()-1)<<checkl; break;
			default: {
				const reexp::exp<P>& v = dynamic_cast<const reexp::exp<P>&>(var);
				print_rel(t, v.rel(), v.state());
				break;
			}
		}
	}



	template <typename P>
	void print_vars(test_tool& t, const reexp::stats<P>& stats) {
		const reexp::lang<P>& lang = stats.data().lang();
		double info = 0;
		for (int i = 0; i < lang.var_count(); i++) {
			const reexp::var_stats<P>& s( stats.var( i ) );
			print_var(t, lang.var(i));
			info += s.information();
			t<<s.freq()<<"/"<<s.n()<<" "<<s.information()<<checkl<<checkl;
		}
		t<<"total info: "<<info<<checkl<<checkl;
	}

	template <typename P>
	void print_rel(test_tool& t, const reexp::rel<P>& rel, int state = -1) {
		int varn = rel.entries().size();
		t<<"[";
		for (int j = 0; j < varn; j++) {
			const reexp::rel_entry<P>& e = rel.entries()[j];
			const reexp::var<P>& var( *e.var_ );
			const reexp::cvec<P>& cv( e.shift_ );
			if (state != -1 && !rel.varState(j, state)) t<<"!";
			print_var_name(t, var);
			t<<"(r";
			for (int k = 0; k < P::DIM; k++) {
				if (cv[k] == 1) {
					t<<" + "<<cvarnames[k];
				} else if (cv[k] != 0) {
					t<<" + "<<cv[k]<<cvarnames[k];
				}
			}
			t<<")";
			if (j + 1 < varn) t<<(", ");
		}
		t<<"] ";
		t<<checkl;
	}

	template <typename P>
	void print_rels(test_tool& t, reexp::stats<P>& st) {
		const reexp::lang<P>& lang = st.data().lang();
		for (int i = 0; i < lang.rel_count(); i++) {
			const reexp::rel<P>& rel( lang.rel(i) );
			const reexp::rel_stats<P>& s( st.rel( i ) );

			t<<i<<" ";
			print_rel(t, rel);
			t<<"    ";
			for (int j = 0; j < rel.varCount(); j++) {
				t<<s.varFreqs()[j]<<"/"<<s.n()<<" ";
			}
			t<<checkl<<"    ";
			for (int j = 0; j < rel.stateCount(); j++) {
				t<<s.stateFreqs()[j]<<" ";
			}
			t<<checkl<<"    ";
			for (int j = 0; j < rel.stateCount(); j++) {
				t<<(s.stateNaiveP(j)*s.n())<<" ";
			}
			t<<"   (expected)"<<checkl<<"    ";
			for (int j = 0; j < rel.stateCount(); j++) {
				t<<s.stateBias(j)<<" ";
			}
			t<<checkl;
		}
	}

	void print_multiline(test_tool& t, int pad, const char* txt) {
		std::string p;
		for (int i = 0; i < pad; i++) p += " ";

		t<<p;
		for (size_t i = 0; txt[i]; i++) {
			if (txt[i] == '\n') {
				t<<checkl<<p;
			} else {
				t<<txt[i];
			}
		}
		t<<checkl;
	}


	void run_exp_stats_test(bool& ok) {
		test_tool t("test/digits/exp_stats", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		t<<"before:"<<checkl;

		print_vars(t, stats);
		print_rels(t, stats);

		t<<"naive information "<<stats.naiveInfo()<<checkl;
		t<<"bit size: "<<data.size()<<checkl;
		lang.add_exp(lang.rel(relid::left_right), 3);

		t<<checkl;
		t<<"exp added."<<checkl;

		stats.update();

		t<<checkl;
		t<<"after:"<<checkl;

		print_vars(t, stats);
		print_rels(t, stats);

		t<<"naive information "<<stats.naiveInfo()<<checkl;
		t<<"bit size: "<<data.size()<<checkl;
	}

	void run_add_exp_test(bool& ok) {
		test_tool t("test/digits/add_exp", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		lang.add_exp(lang.rel(relid::left_right), 3);
		t<<"exp added."<<checkl<<checkl;

		print_pixels(t, data.var(varid::pixel));

		const reexp::exp<p>& exp = lang.exp_back();

		t<<"exp:"<<checkl<<checkl;
		print_pixels(t, data.var(exp.id()));
	}

	void run_gen_rel_test(bool& ok) {
		test_tool t("test/digits/gen_rel", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::pinfo info;
		setup_pinfo(info);
		reexp::lang_info<p> li(info, lang);

		t<<"vars:\n\n"<<li.drawn_vars_tostring(cvarid::x, cvarid::y);

		lang.add_exp(lang.rel(relid::up_down), 3);
		t<<"exp added."<<checkl<<checkl;

		t<<"vars:\n\n"<<li.drawn_vars_tostring(cvarid::x, cvarid::y);
		t<<"rels:\n\n"<<li.drawn_rels_tostring(cvarid::x, cvarid::y);
	}

	void run_print_test(bool& ok) {
		test_tool t("test/digits/print", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		print_pixels(t, data.var(varid::pixel));

	}

	void run_stats_test(bool& ok) {
		test_tool t("test/digits/stats", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		print_vars(t, stats);
		print_rels(t, stats);

		t<<"naive information "<<stats.naiveInfo()<<checkl;
		t<<"bit size: "<<data.size()<<checkl;
	}

	template <typename P>
	void print_scan(test_tool& t, reexp::lang<P>& lang, reexp::stats<P>& stats) {
		reexp::learner<P> learner(lang, stats);
		setup_learner(learner);

		std::priority_queue<reexp::candidate<P> > cands;
		learner.scan(cands);

		t<<"top 10 scan results:"<<checkl<<checkl;

		int n = 10;

		while (!cands.empty() && n--) {
			const reexp::candidate<P>& c( cands.top() );
			t<<c.bias_<<"   ";
			print_rel(t, c.rel_->data().rel_, c.state_);

			cands.pop();
		}
		t<<checkl;
	}

	void run_scan_test(bool& ok) {
		test_tool t("test/digits/scan", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		print_rels(t, stats);

		print_scan(t, lang, stats);
	}

	template <typename P>
	void print_var_deps(test_tool& t, reexp::var<P>& v) {
		t<<"deps for ";
		print_var(t, v);

		for (auto i = v.deps().begin(); i != v.deps().end(); ++i) {
			t<<"  "<<i->shift_<<" ";
			switch (i->var_->id()) {
				case 0: t<<"n"<<checkl; break;
				case 1: t<<"px"<<checkl; break;
				default: t<<"e"<<(i->var_->id()-2)<<checkl; break;
			}
		}
	}


	void run_learning_once_test(bool& ok) {
		test_tool t("test/digits/learning_once", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		print_vars(t, stats);

		reexp::learner<p> learner(lang, stats);

		learner.add_exp();

		t<<"expression added."<<checkl;

		print_vars(t, stats);
		print_rels(t, stats);

		reexp::pinfo info;
		setup_pinfo(info);
		reexp::lang_info<p> li(info, lang);

		t<<"rels:\n\n"<<li.drawn_rels_tostring(cvarid::x, cvarid::y);

		print_pixels(t, data.var(1));
		print_pixels(t, data.var(2));

		print_var_deps(t, lang.var(1));
		print_var_deps(t, lang.var(2));

	}


	void run_learning_test(bool& ok) {
		test_tool t("test/digits/learning", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		print_vars(t, stats);

		reexp::learner<p> learner(lang, stats);

		int rounds = 7;

		for (int i = 0; i < rounds; i++) {
			float before = stats.naiveInfo();
			if (!learner.add_exp()) break;
			t<<"expression added."<<checkl;
			float after = stats.naiveInfo();
			t<<"info "<<before<<" -> "<<after<<checkl<<checkl;
		}

		t<<checkl<<"learning done."<<checkl<<checkl;

		reexp::pinfo info;
		setup_pinfo(info);
		reexp::lang_info<p> li(info, lang);

		t<<"vars:\n\n"<<li.drawn_vars_tostring(cvarid::x, cvarid::y);
		print_vars(t, stats);

		t<<"rels:\n\n"<<li.drawn_rels_tostring(cvarid::x, cvarid::y);

/*		print_rels(t, stats);*/

		t<<"px"<<checkl<<checkl;
		print_pixels(t, data.var(1));

		for (int i = 0; i < rounds; i++) {
			t<<"e"<<i<<checkl<<checkl;
			print_pixels(t, data.var(i+2));
		}

		print_var_deps(t, lang.var(1));
		print_var_deps(t, lang.var(2));

		std::priority_queue<reexp::candidate<p> > cands;
		learner.scan(cands);

		while (!cands.empty()) {
			const reexp::candidate<p>& c( cands.top() );
			t<<c.bias_<<"   ";
			print_rel(t, c.rel_->data().rel_, c.state_);

			const_cast<reexp::rel_stats<p>&>(*c.rel_).update(stats);

			cands.pop();
		}

	}

	void run_bitset_test(bool& ok) {
		test_tool t("test/digits/bitset", ok);

		typedef recognition_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, sample_dim());
		setup_reg<p>(lang, data);
	}
}

template <typename P>
void print_prediction(test_tool& t, reexp::stats<P>& stats) {
	reexp::pred<P> pred(stats);

	print_vars<P>(t, stats);
	print_rels<P>(t, stats);

	for (int i = 0; i < DigitCount; ++i) {
		std::vector<double> pr( pred.p( stats.data(), varid::digit0 + i ) );
		reexp::cvec<P> dim( stats.data().var( varid::digit0 + i ).dim() );

		for (int j = 0; j < SampleCount; ++j) {
			reexp::cvec<P> at(0, 0, j);
			int idx = dim.offset(at);
			double sp = pr[idx];
			printf("p(s%d=='%d') = %d%%\n",j, i, int(100*sp));
		}
		printf("\n");
	}
}

void run_predicting_test(bool& ok) {
	test_tool t("test/digits/predicting", ok);

	typedef recognition_problem p;

	reexp::lang<p> lang;
	reexp::data<p> data(lang, sample_dim());
	setup_reg<p>(lang, data);

	reexp::stats<p> stats(data);

	print_prediction(t, stats);
}

void run_learn_picked_predict_test(bool& ok) {
	test_tool t("test/digits/learn_picked_predict", ok);

	typedef recognition_problem p;

	reexp::lang<p> lang;
	reexp::data<p> data(lang, sample_dim());
	setup_reg<p>(lang, data);

	reexp::stats<p> stats(data);


	float before = stats.naiveInfo();
	lang.add_exp(lang.rel(relid::up_down), 3);
	float after = stats.naiveInfo();
	t<<"info "<<before<<" -> "<<after<<checkl<<checkl;
	before = after;
	lang.add_exp(lang.rel(relid::left_right), 3);
	after = stats.naiveInfo();
	t<<"info "<<before<<" -> "<<after<<checkl<<checkl;

	t<<checkl<<"learning done."<<checkl<<checkl;

	print_scan(t, lang, stats);

	reexp::pinfo info;
	setup_pinfo(info);
	reexp::lang_info<p> li(info, lang);

	t<<"vars:\n\n"<<li.drawn_vars_tostring(cvarid::x, cvarid::y);

	t<<"px"<<checkl<<checkl;
	print_pixels(t, data.var(varid::pixel));
	for (int i = 0; i < 1; i++) {
		t<<"exp"<<i<<checkl<<checkl;
		print_pixels(t, data.var(varid::digit9 + i + 1));
	}

	print_prediction<p>(t, stats);
}


void run_learn_predict_test(bool& ok) {
	test_tool t("test/digits/learn_predict", ok);

	typedef recognition_problem p;

	reexp::lang<p> lang;
	reexp::data<p> data(lang, sample_dim());
	setup_reg<p>(lang, data);

	reexp::stats<p> stats(data);

	reexp::learner<p> learner(lang, stats);
	setup_learner(learner);

	int rounds = 3;

	for (int i = 0; i < rounds; i++) {
		float before = stats.naiveInfo();
		if (!learner.add_exp()) break;
		t<<"expression added."<<checkl;
		float after = stats.naiveInfo();
		t<<"info "<<before<<" -> "<<after<<checkl<<checkl;
	}

	t<<checkl<<"learning done."<<checkl<<checkl;

	print_scan(t, lang, stats);

	reexp::pinfo info;
	setup_pinfo(info);
	reexp::lang_info<p> li(info, lang);

	t<<"vars:\n\n"<<li.drawn_vars_tostring(cvarid::x, cvarid::y);

	t<<"px"<<checkl<<checkl;
	print_pixels(t, data.var(varid::pixel));
	for (int i = 0; i < rounds; i++) {
		t<<"exp"<<i<<checkl<<checkl;
		print_pixels(t, data.var(varid::digit9 + i + 1));
	}

	print_prediction<p>(t, stats);
}

void run_exp_predicting_test(bool& ok) {
	test_tool t("test/digits/exp_predicting", ok);

	typedef recognition_problem p;

	reexp::lang<p> lang;
	reexp::data<p> data(lang, sample_dim());
	setup_reg<p>(lang, data);

	reexp::stats<p> stats(data);

	lang.add_exp(lang.rel(relid::up_down), 3);
	t<<"exp added."<<checkl<<checkl;
	lang.add_exp(lang.rel(relid::left_right), 3);
	t<<"exp added."<<checkl<<checkl;

/*	lang.add_exp(lang.rel(relid::left_right), 3);
	t<<"exp added."<<checkl<<checkl;*/

	print_prediction<p>(t, stats);
}


void digitstest(bool& ok) {
	run_learn_picked_predict_test(ok);
	run_learn_predict_test(ok);
	run_exp_predicting_test(ok);
	run_predicting_test(ok);
	run_learning_test(ok);
	run_learning_once_test(ok);
	run_scan_test(ok);
	run_gen_rel_test(ok);
	run_exp_stats_test(ok);
	run_add_exp_test(ok);
	run_stats_test(ok);
	run_print_test(ok);
	run_setup_test(ok);
	run_ctx_test(ok);
}
