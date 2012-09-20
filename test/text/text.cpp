/*
 * rooms.cpp
 *;
 *  Created on: Oct 30, 2011
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"

namespace {

	struct text_problem {
		static const int DIM = 2;
		static const int MAX_REL_VARS = 2; // two max relation variables
	};

	// context variables
	namespace cvarid {
		static const int index 		= 0;
		static const int sample		= 1;
	}
	// bit variables (used to group bits)
	namespace varid {
		static const int ws  = 0;
		static const int a 	= 1;
		static const int z = 1 + 'z' - 'a';
		static const int firstclass = z + 1;
	}
	static const int VarCount = 4;

	// relations (used to organize bits)
	namespace relid {
		static const int next_id 		= 0;
	}
	const char* cvarnames[] =  {
		"i",
		"s"
	};

	struct sample {
		int class_;
		const char* text_;

	};

	template <typename P>
	void populate(explib::data<P>& data, const sample* samples) {
		explib::cvec<P> at;
		for (int i = 0; samples[i].text_; ++i) {
			const sample& s = samples[i];
			at[cvarid::sample] = i;
			*(data.var(varid::firstclass + s.class_)[at]) = true;

			for (int c = 0; s.text_[c]; ++c) {
				at[cvarid::index] = c;
				int ch = tolower(s.text_[c]);
				int vid = -1;
				if (ch >= 'a' && c <= 'z') {
					vid = ch - 'a' + varid::a;
				} else {
					vid = varid::ws;
				}
				*(data.var(vid)[at]) = true;
			}

		}
	}

	template <typename P>
	void setup_forward_rels(explib::lang<P>& lang, explib::var<P>& var1, explib::var<P>& var2) {
		explib::ctx<P> var_ctx(explib::cvec<P>(0, 0));
		explib::rel<P>& rl(lang.alloc_rel(var_ctx)); // prev next
		rl.add_var(explib::cvec<P>(0, 0), var1);
		rl.add_var(explib::cvec<P>(1, 0), var2);
		lang.rel_done();
	}

	template <typename P>
	void setup_class_rels(explib::lang<P>& lang, explib::var<P>& charvar, explib::var<P>& classvar) {
		explib::ctx<P> var_ctx(explib::cvec<P>(0, 0));
		explib::rel<P>& rl(lang.alloc_rel(var_ctx));
		rl.add_var(explib::cvec<P>(0, 0), charvar);
		rl.add_var(explib::cvec<P>(0, 0), classvar);
		lang.rel_done();
	}

	template <typename P>
	void setup_lang(explib::lang<P>& lang, int classes) {
		explib::ctx<P> charctx(explib::cvec<P>(0, 0));
		lang.add_orig(explib::orig<P>(charctx)); // ws
		for (int c = 'a'; c <= 'z'; c++) {
			lang.add_orig(explib::orig<P>(charctx));
		}
		explib::ctx<P> classctx(explib::cvec<P>(-1, 0));
		for (int i = 0; i < classes; ++i) {
			lang.add_orig(explib::orig<P>(classctx));
		}
		for (int i = varid::ws; i <= varid::z; ++i) {
			for (int j = varid::ws; j <= varid::z; ++j) {
				setup_forward_rels<P>(lang, lang.var(i), lang.var(j));
			}
			for (int j = 0; j < classes; ++j) {
				setup_class_rels<P>(lang,
									lang.var(i),
									lang.var(j+varid::firstclass));
			}
		}
	}

	template <typename P>
	void setup_reg(explib::lang<P>& lang, explib::data<P>& data) {
		setup_lang(lang);
		populate(data);
	}

	template <typename P>
	void setup_names(explib::pinfo& info, const char** classNames) {
		info.vnames_.push_back(" "); // whitespace
		for (char i = 'a'; i <= 'z'; ++i) {
			char buf[2] = {i, '\0'};
			info.vnames_.push_back(buf);
		}
		for (int i = 0; classNames[i]; ++i) {
			info.vnames_.push_back(classNames[i]);
		}
		for (int i = 0; i < P::DIM; ++i) {
			info.cvnames_.push_back(cvarnames[i]);
		}
	}

	explib::cvec<text_problem> text_dim(const sample* samples) {
		explib::cvec<text_problem> v = {0, 0};
		for (int i = 0; samples[i].text_; ++i) {
			v[cvarid::index] = std::max(size_t(v[cvarid::index]),
										strlen(samples[i].text_));
			v[cvarid::sample]++;
		}
		return v;
	}

	template <typename Out>
	void write_var(Out& out,
				   const explib::lang<text_problem>& l,
				   explib::pinfo& info,
				   int vi) {
		const explib::var<text_problem>& var = l.var(vi);
		const explib::exp<text_problem>* exp = dynamic_cast<const explib::exp<text_problem>*>(&var);
		if (exp) {
			const std::vector<explib::rel_entry<text_problem> >& e = exp->rel().entries();
			write_var(out, l, info, e[0].var_->id());
			write_var(out, l, info, e[1].var_->id());
		} else {
			out<<info.vnames_[vi];
		}
	}

	void print_text(TestTool& t,
					const explib::data<text_problem>& d,
					explib::pinfo& info,
					const explib::lang_info<text_problem>& li,
					int classes) {
		typedef text_problem p;

		explib::cvec<p> dim = d.dim();

		for (int s = 0; s < dim[cvarid::sample]; ++s) {
			explib::cvec<p> at;
			at[cvarid::sample] = s;
			int c = -1;
			for (int i = 0; i < classes; ++i) {
				const explib::data_var<p>& cv = d.var(varid::firstclass + i);
				if (cv[at] && *cv[at]) { c = i; break; }
			}
			t<<s;
			if (c >= 0) t<<" "<<c;
			t<<":  ";
			for (size_t i = 0; i < size_t(dim[cvarid::index]); ++i) {
				at[cvarid::index] = i;
				int v = -1;
				for (size_t j = 0; j < d.var_count(); ++j) {
					const explib::data_var<p>& dv = d.var(j);
					if (dv.var().ctx().v_[0]>=0&&i<size_t(dv.dim()[0])&&dv[at]&&*dv[at]) {
						v = j;
						break;
					}
				}
				if (v >= 0) {
					write_var(t, d.lang(), info, v);
					t<<"|";
				}
			}
			t<<"\n";
		}
	}

	sample piratevsbear[] = {
		{0, "yaar yar"},
		{1, "mrr mrr "},
		{0, "yarr yaarr mate"},
		{1, "mrr mrr"},
		{0, "yaaar"},
		{0, "yar yaaaarr "},
		{0, "aar mate"},
		{1, "mrrrrr"},
		{1, "mr mrr mr"},
		{0, "mate"},

		{-1, NULL}
	};

	const char* piratevsbear_classes[] = {
		"pirate",
		"bear",
		0
	};

	void setup_test(TestTool& t) {
		typedef text_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, text_dim(piratevsbear));

		setup_lang(lang, 2);
		populate(data, piratevsbear);

		explib::stats<p> stats(data);

		explib::pinfo names;
		setup_names<p>(names, piratevsbear_classes);

		explib::stats_info<p> si(names, stats);

		t<<"text:\n";
		print_text(t, data, names, si.lang_info(), 2);
		t<<"\n\n";
		t<<si.vars_tostring();
	}

	void learning_test(TestTool& t) {
		typedef text_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, text_dim(piratevsbear));

		setup_lang(lang, 2);
		populate(data, piratevsbear);

		explib::stats<p> stats(data);

		explib::pinfo names;
		setup_names<p>(names, piratevsbear_classes);
		explib::stats_info<p> si(names, stats);

		explib::learner<p> learner(lang, stats, 7);
		int exps = learner.reexpress(true);
		t<<exps<<" exps added\n\n";
		explib::pinfo i;

		print_text(t, data, names, si.lang_info(), 2);

		t<<"\nvars:\n";
		for (int i = 0; i < lang.var_count(); ++i) {
			write_var(t, lang, names, i);
			t<<"\n";
		}
	}

	void predicting_test(TestTool& t) {
		typedef text_problem p;

		explib::lang<p> lang;
		explib::data<p> data(lang, text_dim(piratevsbear));

		setup_lang(lang, 2);
		populate(data, piratevsbear);

		explib::stats<p> stats(data);

		explib::pinfo names;
		setup_names<p>(names, piratevsbear_classes);
		explib::stats_info<p> si(names, stats);

		t<<si.preds_tostring(varid::firstclass, varid::firstclass+2, true)<<"\n";
		t<<si.entropy_tostring(varid::firstclass)<<"\n";
		t<<si.entropy_tostring(varid::firstclass+1)<<"\n";

		explib::learner<p> learner(lang, stats, 7);
		int exps = learner.reexpress(true);
		t<<exps<<" exps added\n\n";

		t<<si.preds_tostring(varid::firstclass, varid::firstclass+2, true)<<"\n";
		t<<si.entropy_tostring(varid::firstclass)<<"\n";
		t<<si.entropy_tostring(varid::firstclass+1)<<"\n";
	}

	const char* generated_classes[] = {
		"moreA",
		"moreB",
		"ab",
		"ba",
		0
	};


	void generate(int clazz, std::string& buf, int len) {
		buf.resize(len+1);
		for (int i = 0; i < len; ++i) {
			if (clazz == 0) {
				buf[i] = (rand() % 3) ? 'a' : 'b';
			} else {
				buf[i] = (rand() % 3) ? 'b' : 'a';
			}
		}
		buf[len] = '\0';
	}

	void genpred_test(TestTool& t) {
		srand(0);

		typedef text_problem p;

		int textlen = 8;
		int n = 20;
		int classes = 2;
		std::vector<std::string> texts;

		std::vector<sample> samples;
		t<<"original:\n";
		for (int i = 0; i < n; ++i) {
			std::string text;
			int clazz = i%classes;
			generate(clazz, text, textlen);
			texts.push_back(std::move(text));
			samples.push_back({clazz, texts.back().c_str()});
			t<<texts.back().c_str()<<"\n";
		}
		t<<"\n";
		samples.push_back({-1, NULL});

		explib::cvec<p> dim = {textlen, n};
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);

		setup_lang(lang, classes);
		populate(data, samples.data());

		explib::stats<p> stats(data);

		explib::pinfo names;
		setup_names<p>(names, generated_classes);
		explib::stats_info<p> si(names, stats);

		t<<"text:\n";
		print_text(t, data, names, si.lang_info(), classes);
		t<<"\n\n";

		t<<si.vars_tostring();
		t<<"\n\n";

		for (int i = 0; i < classes; ++i) {
			t<<si.pred_tostring(varid::firstclass+i, true)<<"\n";
			t<<si.entropy_tostring(varid::firstclass+i)<<"\n";
			t<<si.row_logdep_tostring(varid::firstclass+i);
		}


/*		explib::learner<p> learner(lang, stats, 7);
		int exps = learner.reexpress(true);
		t<<exps<<" exps added\n\n";

		t<<si.preds_tostring(varid::firstclass, varid::firstclass+2, true)<<"\n";
		t<<si.entropy_tostring(varid::firstclass, true)<<"\n";
		t<<si.entropy_tostring(varid::firstclass+1, true)<<"\n";*/
	}

	void eval_genpredwith(TestTool& t, int n, int tn) {

		typedef text_problem p;

		int textlen = 5;
		int classes = 2;
		std::vector<std::string> texts;
//		t<<"test generated text with "<<n<<" samples:\n\n";

		explib::cvec<p> dim = {textlen, n};
		explib::lang<p> lang;
		explib::data<p> data(lang, dim);

		setup_lang(lang, classes);
		{
			std::vector<sample> samples;
			for (int i = 0; i < n; ++i) {
				std::string text;
				int clazz = i%classes;
				generate(clazz, text, textlen);
				texts.push_back(std::move(text));
				samples.push_back({clazz, texts.back().c_str()});
			}
			samples.push_back({-1, NULL});

			populate(data, samples.data());
		}

		explib::stats<p> stats(data);

		explib::cvec<p> tdim = {textlen, tn};
		explib::data<p> tdata(lang, tdim);
		{
			std::vector<sample> samples;
			for (int i = 0; i < tn; ++i) {
				std::string text;
				int clazz = i%classes;
				generate(clazz, text, textlen);
				texts.push_back(std::move(text));
				samples.push_back({clazz, texts.back().c_str()});
			}
			samples.push_back({-1, NULL});

			populate(tdata, samples.data());
		}
		tdata.apply_exps();

		explib::pinfo names;
		setup_names<p>(names, generated_classes);
		explib::stats_info<p> si(names, stats);

		std::ostringstream nlabel;
		nlabel<<"n:"<<n;

/*		t<<"for "<<n<<" samples: \n";*/

		explib::pred<p> pred(stats);


		for (int i = 0; i < classes; ++i) {
			double totalinfo = 0;
			double entryinfo = 0;
			int clz = varid::firstclass+i;
			t<<si.row_logdep_tostring(clz, varid::a);
			t<<si.row_logdep_tostring(clz, varid::a+1);

//			t<<names.vnames_[clz]<<":\n";
			{
//				pred.rowP<TestTool&>(data, clz, t);
				pred.info(data, clz, totalinfo, entryinfo);
				t<<"teach total: "<<totalinfo<<ignorel;
				t<<"teach entry: "<<entryinfo<<ignorel;
				t.record({names.vnames_[clz], "sample:teach", "info", nlabel.str()}, totalinfo);
				t.record({names.vnames_[clz], "sample:teach", "entryinfo", nlabel.str()}, entryinfo);
			}
			{
//				pred.rowP<TestTool&>(tdata, clz, t);
				pred.info(tdata, clz, totalinfo, entryinfo);

				t<<"test total: "<<totalinfo<<ignorel;
				t<<"test entry: "<<entryinfo<<ignorel;
				t.record({names.vnames_[clz], "sample:test", "info", nlabel.str()}, totalinfo);
				t.record({names.vnames_[clz], "sample:test", "entryinfo", nlabel.str()}, entryinfo);

				t<<si.pred_tostring(tdata, clz, true);
			}
		}
		t<<".";
	}

	void genpredeval_test(TestTool& t) {
		srand(0);

		t<<"running genpred eval:\n\n";

		for (int i = 0; i < 16; ++i) {
			eval_genpredwith(t, 0,  16);
			eval_genpredwith(t, 1,  16);
			eval_genpredwith(t, 2,  16);
			eval_genpredwith(t, 4,  16);
			eval_genpredwith(t, 8,  16);
			eval_genpredwith(t, 16, 16);
			eval_genpredwith(t, 32, 16);
			eval_genpredwith(t, 64, 16);
			eval_genpredwith(t, 128, 16);
			t<<"\n";
		}

		Table table(
			t.report(ToTable<Average>({"entryinfo"}, "sample:", "n:")));
		std::ostringstream buf;
		table>>buf;
		t<<"\nresults:\n\n"<<buf.str()<<"\n";

/*		plot("test/text/genpredeval.plot",
			 "entryinfo",
			 "n",
			 "sample");*/
	}



	void setup_text_lang(explib::lang<text_problem>& l,
						 std::vector<int>& v,
						 int& len,
						 std::istream& in) {
		v.resize(256);
		for (size_t i = 0; i < v.size(); ++i) {
			v[i] = -1;
		}
		len = 0;
		while (in) {
			int c = in.get();
			if (c < 0) break;
			len++;
			if (v[c] < 0) {
				v[c] = l.var_count();
				l.add_orig(explib::orig<text_problem>(explib::cvec<text_problem>(0, 0)));
			}
		}
		for (int i = 0; i < l.var_count(); ++i) {
			for (int j = 0; j < l.var_count(); ++j) {
				setup_forward_rels<text_problem>(l, l.var(i), l.var(j));
			}
		}
	}

	void populate_text(explib::data<text_problem>& d, std::vector<int>& v, int len, std::istream& in) {
		int at = 0;
		for (size_t i = 0; i < d.var_count(); ++i) {
			d.var(i).defined().fill(true);
		}
		while (in && at < len) {
			int c = in.get();
			explib::data_var<text_problem>& dv = d.var(v[c]);
			*dv[explib::cvec<text_problem>(at, 0)] = true;
#if 0
			for (int i = 0; i < v[c]; ++i) {
				explib::data_var<text_problem>& dv = d.var(i);
				dv[explib::cvec<text_problem>(at, 0)] = false;
			}
#endif
			at++;
		}
	}

	void setup_names(explib::pinfo& info, explib::lang<text_problem>& l, std::vector<int>& v) {
		info.vnames_.resize(l.var_count());
		for (size_t i = 0; i < v.size(); ++i) {
			if (v[i] >= 0) {
				info.vnames_[v[i]] = sup()<<char(i);
			}
		}
		for (int i = 0; i < text_problem::DIM; ++i) {
			info.cvnames_.push_back(cvarnames[i]);
		}

	}

	void code_test(TestTool& t) {
		typedef text_problem p;
		explib::lang<p> lang;

		std::vector<int> char_vars;
		int len;
		{
			std::ifstream file("src/reexp/lang.h");
			setup_text_lang(lang, char_vars, len, file);
		}
		explib::data<p> data(lang, explib::cvec<p>(len, 1));
		{
			std::ifstream file("src/reexp/lang.h");
			populate_text(data, char_vars, len, file);
		}

		explib::stats<p> s(data);
		double th = 100;
		explib::learner<p> l(lang, s, th, 0.4*th);
		explib::pinfo names;
		setup_names(names, lang, char_vars);
		explib::stats_info<p> si(names, s);

		t<<si.vars_tostring()<<"\n";

		t<<"scan:\n"<<si.scan_tostring()<<"\n";

		double naiveInfo = s.naiveInfo();
		int exps = l.reexpress(true);
		double naiveInfo2 = s.naiveInfo();

		t<<exps<<" exps added.\n\n";

		t<<"naive information dropped "<<(naiveInfo/8)<<"B -> "<<(naiveInfo2/8)<<"B for "<<len<<"B of text.\n\n";

		//print_text(t, data, names, si.lang_info(), 2);


		t<<"\nvars:\n";

		std::vector<std::pair<double, int> > ordered;

		for (int i = lang.orig_count(); i < lang.var_count(); ++i) {
			ordered.push_back(std::pair<double, int>(s.var(i).freq(), i));
		}
		std::sort(ordered.begin(), ordered.end());
		std::reverse(ordered.begin(), ordered.end());

		for (std::pair<double, int>& o : ordered) {
			if (o.first >= 5) {

				std::ostringstream buf;
				write_var(buf, lang, names, o.second);
				std::string str = buf.str();
				t<<"\""<<str<<"\"";
				for (int i = 2 + str.length(); i < 30; ++i) t<<" ";
				t<<"["<<o.first<<"]\n";
			}
		}

	}

}

void addtexttest(TestRunner& runner) {
	runner.add("text/setup", 		{"func"}, 	&setup_test);
	runner.add("text/learning", 	{"func"}, 	&learning_test);
	runner.add("text/predicting", 	{"func"},  	&predicting_test);
	runner.add("text/genpred", 		{"func"}, 	&genpred_test);
	runner.add("text/genpredeval", 	{"func"}, 	&genpredeval_test);
	runner.add("text/code", 	    {"func"}, 	&code_test);
}
