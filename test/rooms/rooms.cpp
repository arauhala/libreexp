/*
 * rooms.cpp
 *
 *  Created on: Oct 30, 2011
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"
#include "evaluation/evaluation.h"

namespace {

	struct rooms_problem {
		static const int DIM = 4; // 4 context variables
		static const int MAX_REL_VARS = 2; // two max relation variables
	};

	static const int Width = 5;

	static const int Height = 5;

	static const int MaxTime = 20;

	static const int Games = 200;

	static const int VarCount = 13;

	// context variables
	namespace cvarid {
		static const int x 			= 0;
		static const int y			= 1;
		static const int turn 		= 2;
		static const int game 		= 3;
	}
	// bit variables (used to group bits)
	namespace varid {
		static const int iswall 	= 0;
		static const int isdark     = 1;
		static const int isbright   = 2;
		static const int isexit 	= 3;
		static const int victory 	= 4;
		static const int escape 	= 5;

		static const int goleft 	= 6;
		static const int goright 	= 7;
		static const int goup 		= 8;
		static const int godown 	= 9;
		static const int idle	 	= 10;

		static const int underplayer= 11; // meta variable
		static const int nearplayer	= 12; // meta variable

	}
	// relations (used to organize bits)
	namespace relid {
		static const int left_right 		= 0;
		static const int up_down 			= 1;
		static const int past_future 		= 2;
		static const int past_future_left 	= 3;
		static const int past_future_right 	= 4;
		static const int past_future_up 	= 5;
		static const int past_future_down 	= 6;
		// 25 relations from
	}
	const char* varnames[] =  {
		"iswall",   // 0
		"isdark",   // 1
		"isbright", // 2
		"isexit",
		"victory",  // 4
		"escape",   // 5
		"goleft",
		"goright",
		"goup ",
		"godown",
		"idle",    	// 10
		"underplayer", // 11
		"nearplayer"
	};
	const char* cvarnames[] =  {
		"x",
		"y",
		"turn",
		"game"
	};

	template <typename P>
	void populate(reexp::data<P>& data) {
		std::ifstream in("test/rooms/game.txt");

		reexp::data_var<P>& iswall      = data.var(varid::iswall);
		reexp::data_var<P>& isdark      = data.var(varid::isdark);
		reexp::data_var<P>& isbright    = data.var(varid::isbright);
		reexp::data_var<P>& isexit      = data.var(varid::isexit);
		reexp::data_var<P>& victory     = data.var(varid::victory);
		reexp::data_var<P>& escape 	 	= data.var(varid::escape);

		reexp::data_var<P>& goleft 	 	= data.var(varid::goleft);
		reexp::data_var<P>& goright 	= data.var(varid::goright);
		reexp::data_var<P>& godown 	 	= data.var(varid::godown);
		reexp::data_var<P>& goup 		= data.var(varid::goup);
		reexp::data_var<P>& idle 		= data.var(varid::idle);

		reexp::data_var<P>& underplayer = data.var(varid::underplayer);
		reexp::data_var<P>& nearplayer  = data.var(varid::nearplayer);

		iswall.defined().fill(false);
		isdark.defined().fill(false);
		isbright.defined().fill(false);
		isexit.defined().fill(false);
		victory.defined().fill(false);
		escape.defined().fill(false);

		goleft.defined().fill(false);
		goright.defined().fill(false);
		godown.defined().fill(false);
		goup.defined().fill(false);
		idle.defined().fill(false);

		int game = 0;
		int turn = 0;
		int l = 0;
		reexp::cvec<P> at;
		while (in && game < Games) {
			std::string line;
			std::getline(in, line); ++l;

			if (line == "reset") {
				at[cvarid::game] = game++;
				at[cvarid::turn] = turn = 0;
			} else {
				escape[at] = true;
				bool running = false;
				if (line == "running.") {
					running = true;
				} else if (line == "escaped.") {
					*escape[at] = true;
					victory[at] = true;
					*victory[at] = true;
				} else if (line == "time out.") {
					victory[at] = true;
					*victory[at] = false;
				} else {
					std::ostringstream buf;
					buf<<"expected game state, not "<<l<<": '"<<line<<"'";
					throw std::runtime_error(buf.str());
				}

				for (int y = 0; y < 5; ++y) {
					at[cvarid::y] = y;
					std::getline(in, line); ++l;
					for (int x = 0; x < 5; ++x) {
						at[cvarid::x] = x;
						char c = line[x];
						underplayer[at] = true;
						nearplayer[at] = true;
						switch (c) {
							case ' ':
								isbright[at] = true;
							case '>':
							case '.':
								isdark[at] = true;
							case '#':
								iswall[at] = true;
							case '<':
								isexit[at] = true;
						}
						switch (c) {
							case ' ': *isbright[at] = true; break;
							case '>':
							case '.': *isdark[at] = true; break;
							case '#': *iswall[at] = true; break;
							case '<': *isexit[at] = true; *isbright[at] = true; break;
						}
						if (x == 2 && y == 2) {
							*underplayer[at] = true;
						}
						if (x >= 1 && x < 4 && y >= 1 && y < 4) {
							*nearplayer[at] = true;
						}
					}
				}
				if (running) {
					int decision = -1;
					std::getline(in, line); ++l;
					if (line == "go left.") {
						decision = 0;
					} else if (line == "go right.") {
						decision = 1;
					} else if (line == "go up.") {
						decision = 2;
					} else if (line == "go down.") {
						decision = 3;
					} else if (line == "idle.") {
						decision = 4;
					} else {
						std::ostringstream buf;
						buf<<"expected decision, not "<<l<<": '"<<line<<"'";
						throw std::runtime_error(buf.str());
					}

					if (running) {
						for (int y = 0; y < 5; ++y) {
							at[cvarid::y] = y;
							for (int x = 0; x < 5; ++x) {
								at[cvarid::x] = x;
								idle[at] 	= true;
								godown[at] 	= true;
								goup[at] 	= true;
								goright[at] = true;
								goleft[at] 	= true;

							  /*switch (decision) {
									case 4: idle[at] 	= true;
									case 3: godown[at] 	= true;
									case 2: goup[at] 	= true;
									case 1: goright[at] = true;
									case 0: goleft[at] 	= true;
								}*/
								switch (decision) {
									case 4: *idle[at] = true; 		break;
									case 3: *godown[at] = true; 	break;
									case 2: *goup[at] = true; 		break;
									case 1: *goright[at] = true;	break;
									case 0: *goleft[at] = true; 	break;
								}
							}
						}
					}
				}
				at[cvarid::turn] = ++turn;
			}
		}
	}

	reexp::cvec<rooms_problem> rooms_dim() {
		reexp::cvec<rooms_problem> v;
		v[cvarid::x] = Width;
		v[cvarid::y] = Height;
		v[cvarid::turn] = MaxTime;
		v[cvarid::game] = Games;
		return v;
	}

	template <typename P>
	void setup_heuristic_rels(reexp::lang<P>& lang, reexp::var<P>& var) {
		reexp::rel<P>& rel(lang.alloc_rel(var.ctx()));
		rel.add_var(reexp::cvec<P>(0, 0, 0, 0), var);
		rel.add_var(reexp::cvec<P>(0, 0, 0, 0), lang.var(varid::victory));
		lang.rel_done();
	}

	template <typename P>
	void setup_map_rels(reexp::lang<P>& lang, reexp::var<P>& var, reexp::var<P>& var2) {
		reexp::ctx<P> map_ctx(reexp::cvec<P>(0, 0, 0, 0));

		reexp::rel<P>& rl(lang.alloc_rel(map_ctx)); // right left
		rl.add_var(reexp::cvec<P>(0, 0, 0, 0), var);
		rl.add_var(reexp::cvec<P>(1, 0, 0, 0), var2);
		lang.rel_done();

		reexp::rel<P>& ud(lang.alloc_rel(map_ctx)); // up down
		ud.add_var(reexp::cvec<P>(0, 0, 0, 0), var);
		ud.add_var(reexp::cvec<P>(0, 1, 0, 0), var2);
		lang.rel_done();

		reexp::rel<P>& pf(lang.alloc_rel(map_ctx)); // past future
		pf.add_var(reexp::cvec<P>(0, 0, 0, 0), var);
		pf.add_var(reexp::cvec<P>(0, 0, 1, 0), var2);
		lang.rel_done();

		reexp::rel<P>& pfl(lang.alloc_rel(map_ctx)); // past future left
		pfl.add_var(reexp::cvec<P>(1, 0, 0, 0), var);
		pfl.add_var(reexp::cvec<P>(0, 0, 1, 0), var2);
		lang.rel_done();

		reexp::rel<P>& pfr(lang.alloc_rel(map_ctx)); // past future right
		pfr.add_var(reexp::cvec<P>(0, 0, 0, 0), var);
		pfr.add_var(reexp::cvec<P>(1, 0, 1, 0), var2);
		lang.rel_done();

		reexp::rel<P>& pfu(lang.alloc_rel(map_ctx)); // past future up
		pfu.add_var(reexp::cvec<P>(0, 1, 0, 0), var);
		pfu.add_var(reexp::cvec<P>(0, 0, 1, 0), var2);
		lang.rel_done();

		reexp::rel<P>& pfd(lang.alloc_rel(map_ctx)); // past future down
		pfd.add_var(reexp::cvec<P>(0, 0, 0, 0), var);
		pfd.add_var(reexp::cvec<P>(0, 1, 1, 0), var2);
		lang.rel_done();


	}



	template <typename P>
	void setup_lang(reexp::lang<P>& lang) {

		reexp::ctx<P> map_ctx(reexp::cvec<P>(0, 0, 0, 0));
		reexp::ctx<P> turn_ctx(reexp::cvec<P>(-1, -1, 0, 0));
		reexp::ctx<P> game_ctx(reexp::cvec<P>(-1, -1, -1, 0));
		lang.add_orig(reexp::orig<P>(map_ctx)); // iswall
		lang.add_orig(reexp::orig<P>(map_ctx)); // isdark
		lang.add_orig(reexp::orig<P>(map_ctx)); // isbright
		lang.add_orig(reexp::orig<P>(map_ctx)); // isexit
		lang.add_orig(reexp::orig<P>(game_ctx)); // victory
		lang.add_orig(reexp::orig<P>(turn_ctx)); // escape

		lang.add_orig(reexp::orig<P>(map_ctx)); // go left
		lang.add_orig(reexp::orig<P>(map_ctx)); // go right
		lang.add_orig(reexp::orig<P>(map_ctx)); // go up
		lang.add_orig(reexp::orig<P>(map_ctx)); // go down
		lang.add_orig(reexp::orig<P>(map_ctx)); // idle

		lang.add_orig(reexp::orig<P>(map_ctx)); // under player
		lang.add_orig(reexp::orig<P>(map_ctx)); // near player

		for (int i = varid::iswall; i <= varid::escape; ++i) {
			if (i != varid::victory) {
				setup_heuristic_rels(lang, lang.var(i));
			}
		}
		for (int i = varid::iswall; i <= varid::isexit; ++i) {
			for (int j = varid::iswall; j <= i; j++) {
				setup_map_rels(lang, lang.var(i), lang.var(j));
			}
/*			for (int j = varid::underplayer; j <= varid::nearplayer; ++j) {
				explib::rel<P>& rel(lang.alloc_rel(map_ctx));
				rel.add_var(explib::cvec<P>(0, 0, 0, 0), lang.var(i));
				rel.add_var(explib::cvec<P>(0, 0, 0, 0), lang.var(j));
				lang.rel_done();
			}*/
			for (int j = varid::escape; j <= varid::nearplayer; ++j) {
				reexp::rel<P>& rel(lang.alloc_rel(map_ctx));
				rel.add_var(reexp::cvec<P>(0, 0, 0, 0), lang.var(i));
				rel.add_var(reexp::cvec<P>(0, 0, 0, 0), lang.var(j));
				lang.rel_done();
			}
		}


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
	void draw_map(test_tool& t, const reexp::data<P>& data, int turn, int game) {
		typedef rooms_problem p;

		const reexp::data_var<P>& iswall     = data.var(varid::iswall);
		const reexp::data_var<P>& isdark     = data.var(varid::isdark);
		const reexp::data_var<P>& isbright   = data.var(varid::isbright);
		const reexp::data_var<P>& isexit     = data.var(varid::isexit);

		int w = data.dim()[cvarid::x];
		int h = data.dim()[cvarid::y];
		reexp::cvec<P> at;
		at[cvarid::turn] = turn;
		at[cvarid::game] = game;
		for (int y = 0; y < h; y++) {
			at[cvarid::y] = y;
			for (int x = 0; x < w; x++) {
				at[cvarid::x] = x;
				if (*iswall[at]) {
					t<<"#";
				} else if (*isexit[at]){
					t<<"<";
				} else if (*isbright[at]) {
					t<<" ";
				} else if (*isdark[at]) {
					t<<".";
				} else {
					t<<"?";
				}
			}
			t<<"\n";
		}
	}

	void setup_test(test_tool& t) {
		typedef rooms_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, rooms_dim());

		setup_reg<p>(lang, data);

		reexp::pinfo info;
		setup_names<p>(info);

		reexp::lang_info<p> li(info, lang);
		t<<"vars:\n"<<li.vars_tostring();
		t<<"rels:\n"<<li.rels_tostring();

		const reexp::data_var<p>& iswall     = data.var(varid::iswall);
		const reexp::data_var<p>& victory     = data.var(varid::victory);
		reexp::cvec<p> at = {0, 0, 0, 0};
		for (int i = 0; i < 10; ++i) {
			at[cvarid::game] = i;

			t<<"game "<<i<<". ";
			if (*victory[at]) {
				t<<"victory.\n\n";
			} else {
				t<<"defeat.\n\n";
			}

			for (int j = 0; j < MaxTime; ++j) {
				at[cvarid::turn] = j;
				if (iswall[at]) { // is this defined?
					draw_map<p>(t, data, j, i);
					t<<"\n";
				}
			}
		}
	}

	void scan_test(test_tool& t) {
		typedef rooms_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, rooms_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);

		t<<si.vars_tostring();
		t<<"scans:\n"<<si.scan_tostring(50, 1);

		t<<"\ngoright deps:\n"<<si.var_deps_tostring(varid::goright, 20);
		t<<"\ngoleft deps:\n"<<si.var_deps_tostring(varid::goleft, 20);
		t<<"\ngoup deps:\n"<<si.var_deps_tostring(varid::goup, 20);
		t<<"\ngodown deps:\n"<<si.var_deps_tostring(varid::godown, 20);
		t<<"\nidle deps:\n"<<si.var_deps_tostring(varid::idle, 20);

		t<<"\nescape deps:\n"<<si.var_deps_tostring(varid::escape, 20);
		t<<"\nvictory deps:\n"<<si.var_deps_tostring(varid::victory, 20);
	}

	void learning_test(test_tool& t) {
		typedef rooms_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, rooms_dim());

		setup_reg<p>(lang, data);

		reexp::stats<p> stats(data);

		reexp::learner<p> learner(lang, stats, 50, 25);
		learner.exclude(varid::victory);
		learner.exclude(varid::escape);
		int exps = learner.reexpress(true, 20);
		t<<exps<<" exps added\n\n";
		reexp::pinfo i;
		setup_names<p>(i);
		reexp::stats_info<p> si(i, stats);

		t<<si.vars_tostring();
		t<<si.lang_info().drawn_vars_tostring(cvarid::x, cvarid::y);

		t<<"scans:\n"<<si.scan_tostring(50, 1);

		t<<"\ngoright deps:\n"<<si.var_deps_tostring(varid::goright, 20);
		t<<"\ngoleft deps:\n"<<si.var_deps_tostring(varid::goleft, 20);
		t<<"\ngoup deps:\n"<<si.var_deps_tostring(varid::goup, 20);
		t<<"\ngodown deps:\n"<<si.var_deps_tostring(varid::godown, 20);
		t<<"\nidle deps:\n"<<si.var_deps_tostring(varid::idle, 20);

		t<<"\nescape deps:\n"<<si.var_deps_tostring(varid::escape, 20);
		t<<"\nvictory deps:\n"<<si.var_deps_tostring(varid::victory, 20);
	}

	using namespace evaluation;

	void setup_predproblem(pred_problem<rooms_problem>& pr) {
		pr.predvars_.resize(varid::victory+1);
		pr.predvars_[varid::victory] = true;
		setup_reg<rooms_problem>(pr.lang_, pr.data_);
		setup_names<rooms_problem>(pr.names_);
		pr.costs_[true_positive] = 0;
		pr.costs_[true_negative] = 0;
		pr.costs_[false_positive] = 1;
		pr.costs_[false_negative] = 1;
		pr.samplecvar_ = cvarid::game;
	}

	void heuristics_byfilter_test(test_tool& t) {
		typedef rooms_problem p;

		pred_problem<p> pr(rooms_dim());
		setup_predproblem(pr);

		for (int i = 0; i < 8; ++i) {
			pred_args args(1024, 1024, double(i));
			crossvalidate_run(t, pr, args, 3);
		}

		t.ignored()<<"expression & relations:\n"<<t.report(to_table<average>({"run:out"}, "exp:", "predfilter:"))<<"\n";

		t<<"results (train):\n"<<t.report(to_table<average>({"run:out", "data:train"}, "prop:", "predfilter:"))<<"\n";

		t<<"results (test):\n"<<t.report(to_table<average>({"run:out", "data:test"}, "prop:", "predfilter:"))<<"\n";

		t<<"entropy by filter:\n"
		 <<t.report(to_table<average>({"run:out", "prop:entropy"}, "data:", "predfilter:")).to_plot(3, 20, 0.5)<<"\n";
	}

	void heuristics_byexps_test(test_tool& t) {
		typedef rooms_problem p;

		pred_problem<p> pr(rooms_dim());
		setup_predproblem(pr);

		double predrelfilter = 0;

		int th = 8;
		for (int i = 0; i < 7; ++i) {
			pred_args args( double(th), th*0.2, predrelfilter);
			crossvalidate_run(t, pr, args, 3);
			th<<=1;
		}

		t.ignored()<<"expression & relations:\n"<<t.report(to_table<average>({"run:out"}, "exp:", "threshold:"))<<"\n";

		t<<"results (train):\n"<<t.report(to_table<average>({"run:out", "data:train"}, "prop:", "threshold:"))<<"\n";

		t<<"results (test):\n"<<t.report(to_table<average>({"run:out", "data:test"}, "prop:", "threshold:"))<<"\n";

		t<<"entropy by threshold:\n"
		 <<t.report(to_table<average>({"run:out", "prop:entropy"}, "data:", "threshold:")).to_plot(3, 20, 0.5)<<"\n";
	}

}
void addroomstest(TestRunner& runner) {
	runner.add("rooms/setup", 	 {"func"}, &setup_test);
	runner.add("rooms/scan", {"func"}, &scan_test);
	runner.add("rooms/learning", {"func"}, &learning_test);
	runner.add("rooms/heuristics_byfilter", {"func"}, &heuristics_byfilter_test);
	runner.add("rooms/heuristics_byexps", {"func"}, &heuristics_byexps_test);
}
