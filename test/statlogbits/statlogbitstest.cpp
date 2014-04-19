/*
 * statlogbitstest.cpp
 *
 *  Created on: Jun 23, 2012
 *      Author: arau
 */
#include "reexp/all.h"

#include "exptesttools.h"
#include <stdexcept>
#include <iostream>
#include "evaluation/evaluation.h"

namespace {

	using namespace evaluation;

	typedef reexp::traits1d statlogbits_problem;


	enum undef_t {
		undef_nothing		= 0,
		undef_excluded 		= 1
	};

	void setup_statlogbits(pred_problem<statlogbits_problem>& pr,
						   const char* file,
						   undef_t undefs = undef_nothing) {
		typedef statlogbits_problem p;
		std::ifstream in(std::string("test/statlogbits/") +file + "_bits.txt");
		std::string type;
		int num;
		in>>type>>num;
		int cat; in>>cat;
		std::vector<int> catbits;
		std::vector<int> bitcat;
		int bits = 0;
		for (int i = 0; i < cat; ++i) {
			catbits.push_back(0);
			in>>catbits.back();
			std::string line; std::getline(in, line);
			for (int j = 0; j < catbits.back(); ++j) {
				pr.names_.vnames_.push_back("");
				std::getline(in, pr.names_.vnames_.back());
				bits++;
				bitcat.push_back(i);
			}
		}
		int samples = 0;
		in>>samples;
	    pr.data_.set_dim(reexp::cvec<p>(samples));

		reexp::ctx<p> varctx(reexp::cvec<p>(0));
		for (int i = 0; i < bits; ++i) {
			pr.lang_.add_orig(reexp::orig<p>(varctx));
		}

		std::string b;
		for (int s = 0; s < samples; ++s) {
			in>>b;
			for (int i = 0; i < bits; ++i) {
				pr.data_.var(i)[s] = true;
				bool v = (b[i] == '1');
				*pr.data_.var(i)[s] = v;
				if (v && undefs == undef_excluded) {
					for (int j = i-1; j >= 0 && bitcat[j] == bitcat[i]; j--) {
						pr.data_.var(j)[s] = false; // undefine
					}
				}
			}
		}

		for (int i = 0; i < bits; ++i) {
			int f = pr.data_.var(i).states().popcount();
			int n = pr.data_.var(i).defined().popcount();
			if (f && f != n) { // is not constant
				for (int j = 0; j < i; ++j) {
					int f2 = pr.data_.var(j).states().popcount();
					int n2 = pr.data_.var(j).defined().popcount();
					if (f2 && f2 != n2) { // is not constant
						reexp::rel<p>& rl(pr.lang_.alloc_rel(varctx));
						rl.add_var(reexp::cvec<p>(0), pr.lang_.var(i));
						rl.add_var(reexp::cvec<p>(0), pr.lang_.var(j));
						pr.lang_.rel_done();
					}
				}
			}
		}
/*		for (int i = 0, b = 0; i < cat; ++i) {
			if (catbits[i] == 1) { // treat as binary variable
				pr.lang_.var(b++).setPrioriP(0.5);
			} else {
				for (int j = 0; j < catbits[i]; ++j) {
					pr.lang_.var(b++).setPrioriP(1./(catbits[i]+1));
				}
			}
		}*/

/*  	for (int i = 0; i < bits; ++i) {
			int f = pr.data_.var(i).states().popcount();
			int n = pr.data_.var(i).defined().popcount();
			pr.lang_.var(i).setPrioriP((double)(f+1)/(double)(n+2));
		}*/

		int v = reexp::util::next_version();
		for (size_t i = 0; i < pr.data_.var_count(); ++i) {
			pr.data_.var(i).version_ = v;
		}
	}


	void generic_dist_test(test_tool& t, pred_problem<statlogbits_problem>& pr,
						    pred_args& args, int groups, const std::string& tag) {
		args.dist_record_groups_ = groups;
		args.tags_.insert(tag);
		crossvalidate_run(t, pr, args, 10);

		table dist(
			t.report(to_table<average>({"run:out", tag}, "prop:", "rank:")));

		t<<"distribution: "<<tag<<"\n"<<dist<<"\n";

		std::ofstream f(t.file_path(tag + ".tex"));
		f<<dist.to_latex_pgf_doc("probability estimate");

	}

	void generic_dists_test(test_tool& t, pred_problem<statlogbits_problem>& pr,
						    pred_args& base_args, double th, size_t groups) {
		{
			pred_args args(th, 0.2*th, base_args.predrelfilter_, base_args.prioriweight_, base_args.verbose_, base_args.useLogDepB_, -1);
			generic_dist_test(t, pr, args, groups, "dist:reexp");
		}

		{
			pred_args args(th, 0.2*th, base_args.predrelfilter_, base_args.prioriweight_, base_args.verbose_, base_args.useLogDepB_, base_args.scaling_group_sz_);
			generic_dist_test(t, pr, args, groups, "dist:scaled reexp");
		}

		{
			pred_args args(10000, 10000, base_args.predrelfilter_, base_args.prioriweight_, base_args.verbose_, base_args.useLogDepB_, -1);
			generic_dist_test(t, pr, args, groups, "dist:naive");
		}

		{
			pred_args args(10000, 10000, base_args.predrelfilter_, base_args.prioriweight_, base_args.verbose_, base_args.useLogDepB_, base_args.scaling_group_sz_);
			generic_dist_test(t, pr, args, groups, "dist:scaled naive");
		}

		table dist(
			t.report(to_table<average>({"run:out", "prop:e_p"}, "dist:", "rank:")));
		t<<"e_p: \n"<<dist<<"\n";
		std::ofstream f(t.file_path("e_p.tex"));
		f<<dist.to_latex_pgf_doc("probability estimate");
	}

	void setup_heart(pred_problem<statlogbits_problem>& pr, undef_t undefs = undef_excluded) {
		setup_statlogbits(pr, "heart", undefs);
		pr.predvars_.resize(pr.lang_.var_count());
		pr.predvars_[pr.lang_.var_count()-1] = true;
		pr.costs_[true_positive] = 0;
		pr.costs_[true_negative] = 0;
		pr.costs_[false_positive] = 5;
		pr.costs_[false_negative] = 1;
	}

	void setup_heart_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr);

		reexp::stats<statlogbits_problem> stats(pr.data_);
		reexp::stats_info<statlogbits_problem> si(pr.names_, stats);

		t<<si.vars_tostring();
	}

	void heart_exp_priories_test(test_tool& t) {
		typedef statlogbits_problem p;
		pred_problem<p> pr;
		setup_heart(pr);

		reexp::stats<p> stats(pr.data_);

		reexp::learner<p> learner(pr.lang_, stats, 25.,25*0.2,0);
		learner.reexpress(true);

		reexp::stats_info<p> si(pr.names_, stats);

		double totalPrioriInfo = 0;
		double totalAverInfo = 0;
		std::ostringstream buf;
		buf.precision(3);
		for (int i = 0; i < pr.lang_.var_count(); ++i) {
			const reexp::var<p>& v = pr.lang_.var(i);
			const reexp::var_stats<p>& vs = stats.var(i);
			double prioriInfo = reexp::estimateEntropy(vs.p(), v.prioriP());
			double averInfo = reexp::entropy(vs.p());
			buf<<si.var_tostring(i);
			buf<<"priori:     "; buf.width(6); buf<<v.prioriP()<<"\t aver:     "; buf.width(6); buf<<vs.p()<<"\n";
			buf<<"prioriInfo: "; buf.width(6); buf<<prioriInfo<< "\t averInfo: "; buf.width(6); buf<<averInfo<<"\n";
			totalPrioriInfo += prioriInfo;
			totalAverInfo += averInfo;
		}
		buf<<"\ninfo: "<<totalAverInfo<<"/"<<totalPrioriInfo<<"\n";
		t<<buf.str();
	}

	void run_heart_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr);

		double th = 25;
		pred_args args = {th, th*0.2, 0., 2., true};
		crossvalidate_run(t, pr, args, 10);
	}

	void run_heart_scan_test(test_tool& t) {
		typedef statlogbits_problem p;
		pred_problem<p> pr;
		setup_heart(pr, undef_excluded);

		reexp::var<p>& v = pr.lang_.var(pr.lang_.var_count()-1);

		reexp::stats<statlogbits_problem> s(pr.data_);

		reexp::learner<statlogbits_problem> l(pr.lang_, s, 9, 3);
		l.exclude(pr.predvars_);
		l.reexpress();


		reexp::stats_info<statlogbits_problem> si(pr.names_, s);
		t<<"deps:\n";
		t<<si.var_deps_tostring(v, 30);
	}
	void print_data_by(test_tool& t, const char* tag, bool includeCost, bool crossCompressed=false) {
		table exptable(
			t.report(to_table<average>({"run:out"}, "exp:", tag)));
		t.ignored()<<"expression & relations:\n"<<exptable<<"\n";

		table traintable(
			t.report(to_table<average>({"run:out", "data:train"}, "prop:", tag)));
		t<<"train:\n"<<traintable<<"\n";

		table testtable(
			t.report(to_table<average>({"run:out", "data:test"}, "prop:", tag)));
		t<<"test:\n"<<testtable<<"\n";

		t<<"entropy:\n";

		table table2(
			t.report(to_table<average>({"run:out", "prop:entropy"}, "data:", tag)));

		t<<table2.to_plot(2, 20, 0.7)<<"\n";

		{
			std::ofstream entropy(t.file_path("entropy.tex"));
			entropy<<table2.to_latex_pgf_doc("crossentropy");
		}

		t<<"err%:\n";

		table table3(
			t.report(to_table<average>({"run:out", "prop:err"}, "data:", tag)));

		t<<table3.to_plot(2, 20, 0.2)<<"\n";

		{
			std::ofstream entropy(t.file_path("error.tex"));
			entropy<<table3.to_latex_pgf_doc("error%");
		}


		if (includeCost) {
			t<<"cost:\n";
			table table4(
				t.report(to_table<average>({"run:out", "prop:cost"}, "data:", tag)));
			t<<table4.to_plot(2, 20, 0.5)<<"\n";
			{
				std::ofstream entropy(t.file_path("cost.tex"));
				entropy<<table4.to_latex_pgf_doc("cost");
			}
		}


		if (crossCompressed) {
			t<<"cross compression byte sizes:\n";
			{
				table table4(
					t.report(to_table<average>({"run:out", "prop:encodedB"}, "data:", tag)));
				t<<table4.to_plot(2, 20, 0.5)<<"\n";
				std::ofstream entropy(t.file_path("encodedB.tex"));
				entropy<<table4.to_latex_pgf_doc("bytes");
			}
			{
				table table4(
					t.report(to_table<average>({"run:out", "prop:exndl"}, "data:", tag)));
				t<<table4.to_plot(2, 20, 0.5)<<"\n";
				std::ofstream entropy(t.file_path("exndl.tex"));
				entropy<<table4.to_latex_pgf_doc("bytes");
			}
			{
				table table4(
					t.report(to_table<average>({"run:out", "prop:xndl"}, "data:", tag)));
				t<<table4.to_plot(2, 20, 0.5)<<"\n";
				std::ofstream entropy(t.file_path("xndl.tex"));
				entropy<<table4.to_latex_pgf_doc("bytes");
			}
			{
				table table4(
					t.report(to_table<average>({"run:out", "prop:ndl"}, "data:", tag)));
				t<<table4.to_plot(2, 20, 0.5)<<"\n";
				std::ofstream entropy(t.file_path("ndl.tex"));
				entropy<<table4.to_latex_pgf_doc("bytesB");
			}

/*			table testtable(
				t.report(to_table<average>({"run:out", "data:test"}, "prop:", tag)));
			t<<"test:\n"<<testtable<<"\n";*/

/*			table table5(
				t.report(to_table<average>({"run:out", "data:test"}, "prop:entropy", "prop:encodedB")));
			t<<"mxe by encodedB:\n"<<table5<<"\n";*/
			{
				std::ofstream entropy(t.file_path("test_enc_vs_mxe.tex"));
				entropy<<testtable.to_scatter_latex_pgf_doc("encodedB", {"entropy"}, "mean cross entropy");
			}
			{
				std::ofstream entropy(t.file_path("train_enc_vs_mxe.tex"));
				entropy<<traintable.to_scatter_latex_pgf_doc("encodedB", {"entropy"}, "mean cross entropy");
			}
			table traintable(
				t.report(to_table<average>({"run:out", "data:train"}, "perf:", tag)));
			t.ignored()<<"train performance:\n"<<traintable<<"\n";
			{
				std::ofstream entropy(t.file_path("train_perf.tex"));
				entropy<<traintable.to_latex_pgf_doc("perf");
			}
			table testtable(
				t.report(to_table<average>({"run:out", "data:test"}, "perf:", tag)));
			t.ignored()<<"test performance:\n"<<testtable<<"\n";
			{
				std::ofstream entropy(t.file_path("test_perf.tex"));
				entropy<<testtable.to_latex_pgf_doc("perf");
			}
		} else {
			table table5(
				t.report(to_table<average>({"run:out", "perf:ns"}, "data:", tag)));

			t<<"\nperformance (ns / entry)\n";
			t.ignored()<<table5;
		}

		t<<"\nns\n";

		t.ignored()<<
			table(t.report(to_table<average>({"run:out", "perf:ns"}, "data:", tag)))
				.to_plot(2, 20, 0)<<"\n";


	}

	void run_heart_filter_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr);
		int filtersteps = 10;
		for (int i = 0; i < filtersteps; ++i) {
			pred_args args = {1000., 1000., double(i), 3, false};
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "predfilter:", true);
	}

	void run_heart_crosscompress_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr, undef_nothing);
		int steps = 26;
		for (int i = 0; i < steps; ++i) {
			double th = 5 + i*5;
			pred_args args(th, .4*th, 0, 2., 0, true);
			args.crosscompress_ = true;
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "threshold:", true, true);
	}

	void run_heart_narrowexps_crosscompress_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr, undef_nothing);
		int steps = 21;
		for (int i = 0; i < steps; ++i) {
			double th = 4 + i;
			pred_args args(th, .4*th, 0, 2., 0, true);
			args.crosscompress_ = true;
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "threshold:", true, true);
	}


	void run_heart_prioriw_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr, undef_excluded);
		int filtersteps = 15;
		for (int i = 0; i < filtersteps; ++i) {
			pred_args args = {1000., 1000., 0, double(i), false};
			crossvalidate_run(t, pr, args, 9);
		}

		print_data_by(t, "prioriw:", true);
	}

	void generic_run_heart_exps_test(test_tool& t, undef_t undefs, bool logDepB = false, int scale_group_sz = -1, int thstep = 5) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr, undefs);
		int steps = 26;
		for (int i = 0; i < steps ; ++i) {
			double th = 5 + i * thstep;
//			double th = 4 + i * thstep;
			pred_args args(th, 0.2*th, 0, 3., 0, logDepB, scale_group_sz);
//			args.crosscompress_ = true;
			crossvalidate_run(t, pr, args, 9);
		}

		print_data_by(t, "threshold:", true);
	}

	void run_heart_exps_test_undef(test_tool& t) {
		generic_run_heart_exps_test(t, undef_excluded, false, -1, 2);
	}

	void run_heart_exps_b_test_undef(test_tool& t) {
		generic_run_heart_exps_test(t, undef_excluded, true, -1,  2);
	}

	void run_heart_exps_b_scaled_test_undef(test_tool& t) {
		generic_run_heart_exps_test(t, undef_excluded, true, 50,  2);
	}

	void run_heart_narrowexps_b_test_undef(test_tool& t) {
		generic_run_heart_exps_test(t, undef_excluded, true, -1, 1);
	}

	void run_heart_narrowexps_b_scaled_test_undef(test_tool& t) {
		generic_run_heart_exps_test(t, undef_excluded, true, 50, 1);
	}

	void run_heart_exps_test_noundefs(test_tool& t) {
		generic_run_heart_exps_test(t, undef_nothing);
	}

	void run_heart_exps_b_test_noundefs(test_tool& t) {
		generic_run_heart_exps_test(t, undef_nothing, true);
	}

	void run_heart_exps_b_scaled_test_noundefs(test_tool& t) {
		generic_run_heart_exps_test(t, undef_nothing, true, 50);
	}

	void run_heart_narrowexps_b_test_noundefs(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr, undef_nothing);
		int steps = 20;
		for (int i = 0; i < steps ; ++i) {
//			double th = 5 + i;
			double th = 4 + i;
			pred_args args(th, 0.2*th, -1000., 3., 0, true);
			crossvalidate_run(t, pr, args, 9);
		}
		print_data_by(t, "threshold:", true);
	}

	void run_heart_dist_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_heart(pr, undef_nothing);
		pred_args args(-1, -1, 0, 3., false, true, 50);
		generic_dists_test(t, pr, args, 8, 12);
	}


	void setup_statlog_shuttle(pred_problem<statlogbits_problem>& pr, undef_t undefs = undef_excluded) {
		setup_statlogbits(pr, "statlog_shuttle", undefs);
		pr.predvars_.resize(pr.lang_.var_count());
		for (int i = 0; i < 7; ++i) {
			int v = pr.lang_.var_count()-1-i;
			pr.data_.var(v).defined().fill(true); // mark all these undefined
			pr.predvars_[pr.lang_.var_count()-1-i] = true;
		}
		pr.type_ = classification_problem;
	}

	void run_setup_statlog_shuttle_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_statlog_shuttle(pr, undef_nothing );

		reexp::stats<statlogbits_problem> stats(pr.data_);
		reexp::stats_info<statlogbits_problem> si(pr.names_, stats);

		t<<si.vars_tostring();

	}

	void run_statlog_shuttle_singlerun_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_statlog_shuttle(pr, undef_nothing );
		double th = 100;

		pred_args args(th, th*0.4, 0, 1.);
		separate_train_test_datas_run(t, pr, args, 43500);

		table exptable(
			t.report(to_table<average>({"run:out"}, "exp:", "predfilter:")));
		t.ignored()<<"expression & relations:\n"<<exptable<<"\n";

		table props(
			t.report(to_table<average>({"run:out"}, "prop:", "data:")));
		t<<"properties:\n"<<props.formatted(4)<<"\n";
	}

	void run_generic_statlog_shuttle_filter_test(test_tool& t, undef_t undefs) {
		pred_problem<statlogbits_problem> pr;
		setup_statlog_shuttle(pr, undefs);

		int filtersteps = 10;
		for (int i = 0; i < filtersteps; ++i) {
			pred_args args(100000., 100000.,double(i), 1.);
			separate_train_test_datas_run(t, pr, args, 43500);
		}

		print_data_by(t, "predfilter:", true);

/*		table exptable(
			t.report(to_table<average>({"run:out"}, "exp:", "predfilter:")));
		t.ignored()<<"expression & relations:\n"<<exptable<<"\n";

		table traintable(
			t.report(to_table<average>({"run:out", "data:train"}, "prop:", "predfilter:")));
		t<<"train:\n"<<traintable<<"\n";

		table testtable(
			t.report(to_table<average>({"run:out", "data:test"}, "prop:", "predfilter:")));
		t<<"test:\n"<<testtable<<"\n";

		t<<"entropy:\n";

		table table2(
			t.report(to_table<average>({"run:out", "prop:entropy"}, "data:", "predfilter:")));

		t<<table2.to_plot(2, 20, 0.7)<<"\n";

		t<<"err%:\n";

		table table3(
			t.report(to_table<multiplied<average, 100> >({"run:out", "prop:err"}, "data:", "predfilter:")));

		t<<table3.to_plot(2, 20, 0.2)<<"\n";

		table table5(
			t.report(to_table<average>({"run:out", "perf:ns"}, "data:", "predfilter:")));

		t<<"\nperformance (ns / entry)\n";
		t.ignored()<<table5;

		t<<"\nns by filter\n";

		t.ignored()<<
			table(t.report(to_table<average>({"run:out", "perf:ns"}, "data:", "predfilter:")))
				.to_plot(2, 20, 0)<<"\n";*/

	}
	void run_statlog_shuttle_filter_test(test_tool& t) {
		run_generic_statlog_shuttle_filter_test(t, undef_excluded);
	}
	void run_statlog_shuttle_filter_noundef_test(test_tool& t) {
		run_generic_statlog_shuttle_filter_test(t, undef_nothing);
	}


	void run_generic_statlog_shuttle_exps_test(test_tool& t, undef_t undefs, bool logDepB = false, int scale_group_sz = -1) {
		pred_problem<statlogbits_problem> pr;
		setup_statlog_shuttle(pr, undefs);

		int steps = 9;
		double th = 100;
		for (int i = 0; i < steps ; ++i) {
			pred_args args(th, 0.4*th, 0, 1., 0, logDepB, scale_group_sz);
			separate_train_test_datas_run(t, pr, args, 43500);
			th*=2;
		}

		print_data_by(t, "threshold:", true);

/*		table exptable(
			t.report(to_table<average>({"run:out"}, "exp:", "threshold:")));
		t.ignored()<<"expression & relations:\n"<<exptable<<"\n";

		table traintable(
			t.report(to_table<average>({"run:out", "data:train"}, "prop:", "threshold:")));
		t<<"train:\n"<<traintable.formatted(4)<<"\n";

		table testtable(
			t.report(to_table<average>({"run:out", "data:test"}, "prop:", "threshold:")));
		t<<"test:\n"<<testtable.formatted(4)<<"\n";

		t<<"entropy:\n";

		table table2(
			t.report(to_table<average>({"run:out", "prop:entropy"}, "data:", "threshold:")));

		t<<table2.to_plot(3, 30, 0)<<"\n";

		t<<"err%:\n";

		table table3(
			t.report(to_table<multiplied<average, 100>> ({"run:out", "prop:err"}, "data:", "threshold:")));

		t<<table3.to_plot(3, 30, 0)<<"\n";

		table table5(
			t.report(to_table<average>({"run:out", "perf:ns"}, "data:", "threshold:")));

		t<<"\nperformance (ns / entry)\n";
		t.ignored()<<table5;

		t<<"\nns by filter\n";

		t.ignored()<<
			table(t.report(to_table<average>({"run:out", "perf:ns"}, "data:", "threshold:")))
				.to_plot(3, 30, 0)<<"\n";*/
	}

	void run_statlog_shuttle_exps_test(test_tool& t) {
		run_generic_statlog_shuttle_exps_test(t, undef_excluded);
	}

	void run_statlog_shuttle_exps_noundef_test(test_tool& t) {
		run_generic_statlog_shuttle_exps_test(t, undef_nothing);
	}

	void run_statlog_shuttle_exps_noundef_b_test(test_tool& t) {
		run_generic_statlog_shuttle_exps_test(t, undef_nothing, true);
	}

	void run_statlog_shuttle_exps_noundef_b_scaled_test(test_tool& t) {
		run_generic_statlog_shuttle_exps_test(t, undef_nothing, true, 128);
	}

	void setup_australian(pred_problem<statlogbits_problem>& pr, undef_t undefs = undef_excluded) {
		setup_statlogbits(pr, "australian_credit", undefs);
		pr.predvars_.resize(pr.lang_.var_count());
		for (int i = 0; i < 2; ++i) {
			int v = pr.lang_.var_count()-1-i;
			pr.data_.var(v).defined().fill(true); // mark all these undefined
			pr.predvars_[v] = true;
		}
		pr.type_ = classification_problem;
	}

	void run_australian_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr);

		double th = 25;
		pred_args args = {th, th*0.2, 0., 2, true};
		crossvalidate_run(t, pr, args, 10);
	}

	void run_australian_filter_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr);

		int filtersteps = 10;
		for (int i = 0; i < filtersteps; ++i) {
			pred_args args = {1000., 1000., double(i), 3, false};
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "predfilter:", false);
	}

	void run_australian_exps_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr, undef_nothing);
		int steps = 22;
		for (int i = 0; i < steps ; ++i) {
			double th = 30 + i * 15;
			pred_args args(th, 0.2*th, 0);
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "threshold:", false);
	}

	void run_australian_exps_b_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr, undef_nothing);
		int steps = 22;
		for (int i = 0; i < steps ; ++i) {
			double th = 30 + i * 15;
			pred_args args(th, 0.2*th, 0, 2., false, true);
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "threshold:", false);
	}

	void run_australian_exps_b_scaled_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr, undef_nothing);
		int steps = 22;
		for (int i = 0; i < steps ; ++i) {
			double th = 30 + i * 15;
			pred_args args(th, 0.2*th, 0, 2., false, true, 115);
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "threshold:", false);
	}

	void run_australian_prioriw_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr);

		int filtersteps = 20;
		for (int i = 0; i < filtersteps; ++i) {
			pred_args args = {1000., 1000., double(0), double(i), false};
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "prioriw:", false);
	}

	void run_australian_dist_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_australian(pr, undef_nothing);
		pred_args args(-1, -1, 0, 2., false, true, 115);
		generic_dists_test(t, pr, args, 240, 12);
	}


	void run_heart_exps_segfault_test(test_tool& t) {
		typedef statlogbits_problem p;
		pred_problem<p> pr;
		setup_heart(pr, undef_nothing);
		double th = 5;
		reexp::stats<p> stats(pr.data_);
		reexp::learner<p> learner(pr.lang_, stats, th, 0.2*th, 0.);
		reexp::stats_info<p> si(pr.names_, stats);
		for (int i = 0; i < 60; ++i) {
			if (!learner.add_exp()) break;
		}
		pr.lang_.unset_obs();
		t<<si.vars_tostring();
		learner.add_exp();
		int v = pr.lang_.var_count()-1;
		t<<si.lang_info().var_tostring(v)<<"\n";
		t<<si.lang_info().drawn_var_implmasks_to_string(pr.lang_.var(v));
	}

	void setup_germancredit(pred_problem<statlogbits_problem>& pr, undef_t undefs = undef_excluded) {
		setup_statlogbits(pr, "german_credit", undefs);
		pr.predvars_.resize(pr.lang_.var_count());
		pr.predvars_[pr.lang_.var_count()-1] = true;
		pr.costs_[true_positive] = 0;
		pr.costs_[true_negative] = 0;
		pr.costs_[false_positive] = 5;
		pr.costs_[false_negative] = 1;
	}

	void run_setup_germancredit_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_germancredit(pr, undef_nothing);

		reexp::stats<statlogbits_problem> s(pr.data_);
		reexp::stats_info<statlogbits_problem> si(pr.names_, s);

		t<<si.vars_tostring();
		t<<si.var_deps_tostring(pr.lang_.var_count()-1, 10);
	}

	void run_generic_germancredit_exps_test(test_tool& t, undef_t undefs, bool predLogDepB = false, size_t scaling_group_sz = -1) {
		pred_problem<statlogbits_problem> pr;
		setup_germancredit(pr, undefs);
		int steps = 33;
		for (int i = 0; i < steps ; ++i) {
			double th = 40 + i * 15;
			pred_args args(th, 0.2*th, 0, 2., false, predLogDepB, scaling_group_sz);
			crossvalidate_run(t, pr, args, 10);
		}
		print_data_by(t, "threshold:", true);
	}

	void run_germancredit_exps_test(test_tool& t) {
		run_generic_germancredit_exps_test(t, undef_excluded);
	}

	void run_germancredit_exps_noundef_test(test_tool& t) {
		run_generic_germancredit_exps_test(t, undef_nothing);
	}

	void run_germancredit_exps_noundef_b_test(test_tool& t) {
		run_generic_germancredit_exps_test(t, undef_nothing, true);
	}

	void run_germancredit_dist_test(test_tool& t) {
		pred_problem<statlogbits_problem> pr;
		setup_germancredit(pr, undef_nothing);
		pred_args args(-1, -1, 0, 2., false, true, 128);
		generic_dists_test(t, pr, args, 85, 12);
	}


	void run_germancredit_exps_noundef_b_scaled_test(test_tool& t) {
		run_generic_germancredit_exps_test(t, undef_nothing, true, 128);
	}

}


void addstatlogbitstest(TestRunner& runner) {
	runner.add("statlogbits/setup_heart", 			{"func"}, &setup_heart_test);
	runner.add("statlogbits/heart_exp_priories",    {"func"}, &heart_exp_priories_test);
	runner.add("statlogbits/heart", 				{"func"}, &run_heart_test);
	runner.add("statlogbits/heart_scan", 			{"func"}, &run_heart_scan_test);
	runner.add("statlogbits/heart_filter", 			{"func"}, &run_heart_filter_test);
	runner.add("statlogbits/heart_prioriw", 		{"func"}, &run_heart_prioriw_test);
	runner.add("statlogbits/heart_crosscompress", 	{"func"}, &run_heart_crosscompress_test);
	runner.add("statlogbits/heart_narrowexps_crosscompress", 	{"func"}, &run_heart_narrowexps_crosscompress_test);
	runner.add("statlogbits/heart_exps", 			{"func"}, &run_heart_exps_test_undef);
	runner.add("statlogbits/heart_exps_b", 			{"func"}, &run_heart_exps_b_test_undef);
	runner.add("statlogbits/heart_exps_b_scaled", 			{"func"}, &run_heart_exps_b_scaled_test_undef);
	runner.add("statlogbits/heart_narrowexps_b", 	{"func"}, &run_heart_narrowexps_b_test_undef);
	runner.add("statlogbits/heart_narrowexps_b_scaled", 	{"func"}, &run_heart_narrowexps_b_scaled_test_undef);
 	runner.add("statlogbits/heart_exps_noundefs", 	{"func"}, &run_heart_exps_test_noundefs);
 	runner.add("statlogbits/heart_exps_noundefs_b", 	{"func"}, &run_heart_exps_b_test_noundefs);
 	runner.add("statlogbits/heart_exps_noundefs_b_scaled", 	{"func"}, &run_heart_exps_b_scaled_test_noundefs);
 	runner.add("statlogbits/heart_narrowexps_noundefs_b", 	{"func"}, &run_heart_narrowexps_b_test_noundefs);
	runner.add("statlogbits/heart_exps_segfault", 	{"func"}, &run_heart_exps_segfault_test);
	runner.add("statlogbits/heart_dist", 			{"func"}, &run_heart_dist_test);
	runner.add("statlogbits/setup_statlog_shuttle", {"func"}, &run_setup_statlog_shuttle_test);
	runner.add("statlogbits/statlog_shuttle_singlerun",{"func"}, &run_statlog_shuttle_singlerun_test);
	runner.add("statlogbits/statlog_shuttle_filter",{"func"}, &run_statlog_shuttle_filter_test);
	runner.add("statlogbits/statlog_shuttle_filter_noundef",{"func"}, &run_statlog_shuttle_filter_noundef_test);
	runner.add("statlogbits/statlog_shuttle_exps",  {"func"}, &run_statlog_shuttle_exps_test);
	runner.add("statlogbits/statlog_shuttle_exps_noundef",  {"func"}, &run_statlog_shuttle_exps_noundef_test);
	runner.add("statlogbits/statlog_shuttle_exps_noundef_b",  {"func"}, &run_statlog_shuttle_exps_noundef_b_test);
	runner.add("statlogbits/statlog_shuttle_exps_noundef_b_scaled",  {"func"}, &run_statlog_shuttle_exps_noundef_b_scaled_test);
	runner.add("statlogbits/australian", 			{"func"}, &run_australian_test);
	runner.add("statlogbits/australian_filter", 	{"func"}, &run_australian_filter_test);
	runner.add("statlogbits/australian_exps", 	    {"func"}, &run_australian_exps_test);
	runner.add("statlogbits/australian_exps_b", 	    {"func"}, &run_australian_exps_b_test);
	runner.add("statlogbits/australian_exps_b_scaled",{"func"}, &run_australian_exps_b_scaled_test);
	runner.add("statlogbits/australian_prioriw", 	{"func"}, &run_australian_prioriw_test);
	runner.add("statlogbits/australian_dist", 		{"func"}, &run_australian_dist_test);
	runner.add("statlogbits/setup_germancredit", 	{"func"}, &run_setup_germancredit_test);
	runner.add("statlogbits/germancredit_exps", 	{"func"}, &run_germancredit_exps_test);
	runner.add("statlogbits/germancredit_exps_noundef", 	{"func"}, &run_germancredit_exps_noundef_test);
	runner.add("statlogbits/germancredit_exps_noundef_b", 	{"func"}, &run_germancredit_exps_noundef_b_test);
	runner.add("statlogbits/germancredit_exps_noundef_b_scaled", 	{"func"}, &run_germancredit_exps_noundef_b_scaled_test);
	runner.add("statlogbits/germancredit_dist", 	{"func"}, &run_germancredit_dist_test);
}
