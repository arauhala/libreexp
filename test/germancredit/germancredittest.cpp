/*
 * germancredittest.cpp
 *
 *  Created on: Apr 14, 2012
 *      Author: arau
 */

#include "reexp/all.h"

#include "exptesttools.h"
#include <stdexcept>
#include <iostream>
#include "evaluation/evaluation.h"

namespace {

	using namespace evaluation;

	static const size_t SampleCount = 1000;

	struct germancredit_problem {
		static const int DIM = 1;
		static const int MAX_REL_VARS = 2; // two max relation variables
	};

	enum undef_t {
		undef_nothing 					= 0,
		undef_excluded_discrete_vars 	= 0x1,
		undef_excluded_numeric_vars 	= 0x2,
		undef_all_excluded_vars 		= 0x3
	};

	const std::vector<std::pair<std::string, std::string>> discrete_variables = {
		{"A11", "accountEmpty"},
		{"A12", "accountUnder200DM"},
		{"A13", "accountOver200DM"},
		{"A14", "noAccount"},

		{"A30", "noCreditTakenOrPaidDuly"},
		{"A31", "creditForThisBankPaidDuly"},
		{"A32", "creditForThisBankPaidDulyTillNow"},
		{"A33", "creditDelays"},
		{"A34", "criticalAccount"},

		{"A40", "forNewCar"},
		{"A41", "forUsedCar"},
		{"A42", "forFurnitudeOrEquipment"},
		{"A43", "forRadioTelevision"},
		{"A44", "forDomesticAppliance"},
		{"A45", "forRepairs"},
		{"A46", "forEducation"},
		{"A47", "forVacation"},
		{"A48", "forRetraining"},
		{"A49", "forBusiness"},
		{"A410", "forOthers"},

		{"A61", "savings(<100)"},
		{"A62", "savings(100-500)"},
		{"A63", "savings(500-1k)"},
		{"A64", "savings(>1k)"},
		{"A65", "savingsUnknown"},

		{"A71", "unemployed"},
		{"A72", "employed(<1y)"},
		{"A73", "employed(1y-4y)"},
		{"A74", "employed(4y-7y)"},
		{"A75", "employed(>7y)"},

		{"A91", "male&divorced/separated"},
		{"A92", "femaled&divorced/separated/married"},
		{"A93", "male&single"},
		{"A94", "male&married/widowed"},
		{"A95", "female&single"},

		{"A101", "singleApplicantNoGuarantor"},
		{"A102", "coapplicant"},
		{"A103", "guarantor"},


		{"A121", "hasRealEstate"},
		{"A122", "hasLifeInsurance"},
		{"A123", "hasCarOrOther"},
		{"A124", "unknownOrNoProperty"},

		{"A141", "hasBankPlans"},
		{"A142", "hasStorePlans"},
		{"A143", "hasNoPlan"},

		{"A151", "home(rent)"},
		{"A152", "home(own)"},
		{"A153", "home(free)"},

		{"A171", "unemployedOrUnskilled"},
		{"A172", "unskilled"},
		{"A173", "skilled"},
		{"A173", "managementOrHighSkillOrSelfEmployed"},

		{"A191", "noTelephone"},
		{"A192", "hasTelephone"},

		{"A201", "foreignWorker"},
		{"A202", "notForeignWorker"}
	};

/*
 *     static const int DURATION = 1;
    static const int DURATION_BITS = 5;
    static const int DURATION_SCALE = 3;

    static const int CREDITS = 4;
    static const int CREDITS_BITS = 5;
    static const int CREDITS_SCALE = 500;

    static const int INSTALLMENT_RATE = 7;
    static const int INSTALLMENT_RATE_BITS = 3;
    static const int INSTALLMENT_RATE_SCALE = 1;

    static const int PRESENT_RESIDENCE_SINCE = 10;
    static const int PRESENT_RESIDENCE_SINCE_BITS = 3;
    static const int PRESENT_RESIDENCE_SINCE_SCALE = 1;

    static const int AGE = 12;
    static const int AGE_BITS = 5;
    static const int AGE_SCALE = 10;

    static const int AGE_MINUS = 10;

    static const int CREDITS_IN_BANK = 15;
    static const int CREDITS_IN_BANK_BITS = 3;
    static const int CREDITS_IN_BANK_SCALE = 1;

    static const int MAINTAINS_PEOPLE = 17;
	static const int MAINTAINS_PEOPLE_BITS = 1;
	static const int MAINTAINS_PEOPLE_SCALE = 1;
 *
 */

    static const int DURATION = 1;
    static const int DURATION_BITS = 4;
    static const double DURATION_SCALE = 3;
    static const double DURATION_OFFSET = 3;

    static const int CREDITS = 4;
    static const int CREDITS_BITS = 4;
    static const double CREDITS_SCALE = 500;

    static const int INSTALLMENT_RATE = 7;
    static const int INSTALLMENT_RATE_BITS = 3;
    static const int INSTALLMENT_RATE_SCALE = 1;

    static const int PRESENT_RESIDENCE_SINCE = 10;
    static const int PRESENT_RESIDENCE_SINCE_BITS = 3;
    static const int PRESENT_RESIDENCE_SINCE_SCALE = 1;

    static const int AGE = 12;
    static const int AGE_BITS = 5;
    static const int AGE_SCALE = 10;
    static const int AGE_MINUS = 10;

    static const int CREDITS_IN_BANK = 15;
    static const int CREDITS_IN_BANK_BITS = 2;
    static const int CREDITS_IN_BANK_SCALE = 1;

    static const int MAINTAINS_PEOPLE = 17;
	static const int MAINTAINS_PEOPLE_BITS = 1;
	static const int MAINTAINS_PEOPLE_SCALE = 1;
	static const int MAINTAINS_PEOPLE_OFFSET = 0;

	size_t good_idx() {
		return discrete_variables.size();
	}

	size_t numeric_variable_bitcount() {
	    return DURATION_BITS
	         + CREDITS_BITS
	    	 + INSTALLMENT_RATE_BITS
	         + PRESENT_RESIDENCE_SINCE_BITS
	         + AGE_BITS
	         + CREDITS_IN_BANK_BITS
		     + MAINTAINS_PEOPLE_BITS;
	}

	size_t variable_bitcount() {
		return discrete_variables.size() + 1 + numeric_variable_bitcount();
	}

	std::vector<int> discrete_variable_groups() {
		// figure out the groups of different attributes
		std::vector<int> groups;
		for (size_t i = 0; i < discrete_variables.size(); ++i) {
			const auto& var = discrete_variables[i];
			// the group is encoded in the variable name itself
			const std::string& label = var.first;
			if (label == "A410") {
				// special case
				groups.push_back(4);
			} else {
				std::string attr;
				for (size_t j = 1; j < label.size()-1; ++j) {
					attr += label[j];
				}
				groups.push_back(atoi(attr.c_str()));
			}
		}
		return groups;

	}
	void add_groups(std::vector<int>& groups, int& lastGroup, int bits) {
		++lastGroup;
		for (int i = 0; i < bits; ++i) {
			groups.push_back(lastGroup);
		}
	}
	std::vector<int> variable_groups() {
		std::vector<int> groups = discrete_variable_groups();
		int lastGroup = groups.back();
		groups.push_back(++lastGroup); // for good
		add_groups(groups, lastGroup, DURATION_BITS);
		add_groups(groups, lastGroup, CREDITS_BITS);
		add_groups(groups, lastGroup, INSTALLMENT_RATE_BITS);
		add_groups(groups, lastGroup, PRESENT_RESIDENCE_SINCE_BITS);
		add_groups(groups, lastGroup, AGE_BITS);
		add_groups(groups, lastGroup, CREDITS_IN_BANK_BITS);
		add_groups(groups, lastGroup, MAINTAINS_PEOPLE_BITS);
		return groups;
	}

	void setup_lang(explib::lang<germancredit_problem>& lang, undef_t undefs = undef_nothing) {
		typedef germancredit_problem p;
		explib::ctx<p> varctx(explib::cvec<p>(0));

		std::vector<int> groups = variable_groups();
		int group = -1;
		int groupBegin = 0;
		size_t i = 0;
		for (; i < groups.size(); ++i) {
			lang.add_orig(explib::orig<p>(varctx));
			if (group != groups[i]) {
				double pLeft = 1.;
				for (int j = i-1; j >= groupBegin; --j) {
					double prioriP = pLeft / (j + 1 - groupBegin);
					lang.var(j).setPrioriP(prioriP);
					bool udef = false;
					if (i < discrete_variables.size()) {
						udef = undefs & undef_excluded_discrete_vars;
					} else {
						udef = undefs & undef_excluded_numeric_vars;
					}
					if (!udef) {
						pLeft -= prioriP;
					}
				}
				group = groups[i];
				groupBegin = i;
			}
		}
		while (i < variable_bitcount()) {
			lang.add_orig(explib::orig<p>(varctx));
			++i;
		}
		lang.var(good_idx()).setPrioriP(0.5); // binary variable, not one-in-many option
		for (int i = 0; i < lang.orig_count(); ++i) {
			if (lang.var(i).prioriP() < 1.0) {
				for (int j = 0; j < i; ++j) {
					if (lang.var(j).prioriP() < 1.0) {
						explib::rel<p>& rl(lang.alloc_rel(varctx));
						rl.add_var(explib::cvec<p>(0), lang.var(i));
						rl.add_var(explib::cvec<p>(0), lang.var(j));
						lang.rel_done();
					}
				}
			}
		}
	}

	static const bool FILL_BITS = false;

/*	void writeIntInLogBits(explib::cond_bits& bits, int& offset, int value, int bitCount, int scale, bool undefineExcluded) {
		value /= scale;
		for (int i = 0; i < bitCount; i++) {
			int v = (value>>(bitCount-i-1));
			if (FILL_BITS) {
				*bits[i+offset] = v>0;
			} else {
				if (v == 1) {
					*bits[i+offset] = true;
				} else if (undefineExcluded && !v) {
					bits[i+offset] = false;
				}
			}
		}
		offset += bitCount;
	}

	void writeIntInDivBits(explib::cond_bits& to, int& offset, int value, int bitCount, int scale, bool undefineExcluded) {
		value /= scale;
		for (int i = 0; i < bitCount; i++) {
			int c = value - (bitCount-i);
			if (FILL_BITS) {
				*to[i + offset] = c > 0;
			} else {
				if (c == 1) {
					*to[i+ offset] = true;
				} else if (undefineExcluded && !c) {
					to[i+offset] = false;
				}
			}
		}
		offset += bitCount;
	}*/

	void writeIntInLogBits(explib::cond_bits& bits,
						  int& offset,
						  double value,
						  int bitCount,
						  double scale,
						  bool undefineExcluded) {
		value /= scale;
		for (int i = 0; i < bitCount; i++) {
			int v = (int(value)>>(bitCount-i-1));
			if (FILL_BITS) {
				*bits[i+offset] = v>0;
			} else {
				if (v == 1 || (i+1 == bitCount && v <= 0)) {
//				if (v == 1) {
					*bits[i+offset] = true;
				} else if (undefineExcluded && v <= 0) {
					bits[i+offset] = false;
				}
			}
		}
		offset += bitCount;
	}

	void writeIntInDivBits(explib::cond_bits& to,
						   int& offset,
						   double value,
						   int bitCount,
						   double scale,
						   bool undefineExcluded) {
		value /= scale;
		for (int i = 0; i < bitCount; i++) {
			int c = value - (bitCount-i-1);
			if (FILL_BITS) {
				*to[i + offset] = c > 0;
			} else {
//				if (c == 1) {
				if (c == 1 || (i+1 == bitCount && c <= 0)) {
					*to[i+ offset] = true;
				} else if (undefineExcluded && c <= 0) {
					to[i+offset] = false;
				}
			}
		}
		offset += bitCount;
	}

	void populate(explib::data<germancredit_problem>& data,
				  std::istream& in,
				  undef_t undefs = undef_nothing) {
		typedef germancredit_problem p;

		explib::cvec<p> at;

		for (int i = 0; i < data.lang().orig_count(); ++i) {
			data.var(i).defined().fill(true);
		}

		std::vector<int> groups = discrete_variable_groups();

		for (size_t i = 0; i < SampleCount; ++i) {
			std::vector<std::string> items( exptest::read_split_line(in) );
			at[0] = i;
			for (size_t j = 0; j < items.size(); ++j) {
				std::string v = items[j];
				for (size_t k = 0; k < discrete_variables.size(); ++k) {
					const auto& var = discrete_variables[k];
					if (v == var.first) {
						*data.var(k)[at] = true;
						if (undefs & undef_excluded_discrete_vars) {
							int group = groups[k];
							for (size_t l = k-1; l >= 0 && groups[l] == group; --l) {
								data.var(l)[at] = false; // mark undefined
							}
						}
						break;
					}
				}
			}
			*data.var(good_idx())[at] = (items.back() == "1");

			bool undefNumeric = undefs & undef_excluded_numeric_vars;

			explib::cond_bits bits;
			bits.resize(numeric_variable_bitcount());
			bits.defined().fill(true);
			int offset = 0;
			writeIntInLogBits(bits,
							  offset,
							  atoi(items[DURATION].c_str()) - DURATION_OFFSET,
							  DURATION_BITS, DURATION_SCALE,
							  undefNumeric);
			writeIntInLogBits(bits,
					          offset,
					          atoi(items[CREDITS].c_str()),
					          CREDITS_BITS,
					          CREDITS_SCALE,
					          undefNumeric);
			writeIntInDivBits(bits,
							  offset,
							  atoi(items[INSTALLMENT_RATE].c_str()),
							  INSTALLMENT_RATE_BITS,
							  INSTALLMENT_RATE_SCALE,
							  undefNumeric);
			writeIntInDivBits(bits,
							  offset,
							  atoi(items[PRESENT_RESIDENCE_SINCE].c_str()),
							  PRESENT_RESIDENCE_SINCE_BITS,
							  PRESENT_RESIDENCE_SINCE_SCALE,
							  undefNumeric);
			writeIntInDivBits(bits,
							  offset,
							  atoi(items[AGE].c_str())-AGE_MINUS,
							  AGE_BITS,
							  AGE_SCALE,
							  undefNumeric);
			writeIntInDivBits(bits,
							  offset,
							  atoi(items[CREDITS_IN_BANK].c_str()),
							  CREDITS_IN_BANK_BITS,
							  CREDITS_IN_BANK_SCALE,
							  undefNumeric);
			writeIntInDivBits(bits,
							  offset,
							  atoi(items[MAINTAINS_PEOPLE].c_str()) - MAINTAINS_PEOPLE_OFFSET,
							  MAINTAINS_PEOPLE_BITS,
							  MAINTAINS_PEOPLE_SCALE,
							  undefNumeric);
			int numeric_begin = discrete_variables.size() + 1; // for good
			for (int j = 0; j < bits.size(); j++) {
				int var = numeric_begin + j;
				data.var(var)[at] = bool(bits[j]);
				*data.var(var)[at] = bool(*bits[j]);
			}
		}
	}

	void add_numericlabels(explib::pinfo& names, const char* title, int bits) {
		for (int i = 0; i < bits; ++i) {
			names.vnames_.push_back(sup()<<title<<i);
		}

	}
	void setup_names(explib::pinfo& names) {
		for (const auto& var : discrete_variables) {
			names.vnames_.push_back(var.second);
		}
		names.vnames_.push_back("good");
		add_numericlabels(names, "duration", DURATION_BITS);
		add_numericlabels(names, "credits", CREDITS_BITS);
		add_numericlabels(names, "installmentRate", INSTALLMENT_RATE_BITS);
		add_numericlabels(names, "presentResidence", PRESENT_RESIDENCE_SINCE_BITS);
		add_numericlabels(names, "age", AGE_BITS);
		add_numericlabels(names, "creditsInBank", CREDITS_IN_BANK_BITS);
		add_numericlabels(names, "maintainsPeople", MAINTAINS_PEOPLE_BITS);
		names.cvnames_.push_back("s"); // sample
	}

	void setup(explib::lang<germancredit_problem>& lang,
			   explib::data<germancredit_problem>& data,
			   explib::pinfo& names,
			   std::istream& in,
			   undef_t undefs = undef_nothing) {
		setup_lang(lang, undefs);
		setup_names(names);
		populate(data, in, undefs);
	}

	struct germancredit {
		typedef germancredit_problem p;
		explib::lang<p> lang_;
		explib::cvec<p> cv_;
		explib::data<p> data_;
		explib::pinfo names_;
		explib::stats<p> stats_;
		explib::stats_info<p> si_;
		germancredit(undef_t undefs = undef_nothing)
		: lang_(), cv_(SampleCount), data_(lang_, cv_), names_(), stats_(data_), si_(names_, stats_) {
			std::ifstream in("test/germancredit/german.data");
			setup(lang_, data_, names_, in, undefs);
			int v = explib::util::next_version();
			for (size_t i = 0; i < data_.var_count(); ++i) {
				data_.var(i).version_ = v;
			}
			stats_.update();
		}
	};

	void vars_test(TestTool& t) {
		std::vector<int> groups = discrete_variable_groups();
		for (size_t i = 0; i < discrete_variables.size(); ++i) {
			t<<discrete_variables[i].first<<" ["<<groups[i]<<"] = "<<discrete_variables[i].second<<"\n";
		}
	}

	void setup_test(TestTool& t) {
		germancredit pr(undef_nothing);

		t<<pr.si_.vars_tostring();
		t<<"\ntop scan:\n";
		t<<pr.si_.scan_tostring(20, 1);
		t<<"\ntop good deps:\n";
		t<<pr.si_.var_deps_tostring(good_idx(), 20);
	}

	void printdata(TestTool& t, const germancredit& cr) {
		std::vector<std::string> labels;
		for (int c = 0; c < 3; ++c) {
			for (size_t i = 0; int(i) < cr.lang_.var_count(); ++i) {
				int n = i;
				for (int d = 2; d > c; --d) {
					n /= 10;
				}
				if (n || c == 2) {
					t<<char('0'+(n%10));
				} else {
					t<<" ";
				}
			}
			t<<"\n";
		}
		for (size_t i = 0; int(i) < cr.lang_.var_count(); ++i) {
			if (int(i) < cr.lang_.orig_count()) {
				t<<"V";
			} else {
				t<<"E";
			}
		}
		t<<"\n";

		explib::cvec<germancredit_problem> at;
		for (int s = 0; s < cr.data_.dim()[0]; ++s) {
			at[0] = s;
			for (int v = 0; v < cr.lang_.var_count(); ++v) {
				const explib::data_var<germancredit_problem>& dv = cr.data_.var(v);
				if (dv[at]) {
					if (*dv[at]) {
						t<<"1";
					} else {
						t<<"0";
					}
				} else {
					t<<".";
				}
			}
			t<<" s"<<s<<"\n";
		}
	}

	void printdata_test(TestTool& t) {
		germancredit pr(undef_nothing);
		printdata(t, pr);

		t<<"\n"<<pr.si_.vars_tostring();
	}

	void printdata_undefs_test(TestTool& t) {
		germancredit pr(undef_all_excluded_vars);
		printdata(t, pr);

		t<<"\n"<<pr.si_.vars_tostring();
	}

	void applynexps_test(TestTool& t, int exps, undef_t undefs = undef_nothing) {
		germancredit pr(undefs);

//		t<<"\n"<<pr.si_.rels_tostring();

		explib::learner<germancredit_problem> learner(pr.lang_, pr.stats_, 70., 3.);
		for (int i = 0; i < exps; ++i) learner.add_exp();
		printdata(t, pr);

		t<<"\n";
		for (int i = 0; i < pr.lang_.exp_count(); ++i) {
			t<<pr.si_.var_tostring(pr.lang_.exp(i).id());
		}

//		t<<"\n"<<pr.si_.rels_tostring();
	}

	void apply1exp_test(TestTool& t) {
		applynexps_test(t, 1);
	}

	void apply2exp_test(TestTool& t) {
		applynexps_test(t, 2);
	}

	void apply3exp_test(TestTool& t) {
		applynexps_test(t, 3);
	}

	void apply4exp_test(TestTool& t) {
		applynexps_test(t, 4);
	}

	void apply10exp_test(TestTool& t) {
		applynexps_test(t, 10);
	}

	void apply1exp_undefs_test(TestTool& t) {
		applynexps_test(t, 1, undef_all_excluded_vars);
	}

	void apply2exp_undefs_test(TestTool& t) {
		applynexps_test(t, 2, undef_all_excluded_vars);
	}

	void apply3exp_undefs_test(TestTool& t) {
		applynexps_test(t, 3, undef_all_excluded_vars);
	}

	void apply4exp_undefs_test(TestTool& t) {
		applynexps_test(t, 4, undef_all_excluded_vars);
	}

	void apply10exp_undefs_test(TestTool& t) {
		applynexps_test(t, 10, undef_all_excluded_vars);
	}

	void setupundefs_test(TestTool& t) {
		germancredit pr(undef_all_excluded_vars);

		t<<pr.si_.vars_tostring();
		t<<"\ntop scan:\n";
		t<<pr.si_.scan_tostring(20, 1);
		t<<"\ntop good deps:\n";
		t<<pr.si_.var_deps_tostring(good_idx(), 20);
	}

	void prioriesnoundef_test(TestTool& t) {
		germancredit pr(undef_nothing);

		std::vector<int> groups = variable_groups();

		for (int i = 0; i < pr.lang_.var_count(); ++i) {
			const explib::var<germancredit_problem>& v( pr.lang_.var(i) );
			t<<"prioriP="<<v.prioriP()<<"  "<<"group:"<<groups[i]<<"  "<<pr.si_.var_tostring(i)<<"\n";
		}
	}

	void priories_test(TestTool& t) {
		germancredit pr(undef_all_excluded_vars);

		std::vector<int> groups = variable_groups();

		for (int i = 0; i < pr.lang_.var_count(); ++i) {
			const explib::var<germancredit_problem>& v( pr.lang_.var(i) );
			t<<"prioriP="<<v.prioriP()<<"  "<<"group:"<<groups[i]<<"  "<<pr.si_.var_tostring(i)<<"\n";
		}
		t<<"naive entropy: "<<pr.stats_.naiveInfo()<<"\n";
	}

	void exppriories_test(TestTool& t) {
		germancredit pr(undef_all_excluded_vars);
		explib::learner<germancredit_problem> learner(pr.lang_, pr.stats_, 70., 3.);
		learner.reexpress(true);

		std::vector<int> groups = variable_groups();
		for (int i = 0; i < pr.lang_.var_count(); ++i) {
			const explib::var<germancredit_problem>& v( pr.lang_.var(i) );
			int group = -1;
			if (size_t(i) < groups.size()) group = groups[i];
			t<<"prioriP="<<v.prioriP()<<"  "<<"group:"<<group<<"  "<<pr.si_.var_tostring(i)<<"\n";
		}
		t<<"naive entropy: "<<pr.stats_.naiveInfo()<<"\n";
	}

	void logdeps_test(TestTool& t) {
		germancredit pr(undef_all_excluded_vars);

		t<<pr.si_.row_logdep_tostring(good_idx());
	}

	void explogdeps_test(TestTool& t) {
		germancredit pr(undef_all_excluded_vars);
		explib::learner<germancredit_problem> learner(pr.lang_, pr.stats_, 70., 3.);
		learner.reexpress(true);

		t<<pr.si_.row_logdep_tostring(good_idx());
	}


	void setup_preproblem(pred_problem<germancredit_problem>& pr,
						  undef_t undefs = undef_all_excluded_vars) {
		pr.predvars_.resize(good_idx()+1);
		pr.predvars_[good_idx()] = true;
		std::ifstream in("test/germancredit/german.data");
		setup(pr.lang_, pr.data_, pr.names_, in, undefs);
		pr.costs_[true_positive] = 0;
		pr.costs_[true_negative] = 0;
		pr.costs_[false_positive] = 5;
		pr.costs_[false_negative] = 1;
	}

	void naivepredict_test(TestTool& t) {
		explib::cvec<germancredit_problem> dim(SampleCount);
		pred_problem<germancredit_problem> p(dim);
		setup_preproblem(p);

		pred_args args(1000, 1000, 0, 2., true);
		crossvalidate_run(t, p, args, 10);
	}

	void verboserun_test(TestTool& t) {
		explib::cvec<germancredit_problem> dim(SampleCount);
		pred_problem<germancredit_problem> p(dim);
		setup_preproblem(p);

		pred_args args(1000, 1000, 0, 2., 2);
		pred_stats stats;

		explib::bits testsamples;
		testsamples.resize(1000);
		for (int j = 0; j < 100; ++j) {
			testsamples[900 + j] = true;
		}
		teach_test_measure(t, p, testsamples, args, stats);

		if (args.verbose_) {
			t<<"exp treshold: "<<args.expthreshold_<<"\n";
			t<<"exp rel filter: "<<args.exprelfilter_<<"\n\n";
			t<<"pred rel filter: "<<args.predrelfilter_<<"\n\n";
			t<<"\nfor train data:\n";
			t<<"entropy: "<<stats.train_.unitEntropy()<<"/"<<stats.train_.unitNaiveEntropy()<<"\n";
			t<<"error rate: "<<stats.train_.errorRate()<<"/"<<stats.train_.naiveErrorRate()<<"\n\n";
			t<<"cost: "<<stats.train_.unitCost()<<"/"<<stats.train_.unitNaiveCost()<<"\n\n";
			t<<"\nfor test data:\n";
			t<<"entropy: "<<stats.test_.unitEntropy()<<"/"<<stats.test_.unitNaiveEntropy()<<"\n";
			t<<"error rate: "<<stats.test_.errorRate()<<"/"<<stats.test_.naiveErrorRate()<<"\n\n";
			t<<"cost: "<<stats.test_.unitCost()<<"/"<<stats.test_.unitNaiveCost()<<"\n\n";
		}
	}

	void singlerun_test(TestTool& t) {
		explib::cvec<germancredit_problem> dim(SampleCount);
		pred_problem<germancredit_problem> p(dim);
		setup_preproblem(p, undef_all_excluded_vars);

		double th = 85;
		pred_args args(th, th*0.2, 5);
		crossvalidate_run(t, p, args, 10);

		Table table(
			t.report(ToTable<Average>({"run:out", "data:train"}, "prop:", "predfilter:")));
		t<<"for train data:\n"<<table<<"\n";

		Table table2(
			t.report(ToTable<Average>({"run:out", "data:test"}, "prop:", "predfilter:")));
		t<<"for test data:\n"<<table2<<"\n";
	}

	void genericfilter_test(TestTool& t, double filterbegin, double filterstep, int filtersteps) {
		explib::cvec<germancredit_problem> dim(SampleCount);
		pred_problem<germancredit_problem> p(dim);
		setup_preproblem(p);

		for (int i = 0; i < filtersteps; ++i) {
			pred_args args(1000, 1000, filterbegin + double(i)*filterstep);
			crossvalidate_run(t, p, args, 10);
		}

		Table exptable(
			t.report(ToTable<Average>({"run:out"}, "exp:", "predfilter:")));
		t.ignored()<<"expression & relations:\n"<<exptable<<"\n";

		Table traintable(
			t.report(ToTable<Average>({"run:out", "data:train"}, "prop:", "predfilter:")));
		t<<"train:\n"<<traintable<<"\n";

		Table testtable(
			t.report(ToTable<Average>({"run:out", "data:test"}, "prop:", "predfilter:")));
		t<<"test:\n"<<testtable<<"\n";

		t<<"entropy:\n";

		Table table2(
			t.report(ToTable<Average>({"run:out", "prop:entropy"}, "data:", "predfilter:")));

		t<<table2.toplot(2, 20, 0.7)<<"\n";

		t<<"err:\n";

		Table table3(
			t.report(ToTable<Average>({"run:out", "prop:err"}, "data:", "predfilter:")));

		t<<table3.toplot(2, 20, 0.2)<<"\n";

		t<<"cost:\n";

		Table table4(
			t.report(ToTable<Average>({"run:out", "prop:cost"}, "data:", "predfilter:")));

		t<<table4.toplot(2, 20, 0.5)<<"\n";

		Table table5(
			t.report(ToTable<Average>({"run:out", "perf:ns"}, "data:", "predfilter:")));

		t<<"\nperformance (ns / entry)\n";
		t.ignored()<<table5;

		t<<"\nns by filter\n";

		t.ignored()<<
			Table(t.report(ToTable<Average>({"run:out", "perf:ns"}, "data:", "predfilter:")))
				.toplot(2, 20, 0)<<"\n";
	}

	void narrowfilter_test(TestTool& t) {
		genericfilter_test(t, 0, 0.5, 30);
	}

	void filter_test(TestTool& t) {
		genericfilter_test(t, 0, 5, 30);
	}

	void generic_exps_test(TestTool& t, double expthbegin, double expthstep, int expthsteps,
										double filterbegin, double filterstep, int filtersteps,
										undef_t undefs = undef_all_excluded_vars) {
		explib::cvec<germancredit_problem> dim(SampleCount);
		pred_problem<germancredit_problem> p(dim);
		setup_preproblem(p, undefs);

		for (int i = 0; i < expthsteps; ++i) {
			double th = expthbegin + double(i*expthstep);
			t<<"threshold: "<<th<<"\n\n";
			for (int j = 0; j < filtersteps; ++j) {
				pred_args args(th, th*0.2, filterbegin+double(j) * filterstep);
				crossvalidate_run(t, p, args, 10);
			}
			std::string thstr = sup()<<"threshold:"<<th;

			Table exptable(
				t.report(ToTable<Average>({thstr, "run:out"}, "exp:", "predfilter:")));
			t.ignored()<<"expression & relations:\n"<<exptable<<"\n";
			Table train(
				t.report(ToTable<Average>({thstr, "data:train", "run:out"}, "prop:", "predfilter:")));
			t<<"train:\n"<<train<<"\n";
			Table test(
				t.report(ToTable<Average>({thstr, "data:test", "run:out"}, "prop:", "predfilter:")));
			t<<"test:\n"<<test<<"\n";
		}
		std::set<std::string> traintags;
		traintags.insert("run:out");
		traintags.insert("data:train");
		std::set<std::string> testtags;
		testtags.insert("run:out");
		testtags.insert("data:test");

		Table traintable(
			t.report(ToTable<Average>(traintags, "prop:", "threshold:")));
		t<<"train data:\n"<<traintable<<"\n";

		Table testtable(
			t.report(ToTable<Average>(testtags, "prop:", "threshold:")));
		t<<"test data:\n"<<testtable<<"\n";

		t<<"\ntrain data entropy:\n";

		Table table2(
			t.report(ToTable<Average>(traintags+"prop:entropy", "threshold:", "predfilter:")));

		t<<table2.toplot(2, 20, 0.74)<<"\n";

		t<<"following statistics are for test data only:";
		t<<"\nentropy:\n";

		t<<t.report(ToTable<Average>(testtags+"prop:entropy", "threshold:", "predfilter:")).toplot(2, 20, 0.74)<<"\n";

		t<<"\nentropy by threshold:\n";

		t<<t.report(ToTable<Average>(testtags+"prop:entropy", "predfilter:", "threshold:"))
			.toplot(2, 20, 0.74)<<"\n";

		t<<"\nerr:\n";

		Table table3(
			t.report(ToTable<Average>(testtags+"prop:err", "threshold:", "predfilter:")));

		t<<table3.toplot(2, 20, 0.25)<<"\n";

		t<<"\nerr by threshold:\n";

		t<<t.report(ToTable<Average>(testtags+"prop:entropy", "predfilter:", "threshold:"))
			.toplot(2, 20, 0.74)<<"\n";

		t<<"\ncost:\n";

		Table table5(
			t.report(ToTable<Average>(testtags+"prop:cost", "threshold:", "predfilter:")));

		t<<table5.toplot(2, 20, 0.5)<<"\n";

		t<<"\ncost by threshold:\n";

		Table table6(
			t.report(ToTable<Average>(testtags+"prop:cost", "predfilter:", "threshold:")));

		t<<table6.toplot(2, 20, 0.5)<<"\n";

		Table table7(
			t.report(ToTable<Average>(traintags+"perf:ns", "threshold:", "predfilter:")));

		t<<"\ntrain performance (ns / entry)\n";
		t.ignored()<<table7;

		Table table8(
			t.report(ToTable<Average>(testtags+"perf:ns", "threshold:", "predfilter:")));

		t<<"\ntest performance (ns / entry)\n";
		t.ignored()<<table8;

		t<<"\nns:\n";

		t.ignored()<<
			Table(t.report(ToTable<Average>(testtags+"perf:ns", "threshold:", "predfilter:")))
				.toplot(2, 20, 0)<<"\n";
	}

	void exps_test(TestTool& t) {
		generic_exps_test(t, 50, 20, 8, 0, 1, 10);
	}

	void exps_noundef_test(TestTool& t) {
		generic_exps_test(t, 50, 20, 8, 0, 1, 10, undef_all_excluded_vars);
	}

	void exps_nonumundef_test(TestTool& t) {
		generic_exps_test(t, 50, 20, 8, 0, 1, 10, undef_excluded_discrete_vars);
	}

	void generic_wideexps_test(TestTool& t, undef_t undefs, double predfilter) {
		explib::cvec<germancredit_problem> dim(SampleCount);
		pred_problem<germancredit_problem> p(dim);
		setup_preproblem(p, undefs);

		for (int i = 0; i < 25; ++i) {
			double th = 50 + double(i*12.5);
			pred_args args(th, 0.2*th, predfilter);
			crossvalidate_run(t, p, args, 10);
		}

		std::set<std::string> tags;
		tags.insert("run:out");

		Table exptable(
			t.report(ToTable<Average>({"run:out"}, "exp:", "threshold:")));
		t.ignored()<<"expression & relations:\n"<<exptable<<"\n";

		Table traintable(
			t.report(ToTable<Average>(tags+"data:train", "prop:", "threshold:")));
		t<<"train data:\n"<<traintable<<"\n";

		Table testtable(
			t.report(ToTable<Average>(tags+"data:test", "prop:", "threshold:")));
		t<<"test data:\n"<<testtable<<"\n";

		t<<"\nentropy by threshold:\n";

		t<<t.report(ToTable<Average>(tags+"prop:entropy", "data:", "threshold:"))
			.toplot(2, 20, 0.74)<<"\n";

		t<<"\nerr by threshold:\n";

		t<<t.report(ToTable<Average>(tags+"prop:err", "data:", "threshold:"))
			.toplot(2, 20, 0.74)<<"\n";

		t<<"\ncost by threshold:\n";

		Table table6(
			t.report(ToTable<Average>(tags+"prop:cost", "data:", "threshold:")));

		t<<table6.toplot(2, 20, 0.5)<<"\n";

		Table table7(
			t.report(ToTable<Average>(tags+"perf:ns", "data:", "threshold:")));

		t<<"\nns / entry\n";
		t.ignored()<<table7;

		Table table8(
			t.report(ToTable<Average>(tags+"perf:ns", "data:", "threshold:")));

		t<<"\nns / entrys:\n";

		t.ignored()<<
			Table(t.report(ToTable<Average>(tags+"perf:ns", "data:", "threshold:")))
				.toplot(2, 20, 0)<<"\n";

	}

	void wideexps_test(TestTool& t) {
		generic_wideexps_test(t, undef_all_excluded_vars, 5);
	}

	void wideexps_noundef_test(TestTool& t) {
		generic_wideexps_test(t, undef_nothing, 0);
	}

	void wideexps_nonumundef_test(TestTool& t) {
		generic_wideexps_test(t, undef_excluded_discrete_vars, 4);
	}

	void peak_test(TestTool& t) {
		generic_exps_test(t, 80, 5, 7, 0, 0.5, 5);
	}

	void entropypeak_test(TestTool& t) {
		generic_exps_test(t, 50, 10, 10, 1, 2, 6);
	}

}

void addgermancredittest(TestRunner& runner) {
	runner.add("germancredit/vars", {"func"}, &vars_test);
	runner.add("germancredit/setup", {"func"}, &setup_test);

	runner.add("germancredit/printdata", {"func", "apply"}, &printdata_test);
	runner.add("germancredit/apply1exp", {"func", "apply"}, &apply1exp_test);
	runner.add("germancredit/apply2exp", {"func", "apply"}, &apply2exp_test);
	runner.add("germancredit/apply3exp", {"func", "apply"}, &apply3exp_test);
	runner.add("germancredit/apply4exp", {"func", "apply"}, &apply4exp_test);
	runner.add("germancredit/apply10exp", {"func", "apply"}, &apply10exp_test);
	runner.add("germancredit/printdata_undefs", {"func", "apply"}, &printdata_undefs_test);
	runner.add("germancredit/apply1exp_undefs", {"func", "apply"}, &apply1exp_undefs_test);
	runner.add("germancredit/apply2exp_undefs", {"func", "apply"}, &apply2exp_undefs_test);
	runner.add("germancredit/apply3exp_undefs", {"func", "apply"}, &apply3exp_undefs_test);
	runner.add("germancredit/apply4exp_undefs", {"func", "apply"}, &apply4exp_undefs_test);
	runner.add("germancredit/apply10exp_undefs", {"func", "apply"}, &apply10exp_undefs_test);

	runner.add("germancredit/prioriesnoundef", {"func"}, &prioriesnoundef_test);
	runner.add("germancredit/priories", {"func"}, &priories_test);
	runner.add("germancredit/exppriories", {"func"}, &exppriories_test);
	runner.add("germancredit/logdeps", {"func"}, &logdeps_test);
	runner.add("germancredit/explogdeps", {"func"}, &explogdeps_test);
	runner.add("germancredit/setupundefs", {"func"}, &setupundefs_test);
	runner.add("germancredit/naivepredict", {"func"}, &naivepredict_test);
	runner.add("germancredit/singlerun", {"func"}, &singlerun_test);
	runner.add("germancredit/verboserun", {"func"}, &verboserun_test);
	runner.add("germancredit/narrowfilter", {"perf"}, &narrowfilter_test);
	runner.add("germancredit/filter", {"perf"}, &filter_test);
	runner.add("germancredit/exps", {"perf"}, &exps_test);
	runner.add("germancredit/wideexps", {"perf"}, &wideexps_test);
	runner.add("germancredit/peak", {"perf"}, &peak_test);
	runner.add("germancredit/entropypeak", {"perf"}, &entropypeak_test);

	runner.add("germancredit/exps_nonumundef", {"perf"}, 	 &exps_nonumundef_test);
	runner.add("germancredit/wideexps_nonumundef", {"perf"}, &wideexps_nonumundef_test);

	runner.add("germancredit/exps_noundef", {"perf"}, &exps_noundef_test);
	runner.add("germancredit/wideexps_noundef", {"perf"}, &wideexps_noundef_test);

}

