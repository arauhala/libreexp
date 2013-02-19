/*
 * varstest.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"

namespace {

	struct vars_problem {
		static const int DIM = 1;
		static const int MAX_REL_VARS = 2; // two max relation variables
	};


	class ivarstestcase {
		public:
			virtual int var_count() const = 0;
			virtual const char* var_name(int i) const = 0;

			virtual bool is_output(int i) const = 0;
			bool is_input(int i) const  {
				return !is_output(i);
			}
			virtual double p(int var, const reexp::cond_bits& prevstates) const = 0;
			bool val(int i, const reexp::cond_bits& states) const {
				return (rand() / double(RAND_MAX)) < p(i, states);
			}
			virtual void populate_states(reexp::cond_bits& states) const{
				for (int i = 0; i < var_count(); ++i) {
					states[i] = true;
					*(states[i]) = val(i, states);
				}
			}
			// turn explib::cond_bits to cond_bits
			reexp::cond_bits states() const {
				reexp::cond_bits states;
				states.resize(var_count());
				populate_states(states);
				return states;
			}
	};

	class xorproblem : public ivarstestcase {
		public:
			int var_count() const {
				return 3;
			}
			const char* var_name(int i) const {
				switch (i) {
					case 0: return "a";
					case 1: return "b";
					case 2: return "x";
				}
				return 0;
			}
			bool is_output(int i) const {
				return i == 2;
			}
			double p(int i, const reexp::cond_bits& states) const {
				if (i < 2) return 0.5;
				else return *states[0] xor *states[1] ? 1. : 0.;
			}
	};

	class andproblem : public ivarstestcase {
		public:
			int var_count() const {
				return 3;
			}
			const char* var_name(int i) const {
				switch (i) {
					case 0: return "a";
					case 1: return "b";
					case 2: return "x";
				}
				return 0;
			}
			bool is_output(int i) const {
				return i == 2;
			}
			double p(int i, const reexp::cond_bits& states) const {
				if (i < 2) return 0.5;
				else return *states[0] && *states[1] ? 1. : 0.;
			}
	};

	class orproblem : public ivarstestcase {
		public:
			int var_count() const {
				return 3;
			}
			const char* var_name(int i) const {
				switch (i) {
					case 0: return "a";
					case 1: return "b";
					case 2: return "x";
				}
				return 0;
			}
			bool is_output(int i) const {
				return i == 2;
			}
			double p(int i, const reexp::cond_bits& states) const {
				if (i < 2) return 0.5;
				else return *states[0] || *states[1] ? 1. : 0.;
			}
	};

	class or4problem : public ivarstestcase {
		public:
			int var_count() const {
				return 5;
			}
			const char* var_name(int i) const {
				switch (i) {
					case 0: return "a";
					case 1: return "b";
					case 2: return "c";
					case 3: return "d";
					case 4: return "x";
				}
				return 0;
			}
			bool is_output(int i) const {
				return i == 4;
			}
			bool val(int i, const reexp::cond_bits& states) const {
				return states[0] || states[1] || states[2] || states[3];
			}
			double p(int i, const reexp::cond_bits& states) const {
				if (is_input(i)) return 0.5;
				else return *states[0] || *states[1] || *states[2] || *states[3] ? 1. : 0.;
			}
	};

	class redundantproblem : public ivarstestcase {
		private:
			double p_;
			std::vector<std::string> constantInputs_;
		public:
			redundantproblem(double p, int constantInputs = 0)
			: p_(p), constantInputs_()  {
				for (int i = 0; i < constantInputs; ++i) {
					constantInputs_.push_back(sup()<<"c"<<i);
				}
			}
			int var_count() const {
				return constantInputs_.size() + 3;
			}
			const char* var_name(int i) const {
				switch (i) {
					case 0: return "a1";
					case 1: return "a2";
					case 2: return "x";
				}
				i-=3;
				if (i < int(constantInputs_.size())) {
					return constantInputs_[i].c_str();
				}
				return 0;
			}
			double p(int i, const reexp::cond_bits& states) const {
				if (i == 0) return 0.5;
				if (i == 1) return *states[0]?1.:0.;
				if (i == 2) return *states[0]?0.75:0.25;
				return 0.; // constant
			}
			bool is_output(int i) const {
				return i == 2;
			}
	};

	/**
	 * uniform distribution for input variables.
	 */
	double rand_dep(double lvar, double p1, double p2) {
		double ldep = 0;
		while (rand() % 20) {
			ldep++;
		}
		ldep = 1 + lvar*(ldep/20.);
		if (rand() % 2) ldep = 1. / ldep;
		ldep = std::min(ldep, 1./p1);
		ldep = std::min(ldep, 1./p2);
		ldep = std::max(ldep, (1-(1-p1)/p2)/p1);
		ldep = std::max(ldep, (1-(1-p2)/p1)/p2);
		return ldep;
	}
	class dep_problem : public ivarstestcase {
		private:
			std::vector<std::string> input_;
			std::vector<double> ps_;
			std::vector<double> deps_;
			std::vector<bool> hidden_;
		public:
			int depidx(int i, int j) const {
				if (i < j) return depidx(j,i);
				return (i*(i+1))/2+j;
			}
			double& dep(int i, int j) {
				return deps_[depidx(i,j)];
			}
			double dep(int i, int j) const {
				return deps_[depidx(i,j)];
			}
			static double depToNeg(double d, double p) {
				return (1.-d*p) / (1.-p);
			}
			double not_dep(int i, int j, double p) const {
				return depToNeg(dep(i, j), p);
			}
			dep_problem(int vars, double ldep = 1, double inputldep = 0, int hidden = 0) : input_(), ps_(), deps_() {
				for (int i = 0; i < vars; ++i) {
					input_.push_back(sup()<<"in"<<i);
					ps_.push_back(rand() / double(RAND_MAX));
					hidden_.push_back(i < hidden);
				}
				ps_.push_back(rand() / double(RAND_MAX));
				deps_.resize(depidx(vars, vars));
				for (int i = 0; i < vars; ++i) {
					for (int j = 0; j < i; ++j) {
						dep(j, i) = rand_dep(inputldep, ps_[i], ps_[j]);
					}
				}
				for (int i = 0; i < vars; ++i) {
					dep(i, vars) = rand_dep(ldep, ps_[i], ps_[vars]);
				}
				hidden_.push_back(false);
			}
			int var_count() const {
				return ps_.size();
			}
			const char* var_name(int i) const {
				if (i < int(input_.size())) {
					return input_[i].c_str();
				}
				return "x";
			}
			bool is_output(int i) const {
				return i == int(input_.size());
			}
			double p(int v, const reexp::cond_bits& states) const {
				double pTrue = ps_[v];
				double pFalse = 1-ps_[v];
				for (int i = 0; i < v; ++i) {
					if (*states[i]) {
						pTrue *= dep(v, i);
						pFalse *= depToNeg(dep(v, i), ps_[v]);
					} else {
						double ndep = not_dep(v, i, ps_[i]);
						pTrue *= ndep;
						pFalse *= depToNeg(ndep, ps_[v]);
						if (pFalse < 0) {
							pFalse = 0;
//							printf("negative pFalse?\n");
						}
					}
				}
				if (pTrue + pFalse == 0) {
//					printf("zero denominator?\n");
					return 0.;
				}
				return pTrue / (pTrue+pFalse);
			}
			double output_entropy(int from, reexp::cond_bits& states) {
				double rv = 0;
				if (is_output(from)) {
					double op = p(from, states);
					rv = reexp::entropy(op);
				} else {
					double varp = p(from, states);
					*states[from] = true;
					rv += varp*output_entropy(from+1, states);
					*states[from] = false;
					rv += (1.-varp)*output_entropy(from+1, states);
				}
				return rv;
			}
			double output_entropy() {
				reexp::cond_bits states;
				states.resize(var_count());
				states.defined().fill(true);
				return output_entropy(0, states);
			}
			double real_p(int var, int from, reexp::cond_bits& states) {
				double rv = 0;
				if (from == var) {
					rv = p(from, states);
				} else {
					double varp = p(from, states);
					*states[from] = true;
					rv += varp*real_p(var, from+1, states);
					*states[from] = false;
					rv += (1.-varp)*real_p(var, from+1, states);
				}
				return rv;
			}
			double real_p(int var) {
				reexp::cond_bits states;
				states.resize(var_count());
				states.defined().fill(true);
				return real_p(var, 0, states);
			}
			double output_naive_entropy() {
				return reexp::entropy(real_p(var_count()-1));
			}
			bool val(int i, const reexp::cond_bits& states) const {
				return (rand() / double(RAND_MAX)) < p(i, states);
			}
			void populate_states(reexp::cond_bits& states) const{
				for (int i = 0; i < var_count(); ++i) {
					states[i] = !hidden_[i];
					*(states[i]) = val(i, states);
				}
			}
			void write(std::ostream& out) const {
				for (int i = 0; i < var_count(); ++i) {
					out<<"p("<<var_name(i)<<")="<<ps_[i]<<"\n";
					for (int j = 0; j < i; ++j) {
						out<<"d("<<var_name(j)<<";"<<var_name(i)<<")="<<dep(i,j)<<"\n";
						out<<"d(!"<<var_name(j)<<";"<<var_name(i)<<")="<<not_dep(i,j, ps_[j])<<"\n";
					}
				}
			}
	};
	class class_problem : public ivarstestcase {
		private:
			std::vector<std::vector<double> > classPs_;
			std::vector<std::string> names_;
			int vars_;
		public:
			class_problem(int classes, int vars) {
				vars_ = vars;
				for (int i = 0; i < classes; ++i) {
					classPs_.push_back(std::vector<double>());
					std::vector<double>& v = classPs_.back();
					for (int j = 0; j < vars+1; ++j) {
						v.push_back(rand() / double(RAND_MAX));
					}
					names_.push_back(sup()<<"class"<<i);
				}
				for (int i = 0; i < vars; ++i) {
					names_.push_back(sup()<<"input"<<i);
				}
				names_.push_back("output");
			}
			int class_count() const {
				return classPs_.size();
			}
			int var_count() const {
				return names_.size();
			}
			const char* var_name(int i) const {
				return names_[i].c_str();
			}
			bool is_output(int i) const {
				return i == int(names_.size()-1);
			}
			double p(int var, const reexp::cond_bits& prevstates) const {
				if (var < int(classPs_.size())) {
					for (int i = 0; i < var; ++i) {
						if (*prevstates[var]) return 0.;
					}
					return 1. / double(classPs_.size()-var);
				}
				int clazz = -1;
				for (int i = 0; i < int(classPs_.size()); ++i) {
					if (*prevstates[i]) clazz = i;
				}
				return classPs_[clazz][var-classPs_.size()];
			}
			double output_entropy() const {
				double rv = 0;
				for (size_t i = 0; i < classPs_.size(); ++i) {
					rv += reexp::entropy(classPs_[i].back()) / double(classPs_.size());
				}
				return rv;
			}

			virtual void populate_states(reexp::cond_bits& states) const{
				int clazz = rand() % classPs_.size();
				for (int i = 0; i < int(classPs_.size()); ++i) {
					*states[i] = (i == clazz);
					states[i] = false; // hide the class
				}
				for (int i = classPs_.size(); i < var_count(); ++i) {
					*states[i] = val(i,states);
					states[i] = true;
				}
			}
			std::string tostring() const {
				std::ostringstream buf;
				for (int i = 0; i < int(classPs_.size()); ++i) {
					buf<<"class"<<i<<"\n";
					for (int j = 0; j < int(classPs_[i].size()); ++j) {
						buf<<"p("<<var_name(j+classPs_.size())<<"|"<<var_name(i)<<")="<<classPs_[i][j]<<", info="<<reexp::entropy(classPs_[i][j])<<"\n";
					}
				}
				return buf.str();
			}
	};

	std::ostream& operator<<(std::ostream& out, const dep_problem& dp) {
		dp.write(out);
		return out;
	}


	template <typename P>
	void populate(reexp::data<P>& data, ivarstestcase& vars, int n) {
		reexp::cvec<P> at;
		for (int i = 0; i < n; ++i) {
			at[0] = i;
			reexp::cond_bits states = vars.states();
			for (int j = 0; j < vars.var_count(); ++j) {
				data.var(j)[at] = bool(states[j]);
				*data.var(j)[at] = bool(*states[j]);
			}
		}
	}
	template <typename P>
	void setup_lang(reexp::lang<P>& lang, ivarstestcase& vars) {
		reexp::ctx<P> varctx(reexp::cvec<P>(0));
		for (int i = 0; i < vars.var_count(); ++i) {
			lang.add_orig(reexp::orig<P>(varctx));
		}
		for (int i = 0; i < vars.var_count(); ++i) {
			for (int j = 0; j < i; ++j) {
				reexp::rel<P>& rl(lang.alloc_rel(varctx));
				rl.add_var(reexp::cvec<P>(0), lang.var(i));
				rl.add_var(reexp::cvec<P>(0), lang.var(j));
				lang.rel_done();
			}
		}
	}

	template <typename P>
	void setup_reg(reexp::lang<P>& lang, reexp::data<P>& data, ivarstestcase& vars, int n) {
		setup_lang(lang, vars);
		populate(data, vars, n);
	}

	void setup_names(reexp::pinfo& names, const ivarstestcase& p) {
		for (int i = 0; i < p.var_count(); ++i) {
			names.vnames_.push_back(p.var_name(i));
		}
		names.cvnames_.push_back("s"); // sample
	}

	void setup_test(test_tool& t) {
		xorproblem pr;

		typedef vars_problem p;

		srand(0);
		int n = 256;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, reexp::cvec<p>(n));
		setup_reg(lang, data, pr, n);

		reexp::stats<p> stats(data);

		reexp::pinfo names;
		setup_names(names, pr);

		reexp::stats_info<p> si(names, stats);
		t<<si.vars_tostring();
		t<<si.rels_tostring();
	}

	double ideal_info(const ivarstestcase& test, const reexp::data<vars_problem>& data, int var) {
		double rv = 0;
		reexp::cond_bits states;
		states.resize(test.var_count());
		states.defined().fill(true);
		int n = data.dim()[0];
		for (int s = 0; s < n; ++s) {
			reexp::cvec<vars_problem> at(s);
			for (int j = 0; j < test.var_count(); ++j) {
				states[j] = bool(data.var(j)[at]);
				*states[j] = bool(*data.var(j)[at]);
			}
			double p = test.p(var, states);
			double info = *states[var]?-log2(p):-log2(1-p);
			if (info != info) { // test for nan
				printf("nan for p=%f state=%d\n", p, int(*states[var]));
			} else {
				rv += info;
			}
		}
		return rv / n;
	}

	double naive_pred_info(const reexp::data<vars_problem>& data, const reexp::stats<vars_problem>& stats, int var) {
		double p = stats.var(var).eP();
		const reexp::data_var<vars_problem>& dv = data.var(var);
		double rv = 0;
		int n = data.dim()[0];
		for (int s = 0; s < n; ++s) {
			rv += *dv[s]?-log2(p):-log2(1-p);
		}
		return rv / n;
	}

	void eval_genpredwith(test_tool& t, const std::set<std::string>& runtags, ivarstestcase& test, int n, int tn, double threshold, bool filter = true, int maxexps = -1, double prioriW = 2., int scaling_groups = -1) {
		std::ostringstream nlabel;
		nlabel<<"n:"<<n;
		std::set<std::string> tags(runtags);
		tags.insert(nlabel.str());

		typedef vars_problem p;

		// setup the teach sample
		reexp::lang<p> lang;
		reexp::data<p> data(lang, reexp::cvec<p>(n));
		size_t ms;
		{
			time_sentry time;
			setup_reg(lang, data, test, n);
			ms = time.ms();
			t.record(tags+"prop:setup ms", ms);
		}
		time_sentry time1;
		reexp::stats<p> stats(data);
		ms = time1.ms();
		t.record(tags+"prop:stats ms", ms);
		reexp::learner<p> learner(lang, stats, threshold);
		for (int i = 0; i < lang.var_count(); ++i) {
			if (test.is_output(i)) {
				learner.exclude(i);
			}
		}
		t.record(tags+"prop:entropy orig", stats.naiveInfo()/n);
		int exps = 0;
		if (threshold >= 0){
			time_sentry time;
			exps = learner.reexpress(filter, maxexps);
			long ms = time.ms();
			t.record(tags+"prop:reexp ms", ms);
		}
		t.record(tags+"prop:entropy reexp", stats.naiveInfo()/n);

		t.record(tags+"prop:exps", double(exps));

		// test sample
		reexp::data<p> tdata(lang, reexp::cvec<p>(tn));
		populate(tdata, test, tn);
		tdata.apply_exps();

		reexp::pinfo names;
		setup_names(names, test);
		reexp::stats_info<p> si(names, stats);

		reexp::pred<p> pred(stats, prioriW);
		for (int i = 0; i < test.var_count(); ++i) {
			if (test.is_output(i)) {
				reexp::group_scaler<p> scaler(pred, i, scaling_groups);

				double totalinfo = 0;
				double entryinfo = 0;
				{
					double idealinfo = ideal_info(test, data, i);
					double naivepredinfo = naive_pred_info(data, stats, i);
					std::vector<double> ps = pred.p(data, i);
					if (scaling_groups > 0) scaler.scale(ps);
					pred.info(ps, data, i, totalinfo, entryinfo);
//					t<<"train total: "<<totalinfo<<ignorel;
//					t<<"train entry: "<<entryinfo<<ignorel;
					t.record(tags+names.vnames_[i]+"sample:train"+"prop:info", totalinfo);
					t.record(tags+names.vnames_[i]+"sample:train"+"prop:entryinfo", entryinfo);
					t.record(tags+names.vnames_[i]+"sample:train"+"prop:naivepredinfo", naivepredinfo);
					t.record(tags+names.vnames_[i]+"sample:train"+"prop:idealinfo", idealinfo);
				}
				{
					double idealinfo = ideal_info(test, tdata, i);
					double naivepredinfo = naive_pred_info(tdata, stats, i);
					std::vector<double> ps = pred.p(tdata, i);
					if (scaling_groups > 0) scaler.scale(ps);
					pred.info(ps, tdata, i, totalinfo, entryinfo);
//					t<<"test total: "<<totalinfo<<ignorel;
//					t<<"test entry: "<<entryinfo<<ignorel;
					t.record(tags+names.vnames_[i]+"sample:test"+"prop:info", totalinfo);
					t.record(tags+names.vnames_[i]+"sample:test"+"prop:entryinfo", entryinfo);
					t.record(tags+names.vnames_[i]+"sample:test"+"prop:naivepredinfo", naivepredinfo);
					t.record(tags+names.vnames_[i]+"sample:test"+"prop:idealinfo", idealinfo);
//					t<<si.pred_tostring(tdata, i, true);
				}
			}
		}
//		t<<".";
	}

	void record_distribution(test_tool& t,
						  	 const std::set<std::string>& runtags,
						  	 ivarstestcase& test, int n, int tn,
						  	 double threshold, bool filter = true,
						  	 int maxexps = -1, double prioriW = 2.,
						  	 int scaling_groups = -1) {
		typedef vars_problem p;

		reexp::lang<p> lang;
		reexp::data<p> data(lang, reexp::cvec<p>(n));
		setup_reg(lang, data, test, n);
		reexp::stats<p> stats(data);
		reexp::learner<p> learner(lang, stats, threshold);
		for (int i = 0; i < lang.var_count(); ++i) {
			if (test.is_output(i)) {
				learner.exclude(i);
			}
		}
		learner.reexpress(filter, maxexps);
		t.record(tags+"prop:entropy reexp", stats.naiveInfo()/n);

		t.record(tags+"prop:exps", double(exps));

		// test sample
		reexp::data<p> tdata(lang, reexp::cvec<p>(tn));
		populate(tdata, test, tn);
		tdata.apply_exps();

		reexp::pinfo names;
		setup_names(names, test);
		reexp::stats_info<p> si(names, stats);

		reexp::pred<p> pred(stats, prioriW);
	}


	void measure_reexp(test_tool& t,
					   const std::set<std::string>& tags,
					   ivarstestcase& test, int n, int tn, double threshold, bool filter = true, int maxexps = -1) {
		std::ostringstream nlabel;
		nlabel<<"n:"<<n;
		std::set<std::string> tags2(tags);
		tags2.insert(nlabel.str());

		typedef vars_problem p;
		// setup the teach sample
		reexp::lang<p> lang;
		reexp::data<p> data(lang, reexp::cvec<p>(n));
		time_sentry times;
		setup_reg(lang, data, test, n);
		long ms = times.ms();
		t.record(tags2+"prop:setup ms", ms);

		time_sentry time1;
		reexp::stats<p> stats(data);
		ms = time1.ms();
		t.record(tags2+"prop:stats ms", ms);
		reexp::learner<p> learner(lang, stats, threshold);
		for (int i = 0; i < lang.var_count(); ++i) {
			if (test.is_output(i)) {
				learner.exclude(i);
			}
		}
		int exps = 0;
		if (threshold >= 0){
			time_sentry time;
			exps = learner.reexpress(filter, maxexps);
			ms = time.ms();
			t.record(tags2+"prop:reexp ms", ms);
		}
		t.record(tags2+"prop:exps", double(exps));

		// test sample
		reexp::data<p> tdata(lang, reexp::cvec<p>(tn));
		time_sentry timep;
		populate(tdata, test, tn);
		ms = timep.ms();
		t.record(tags2+"prop:populate ms", ms);
		time_sentry time2;
		tdata.apply_exps();
		ms = time2.ms();
		t.record(tags2+"prop:apply ms", ms);

		int activeRels = 0;
		for (int i = 0; i < lang.rel_count(); ++i) {
			if (!lang.rel(i).disabled()) activeRels++;
		}
		t.record(tags2+"prop:active rels", activeRels);
		t.record(tags2+"prop:rels", lang.rel_count());

	}


	void vars_test(test_tool& t, ivarstestcase& test) {
		srand(0);

		t<<"example:\n\n";
		for (int s = 0; s < 8; ++s) {
			reexp::cond_bits states = test.states();
			t<<"{";
			for (int v = 0; v < test.var_count(); ++v) {
				if (v) t<<", ";
				t<<(*states[v]?"":"!")<<test.var_name(v);
			}
			t<<"}\n";
		}
		t<<"\n";

		std::set<std::string> naivemode = {"mode:naive"};
		for (int i = 0; i < 8; ++i) {
			eval_genpredwith(t, naivemode, test, 1,  32, -1);
			eval_genpredwith(t, naivemode, test, 2,  32, -1);
			eval_genpredwith(t, naivemode, test, 4,  32, -1);
			eval_genpredwith(t, naivemode, test, 8,  32, -1);
			eval_genpredwith(t, naivemode, test, 16, 32, -1);
			eval_genpredwith(t, naivemode, test, 32, 32, -1);
			eval_genpredwith(t, naivemode, test, 64, 32, -1);
			eval_genpredwith(t, naivemode, test, 128, 32, -1);
			eval_genpredwith(t, naivemode, test, 256, 32, -1);
			eval_genpredwith(t, naivemode, test, 512, 32, -1);
		}

		double threshold = 5;

		std::set<std::string> reexpmode = {"mode:reexp"};
		for (int i = 0; i < 8; ++i) {
			eval_genpredwith(t, reexpmode, test, 1,  32, threshold);
			eval_genpredwith(t, reexpmode, test, 2,  32, threshold);
			eval_genpredwith(t, reexpmode, test, 4,  32, threshold);
			eval_genpredwith(t, reexpmode, test, 8,  32, threshold);
			eval_genpredwith(t, reexpmode, test, 16, 32, threshold);
			eval_genpredwith(t, reexpmode, test, 32, 32, threshold);
			eval_genpredwith(t, reexpmode, test, 64, 32, threshold);
			eval_genpredwith(t, reexpmode, test, 128, 32, threshold);
			eval_genpredwith(t, reexpmode, test, 256, 32, threshold);
			eval_genpredwith(t, reexpmode, test, 512, 32, threshold);
		}

		std::set<std::string> tags;
		tags.insert("run:out");

		table train(
			t.report(to_table<average>(tags+"prop:entryinfo"+"sample:train", "mode:", "n:")));
		table teste(
			t.report(to_table<average>(tags+"prop:entryinfo"+"sample:test", "mode:", "n:")));
		table exps(
			t.report(to_table<average>(tags+"prop:exps", "mode:", "n:")));
		std::ostringstream buf;
		buf<<"for train data:\n\n";
		train>>buf;
		buf<<"\n";
		buf<<"for test data:\n\n";
		teste>>buf;
		buf<<"\n";
		buf<<"exps:\n\n";
		exps>>buf;
		t<<"\nresults:\n\n"<<buf.str()<<"\n";
	}

	void xor_test(test_tool& t) {
		xorproblem pr;
		vars_test(t, pr);
	}

	void and_test(test_tool& t) {
		andproblem pr;
		vars_test(t, pr);
	}

	void or_test(test_tool& t) {
		orproblem pr;
		vars_test(t, pr);

		table teste(
			t.report(to_table<average>({"run:out","prop:entryinfo","sample:test"}, "mode:", "n:")));

		std::ofstream o( t.file_path("test_entropy.tex") );
		o<<teste.to_latex_pgf_doc("crossentropy");
	}

	void or4_test(test_tool& t) {
		or4problem pr;
		vars_test(t, pr);
	}

	void or4_M_test(test_tool& t) {
		or4problem pr;
		measure_reexp(t, {}, pr, 1024*1024, 1024, -1);
		t<<t.report(to_table<average>({"run:out"}, "n:", "prop:"));
	}

	void dists_10var_M_test(test_tool& t) {
		dep_problem pr(10);
		measure_reexp(t, {}, pr, 1024*1024, 1024, -1);
		t<<t.report(to_table<average>({"run:out"}, "n:", "prop:"));
	}

	void dists_20var_M_test(test_tool& t) {
		dep_problem pr(20);
		measure_reexp(t, {}, pr, 1024*1024, 1024, -1);
		t<<t.report(to_table<average>({"run:out"}, "n:", "prop:"));
	}

	void dists_50var_M_test(test_tool& t) {
		dep_problem pr(50);
		measure_reexp(t, {}, pr, 1024*1024, 1024, -1);
		t<<t.report(to_table<average>({"run:out"}, "n:", "prop:"));
	}

	void redundant_test(test_tool& t) {
		redundantproblem pr(0.25);
		vars_test(t, pr);
	}

	void sparseredundant_test(test_tool& t) {
		redundantproblem pr(0.001);
		vars_test(t, pr);
	}

	void redundantperf_test(test_tool& t) {
		int samples = 256*1024;
		int tsamples = 256*1024;
		double threshold = 5;
		int constantInputs = 100;
		{
			redundantproblem test(0.1, constantInputs);
			measure_reexp(t, {"p:0.1"}, test, samples, tsamples, threshold, false, 1);
		}
		{
			redundantproblem test(0.01, constantInputs);
			measure_reexp(t, {"p:0.01"}, test, samples, tsamples, threshold, false, 1);
		}
		{
			redundantproblem test(0.001, constantInputs);
			measure_reexp(t, {"p:0.001"}, test, samples, tsamples, threshold, false, 1);
		}
		{
			redundantproblem test(0.0001, constantInputs);
			measure_reexp(t, {"p:0.0001"}, test, samples, tsamples, threshold, false, 1);
		}
		{
			redundantproblem test(0.00001, constantInputs);
			measure_reexp(t, {"p:0.00001"}, test, samples, tsamples, threshold, false, 1);
		}

		t.ignored()<<t.report(to_table<average>({}, "prop:", "p:"));
	}


	void dists_generic_setup_test(test_tool& t, double indeps, int hidden) {
		srand(0);

		for (int v = 1; v < 5; ++v) {
			dep_problem pr(v, 1., indeps, hidden);

			int n = 10000;

			typedef vars_problem p;

			reexp::lang<p> lang;
			reexp::data<p> data(lang, reexp::cvec<p>(n));
			setup_reg(lang, data, pr, n);

			reexp::stats<p> stats(data);

			reexp::pinfo names;
			setup_names(names, pr);

			t<<"test:\n"<<pr<<"\n";
			reexp::stats_info<p> si(names, stats);
			t<<"stats:\n";
			for (int v = 0; v < lang.var_count(); ++v) {
				t<<"p("<<si.lang_info().var_tostring(v)<<")="<<stats.var(v).p()<<"\n";
			}
			t<<"\n";
			t<<"real:\n";
			for (int i = 0; i < pr.var_count(); ++i) {
				t<<"p("<<si.lang_info().var_tostring(i)<<"): "<<pr.real_p(i)<<"\n";
			}
			t<<"x entropy: "<<pr.output_entropy()<<" / "<<reexp::entropy(pr.real_p(pr.var_count()-1))<<"\n";
			t<<"data info: "<<ideal_info(pr, data, pr.var_count()-1)<<"\n";

			t<<"\n";

			t<<"probabilities for x:\n";
			for (int i = 0; i < (1<<(pr.var_count()-1)); ++i) {
				reexp::cond_bits b;
				b.resize(pr.var_count());
				b.defined().fill(true);
				for (int j = 0; j < pr.var_count()-1; ++j) {
					*b[j] = ((i>>j)&0x1)?true:false;
				}
				t<<"p(x|"<<reexp::vector_todensestring(b.states())<<")="<<pr.p(pr.var_count()-1, b)<<"\n";
			}

			t<<"\nexample samples:\n";

			for (int i = 0; i < 8; ++i) {
				reexp::cond_bits b(pr.states());
				t<<"s:"<<reexp::vector_todensestring(b.states());
				t<<", d:"<<reexp::vector_todensestring(b.defined())<<"\n";
			}
		}

	}

	void dists_setup_test(test_tool& t) {
		dists_generic_setup_test(t, 0., 0);
	}

	void dists_indep_setup_test(test_tool& t) {
		dists_generic_setup_test(t, 1., 0);
	}

	void dists_hidden_setup_test(test_tool& t) {
		dists_generic_setup_test(t, 1., 1);
	}

	void dists_simple_test(test_tool& t) {
		int invars = 8;
		int maxtrainsamples = 512;
		int testsamples = 50;
		int repeats = 30;

		srand(0);
		for (int i = 1; i <= invars; i*=2) {
			for (int r = 0; r < repeats; r++) {
				dep_problem test(invars);
				for (int j = 1; j <= maxtrainsamples ; j*=2) {
					eval_genpredwith(t, {sup()<<"vars:"<<i, sup()<<"n:"<<j}, test, j, testsamples, 10000, true, -1, 2.);
				}
			}
		}

		t<<"this test tests learning with:\n";
		t<<"  1. varying number of independent variables\n";
		t<<"  2. varying number of train data\n\n";

		std::set<std::string> tags = {"sample:test", "run:out"};
		t<<"predicted variable entropy by input variable count and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "vars:", "n:"));

		t<<"\nsame plotted:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "vars:", "n:")).to_plot(3, 20);
	}

	void dists_prioriw_test(test_tool& t) {
		int maxvars = 36;
		int maxprioriw = 8;
		int maxtrainsamples = 512;
		int testsamples = 100;
		int repeats = 50;

		for (int v = 1; v <= maxvars; v*=6) {
			for (double i = 0.5; i <= maxprioriw; i*=2) {
				srand(0);
				for (int r = 0; r < repeats; r++) {
					dep_problem test(v);
					for (int j = 1; j <= maxtrainsamples; j*=2) {
						eval_genpredwith(t, {sup()<<"vars:"<<v, sup()<<"prioriw:"<<i, sup()<<"n:"<<j}, test, j, testsamples, 10000, true, -1, i);
					}
				}
			}
			t<<"\n\nfor "<<v<<" vars:\n\n";
			t<<"this test tests learning with:\n";
			t<<"  1. varying number of priori weight\n";
			t<<"  2. varying number of train data\n\n";

			std::set<std::string> tags = {"sample:test", "run:out", sup()<<"vars:"<<v};
			t<<"predicted variable entropy by prioriw and train data:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "prioriw:", "n:"));

			t<<"\nsame plotted:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "prioriw:", "n:")).to_plot(3, 20);
		}

		std::set<std::string> tags = {"sample:test", "run:out"};
		t<<"predicted variable entropy by prioriw and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "prioriw:", "n:"));

		t<<"\nfor prioriw 1, by vars and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo"+"prioriw:1", "vars:", "n:")).to_plot(3, 20);

		t<<"\nfor prioriw 4, by vars and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo"+"prioriw:4", "vars:", "n:")).to_plot(3, 20);

	}

	void generic_dists_exps_test(test_tool& t, double inputldeps, int hidden = 0) {
		int maxvars = 32;
		int maxthreshold = 1024;
		int minthreshold = 16;
		int thresholdstep = 4;
		int maxtrainsamples = 2048;
		int testsamples = 200;
		int repeats = 50;

		for (int v = 4; v <= maxvars; v*=2, hidden *= 2) {
			for (double i = maxthreshold; i >= minthreshold; i/=thresholdstep) {
				srand(0);
				for (int r = 0; r < repeats; r++) {
					dep_problem test(v, 0.5, inputldeps, hidden);
					for (int j = 1; j <= maxtrainsamples; j*=2) {
						int exps = -1;
						if (i > maxthreshold-1) exps = 0;
						std::string thtag("threshold:");
						if (exps<0) thtag += sup()<<i;
						else thtag += "naive";

						eval_genpredwith(t,
										 {sup()<<"vars:"<<v,
										  thtag,
										  sup()<<"n:"<<j},
										  test, j, testsamples, i, true, exps);
					}
/*					double e = test.outputEntropy();
					double ne = test.outputNaiveEntropy();
					std::set<std::string> tags = {"sample:test", sup()<<"vars:"<<v, sup()<<"threshold:"<<i};
					t.record(tags+"prop:ideal_info", double(e));
					t.record(tags+"prop:naive_info", double(ne));*/
				}
			}
			t<<"\n\nfor "<<v<<" vars:\n\n";
			t<<"this test tests learning with:\n";
			t<<"  1. varying number of threshold\n";
			t<<"  2. varying number of train data\n\n";

			std::set<std::string> tags = {"sample:test", "run:out", sup()<<"vars:"<<v};
			t<<"predicted variable entropy by threshold and train data:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:"));

			t<<"\nsame plotted:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:")).to_plot(3, 30);

			t<<"\nexpression counts by threshold and train data:\n";
			t<<t.report(to_table<average>({"run:out", sup()<<"vars:"<<v, "prop:exps"}, "threshold:", "n:"));
		}

		std::set<std::string> tags = {"sample:test", "run:out"};

		t<<"\naverage ideal variable entropy by threshold and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:idealinfo", "threshold:", "vars:"));

		t<<"\n\npredicted variable entropy by threshold and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:"));

		t<<"\nfor naive bayesian, by vars and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo"+"threshold:naive", "vars:", "n:")).to_plot(3, 20);
		t<<"\nfor threshold "<<minthreshold<<", by vars and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo"+(sup()<<"threshold:"<<minthreshold), "vars:", "n:")).to_plot(3, 20);
	}

	void dists_exps_test(test_tool& t) {
		generic_dists_exps_test(t, 0.15, 0);
	}

	void dists_exps_h_test(test_tool& t) {
		generic_dists_exps_test(t, 0.15, 1);
	}

	void dists_exps_h2_test(test_tool& t) {
		generic_dists_exps_test(t, 0.15, 2);
	}

	void dists_exps_2_test(test_tool& t) {
		generic_dists_exps_test(t, 0.25, 0);
 	}

	void dists_exps_2_h_test(test_tool& t) {
		generic_dists_exps_test(t, 0.25, 1);
 	}

	void dists_exps_3_test(test_tool& t) {
		generic_dists_exps_test(t, 0.50, 0);
 	}
	void dists_exps_3_h_test(test_tool& t) {
		generic_dists_exps_test(t, 0.50, 1);
 	}
	void dists_exps_4_test(test_tool& t) {
		generic_dists_exps_test(t, 1., 0);
 	}
	void dists_exps_4_h_test(test_tool& t) {
		generic_dists_exps_test(t, 1., 1);
 	}

	void generic_dists_ldeps_exps_test(test_tool& t, int vars, int mintrainsamples, int maxtrainsamples, int scaling_groups = -1) {
		int maxthreshold = 1024;
		int minthreshold = 64;
		int thresholdstep = 4;
		int testsamples = 400;
		int repeats = 20;

		for (double l = 0; l < 1.0; l+= 0.2) {
			for (double i = maxthreshold; i >= minthreshold; i/=thresholdstep) {
				srand(0);
				for (int r = 0; r < repeats; r++) {
					dep_problem test(vars, 0.2, l, 0);
					for (int j = mintrainsamples; j <= maxtrainsamples; j*=2) {
						int exps = -1;
						if (i > maxthreshold-1) exps = 0;
						std::string thtag("threshold:");
						if (exps<0) thtag += sup()<<i;
						else thtag += "naive";

						eval_genpredwith(t,
										 {sup()<<"ldeps:"<<l,
										  thtag,
										  sup()<<"n:"<<j},
										  test, j, testsamples, i, true, exps, 2., scaling_groups);
					}
				}
			}
			t<<"\n\nfor "<<l<<" ldeps:\n\n";
			t<<"this test tests learning with:\n";
			t<<"  1. varying number of threshold\n";
			t<<"  2. varying number of train data\n\n";

			std::set<std::string> tags = {"sample:test", "run:out", sup()<<"ldeps:"<<l};
			t<<"ldeps entropy by threshold and train data:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:"));

			t<<"\nsame plotted:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:")).to_plot(3, 30);

			t<<"\nexpression counts by threshold and train data:\n";
			t<<t.report(to_table<average>({"run:out", sup()<<"ldeps:"<<l, "prop:exps"}, "threshold:", "n:"));
		}

		std::set<std::string> tags = {"sample:test", "run:out"};

		table datainfo(t.report(to_table<average>({"run:out", sup()<<"n:"<<maxtrainsamples, "prop:entropy orig"}, "threshold:", "ldeps:")));
 		t<<"\ntrain data naive entropy before re-expression:\n";
		t<<datainfo;

		table datainfo2(t.report(to_table<average>({"run:out", sup()<<"n:"<<maxtrainsamples, "prop:entropy reexp"}, "threshold:", "ldeps:")));
 		t<<"\ntrain data naive entropy after re-expression:\n";
		t<<datainfo2;

		table naive(t.report(to_table<average>(tags+"prop:naivepredinfo"+(sup()<<"n:"<<maxtrainsamples), "threshold:", "ldeps:")));
		t<<"\nentropy based on average p, input ignored:\n";
		t<<naive;

		table ideal(t.report(to_table<average>(tags+"prop:idealinfo", "threshold:", "ldeps:")));
		t<<"\naverage ideal variable entropy by threshold and train data:\n";
		t<<ideal;

		t<<"\n\npredicted variable entropy by threshold and ldeps for "<<maxtrainsamples<<" train samples:\n";
		table measured( t.report(to_table<average>(tags+"prop:entryinfo"+(sup()<<"n:"<<maxtrainsamples), "threshold:", "ldeps:")) );
		t<<measured;
		t<<"\nsame plotted:\n";
		t<<measured.to_plot(4, 30);
		table delta( measured - ideal );
		t<<"\n(measured - ideal) entropy by threshold and ldeps plotted:\n";
		t<<delta.to_plot(4, 30);

		{
			table tldep02(
				t.report(to_table<average>({"sample:test", "run:out", "ldeps:0.2","prop:entryinfo"}, "threshold:", "n:")));

			std::ofstream o( t.file_path("ldeps_0_2_entropy.tex") );
			o<<tldep02.to_latex_pgf_doc("crossentropy");
		}

		{
			table tldep08(
				t.report(to_table<average>({"sample:test", "run:out", "ldeps:0.8","prop:entryinfo"}, "threshold:", "n:")));

			std::ofstream o( t.file_path("ldeps_0_8_entropy.tex") );
			o<<tldep08.to_latex_pgf_doc("crossentropy");
		}
	}

	void dists_ldeps_exps_test(test_tool& t) {
		generic_dists_ldeps_exps_test(t, 64, 1, 4096);
	}

	void dists_ldeps_exps_scaled_test(test_tool& t) {
		generic_dists_ldeps_exps_test(t, 64, 1, 4096, 128);
	}

	void dists_8var_ldeps_exps_test(test_tool& t) {
		generic_dists_ldeps_exps_test(t, 8, 1, 4096);
	}

	void dists_8var_ldeps_exps_scaled_test(test_tool& t) {
		generic_dists_ldeps_exps_test(t, 8, 1, 4096, 128);
	}

	void dists_64var_fast_ldeps_exps_test(test_tool& t) {
		generic_dists_ldeps_exps_test(t, 64, 4096, 4096);
	}


	void classes_setup_test(test_tool& t) {
		srand(0);

		for (int v = 1; v < 5; ++v) {
			class_problem pr(2, v);

			int n = 10000;

			typedef vars_problem p;

			reexp::lang<p> lang;
			reexp::data<p> data(lang, reexp::cvec<p>(n));
			setup_reg(lang, data, pr, n);

			reexp::stats<p> stats(data);

			reexp::pinfo names;
			setup_names(names, pr);

			t<<"test:\n"<<pr.tostring()<<"\n";
			reexp::stats_info<p> si(names, stats);
			t<<"stats:\n";
			for (int v = 0; v < lang.var_count(); ++v) {
				t<<"p("<<si.lang_info().var_tostring(v)<<")="<<stats.var(v).p()<<"\n";
			}
			t<<"\n";
			t<<"data info: "<<ideal_info(pr, data, pr.var_count()-1)<<"\n";
			t<<"\n";
			t<<"\nexample samples:\n";

			for (int i = 0; i < 8; ++i) {
				reexp::cond_bits b(pr.states());
				t<<"s:"<<reexp::vector_todensestring(b.states());
				t<<", d:"<<reexp::vector_todensestring(b.defined())<<"\n";
			}
		}
	}

	void classes_simple_test(test_tool& t) {
		int invars = 8;
		int maxtrainsamples = 2048;
		int testsamples = 50;
		int repeats = 30;

		srand(0);
		for (int i = 1; i <= invars; i*=2) {
			for (int r = 0; r < repeats; r++) {
				class_problem test(2, invars);
				for (int j = 1; j <= maxtrainsamples ; j*=2) {
					eval_genpredwith(t, {sup()<<"vars:"<<i, sup()<<"n:"<<j}, test, j, testsamples, 10000, true, -1, 2.);
				}
			}
		}

		t<<"this test tests learning with:\n";
		t<<"  1. varying number of independent variables\n";
		t<<"  2. varying number of train data\n\n";

		std::set<std::string> tags = {"sample:test", "run:out"};
		t<<"predicted variable entropy by input variable count and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "vars:", "n:"));

		t<<"\nsame plotted:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "vars:", "n:")).to_plot(3, 20);
	}

	void classes_exps_test(test_tool& t) {
		int maxvars = 16;
		int maxthreshold = 1024;
		int minthreshold = 16;
		int thresholdstep = 4;
		int maxtrainsamples = 2048;
		int testsamples = 200;
		int repeats = 50;

		for (int v = 4; v <= maxvars; v*=2) {
			for (double i = maxthreshold; i >= minthreshold; i/=thresholdstep) {
				srand(0);
				for (int r = 0; r < repeats; r++) {
					class_problem test(2, v);
					int exps = -1;
					if (i > maxthreshold-1) exps = 0;
					std::string thtag("threshold:");
					if (exps<0) thtag += sup()<<i;
					else thtag += "naive";
					for (int j = 1; j <= maxtrainsamples; j*=2) {
						eval_genpredwith(t,
										 {sup()<<"vars:"<<v,
										  thtag,
										  sup()<<"n:"<<j},
										  test, j, testsamples, i, true, exps);
					}
				}
			}
			t<<"\n\nfor "<<v<<" vars:\n\n";
			t<<"this test tests learning with:\n";
			t<<"  1. varying number of threshold\n";
			t<<"  2. varying number of train data\n\n";

			std::set<std::string> tags = {"sample:test", "run:out", sup()<<"vars:"<<v};
			t<<"predicted variable entropy by threshold and train data:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:"));

			t<<"\nsame plotted:\n";
			t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:")).to_plot(3, 30);

			t<<"\nexpression counts by threshold and train data:\n";
			t<<t.report(to_table<average>({"run:out", sup()<<"vars:"<<v, "prop:exps"}, "threshold:", "n:"));
		}

		std::set<std::string> tags = {"sample:test", "run:out"};

		t<<"\naverage ideal variable entropy by threshold and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:idealinfo", "threshold:", "vars:"));

		t<<"\n\npredicted variable entropy by threshold and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "threshold:", "n:"));

		t<<"\nfor naive bayesian, by vars and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo"+"threshold:naive", "vars:", "n:")).to_plot(3, 20);
		t<<"\nfor threshold "<<minthreshold<<", by vars and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo"+(sup()<<"threshold:"<<minthreshold), "vars:", "n:")).to_plot(3, 20);
	}

	void classes_single_test(test_tool& t) {
		srand(0);
		class_problem pr(2, 4);

		t<<pr.tostring();

		t<<"output ideal entropy: "<<pr.output_entropy()<<"\n";

		for (int th = 10; th <= 1000; th*=10) {
			for (int j = 1; j <= 8*1024; j*=2) {
				srand(0);
				for (int i = 0; i < 20; ++i) {
					eval_genpredwith(t, {sup()<<"th:"<<th, sup()<<"n:"<<j}, pr, j, 100, th, true);
				}
			}
		}
		{
			typedef vars_problem p;
			// setup the teach sample
			reexp::lang<p> lang;
			int n = 1024;
			reexp::data<p> data(lang, reexp::cvec<p>(n));
			setup_reg(lang, data, pr, n);
			time_sentry time1;
			reexp::stats<p> stats(data);
			reexp::learner<p> learner(lang, stats, 100);
			for (int i = 0; i < lang.var_count(); ++i) {
				if (pr.is_output(i)) {
					learner.exclude(i);
				}
			}
			learner.reexpress(true, -1);

			reexp::pinfo names;
			setup_names(names, pr);
			reexp::stats_info<p> si(names, stats);

			t<<si.vars_tostring();

			t<<si.rels_tostring();
		}

		std::set<std::string> tags = {"run:out"};
		t<<"predicted variable entropy by threshold and train data:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "th:", "n:"));

		t<<"\nsame plotted:\n";
		t<<t.report(to_table<average>(tags+"prop:entryinfo", "th:", "n:")).to_plot(3, 30);

		t<<"\nexpressions:\n";
		t<<t.report(to_table<average>(tags+"prop:exps", "th:", "n:"));

	}

}

void addvarstest(TestRunner& runner) {
	runner.add("vars/setup", 			{"func"},  &setup_test);
	runner.add("vars/xor",   			{"func"},  &xor_test);
	runner.add("vars/and",  			{"func"},  &and_test);
	runner.add("vars/or",    			{"func"},  &or_test);
	runner.add("vars/or4",   			{"func"},  &or4_test);
	runner.add("vars/redundant",		{"func"},  &redundant_test);
	runner.add("vars/sparseredundant",  {"func"},  &sparseredundant_test);
	runner.add("vars/redundantperf",    {"func"},  &redundantperf_test);

	runner.add("vars/or4_M",   			{"perf"},  &or4_M_test);
	runner.add("vars/dists_10var_M",   	{"perf"},  &dists_10var_M_test);
	runner.add("vars/dists_20var_M",   	{"perf"},  &dists_20var_M_test);
	runner.add("vars/dists_50var_M",   	{"perf"},  &dists_50var_M_test);

	runner.add("vars/dists_setup",		{"func"},  &dists_setup_test);
	runner.add("vars/dists_indep_setup",{"func"},  &dists_indep_setup_test);
	runner.add("vars/dists_hidden_setup",{"func"}, &dists_hidden_setup_test);
	runner.add("vars/dists_simple",		{"func"},  &dists_simple_test);
	runner.add("vars/dists_prioriw",	{"func"},  &dists_prioriw_test);
	runner.add("vars/dists_exps",		{"func"},  &dists_exps_test);
	runner.add("vars/dists_exps_h",		{"func"},  &dists_exps_h_test);
	runner.add("vars/dists_exps_h2",	{"func"},  &dists_exps_h2_test);
	runner.add("vars/dists_exps_2",		{"func"},  &dists_exps_2_test);
	runner.add("vars/dists_exps_2_h",	{"func"},  &dists_exps_2_h_test);
	runner.add("vars/dists_exps_3",		{"func"},  &dists_exps_3_test);
	runner.add("vars/dists_exps_3_h",	{"func"},  &dists_exps_3_h_test);
	runner.add("vars/dists_exps_4",		{"func"},  &dists_exps_4_test);
	runner.add("vars/dists_exps_4_h",	{"func"},  &dists_exps_4_h_test);
	runner.add("vars/dists_ldeps_exps",	{"func"},  &dists_ldeps_exps_test);
	runner.add("vars/dists_ldeps_exps_scaled",	{"func"},  &dists_ldeps_exps_scaled_test);
	runner.add("vars/dists_8var_ldeps_exps",	{"func"},  &dists_8var_ldeps_exps_test);
	runner.add("vars/dists_8var_ldeps_exps_scaled",	{"func"},  &dists_8var_ldeps_exps_scaled_test);
	runner.add("vars/dists_64var_fast_ldeps_exps",	{"func"},  &dists_64var_fast_ldeps_exps_test);

	runner.add("vars/classes_setup",	{"func"},  &classes_setup_test);
	runner.add("vars/classes_simple",	{"func"},  &classes_simple_test);
	runner.add("vars/classes_exps",	    {"func"},  &classes_exps_test);
	runner.add("vars/classes_single",	{"func"},  &classes_single_test);

}

