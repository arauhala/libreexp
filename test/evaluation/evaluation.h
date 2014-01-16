/*
 * evaluation.h
 *
 *  Created on: May 8, 2012
 *      Author: arau
 */

#ifndef EVALUATION_H_
#define EVALUATION_H_
#include "tester.h"
#include "reexp/io.h"

namespace evaluation {

	enum pred_problem_type {
		var_pred_problem,
		classification_problem
	};

	template <typename P>
	class pred_problem {
		public:
			reexp::lang<P> lang_;
			reexp::data<P> data_;
			reexp::bits predvars_;
			reexp::pinfo names_;
			int costs_[4];
			int samplecvar_;
			pred_problem_type type_;
			inline pred_problem(const reexp::cvec<P>& dim = reexp::cvec<P>(), pred_problem_type type = var_pred_problem)
			: lang_(), data_(lang_, dim), predvars_(),
			  names_(), costs_(), samplecvar_(), type_(type) {}
	};

	struct pred_args {
		double expthreshold_;
		double exprelfilter_;
		double predrelfilter_;
		double prioriweight_;
		int verbose_;
		bool useLogDepB_;
		int scaling_group_sz_;
		int dist_record_groups_;
		bool crosscompress_;
		std::set<std::string> tags_;
		pred_args(double expthreshold,
				  double exprelfilter,
				  double predrelfilter,
				  double prioriweight = 2.,
				  int verbose = 0,
				  bool useLogDepB = false,
				  int scaling_group_sz = -1,
				  bool crosscompress = false)
		:   expthreshold_(expthreshold),
		    exprelfilter_(exprelfilter),
		    predrelfilter_(predrelfilter),
		    prioriweight_(prioriweight),
		    verbose_(verbose),
		    useLogDepB_(useLogDepB),
		    scaling_group_sz_(scaling_group_sz),
		    dist_record_groups_(0),
		    crosscompress_(crosscompress),
		    tags_() {
			tags_.insert(sup()<<"threshold:"<<expthreshold_);
			tags_.insert(sup()<<"relfilter:"<<exprelfilter_);
			tags_.insert(sup()<<"predfilter:"<<predrelfilter_);
			tags_.insert(sup()<<"prioriw:"<<prioriweight_);

		}
		const std::set<std::string>& tags() const {
			return tags_;
		}
	};
	enum {
		true_positive = 0,
		true_negative,
		false_positive,
		false_negative
	};



	struct pred_data_stats {
		double e_; // entropy
		double naiveE_;
		double us_;
		int errors_;
		int naiveErrors_;
		int costs_;
		int naiveCosts_;
		int n_;
		double ndl_;
		double xndl_; // maybe infinite
		double endl_;
		double exndl_;
		int encodedN_;
		double encodedB_;
		long encodeus_;
		long decodeus_;
		inline pred_data_stats()
		: e_(), naiveE_(), us_(), errors_(), naiveErrors_(),
		  costs_(), naiveCosts_(), n_(), ndl_(), xndl_(), endl_(), exndl_(), encodedN_(), encodedB_(), encodeus_(), decodeus_() {}
		inline double unitEntropy() const {
			return e_ / n_;
		}
		inline double unitNaiveEntropy() const {
			return naiveE_ / n_;
		}
		inline double errorRate() const {
			return errors_ / double(n_);
		}
		inline double naiveErrorRate() const {
			return naiveErrors_ / double(n_);
		}
		inline double unitCost() const {
			return costs_ / double(n_);
		}
		inline double unitNaiveCost() const {
			return naiveCosts_ / double(n_);
		}
		inline double unitNdl() const {
			return ndl_ / double(encodedN_);
		}
		// maybe infinite, because the underlying estimate is plain average
		inline double unitXndl() const {
			return xndl_ / double(encodedN_);
		}
		inline double unitEndl() const {
			return endl_ / double(encodedN_);
		}
		inline double unitExndl() const {
			return exndl_ / double(encodedN_);
		}
		inline double unitEncodedSizeB() const {
			return encodedB_ / double(encodedN_);
		}
		inline double unitEncodeUs() const {
			return encodeus_/ double(encodedN_);
		}
		inline double unitDecodeUs() const {
			return decodeus_/ double(encodedN_);
		}

		inline void record(test_tool& t, const std::set<std::string>& tags) {
			t.record(tags+"prop:entropy",    unitEntropy());
			t.record(tags+"prop:naive e", 	 unitNaiveEntropy());
			t.record(tags+"prop:err",		 errorRate());
			t.record(tags+"prop:naive err"	,naiveErrorRate());
			t.record(tags+"prop:cost", 		 unitCost());
			t.record(tags+"prop:naiveCost",  unitNaiveCost());
			t.record(tags+"perf:ns",	     (1000*us_)/(n_));

			if (encodedN_) {
				t.record(tags+"prop:ndl",  		unitNdl());
				t.record(tags+"prop:xndl",  	unitXndl());
				t.record(tags+"prop:endl",  	unitEndl());
				t.record(tags+"prop:exndl",  	unitExndl());
				t.record(tags+"prop:encodedB",  unitEncodedSizeB());
				t.record(tags+"perf:encodeUs",  unitEncodeUs());
				t.record(tags+"perf:decodeUs",  unitDecodeUs());
			}
	//			t.record({"perf:round us", thstr, filstr}, 	 (stats.us_)/(stats.rounds_));

		}

	};

	struct pred_stats {
		int activerels_;
		int allrels_;
		int exps_;
		double ndl_;
		long reexpus_;
		long applyus_;
		int rounds_;
		pred_data_stats train_;
		pred_data_stats test_;
		inline pred_stats()
		: activerels_(), allrels_(), exps_(), ndl_(), reexpus_(), applyus_(), rounds_(), train_(), test_() {}
		inline double averReexpUs() const {
			return reexpus_ / double(rounds_);
		}
		inline double averApplyUs() const {
			return applyus_ / double(rounds_);
		}
		inline double averExps() const {
			return exps_ / double(rounds_);
		}
		inline double averndl() const {
			return ndl_ / double(rounds_);
		}
		inline double averActiveRels() const {
			return activerels_ / double(rounds_);
		}
		inline double averAllRels() const {
			return allrels_ / double(rounds_);
		}
	};

	template <typename P>
	void measure_var_pred(std::ostream& buf,
						  pred_problem<P>& pr,
						  pred_args& args,
						  pred_data_stats& s,
						  const reexp::stats<P>& trainstats,
						  const reexp::data<P>& predicted,
						  const std::vector<std::vector<double> >& allps) {
		for (int v = 0; v < pr.predvars_.size(); v++) {
			if (pr.predvars_[v]) {
				const std::vector<double>& ps = allps[v];
				const reexp::data_var<P>& datavar = predicted.var(v);
				const reexp::var_stats<P>& vs = trainstats.var(v);
				reexp::cvec<P> at;
				for (size_t i = 0; i < ps.size(); i++) {
					at[pr.samplecvar_] = i;
					bool value = *datavar[at];
					double p = ps[i];
					double naiveP = vs.eP();

					// calculate entropy
					double e = (value?-log2(p):-log2(1-p));
					double naiveE = (value?-log2(naiveP):-log2(1-naiveP));
					s.e_ += e;
					s.naiveE_ += naiveE;

					bool decision = p >= 0.5;
					bool naiveDecision = naiveP >= 0.5;
					bool error = value != decision;
					bool naiveError = value != naiveDecision;
					s.errors_ += error?1:0;
					s.naiveErrors_ += naiveError?1:0;

					bool costDecision;
					{
						double exptruecost = 	p*pr.costs_[true_positive] + (1-p) * pr.costs_[false_positive];
						double expfalsecost = 	(1-p)*pr.costs_[true_negative] + p * pr.costs_[false_negative];
						costDecision = exptruecost < expfalsecost;

						if (costDecision != value) {
							s.costs_ += pr.costs_[costDecision?false_positive:false_negative];
						} else {
							s.costs_ += pr.costs_[costDecision?true_positive:true_negative];
						}
					}
					bool costNaiveDecision;
					{
						double exptruecost = naiveP * true_positive + (1-naiveP) * pr.costs_[false_positive];
						double expfalsecost = (1-naiveP)* true_negative + naiveP * pr.costs_[false_negative];
						costNaiveDecision = exptruecost < expfalsecost;

						if (costNaiveDecision != value) {
							s.naiveCosts_ += pr.costs_[costNaiveDecision?false_positive:false_negative];
						} else {
							s.naiveCosts_ += pr.costs_[costNaiveDecision?true_positive:true_negative];
						}
					}
					s.n_++;
					buf.width(14);
					buf<<pr.names_.vnames_[v];
					buf<<p<<"/"<<naiveP;
					buf<<"   ["<<(value?"true":"false")<<"]\n";
				}
			}
		}
	}

	template <typename P>
	void measure_classification(std::ostream& buf,
								pred_problem<P>& pr,
								pred_args& args,
								pred_data_stats& s,
								const reexp::stats<P>& trainstats,
								const reexp::data<P>& predicted,
  							    const std::vector<std::vector<double> >& allps) {
		size_t samples = predicted.dim()[0];
		for (size_t i = 0; i < samples; i++) {
			reexp::cvec<P> at(i);
			double sum = 0, naiveSum = 0, selectedP = -1, naiveSelectedP= -1;
			int selected = -1, naiveSelected = -1;
			int correct = -1;
			for (int v = 0; v < pr.predvars_.size(); v++) {
				const reexp::data_var<P>& datavar = predicted.var(v);
				const reexp::var_stats<P>& vs = trainstats.var(v);
				if (pr.predvars_[v]) {
					double p = allps[v][i];
					double naiveP = vs.eP();
					if (p > selectedP) {
						selected = v;
						selectedP = p;
					}
					if (naiveP > naiveSelectedP) {
						naiveSelected = v;
						naiveSelectedP = naiveP;
					}
					sum += p;
					naiveSum += naiveP;
					if (*(datavar[at])) {
						correct = v;
					}
				}
			}
			selectedP /= sum;
			naiveSelectedP /= naiveSum;
			// calculate entropy
			double e = -log2(selectedP);
			double naiveE = -log2(naiveSelectedP);
			s.e_ += e;
			s.naiveE_ += naiveE;

			bool error = (selected != correct);
			bool naiveError = (naiveSelected != correct);
			s.errors_ += error?1:0;
			s.naiveErrors_ += naiveError?1:0;

			s.n_++;
			buf.width(14);
			buf<<"selected:"<<pr.names_.vnames_[selected]<<", p="<<selectedP<<". correct:"<<pr.names_.vnames_[correct]<<"  "<<(error?"[error]":"[ok]")<<"\n";
		}
	}

	template <typename P>
	void record_distribution(test_tool& t,
							 const std::set<std::string>& runtags,
							 const reexp::stats<P>& trainstats,
							 const reexp::data<P>& predicted,
  							 const std::vector<std::vector<double> >& allps,
  							 size_t groups) {
		for (size_t i = 0; i < allps.size(); ++i) {
			if (allps[i].size()) {
				const std::vector<double>& ps = allps[i];

				const reexp::data_var<P>& dv = predicted.var(i);

				std::vector<std::pair<double, size_t> > pi;
				for (size_t i = 0; i < ps.size(); ++i) {
					pi.push_back(std::pair<double, size_t>{ps[i], i});
				}
				std::sort(pi.begin(), pi.end());
				std::reverse(pi.begin(), pi.end());
;
				size_t at = 0;
				std::vector<double> e_p_avers;
				std::vector<double> f_avers;
				std::vector<size_t> szs;
				for (size_t i = 0; i < groups; ++i) {
					size_t group_sz = (pi.size() - at) / (groups - i);
					double f_sum = 0;
					double e_p_sum = 0;
		//					double i_sum = 0;
					for (size_t j = 0; j < group_sz; ++j) {
						double p = pi[at].first;
						bool state = dv.states()[pi[at].second];
						e_p_sum += p;
						f_sum += state?1:0;
						at++;
					}
					e_p_avers.push_back(e_p_sum / group_sz);
					f_avers.push_back(f_sum / double(group_sz));
					szs.push_back(group_sz);
				}
				// recording structure
				//     input tags
				//     index
				//     real probability
				//
				size_t group = 0;
				size_t of_group = 0;
				for (size_t i = 0; i < pi.size(); ++i) {
					t.record(runtags + (sup()<<"rank:"<<i) + "prop:e_p", pi[i].first);
					//t.record(runtags + (sup()<<"rank:"<<i) + "prop:state", dv.states()[pi[i].second]);

					if (of_group >= szs[group]) {
						group++; of_group = 0;
					}
					t.record(runtags + (sup()<<"rank:"<<i) + "prop:aver_e_p", 	 e_p_avers[group]);
					t.record(runtags + (sup()<<"rank:"<<i) + "prop:aver_f",      f_avers[group]);

					of_group++;
					//t.record(runtags + (sup()<<"rank:"<<i) + "prop:aver_delta",  e_p_avers[group]-p_avers[group]);
				}
			}
		}
	}


	template <typename P>
	void crosscompress_and_measure(test_tool& t,
								   pred_problem<P>& pr,
								   pred_args& args,
								   pred_data_stats& s,
								   const reexp::stats<P>& trainstats,
								   reexp::data<P>& compressed) {
		reexp::stats<P> compressed_stats(compressed);

		reexp::bits buffer(64*1024*1024); // big enough buffer
		reexp::bit_ostream bout = buffer.ostream(0);
		reexp::arithmetic_bit_ostream out(bout);

		reexp::io<P> io(trainstats);
		long wus;
		{
			time_sentry t;
			io.write_states(out, compressed);
			out.finish();
			wus = t.us();
		}
		s.encodedB_ += bout.pos();
		s.encodeus_ += wus;

		reexp::bit_istream bin = buffer.istream(0);
		reexp::arithmetic_bit_istream in(bin);
		reexp::data<P> data2(compressed.lang(), compressed.dim());
		data2.prepare_exp_vars();
		long rus;
		reexp::index_over_var_bits<P> indexspace(data2);
		std::vector<int> knownAt;
		{
			time_sentry t;
			io.read_states(in, data2, indexspace, knownAt);
			rus = t.us();
		}
		if (bin.pos() != bout.pos()) {
			reexp::lang_info<P> li(pr.names_, compressed.lang());

			t<<li.vars_tostring()<<"\n";
			t<<li.invorder_diff_tostring(compressed, data2, 1, &indexspace, &knownAt);

			throw std::runtime_error(sup()<<"read and wrote different number of bits, "<<bin.pos()<<" != "<<bout.pos());
		}

		s.decodeus_ += rus;
		s.encodedN_ += compressed.dim().volume();
		s.ndl_ 		+= compressed_stats.ndl();
		s.xndl_ 	+= trainstats.xndl(compressed);
		s.endl_ 	+= compressed_stats.endl();
		s.exndl_ 	+= trainstats.exndl(compressed);

//		t<<bin.pos()<<" bits read from encoded blob in "<<rms<<" ms\n";
	}

	template <typename P>
	void predict_and_measure(test_tool& t,
							 pred_problem<P>& pr,
							 pred_args& args,
							 pred_data_stats& s,
							 const reexp::stats<P>& trainstats,
							 const reexp::data<P>& predicted,
							 bool test = false) {
		reexp::pred<P> pred(trainstats, args.prioriweight_, args.useLogDepB_);

		std::ostringstream buf;
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);
		buf.setf(std::ios::left, std::ios::adjustfield);
		std::vector<std::vector<double> > ps;
		for (int v = 0; v < pr.predvars_.size(); v++) {
			ps.push_back(std::vector<double>());
			if (pr.predvars_[v]) {
				reexp::group_scaler<P> scaler(pred, v, args.scaling_group_sz_);
				time_sentry time;
				ps.back() = pred.bitP(predicted, v);
				scaler.scale(ps.back());
				s.us_ += time.us();
			}
		}

		switch (pr.type_) { // should be done with inheritance
		case var_pred_problem:
			measure_var_pred(buf, pr, args, s, trainstats, predicted, ps);
			break;
		case classification_problem:
			measure_classification(buf, pr, args, s, trainstats, predicted, ps);
			break;
		}

		if (args.verbose_)
			t<<buf.str()<<"\n";

		if (test && args.dist_record_groups_ > 0) {
			record_distribution<P>(t, args.tags(), trainstats, predicted, ps, args.dist_record_groups_);
		}
	}

	template <typename P>
	void teach_test_measure(test_tool& t,
							pred_problem<P>& pr,
							reexp::bits& testsamples,
							pred_args& args,
							pred_stats& s) {
		int samplecvar = pr.samplecvar_; // assumption
		int sampleN = pr.data_.dim()[samplecvar];
		int testN = testsamples.popcount();
		int trainN = sampleN - testN;
		reexp::lang<P> lang(pr.lang_);
		reexp::cvec<P> traindim = pr.data_.dim();
		traindim[samplecvar] = trainN;
		reexp::data<P> train(lang, traindim);
		reexp::cvec<P> testdim = pr.data_.dim();
		testdim[samplecvar] = testN;
		reexp::data<P> test(lang, testdim);
		lang.set_obs(train); // notify solely the train data of new exps

		for (int v = 0; v < pr.lang_.var_count(); ++v) {
			const reexp::data_var<P>& datavar = pr.data_.var(v);
			reexp::cvec<P> dataAt;
			reexp::data_var<P>& trainvar = train.var(v);
			reexp::cvec<P> trainAt;
			reexp::data_var<P>& testvar = test.var(v);
			reexp::cvec<P> testAt;
			for (int s = 0; s < sampleN; ++s) {
				auto databit = datavar[dataAt];
				if (testsamples[s]) {
					testvar[testAt] = bool(databit);
					*testvar[testAt] = bool(*databit);
					testAt[samplecvar]++;
				} else {
					trainvar[trainAt] = bool(databit);
					*trainvar[trainAt] = bool(*databit);
					trainAt[samplecvar]++;
				}
				dataAt[samplecvar]++;
			}
		}
		reexp::stats<P> trainstats(train);

		reexp::learner<P> learner(lang, trainstats, args.expthreshold_, args.exprelfilter_, args.predrelfilter_);
		learner.exclude(pr.predvars_);
		int exps;
		{
			time_sentry time;
			exps = learner.reexpress(true);
			s.reexpus_ += time.us();
		}
		learner.unexclude(pr.predvars_);

		// also reexpress the test data with learned expressions
		{
			time_sentry time;
			test.apply_exps();
			s.applyus_ += time.us();
		}

		// count active relations
		int relsactive = 0;
		for (int i = 0; i < lang.rel_count(); ++i) {
			if (!lang.rel(i).disabled()) relsactive++;
		}

		s.ndl_ += trainstats.ndl();
		s.exps_ += exps;
		s.activerels_ += relsactive;
		s.allrels_ += lang.rel_count();
		s.rounds_++;

		if (args.verbose_) {
			t<<"test run, "<<relsactive<<"/"<<lang.rel_count()<<" rels active, "<<exps<<" exps\n";
		}

		// do crosscompression before predicted values are undefined. compression method cannot know,
		// which values were undefined by hand, which breaks the encoding
		if (args.crosscompress_) {
			crosscompress_and_measure(t, pr, args, s.train_, trainstats, train);
			crosscompress_and_measure(t, pr, args, s.test_, trainstats, test);
		}

		// prevent predictor from using predicted variables
		// for predicting other predicted variables
		for (int v = 0; v < pr.predvars_.size(); v++) {
			if (pr.predvars_[v]) {
				train.var(v).defined().fill(false);
				test.var(v).defined().fill(false);
			}
		}
		predict_and_measure(t, pr, args, s.train_, trainstats, train);
		predict_and_measure(t, pr, args, s.test_,  trainstats, test, true);

	}

	inline void report_measurements(test_tool& t, pred_args& args, pred_stats& stats) {
		const std::set<std::string>& tags = args.tags();

		stats.test_.record(t, tags+"data:test");
		stats.train_.record(t, tags+"data:train");


		t.record(tags+"exp:reexp us",   stats.averReexpUs());
		t.record(tags+"exp:apply us",   stats.averApplyUs());
		t.record(tags+"exp:exps", 		stats.averExps());
		t.record(tags+"exp:naive info", stats.averndl());
		t.record(tags+"exp:act rels", 	stats.averActiveRels());
		t.record(tags+"exp:all rels",   stats.averAllRels());
		if (args.verbose_) {
			t<<"exp treshold: "<<args.expthreshold_<<"\n";
			t<<"exp rel filter: "<<args.exprelfilter_<<"\n\n";
			t<<"pred rel filter: "<<args.predrelfilter_<<"\n\n";
			t<<"priori weight: "<<args.prioriweight_<<"\n\n";
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

	template <typename P>
	void crossvalidate_run(test_tool& t, pred_problem<P>& pr, pred_args& args, int split) {
		pred_stats stats;
		int samples = pr.data_.dim()[pr.samplecvar_];
		int chunkSize = samples / split;
		reexp::bits testsamples;
		testsamples.resize(samples);
		for (int i = 0; i < split; ++i) {
			testsamples.fill(false);
			for (int j = 0; j < chunkSize; ++j) {
				testsamples[i*chunkSize + j] = true;
			}
			teach_test_measure(t, pr, testsamples, args, stats);
		}
		report_measurements(t, args, stats);
	}

	template <typename P>
	void random_crossvalidate_run(test_tool& t, pred_problem<P>& pr, pred_args& args) {
		pred_stats stats;
		int samples = pr.data_.dim()[pr.samplecvar_];
		int left = samples/2;
		reexp::bits testsamples;
		testsamples.resize(samples);
		for (int i = 0; i < samples; ++i) {
			double p = double(left)/(samples-i);
			if (rand() < p*RAND_MAX) {
				testsamples[i] = true;
				left--;
			}
		}
		teach_test_measure(t, pr, testsamples, args, stats);
		for (int i = 0; i < samples; ++i) { // invert
			testsamples[i] = !testsamples[i];
		}
		teach_test_measure(t, pr, testsamples, args, stats);

		report_measurements(t, args, stats);
	}

	template <typename P>
	void separate_train_test_datas_run(test_tool& t, pred_problem<P>& pr, pred_args& args, int train) {
		pred_stats stats;
		int samples = pr.data_.dim()[pr.samplecvar_];
		reexp::bits testsamples;
		testsamples.resize(samples);
		testsamples.fill(false);
		for (int i = train; i < samples ; ++i) {
			testsamples[i] = true;
		}
		teach_test_measure(t, pr, testsamples, args, stats);

		report_measurements(t, args, stats);
	}

}

#endif /* EVALUATION_H_ */
