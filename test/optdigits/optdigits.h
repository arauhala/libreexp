/*
 * optdigits.h
 *
 *  Created on: Sep 13, 2012
 *      Author: arau
 */

#ifndef OPTDIGITS_H_
#define OPTDIGITS_H_

namespace optdigits {
	static const int OriginalWidth = 32;

	static const int OriginalHeight = 32;

	extern const char* TRA_DATA_FILE;
	extern int TRA_DATA_FILE_SAMPLES;
	extern const char* CV_DATA_FILE;
	extern int CV_DATA_FILE_SAMPLES;
	extern const char* WDEP_DATA_FILE;
	extern int WDEP_DATA_FILE_SAMPLES;
	// use this one testing
	extern const char* WINDEP_DATA_FILE;
	extern int WINDEP_DATA_FILE_SAMPLES;

	static int DATA_FILE_HEADER_LINES = 21;

	namespace varid {
		static const int digit0 = 0;
		static const int digit1 = 1;
		static const int digit2 = 2;
		static const int digit3 = 3;
		static const int digit4 = 4;
		static const int digit5 = 5;
		static const int digit6 = 6;
		static const int digit7 = 7;
		static const int digit8 = 8;
		static const int digit9 = 9;
	}

	template <typename P>
	void evaluate(TestTool& t,
				  const explib::data<P>& data,
				  const explib::stats<P>& stats,
				  const std::set<std::string>& tags,
				  int sample_cvarid) {
		explib::pred<P> pred(stats);
		std::vector<double> ps[10];

		TimeSentry timer;
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
		const explib::data_var<P>* vars[10];

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
		explib::cvec<P> at;
		int ok = 0;
		double info = 0;
		for (int s = 0; s < samples; s++) {
			int digit = -1;
			// find out the right answer
			at[sample_cvarid] = s;
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

	template <typename P>
	void do_evaluate(TestTool& t,
					 explib::lang<P>& lang,
					 explib::data<P>& data,
					 explib::data<P>& tdata,
					 explib::stats<P>& stats,
					 explib::learner<P>& learner,
					 int sample_cvarid) {
		tdata.apply_exps();
		std::set<std::string> tags = {sup()<<"exps: "<<lang.exp_count()};
		t.record(tags+"gen:info", stats.naiveInfo());
		evaluate(t, data, stats, tags+"data:train", sample_cvarid);
		evaluate(t, tdata, stats, tags+"data:test", sample_cvarid);
	}
}

#endif /* OPTDIGITS_H_ */
