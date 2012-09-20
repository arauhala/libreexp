/*
 * pred.h
 *
 *  Created on: Aug 16, 2011
 *      Author: arau
 */

#ifndef PRED_H_
#define PRED_H_


#include "stats.h"

namespace explib {

	inline int pow3i(int v) {
		int rv = 1;
		while (v-- > 0) rv *= 3;
		return rv;
	}
	template <typename P>
	struct rel_inputvars {
		const explib::rel<P>& r_;
		const int& predicted_;
		rel_inputvars(const explib::rel<P>& r,
				      const int& predicted)
		: 	r_(r), predicted_(predicted) {}
		inline int stateCount() const {
			return pow3i(r_.varCount()-1);
		}
		enum DefState {
			FALSE = 0,
			TRUE = 1,
			UNDEFINED = 2
		};
		inline int maskState(int dmask, int smask) const {
			int rv = 0; // here we return the riv state
			int exp = 1, mask = 1; // exp is used for writing riv, mask is used for reading mask
			for (int i = 0; i < r_.varCount(); ++i) {
				if (i != predicted_) {
					if (!(dmask&mask)) { // if undefined
						rv += UNDEFINED*exp;
					} else if (smask&mask) { // if true
						rv += TRUE*exp;
					}
					exp*=3;
				}
				mask<<=1;
			}
			return rv;
		}
		inline int varDefState(int var, int state) const {
			if (var > predicted_) var--;
			return (state / pow3i(var)) % 3;
		}
		inline bool includesRelState(int istate, int rstate) const {
			int i = predicted_;
			while (istate || rstate || i > 0) {
				if (!i--) {
					rstate /= 2;
				} else {
					int is = istate % 3;
					int rs = rstate % 2;
					if (is != UNDEFINED && is != rs) {
						return false;
					}
					istate /= 3;
					rstate /= 2;
				}
			}
			return true;
		}
	};

	struct inputstate_logdep {
		inline inputstate_logdep(int max_rels)
		: d_(2*size_t(pow3i(max_rels-1))) {}
		inline double& depStateV(int state) {
			return d_[state*2];
		}
		inline double& depStateNotV(int state) {
			return d_[(state*2)+1];
		}
		std::vector<double> d_;
	};

	struct fake_output {
		template <typename T>
		inline fake_output& operator<<(T) {
			return *this;
		}
	};

	template <typename P>
	class pred {
		private:
			const stats<P>& stats_;
			double prioriWeight_;
			bool useLogDepB_;

		public:
			pred(const stats<P>& s, double prioriWeight = 2., bool useLogDepB = false)
			: stats_(s), prioriWeight_(prioriWeight), useLogDepB_(useLogDepB) {}

#if 0

			std::vector<double> bitP(const data<P>& d, int varid) {
				const var_stats<P>& vs = stats_.var(varid);
				const data_var<P>& dv = d.var(varid);
				const var<P>& v = dv.var_;
				const cvec<P>& dim = dv.dim();

				int volume = dim.volume();
				std::vector<double> rv(volume);
				std::vector<long> wTrue(volume);
				std::vector<long> wFalse(volume);

				for (int i = 0; i < volume; ++i) {
					wTrue[i] = vs.freq() + 1;
					wFalse[i] = vs.n() - vs.freq() + 1;
				}

				// use relations for prediction work
				const std::vector<rel<P>*>& rels = v.rels();
				for (auto i = rels.begin(); i != rels.end(); ++i) {
					const rel<P>& r = **i;
					if (!r.disabled()) {
						const rel_stats<P>& rs = stats_.rel(r.id());
						const data_rel<P>& dr = d.rel(r.id());
						auto re = r.entries();


						for (dim_iterator<P> j(dr.dim()); j; ++j) {
							const cvec<P>& at = *j;

							int vidx = 0; // this is the predicted one
							int maskd = 0; // which vars are defined
							int masks = 0; // which defined vars are true
							for (int j = 0; j < r.varCount(); ++j) {
								if (re[j].var_ == &v) {
									vidx = j;
								} else {
									auto b = dr.var(j)[at];
									if (b) {
										maskd |= 1<<j;
										masks |= (*b != 0)<<j;
									}
								}
							}
							// So we have a relation of for r(y, x1, x2, ...), where y
							// is interesting variable. In here, we are interesting of
							// p(x1&x2&...|y) and p(x1&x2&...|!y) values, which are
							// of form p(x1&x2&...|y) = p(x1&x2&...&y) / p(y)

							long freqTrue = 1;
							long freqFalse = 1;

							for (int s = 0; s < r.stateCount(); ++s) {
								if ((s&maskd) == masks) {
									// a valid interesting state, based on data
									int sf = rs.stateFreqs()[s];
									bool vs = r.varState(vidx, s);
									if (vs) {
										freqTrue += sf;
									} else {
										freqFalse += sf;
									}
								}
							}

							int sidx = dim.offset(at);
							long& wt = wTrue[sidx];
							long& wf = wFalse[sidx];
							wt *= freqTrue * (rs.n() - rs.varFreqs()[vidx] + 1);
							wf *= freqFalse * (1+rs.varFreqs()[vidx]);
							const long bigNumber = 0x10000000;
							while (wt > bigNumber || wf > bigNumber) {
								wt >>= 1;
								wf >>= 1;
							}
						}
					}
				}

				for (int i = 0; i < volume; ++i) {
					rv[i] = wTrue[i] / double(wTrue[i] + wFalse[i]);
				}

				return rv;
			}
#endif

			template <typename Out>
			inline void getInputStateLogDepA(const explib::rel_stats<P>& rs,
											 const rel_inputvars<P>& riv,
											 inputstate_logdep& inputStateLogDep,
											 Out& explain) {
				const explib::rel<P>& r = rs.data().rel();
				int vidx = riv.predicted_;
				double n = rs.n();
				double f = rs.varFreqs()[vidx];
				double varPrioriP = r.entries()[vidx].var_->prioriP();
//				double p = (rs.varFreqs()[vidx]+varPriori) / double(n + 2);
//				double p = (rs.varFreqs()[vidx]+1) / double(n + 2);
				// cache & prepopulate dependency values
				for (int i = 0; i < riv.stateCount(); ++i) {
/*					double predVarTrue = 2*p;
					double predVarFalse = 2*(1-p);*/
					double predVarTrue = 0;
					double stateFreq = 0;
					double trueStatePrioriP = 0;
					double falseStatePrioriP = 0;
					for (int s = 0; s < r.stateCount(); ++s) {
						if (riv.includesRelState(i, s)) {
							int sfreq = rs.stateFreqs()[s];
							bool predVarState = r.varState(vidx, s);
							if (predVarState) {
								predVarTrue += sfreq;
								trueStatePrioriP += r.statePrioriP(s);
							} else {
								falseStatePrioriP += r.statePrioriP(s);
							}
							stateFreq += sfreq;
						}
					}
					inputStateLogDep.depStateV(i) = log2((predVarTrue + trueStatePrioriP*prioriWeight_)/(f+varPrioriP*prioriWeight_));
					explain<<"depStateV("<<i<<")=log2(("<<predVarTrue<<"+"<<trueStatePrioriP*prioriWeight_<<")/("<<f<<"+"<<varPrioriP*prioriWeight_<<"))="<<inputStateLogDep.depStateV(i)<<"\n";
					inputStateLogDep.depStateNotV(i) = log2(((stateFreq-predVarTrue)+falseStatePrioriP*prioriWeight_)
														    /(n - f + prioriWeight_*(1-varPrioriP)));
					explain<<"depStateNotV("<<i<<")=log2(("<<(stateFreq-predVarTrue)<<"+"<<falseStatePrioriP*prioriWeight_
														  <<")/("<<(n - f)<<"+"<<prioriWeight_*(1-varPrioriP)<<"))="<<inputStateLogDep.depStateNotV(i)<<"\n";
				}
			}

			template <typename Out>
			inline void getInputStateLogDepB(const explib::rel_stats<P>& rs,
											 const rel_inputvars<P>& riv,
											 inputstate_logdep& inputStateLogDep,
											 Out& explain) {
				const explib::rel<P>& r = rs.data().rel();
				int vidx = riv.predicted_;
				double n = rs.n();
				double f = rs.varFreqs()[vidx];
				for (int i = 0; i < riv.stateCount(); ++i) {
					double predVarTrue = 0;
					double stateFreq = 0;
					for (int s = 0; s < r.stateCount(); ++s) {
						if (riv.includesRelState(i, s)) {
							int sfreq = rs.stateFreqs()[s];
							bool predVarState = r.varState(vidx, s);
							if (predVarState) {
								predVarTrue += sfreq;
							}
							stateFreq += sfreq;
						}
					}
					inputStateLogDep.depStateV(i) = log2(predVarTrue + prioriWeight_) - log2(f + prioriWeight_);
					explain<<"depStateV("<<i<<")=log2("<<predVarTrue<<"+"<<prioriWeight_<<") - log2("<<f<<"+"<<prioriWeight_<<")="<<inputStateLogDep.depStateV(i)<<"\n";
					inputStateLogDep.depStateNotV(i) = log2((stateFreq-predVarTrue)+prioriWeight_)  - log2(n - f + prioriWeight_);
					explain<<"depStateNotV("<<i<<")=log2("<<(stateFreq-predVarTrue)<<"+"<<prioriWeight_<<") - log2("<<(n-f)<<"+"<<prioriWeight_<<")="<<inputStateLogDep.depStateNotV(i)<<"\n";
				}
			}
			template <typename Out>
			inline void getInputStateLogDep(const explib::rel_stats<P>& rs,
											const rel_inputvars<P>& riv,
											inputstate_logdep& inputStateLogDep,
											Out& explain) {
				if (!useLogDepB_) {
					getInputStateLogDepA(rs, riv, inputStateLogDep, explain);
				} else {
					getInputStateLogDepB(rs, riv, inputStateLogDep, explain);
				}
			}



			/**
			 * Bit by bit predicting
			 */
			template <typename Out>
			std::vector<double> bitP(const data<P>& d,
									 int varid,
									 Out& explain) {
				const data_var<P>& dv = d.var(varid);
				const var<P>& v = dv.var_;
				const cvec<P>& dim = dv.dim();

				int volume = dim.volume();
				std::vector<double> rv(volume);
				std::vector<double> wTrue(volume);
				std::vector<double> wFalse(volume);
				inputstate_logdep inputStateLogDep(P::MAX_REL_VARS);
				{
					const var_stats<P>& vs = stats_.var(varid);
					double logPrioriTrue = log2(vs.freq() + 1);
					double logPrioriFalse = log2(vs.n() - vs.freq() + 1);
					for (int i = 0; i < volume; ++i) {
						wTrue[i] = logPrioriTrue;
						wFalse[i] = logPrioriFalse;
					}
					explain<<"log priori true:  "<<logPrioriTrue<<"\n";
					explain<<"log priori false: "<<logPrioriFalse<<"\n";
				}

				// predict using target variable relations
				const std::vector<rel<P>*>& rels = v.rels();
				for (auto i = rels.begin(); i != rels.end(); ++i) {
					const rel<P>& r = **i;
					if (!r.disabled()) {
						const rel_stats<P>& rs = stats_.rel(r.id());
						const data_rel<P>& dr = d.rel(r.id());
						auto re = r.entries();

						explain<<"rel "<<r.id()<<"\n";

						int vidx = 0; // this is the predicted one
						for (int j = 0; j < r.varCount(); ++j) {
							if (re[j].var_ == &v) {
								vidx = j;
							}
						}

						rel_inputvars<P> riv(r, vidx); // helper class
						getInputStateLogDep(rs, riv, inputStateLogDep, explain);

						explain<<"states:\n";
						for (int i = 0; i < riv.stateCount(); ++i) {
							explain<<i<<"   "<<inputStateLogDep.depStateV(i)<<"   "<<inputStateLogDep.depStateNotV(i)<<"\n";
						}

						for (dim_iterator<P> j(dr.dim()); j; ++j) {
							const cvec<P>& at = *j;

							int maskd = 0; // which vars are defined
							int masks = 0; // which defined vars are true
							for (int j = 0; j < r.varCount(); ++j) {
								if (re[j].var_ != &v) {
									auto b = dr.var(j)[at];
									if (b) {
										maskd |= 1<<j;
										masks |= (*b != 0)<<j;
									}
								}
							}
							int rivs = riv.maskState(maskd, masks);

							int sidx = dim.offset(at);

							wTrue[sidx]  += inputStateLogDep.depStateV(rivs);
							wFalse[sidx] += inputStateLogDep.depStateNotV(rivs);
						}
					}
				}

				for (int i = 0; i < volume; ++i) {
					rv[i] = exp2(wTrue[i]) / (exp2(wTrue[i]) + exp2(wFalse[i]));
					explain<<"var "<<varid<<" at "<<i<<":\n";
					explain<<rv[i]<<" = exp2("<<wTrue[i]<<")/(exp2("<<wTrue[i]<<")+exp2("<<wFalse[i]<<"))\n";
					explain<<rv[i]<<" = "<<exp2(wTrue[i])<<"/("<<exp2(wTrue[i])<<"+"<<exp2(wFalse[i])<<")\n";
				}

				return rv;
			}

			inline std::vector<double> bitP(const data<P>& d, int varid) {
				fake_output out;
				return bitP(d, varid, out);
			}

			/**
			 * Row predicting. This method of predicting is especially optimized
			 * for situation, where relation has many-to-single bits mapping.
			 */
			template <typename Out = fake_output>
			std::vector<double> rowP(const data<P>& d, int varid, Out explain = Out()) {
				const data_var<P>& dv = d.var(varid);
				const var<P>& v = dv.var_;
				const cvec<P>& dim = dv.dim();

				int volume = dim.volume();
				std::vector<double> rv(volume);
				std::vector<double> wTrue(volume);
				std::vector<double> wFalse(volume);
				inputstate_logdep inputStateLogDep(P::MAX_REL_VARS);

				{
					const var_stats<P>& vs = stats_.var(varid);
					double logPrioriTrue = log2(vs.freq() + 1);
					double logPrioriFalse = log2(vs.n() - vs.freq() + 1);
					for (int i = 0; i < volume; ++i) {
						wTrue[i] = logPrioriTrue;
						wFalse[i] = logPrioriFalse;
					}
					explain<<"log priori true:  "<<logPrioriTrue<<"\n";
					explain<<"log priori false: "<<logPrioriFalse<<"\n";
				}

				bits defrows[P::MAX_REL_VARS]; // which vars are defined
				bits staterows[P::MAX_REL_VARS]; // which defined vars are true
				bits relstates;
				bits tmpbits;

				const std::vector<rel<P>*>& rels = v.rels(); // use these for prediction work
				for (auto i = rels.begin(); i != rels.end(); ++i) {
					const rel<P>& r = **i;
					if (!r.disabled()) {
						const rel_stats<P>& rs = stats_.rel(r.id());
						const data_rel<P>& dr = d.rel(r.id());
						auto re = r.entries();

						explain<<"rel "<<r.id()<<"\n";

						int vidx = 0; // this is the predicted one
						for (int j = 0; j < r.varCount(); ++j) {
							if (re[j].var_ == &v) {
								vidx = j;
							}
						}

						// helper class
						rel_inputvars<P> riv(r, vidx);

						int rowdim = 0;
						const explib::ctx<P>& c = r.ctx();
						for (; c.v_[rowdim] < 0; ++rowdim) {}
						dim_row_iterator<P> j(dr.dim(), rowdim);

						getInputStateLogDep(rs, riv, inputStateLogDep, explain);

						for (int i = 0; i < P::MAX_REL_VARS; ++i) {
							staterows[i].resize(j.length());
							defrows[i].resize(j.length());
						}
						relstates.resize(j.length());
						tmpbits.resize(j.length());

						for (; j; ++j) {
							const cvec<P>& at = j.begin();

							// first populate state rows
							for (int k = 0; k < r.varCount(); ++k) {
								if (k != vidx) {
									shifted_data_var<P> v = dr.var(k);
									cond_bits_ref row = v.bitrow(at, j.length());
									defrows[k] = row; // defines
									staterows[k] = *row; // states
								}
							}
							// So we have a relation of for r(y, x1, x2, ...), where y
							// is interesting variable. In here, we are interested of
							// p(x1&x2&...|y) and p(x1&x2&...|!y) values, which are
							// of form p(x1&x2&...|y) = p(x1&x2&...&y) / p(y)

							// the predicted variable for this row
							int sidx = dim.offset(at);

							explain<<"var "<<varid<<" at "<<at<<"("<<sidx<<"):\n";
							for (int s = 0; s < riv.stateCount(); ++s) {
								if (inputStateLogDep.depStateV(s) != 0) {
									relstates.fill(true);
									for (int k = 0; k < r.varCount(); ++k) {
										if (k != vidx) {
											int defstate = riv.varDefState(k, s);

											if (defstate == rel_inputvars<P>::UNDEFINED) {
												relstates.andNeg(defrows[k]);
											} else {
												relstates &= defrows[k];
												if (defstate == rel_inputvars<P>::FALSE) {
													relstates.andNeg(staterows[k]);
												} else {
													relstates &= staterows[k];
												}
											}
										}
									}
									int statecount = relstates.popcount();

									wTrue[sidx]  += statecount * inputStateLogDep.depStateV(s);
									wFalse[sidx] += statecount * inputStateLogDep.depStateNotV(s);

									explain<<"state "<<s<<", freq="<<statecount<<"\n";
									explain<<"logTrue += "<<statecount<<" * "<<inputStateLogDep.depStateV(s)<<"\n";
									explain<<"logTrue -> "<<wTrue[sidx]<<"\n";
									explain<<"logFalse += "<<statecount<<" * "<<inputStateLogDep.depStateNotV(s)<<"\n";
									explain<<"logFalse -> "<<wFalse[sidx]<<"\n";
								}
							}
						}
					}
				}

				for (int i = 0; i < volume; ++i) {
					rv[i] = exp2(wTrue[i]) / (exp2(wTrue[i]) + exp2(wFalse[i]));

					explain<<"var "<<varid<<" at "<<i<<":\n";
					explain<<rv[i]<<" = exp2("<<wTrue[i]<<")/(exp2("<<wTrue[i]<<")+exp2("<<wFalse[i]<<"))\n";
					explain<<rv[i]<<" = "<<exp2(wTrue[i])<<"/("<<exp2(wTrue[i])<<"+"<<exp2(wFalse[i])<<")\n";
				}

				return rv;
			}

			std::vector<double> p(const data<P>& d, int varid) {
				const data_var<P>& dv = d.var(varid);
				const var<P>& v = dv.var_;
				if (v.ctx().v_[0] < 0) {
					return rowP(d, varid);
				} else {
					return bitP(d, varid);
				}
			}

			void info(const data<P>&d,
					  int varid,
					  double& totalinfo,
					  double& entryinfo) {
				std::vector<double> ps = p(d, varid);
				const explib::data_var<P>& dv = d.var(varid);
				explib::cvec<P> dim = dv.dim();
				double info = 0;
				int defined = 0;
				for (size_t i = 0; i < ps.size(); ++i) {
					if (dv[dim.at(i)]) {
						if (*dv[dim.at(i)]) {
							info += -log2(ps[i]);
						} else {
							info += -log2(1-ps[i]);
						}
						defined++;
					}
				}
				totalinfo = info;
				entryinfo = info / defined;
			}



	};

}

#endif /* PRED_H_ */
