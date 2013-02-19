/*
 * pred.h
 *
 *  Created on: Aug 16, 2011
 *      Author: arau
 */

#ifndef REEXP_PRED_H_
#define REEXP_PRED_H_


#include "stats.h"

namespace reexp {

	inline int pow3i(int v) {
		int rv = 1;
		while (v-- > 0) rv *= 3;
		return rv;
	}
	template <typename P>
	struct rel_inputvars {
		const reexp::rel<P>& r_;
		int predicted_;
		rel_inputvars(const reexp::rel<P>& r,
					  int predicted)
		: 	r_(r), predicted_(predicted) {}
		inline int stateCount() const {
			return pow3i(r_.varCount()-1);
		}
		enum DefState {
			FALSE = 0,
			TRUE = 1,
			UNDEFINED = 2
		};
		/**
		 * Transforms relation state into rel_inputvars state.
		 * dmask contains the definition bits, while smask contains
		 * the state bits.
		 */
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
			const reexp::stats<P>& stats_;
			double prioriWeight_;
			bool useLogDepB_;

		public:
			const reexp::stats<P>& stats() const { return stats_; }

			pred(const reexp::stats<P>& s, double prioriWeight = 2., bool useLogDepB = false)
			: stats_(s), prioriWeight_(prioriWeight), useLogDepB_(useLogDepB) {}


#if 0
			/**
			 * Probabilistic re-expression function. Input data
			 * must be in original language format.
			 */
			std::vector<double> reexpress(const reexp::data<P>& d) {
				std::vector<double> states2;
				std::vector<double> states;
				std::vector<double> ns;
				const lang<P>& l = d.lang();
				const lang<P>& tl = stats_.data().lang();

				int volume = 0;
				for (int i = 0; i < tl.var_count(); ++i) {
					volume += tl.var(i).ctx().dim(d.dim()).volume();
				}
				states.resize(volume);
				states2.resize(volume);
				ns.resize(volume);
				int w = 0;

				std::vector<int> voffsets;
				for (int i = 0; i < tl.var_count(); ++i) {
					voffsets.push_back(w);
					const var<P>& v = l.var(i);
					const exp<P>* e = dynamic_cast<const exp<P>*>(&v);
					if (e) {
						const rel<P>& er = e->rel();
						const exp_rel_stats<P>& ers = stats_.exp_rel(e->id());
						size_t estate = e->state();
						cvec<P> edim = e->ctx().dim(d.dim());
						size_t evol = edim.volume();
						std::vector<double> notfirst;
						notfirst.resize(evol);

						for (double& p : notfirst) p = 1.;

						for (size_t i = 0; i < evol; ++i) {
							states[w + i] = 1.;
						}

						//
						// There are basically three steps
						//
						// 1. Adjust the priori probabilities based on relation
						// 2. Estimate probabilities for expressions
						// 3. Make redundancy elimination on probabilities
						//

						//
						// 1. Adjust priori probabilities
						//
						// 1.1. first copy the priories into the lookup
						//
						states2 = states;

						//
						// 1.2. The basic idea is to adjust the estimate for variable based on
						// 		each variable's state.
						//
						const std::vector<rel_entry<P> >& eres = er.entries();
						for (size_t vi = 0; vi < eres.size(); ++vi) {
							const rel_entry<P>& ere = eres[vi];
							const var<P>& erv = *ere.var_;
							const cvec<P>& shift = ere.shift_;
							const data_var<P>& erdv = d.var(erv.id());
							cvec<P> ervdim = erdv.dim();
							size_t off = voffsets[erv.id()] + ervdim.offset(shift);


							rel_inputvars<P> riv(er, vi);
							inputstate_logdep logdep(P::MAX_REL_VARS);
							fake_output out;
							getInputStateLogDepB(ers.n(), ers.varFreqs(), ers.stateFreqs(), riv, logdep, out);

							// this may be good to turn into row form
							dim_iterator<P> it(edim);
							while (it) {
								int maskd = 0; // which vars are defined
								int masks = 0; // which defined vars are true

								for (int j = 0; j < ere.size(); ++j) {
									if (j != vi) {
										cvec<P> jpos = *it + ere[j].shift_;
										int jid = ere[j].var_->id();
										const data_var<P>& jdv = d.var(jid);
										int jidx = voffsets[j] + jdv.offset(jpos);
										double f = states[jidx];
										double n = ns[jidx];
										(double)f /(double)n;
									}
								}
								int rivs = riv.maskState(maskd, masks);
								++it;

								cvec<P> ervpos = *it + shift;
								int ervidx = w + edim.offset(ervpos);
								double pPos = states2[ervidx];
								double pNeg = 1 - pPos;
								pPos *= exp2(logdep.depStateV(rivs));
								pNeg *= exp2(logdep.depStateNotV(rivs));
								states[ervidx] = pPos / (pPos + pNeg);
							}
						}

						//
						// 2. Estimate expression probabilities
						//
						for (size_t vi = 0; vi < er.entries().size(); ++vi) {
							const rel_entry<P>& ere = er.entries()[vi];
							const var<P>& erv = *ere.var_;
							const cvec<P>& shift = ere.shift_;
							const data_var<P>& erdv = d.var(erv.id());
							cvec<P> ervdim = erdv.dim();
							size_t off = voffsets[erv.id()] + ervdim.offset(shift);

							size_t v = 0;
							dim_row_iterator<P> rows(ervdim);
							size_t rlen = edim[rows.rowdim_];
							while (rows) {
								int rowbegin = off + ervdim.offset(rows.begin());
								for (size_t i = 0; i < rlen; ++i) {
									double d;

									if (er.varState(vi, estate)) {
										d = states[rowbegin + i];
									} else {
										d = (1-states[rowbegin + i]);
									}
									if (vi) notfirst[v] *= d;
									states[w + v++] *= d;
								}
								++rows;
							}
						}

						//
						// 3. make redundancy elimination
						//
						for (size_t vi = 0; vi < er.entries().size(); ++vi) {
							const rel_entry<P>& ere = er.entries()[vi];
							const var<P>& erv = *ere.var_;
							const cvec<P>& shift = ere.shift_;
							const data_var<P>& erdv = d.var(erv.id());
							cvec<P> ervdim = erdv.dim();
							int off = voffsets[erv.id()] + ervdim.offset(shift);

							size_t v = 0;
							dim_row_iterator<P> rows(ervdim);
							size_t rlen = edim[rows.rowdim_];
							while (rows) {
								int rowbegin = off + ervdim.offset(rows.begin());
								for (size_t i = 0; i < rlen; ++i) {
									double expP = states[w + v];
									if (er.varState(vi, estate)) {
										states[rowbegin + i] -= expP;
									}
									ns[rowbegin + i] -= expP;
									if (vi == 0) {
										double d = notfirst[v] - expP;
										if (er.varState(vi, estate)) {
											states[rowbegin + i] -= d;
										}
										ns[rowbegin + i] -= d;
									}
									v++;
								}
								++rows;
							}
						}
						w+= evol;
					} else {
						const data_var<P>& dv = d.var(i);
						double defP = stats_.var(i).eP();
						int v = dv.dim().volume();
						const cond_bits& b = dv.bits();
						for (int i = 0; i < v; ++i) {
							double p;
							if (b[i]) {
								p = *b[i] ? 1. : 0.;
							} else {
								p = defP;
							}
							states[w++] = p;
						}
					}
				}
				return states;
			}
#endif


			template <typename Out>
			void getInputStateLogDepA(const rel_stats<P>& rs,
									  const rel_inputvars<P>& riv,
									  inputstate_logdep& inputStateLogDep,
									  Out& explain) const {
				const reexp::rel<P>& r = rs.data().rel();
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
			void getInputStateLogDepB(const rel_stats<P>& rs,
									  const rel_inputvars<P>& riv,
									  inputstate_logdep& inputStateLogDep,
									  Out& explain) const {
				getInputStateLogDepB<Out>(rs.n(),
								     	  rs.varFreqs(),
								     	  rs.stateFreqs(),
								     	  riv,
								     	  inputStateLogDep,
								     	  explain);
			}


			template <typename Out>
			void getInputStateLogDepB(int n,
									  const std::vector<int>& varFreqs,
									  const std::vector<int>& stateFreqs,
									  const rel_inputvars<P>& riv,
									  inputstate_logdep& inputStateLogDep,
									  Out& explain) const {
//				const reexp::rel<P>& r = rs.data().rel();
				int vidx = riv.predicted_;
				double f = varFreqs[vidx];
				for (int i = 0; i < riv.stateCount(); ++i) {
					double predVarTrue = 0;
					double stateFreq = 0;
					for (size_t s = 0; s < stateFreqs.size(); ++s) {
						if (riv.includesRelState(i, s)) {
							int sfreq = stateFreqs[s];
							bool predVarState = rel<P>::varState(vidx, s);
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
			inline void getInputStateLogDep(const reexp::rel_stats<P>& rs,
											const rel_inputvars<P>& riv,
											inputstate_logdep& inputStateLogDep,
											Out& explain) const {
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
									 Out& explain) const {
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

			inline std::vector<double> bitP(const data<P>& d, int varid) const {
				fake_output out;
				return bitP(d, varid, out);
			}

			/**
			 * Row predicting. This method of predicting is especially optimized
			 * for situation, where relation has many-to-single bits mapping.
			 */
			template <typename Out = fake_output>
			std::vector<double> rowP(const data<P>& d, int varid, Out explain = Out()) const {
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
						const reexp::ctx<P>& c = r.ctx();
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

			std::vector<double> p(const data<P>& d, int varid) const {
				const data_var<P>& dv = d.var(varid);
				const var<P>& v = dv.var_;
				if (v.ctx().v_[0] < 0) {
					return rowP(d, varid);
				} else {
					return bitP(d, varid);
				}
			}

			void info(const std::vector<double>& ps,
					  const data<P>& d,
					  int varid,
					  double& totalinfo,
					  double& entryinfo) const {
				const reexp::data_var<P>& dv = d.var(varid);
				reexp::cvec<P> dim = dv.dim();
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

			void info(const data<P>&d,
					  int varid,
					  double& totalinfo,
					  double& entryinfo) const {
				info(p(d, varid), d, varid, totalinfo, entryinfo);
			}

	};

	template <typename P>
	class group_scaler {
		private:
			std::vector<double> scales_;
			std::vector<double> avers_;

		public:
			group_scaler(const pred<P>& pr, int var, int target_group_sz)
			: scales_(), avers_() {
				// try to fit at least 128 items in each group
				int groups = pr.stats().data().var(var).states().size() / target_group_sz;

				if (groups > 3) {
					std::vector<double> ps = pr.p(pr.stats().data(), var);

					std::vector<std::pair<double, int> > pi;
					for (size_t i = 0; i < ps.size(); ++i) {
						pi.push_back({ps[i], i});
					}
					std::sort(pi.begin(), pi.end());
					std::reverse(pi.begin(), pi.end());

					size_t at = 0;
					const data_var<P>& dv = pr.stats().data().var(var);
					double e_p = pr.stats().var(var).eP();
					for (int i = 0; i < groups; ++i) {
						size_t group_sz = (ps.size() - at) / (groups - i);
						double sum_p = 0;
						size_t f = 0;
						for (size_t j = 0; j < group_sz; ++j) {
							if (dv.states()[pi[at].second]) f++;
							sum_p += pi[at].first;
							at++;
						}
						double p_aver = (sum_p + e_p) / double(group_sz + 2);
						double group_aver = (f + e_p) / double(group_sz + 2);
						avers_.push_back(p_aver);
						scales_.push_back(group_aver / p_aver);
					}
				}
			}

			void scale(std::vector<double>& ps) {
				if (avers_.size()) {
					for (double& p : ps) {
						if (p >= avers_[0]) {
							p *= scales_[0];
						} else if (p <= avers_.back()) {
							p *= scales_.back();
						} else {
							size_t i = 0;
							while (p < avers_[i+1]) ++i;
							double w = avers_[i] - p;
							double w2 = p - avers_[i+1];
							double s = ((scales_[i] * w2 + scales_[i+1] * w) / (w+w2));
							p *= s;
						}
						p = std::min(p, 1.);
					}
				}
			}

	};

}

#endif /* REEXP_PRED_H_ */
