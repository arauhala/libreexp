/*
 * printer.i.h
 *
 *  Created on: Dec 1, 2013
 *      Author: arau
 */

#ifndef REEXP_PRINTER_I_H_
#define REEXP_PRINTER_I_H_

#include "printer.h"
#include "pred.i.h"

namespace reexp {

	char char_symbol_for(int i);

	template <typename P>
	int count_diff_in_datas(const reexp::data<P>& d1,
			   	   	   	   	const reexp::data<P>& d2,
			   	   	   	   	int& ddiff,
			   	   	   	   	int& sdiff) {
		ddiff = 0;
		sdiff = 0;
		for (int i = d1.var_count() - 1; i >= 0; --i) {
			const reexp::data_var<P>& dv1 = d1.var(i);
			const reexp::data_var<P>& dv2 = d2.var(i);
			reexp::bits b;
			b = dv1.defined();
			b ^= dv2.defined();
			int dd = b.popcount();
			b = dv1.states();
			b ^= dv2.states();
			int sd = b.popcount();
			ddiff += dd;
			sdiff += sd;
		}
		return ddiff + sdiff;
	}

	template <typename P>
	lang_info<P>::lang_info(const pinfo& i, const reexp::lang<P>& l)
	: info_(i), lang_(l) {}

	template <typename P>
	std::string lang_info<P>::varid_tostring(size_t i) const {
		if (i < info_.vnames_.size()) {
			return info_.vnames_[i];
		} else {
			const reexp::var<P>& v = lang_.var(i);
			const reexp::exp<P>* e = dynamic_cast<const exp<P>*>(&v);
			std::ostringstream buf;
			if (e) {
				buf<<"e"<<i;
			} else {
				buf<<"v"<<i;
			}
			return buf.str();
		}
	}


	template <typename P>
	std::string lang_info<P>::shift_tostring(const cvec<P>& shift) const {
		std::ostringstream buf;
		buf<<"(r";
		for (int i = 0; i < P::DIM; ++i) {
			int s = shift[i];
			if (s > 0) {
				buf<<" + ";
			} else if (s < 0) {
				buf<<" - ";
				s *= -1;
			} else {
				continue;
			}
			if (s > 1) buf<<s;
			buf<<info_.cvnames_[i];
		}
		buf<<")";
		return buf.str();
	}

	template <typename P>
	std::string lang_info<P>::rel_tostring(const rel<P>& r, int state, const cvec<P>& shift) const {
		std::ostringstream buf;
		const std::vector<rel_entry<P> >& e = r.entries();
		for (size_t i = 0; i < e.size(); ++i) {
			if (!r.varState(i, state)) {
				buf<<"!";
			}
			buf<<var_tostring(*e[i].var_, shift+e[i].shift_);
			if (i + 1 < e.size()) {
				buf<<", ";
			}
		}
/*				if (r.disabled()) {
			buf<<" (DISABLED)";
		}*/
//				buf<<"@"<<r.hash();
		return buf.str();
	}

	template <typename P>
	std::string lang_info<P>::rel_tostring(int rid, int state, const cvec<P>& shift) const {
		return rel_tostring(lang_.rel(rid), state, shift);
	}

	template <typename P>
	std::string lang_info<P>::var_tostring(const var<P>& v, const cvec<P>& shift) const {
		std::ostringstream buf;
		if (v.id() < info_.vnames_.size()) {
			buf<<info_.vnames_[v.id()];
			buf<<shift_tostring(shift);
		} else {
			const exp<P>* e = dynamic_cast<const exp<P>*>(&v);
			if (e) {
				buf<<"<";
				buf<<rel_tostring(e->rel(), e->state(), shift);
				buf<<">";
			} else {
				buf<<"v"<<v.id();
				buf<<shift_tostring(shift);
			}
		}
		return buf.str();
	}
	template <typename P>
	std::string lang_info<P>::var_tostring(int i, const cvec<P>& shift) const {
		return var_tostring(lang_.var(i), shift);
	}
	template <typename P>
	std::string lang_info<P>::vars_tostring() const {
		std::ostringstream buf;
		for (int i = 0; i < lang_.var_count(); ++i) {
			buf<<var_tostring(i)<<std::endl;
		}
		return buf.str();
	}
	template <typename P>
	std::string lang_info<P>::rels_tostring() const {
		std::ostringstream buf;
		for (int i = 0; i < lang_.rel_count(); ++i) {
			buf<<rel_tostring(i)<<std::endl;
		}
		return buf.str();
	}

	template <typename P>
	void lang_info<P>::draw_rel(var_graphics& g, const reexp::rel<P>& r, size_t state) {
		for (size_t i = 0; i < r.entries().size(); i++) {
			const reexp::rel_entry<P>& re = r.entries()[i];
			translation_sentry s(g, re.shift_[g.xcvar_], re.shift_[g.ycvar_]);
			draw_var(g, *re.var_, bool((state>>i)&1));
		}
	}

	template <typename P>
	void lang_info<P>::draw_var(var_graphics& g, const reexp::var<P>& v, bool s) {
		const reexp::exp<P>* e = dynamic_cast<const reexp::exp<P>*>(&v);
		if (e) {
			draw_rel(g, e->rel(), s ? e->state() : 0);
		} else {
			g.set(0, 0, 1, s);
		}
	}

	template <typename P>
	std::string lang_info<P>::drawn_implmask_to_string(const implmask<P>& o) const {
		std::ostringstream buf;
		const reexp::var<P>& v = *o.var_;
		buf<<var_tostring(v.id())<<":\n";
		const cvec<P>& dim = o.mask_.ndim_.dim_;
		dim_row_iterator<P> i(dim, dim.rowdim());
		int length = i.length();
		while (i) {
			int offset = dim.offset(i.begin());
			for (int j = 0; j < length; ++j) {
				if (o.mask_.bits_[offset+j]) {
					buf<<'X';
				} else {
					buf<<'.';
				}
			}
			buf<<"    "<<(i.begin()+o.mask_.ndim_.shift_)<<"\n";
			++i;
		}
		return buf.str();
	}

	template <typename P>
	std::string lang_info<P>::drawn_var_rootmasks_to_string(const reexp::var<P>& var) const {
		std::ostringstream buf;
		for (const implmask<P>& o : var.rootmasks()) {
			buf<<drawn_implmask_to_string(o);
		}
		return buf.str();
	}

	template <typename P>
	std::string lang_info<P>::drawn_var_expmasks_to_string(const reexp::var<P>& var) const {
		std::ostringstream buf;
		for (const implmask<P>& o : var.expmasks()) {
			buf<<drawn_implmask_to_string(o);
		}
		return buf.str();
	}

	template <typename P>
	std::string lang_info<P>::drawn_var_implmasks_to_string(const reexp::var<P>& var) const {
		std::ostringstream buf;
		buf<<"rootmasks:\n";
		buf<<drawn_var_rootmasks_to_string(var);
		buf<<"\nexclmasks:\n";
		buf<<drawn_var_expmasks_to_string(var);
		return buf.str();
	}

	template <typename P>
	std::string lang_info<P>::drawn_var_tostring(int xcvar, int ycvar, const reexp::var<P>& v) const {
		bitmap bm;
		bitmap bm2;
		var_graphics g(bm, bm2, xcvar, ycvar);
		draw_var(g, v, true);
		return g.to_string(char_symbol_for(v.id()));
	}

	template <typename P>
	std::string lang_info<P>::drawn_rel_tostring(int xcvar, int ycvar, const reexp::rel<P>& r) const {
		bitmap bm;
		bitmap bm2;
		var_graphics g(bm, bm2, xcvar, ycvar);
		draw_rel(g, r);
		return g.to_string('X');
	}


	template <typename P>
	std::string lang_info<P>::drawn_vars_tostring(int xcvar, int ycvar) const {
		std::ostringstream buf;
		for (int i = 0; i < lang_.var_count(); ++i) {
			buf<<var_tostring(i)<<"\n\n";
			buf<<drawn_var_tostring(xcvar, ycvar, lang_.var(i))<<"\n\n";
		}
		return buf.str();
	}


	template <typename P>
	std::string lang_info<P>::drawn_rels_tostring(int xcvar, int ycvar) const {
		std::ostringstream buf;
		for (int i = 0; i < lang_.rel_count(); ++i) {
			buf<<rel_tostring(i)<<"\n\n";
			buf<<drawn_rel_tostring(xcvar, ycvar, lang_.rel(i))<<"\n\n";
		}
		return buf.str();
	}


	template <typename P>
	std::string lang_info<P>::invorder_diff_tostring(const reexp::data<P>& d1,
													 const reexp::data<P>& d2,
													 int maxprints,
													 const reexp::index_over_var_bits<P>* indexspace,
													 const std::vector<int>* d2KnownAt) const {
		std::ostringstream buf;
		int errors = 0;
		int ddiff = 0;
		int sdiff = 0;
		for (int i = d1.var_count() - 1; i >= 0; --i) {
			const reexp::data_var<P>& dv1 = d1.var(i);
			const reexp::data_var<P>& dv2 = d2.var(i);
			reexp::bits b;
			b = dv1.defined();
			b ^= dv2.defined();
			int dd = b.popcount();
			b = dv1.states();
			b ^= dv2.states();
			int sd = b.popcount();
			ddiff += dd;
			sdiff += sd;

			if (dd || sd) {
				buf<<"variable "<<i<<" "<<var_tostring(i);
				if (sd || dd) {
					buf<<" differs.\n";
				} else {
					buf<<" ok.\n";
				}
				buf<<dd<<" diffs in defined.\n";
				buf<<sd<<" diffs in states.\n";
				cvec<P> dim = dv1.dim();
				for (int i = dv1.defined().size()-1; i >= 0 && errors < maxprints; --i) {
					if (dv1.defined()[i] != dv2.defined()[i]) {
						buf<<"defined "<<dv1.defined()[i]<<" != "<<dv2.defined()[i]<<" at "<<dim.at(i)<<"\n";
						errors++;
						if (d2KnownAt && !dv2.defined()[i]) {
							int index = (*d2KnownAt)[indexspace->from_var_bit(dv2.var_.id(), i)];
							int detvarid = 0, detvarstate = 0;
							if (indexspace->to_var_bit(index, detvarid, detvarstate)) {
								const data_var<P>& detvar = d2.var(detvarid);
								cvec<P> detvarat = detvar.dim().at(detvarstate);
								buf<<"d2 determined by "<<detvarid<<" "<<var_tostring(detvarid)<<" state "<<detvar.states()[detvarstate]<<" at "<<detvarat<<"\n";
								const reexp::exp<P>* detexp = dynamic_cast<const reexp::exp<P>*>(&detvar.var());
								if (detexp) {
									for (const rel_entry<P>& e : detexp->rel().entries()) {
										if (e.var_ != &dv2.var()) {
											// print the other variable known at
											cvec<P> rvarat = detvarat + e.shift_;
											const data_var<P>& rvar = d2.var(e.var_->id());
											int rvidx = rvar.index(rvarat);
											int index = indexspace->from_var_bit(e.var_->id(), rvidx);
											buf<<"   " <<e.var_->id()<<" "<<var_tostring(*e.var_)<<":\n";
											buf<<"   defined " << rvar.defined()[rvidx]<<"\n";
											buf<<"   state "<< rvar.states()[rvidx] << " at " <<rvarat<<"\n";

											int kat = (*d2KnownAt)[index];
											int detvarid = 0, detvarstate = 0;
											if (indexspace->to_var_bit(kat, detvarid, detvarstate)) {
												const data_var<P>& detvar = d2.var(detvarid);
												cvec<P> detvarat = detvar.dim().at(detvarstate);
												buf<<"   determined by "<<detvarid<<" " <<var_tostring(detvarid)<<" state "<<*detvar[detvarat]<<" at "<<detvarat<<"\n";
											}
										}
									}
								}
							}
						}
					}
					if (dv1.states()[i] != dv2.states()[i]) {
						buf<<"state "<<dv1.states()[i]<<" != "<<dv2.states()[i]<<" at "<<dim.at(i)<<"\n";
						errors++;
					}
				}
				buf<<"\n";
			}
		}
		if (errors >= maxprints) {
			buf<<"no more differences printed as maximum ("<<maxprints<<") reached.\n";
		}
		if (ddiff || sdiff) {
			buf<<ddiff<<" diffs in defined.\n";
			buf<<sdiff<<" diffs in states.\n";
		} else {
			buf<<"parts are identical\n";
		}
		return buf.str();
	}

	template <typename P>
	const reexp::lang<P>& lang_info<P>::lang() const {
		return lang_;
	}

	template <typename P>
	stats_info<P>::stats_info(const pinfo& i, const reexp::stats<P>& s)
	: info_(i), lang_(i, s.data().lang()), stats_(s) {}

	template <typename P>
	const pinfo& stats_info<P>::info() const {
		return info_;
	}

	template <typename P>
	std::string stats_info<P>::drawn_data_var_tostring(int xvarid,
													   int yvarid,
													   const reexp::data_var<P>& pixels,
													   const reexp::cvec<P>& at) const {
		std::ostringstream buf;
		reexp::cvec<P> d = pixels.dim();
		reexp::cvec<P> a( at );
		for (int i = at[yvarid]; i < d[yvarid]; i++) {
			a[yvarid] = i;
			for (int j = at[xvarid]; j < d[xvarid]; j++) {
				a[xvarid] = j;
				auto b = pixels[a];
				if (b) {
					if (*b) {
						buf<<"X";
					} else {
						buf<<".";
					}
				} else {
					if (*b) buf<<"x";
					else buf<<" ";
				}
			}
			buf<<"\n";
		}
		return buf.str();
	}
	template <typename P>
	std::string stats_info<P>::drawn_data_tostring(cvec<P> at, int xcvar, int ycvar) const {
		const reexp::data<P>& d = stats_.data();
		return drawn_data_tostring(at, xcvar, ycvar, d);
	}

	template <typename P>
	std::string stats_info<P>::drawn_data_tostring(cvec<P> at, int xcvar, int ycvar, const reexp::data<P>& d) const {
		int w = d.dim()[xcvar];
		int h = d.dim()[ycvar];
		std::string chars;
		chars.resize((w+1)*h); // leave \n marks
		for (size_t i = 0; i < chars.size(); i++) {
			if (int(i % (w+1)) == w) {
				chars[i] = '\n';
			} else {
				chars[i] = '.';
			}
		}

		for (int i = 0; i < d.lang().var_count(); ++i) {
			char c = char_symbol_for(i);
			const data_var<P>& dv = d.var(i);
			for (int x = 0; x < dv.dim()[xcvar]; x++) {
				at[xcvar] = x;
				for (int y = 0; y < dv.dim()[ycvar]; y++) {
					at[ycvar] = y;
					const_cond_bit_ref cbit = dv[at];
					if (cbit && *cbit) {
						bitmap bm;
						bitmap bm2;
						var_graphics g(bm, bm2, xcvar, ycvar);
						lang_.draw_var(g, dv.var(), true);
						for (int x2 = 0; x2 < bm.w_; ++x2) {
							for (int y2 = 0; y2 < bm.h_; ++y2) {
								if (bm.get(x2, y2) && bm2.get(x2, y2)) {
									chars[(x+x2)+(y+y2)*(w+1)] = c;
								}
							}
						}
					}
				}
			}
		}
		return chars;
	}

	template <typename P>
	std::string stats_info<P>::drawn_data_tostring(cvec<P> at, int xcvar, int ycvar, int itercvar, int n) const {
		std::ostringstream buf;
		for (int i = 0; i < n; ++i) {
			buf<<i<<":\n";
			buf<<drawn_data_tostring(at, xcvar, ycvar)<<"\n";
			at[itercvar] = at[itercvar]+1;
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::var_tostring(const var_stats<P>& vs, const cvec<P>& shift) const {
		std::ostringstream buf;
		buf<<lang_.var_tostring(vs.data().var(), shift)<<"  "<<vs.freq()<<"/"<<vs.n()<<" ["<<vs.dl()<<"]\n";
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::var_tostring(int i, const cvec<P>& shift) const {
		return var_tostring(stats_.var(i), shift);
	}

	template <typename P>
	std::string stats_info<P>::vars_tostring() const {
		std::ostringstream buf;
		for (int i = 0; i < stats_.data().lang().var_count(); ++i) {
			buf<<var_tostring(i)<<std::endl;
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::detailed_vars_tostring(int xcvar, int ycvar) const {
		std::ostringstream buf;
		const lang<P>& lang = stats_.data().lang();
		for (int i = 0; i < lang.var_count(); ++i) {
			buf<<lang_.var_tostring(i)<<"\n\n";
			buf<<lang_.drawn_var_tostring(xcvar, ycvar, lang.var(i))<<"\n";
			buf<<var_deps_tostring(i)<<"\n\n";
		}
		return buf.str();
	}


	template <typename P>
	std::string stats_info<P>::rel_tostring(const rel_stats<P>& rs, int state, const cvec<P>& shift) const {
		std::ostringstream buf;
		const rel<P>& r = rs.data().rel();
		buf<<lang_.rel_tostring(r, state, shift)<<"\n    ";
		if (!r.disabled()) {
			for (int i = 0; i < rs.stateCount(); ++i) {
				buf<<rs.stateFreqs()[i]<<" ";
			}
			buf<<"/ "<<rs.n()<<"\n    ";
			for (int i = 0; i < rs.stateCount(); ++i) {
				buf<<rs.stateEntropy(i)<<" ";
			}
			buf<<"/ "<<rs.entropy()<<"\n    ";
			double sum = 0;
			for (int i = 0; i < rs.stateCount(); ++i) {
				double ne = rs.eStateP(i) * rs.eStateDl(i);
				buf<<ne<<" ";
				sum += ne;
			}
			buf<<"/ "<<sum<<"\n    ";
			for (int i = 0; i < rs.stateCount(); ++i) {
				buf<<rs.eStateBias(i)<<" ";
			}
			buf<<"\n";
		} else {
			buf<<"DISABLED\n";
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::rel_tostring(int i, int state, const cvec<P>& shift) const {
		return rel_tostring(stats_.rel(i), state, shift);
	}

	template <typename P>
	std::string stats_info<P>::rels_tostring() const {
		std::ostringstream buf;
		for (int i = 0; i < stats_.data().lang().rel_count(); ++i) {
			buf<<rel_tostring(i)<<std::endl;
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::scan_tostring(size_t max, int details) const {
		std::ostringstream buf;
		reexp::learner<P> learner(const_cast<reexp::lang<P>&>(stats_.data().lang()),
								   const_cast<reexp::stats<P>&>(stats_));

		std::priority_queue<reexp::candidate<P> > cands;
		learner.scan(cands);

		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);
		buf.setf(std::ios::left, std::ios::adjustfield);

		size_t n = std::min(max, cands.size());
		for (size_t i = 0; i < n; ++i) {
			const reexp::candidate<P>& c( cands.top() );
			buf<<c.bias_<<"   ";
			buf.width(50);
			buf<<lang_.rel_tostring(c.rel_->data().rel_, c.state_);
			if (details & 1) {
				const rel_stats<P>& rs = *c.rel_;
				buf<<"    [f(s)="<<rs.stateFreqs()[c.state_];
				for (int i = 0; i < rs.varCount(); ++i) {
					bool s = rs.data().rel().varState(i, c.state_);
					buf<<", f(";
					if (!s) buf<<"!";
					buf<<lang_.var_tostring(rs.data().rel().entries()[i].var_->id());
					int f = s?rs.varFreqs()[i]:rs.n()-rs.varFreqs()[i];
					buf<<")="<<f;
				}
				buf<<"]";
			}
			if (details & 2) {
				const rel_stats<P>& rs = *c.rel_;
				buf<<"    [n(s)="<<rs.stateFreqs()[c.state_];
				buf<<",ndl(s)="<<rs.stateDl(c.state_);
				buf<<",info(s)="<<rs.stateInfo(c.state_);
				buf<<",endl(s)="<<rs.eStateDl(c.state_);
				buf<<",eInfo(s)="<<rs.eStateInfo(c.state_);
				buf<<"]";
			}
			buf<<"\n";
			cands.pop();
		}
		return buf.str();
	}

	template <typename P>
	void stats_info<P>::var_scan(const reexp::var<P>& var,
				  std::priority_queue<candidate<P>>& cands) const {
		const std::vector<rel<P>*>& rels = var.rels();
		for (size_t i = 0; i < rels.size(); i++) {
			const rel<P>& r( *rels[i] );
			const rel_stats<P>& rs( stats_.rel(r.id()) );
			int stateCount = r.stateCount();
			for (int j = 0; j < stateCount; j++) {
				double b = rs.eStateBias(j);
				if (b > 0) {
					candidate<P> cand(rs, j, b);
					cands.push(cand);
				}
			}
		}
	}

	template <typename P>
	double stats_info<P>::calc_influence(int from, const reexp::bits& to) const {
		const reexp::stats<P>& s = stats();
		const reexp::lang<P>& l = s.data().lang();
		double rv = 0;

		for (int i = 0; i < l.rel_count(); ++i) {
			const reexp::rel<P>& r = l.rel(i);
			bool inflFound = 0;
			bool efound = 0;
			for (size_t j = 0; j < r.entries().size(); ++j) {
				int id = r.entries()[j].var_->id();
				if (id == from) inflFound = true;
				if (id < to.size() && to[id]) efound = true;
			}
			if (inflFound && efound) {
				const reexp::rel_stats<P>& rs = s.rel(i);
				rv += rs.n() * (rs.eNaiveEntropy() - rs.eEntropy());
			}
		}
		return rv;
	}

	template <typename P>
	std::string stats_info<P>::top_influence_tostring(const reexp::bits& from,
									   	   	   	      const reexp::bits& to,
									   	   	   	      int cap) const {
		std::priority_queue< std::pair<double, int> > cands;

		for (int i = 0; i < from.size(); ++i) {
			if (from[i]) {
				cands.push(std::pair<double, int>(calc_influence(i, to), i) );
			}
		}

		std::ostringstream buf;
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);
		buf.setf(std::ios::left, std::ios::adjustfield);
		while (!cands.empty() && cap--) {
			const std::pair<double, int>& c = cands.top();
			buf.width(6);
			buf<<c.first<<" "<<lang_.var_tostring(c.second)<<"\n";

			cands.pop();
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::var_deps_tostring(const reexp::var<P>& var, int cap) const {
		std::priority_queue<candidate<P>> cands;
		var_scan(var, cands);

		std::ostringstream buf;
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);
		buf.setf(std::ios::left, std::ios::adjustfield);

		while (!cands.empty() && cap--) {
			const candidate<P>& c = cands.top();
			const rel_stats<P>& rs = *c.rel_;
			const rel<P>& r = rs.data().rel();
			buf<<c.bias_<<" ";
			buf.width(50);
			buf<<lang_.rel_tostring(c.rel_->data().rel(), c.state_);
			buf<<"    [f(s)="<<rs.stateFreqs()[c.state_]<<", n="<<rs.n();
			for (int i = 0; i < rs.varCount(); ++i) {
				bool s = r.varState(i, c.state_);
				buf<<", f(";
				if (!s) buf<<"!";
				buf<<lang_.var_tostring(r.entries()[i].var_->id());
				int f = s?rs.varFreqs()[i]:rs.n()-rs.varFreqs()[i];
				buf<<")="<<f;
			}
			buf<<"]\n";
			cands.pop();
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::var_deps_tostring(int var, int cap) const {
		return var_deps_tostring(stats_.data().lang().var(var), cap);
	}

	template <typename P>
	std::string stats_info<P>::pred_tostring(const reexp::data<P>& d, int var, bool showbit) const {
		std::ostringstream buf;
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);

		reexp::pred<P> pred(stats_);
		std::vector<double> p = pred.p(d, var);
		const reexp::data_var<P>& dv = d.var(var);
		reexp::cvec<P> dim = dv.dim();
		for (size_t i = 0; i < p.size(); ++i) {
			buf<<"p("<<lang_.varid_tostring(var)<<"|"<<dim.at(i)<<") = "<<p[i]<<"%";
			if (showbit) {
				if (*dv[dim.at(i)]) {
					buf<<" [true]\n";
				} else {
					buf<<" [false]\n";
				}
			} else {
				buf<<"\n";
			}

		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::pred_tostring(int var, bool showbit) const {
		return pred_tostring(stats_.data(), var, showbit);
	}

	template <typename P>
	std::string stats_info<P>::row_logdep_tostring(int var, int var2) const {
		std::ostringstream buf;
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);

		reexp::pred<P> pred(stats_);
		const reexp::var<P>& v = stats_.data().lang().var(var);
		const std::vector<rel<P>*>& rels = v.rels(); // use these for prediction work
		inputstate_logdep logdep(P::MAX_REL_VARS);
		for (auto i = rels.begin(); i != rels.end(); ++i) {
			const rel<P>& r = **i;
			const rel_stats<P>& rs = stats_.rel(r.id());
			auto re = r.entries();
			int vidx = -1, v2idx = -1;
			for (int j = 0; j < r.varCount(); ++j) {
				if (re[j].var_ == &v) vidx = j;
				if (re[j].var_->id() == size_t(var2)) v2idx = j;
			}
			if (var2 >= 0 && v2idx == -1) continue;
			rel_inputvars<P> riv(r, vidx);

			buf<<"predicted: "<<lang_.var_tostring(v)<<"\n";
			for (int s = 0; s < riv.stateCount(); ++s) {
				std::ostringstream state;
				int istates = 0;
				for (int j = 0; j < rs.varCount(); ++j) {
					if (j != vidx) {
						int ds = riv.varDefState(j, s);
						if (ds == 0) {
							if (istates++) buf<<", ";
							state<<"!"<<lang_.var_tostring(*re[j].var_);
						} else if (ds == 1) {
							if (istates++) buf<<", ";
							state<<lang_.var_tostring(*re[j].var_);
						}
					}
				}
				buf<<"state "<<s<<": "<<state.str()<<"\n";
			}
			pred.getInputStateLogDep(rs, riv, logdep, static_cast<std::ostream&>(buf));
		}
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::entropy_tostring(int var) const {
		std::ostringstream buf;
		buf.setf(std::ios::fixed,std::ios::floatfield);
		buf.precision(2);

		double info = 0, entryinfo = 0;
		reexp::pred<P> pred(stats_);
		pred.info(stats_.data(), var, info, entryinfo);

		buf<<"total info: "<<info<<"\n";
		buf<<"entry info: "<<entryinfo<<"\n";
		return buf.str();
	}

	template <typename P>
	std::string stats_info<P>::preds_tostring(int begin, int end, bool showbit) const {
		std::ostringstream buf;
		for (int i = begin; i < end; ++i) {
			buf<<pred_tostring(i, showbit);
		}
		return buf.str();
	}

	template <typename P>
	const reexp::lang_info<P>& stats_info<P>::lang_info() const {
		return lang_;
	}

	template <typename P>
	const reexp::stats<P>& stats_info<P>::stats() const {
		return stats_;
	}

	template <typename P>
	double stats_info<P>::var_naive_p(const std::vector<double>& origvarPs, int v) const {
		const reexp::exp<P>* e = dynamic_cast<const reexp::exp<P>* >(&lang_.lang().var(v));
		double rv = 1;
		if (e) {
			for (auto entry : e->rel().entries()) {
				rv *= var_naive_p(origvarPs, entry.var_->id());
			}
		} else {
			rv *= origvarPs[v];
		}
		return rv;
	}

	template <typename P>
	std::string stats_info<P>::drawn_vars_tostring_byexpbias(const std::vector<double>& origvarPs, int xcvar, int ycvar) const {
		std::ostringstream buf;
		std::vector<std::pair<double, int> > ordered;
		for (int i = 0; i < lang_.lang().var_count(); ++i) {
			// simplistic math to get the most significant
			// expressions visible.
			//
			double p = var_naive_p(origvarPs, i);
			double cvalue = -log(p) * stats_.var(i).freq();
			ordered.push_back(std::pair<double, int>(cvalue, i));
		}
		std::sort(ordered.begin(),
				  ordered.end(),
				  [] (const std::pair<double, int>& e1,
					  const std::pair<double, int>& e2) {
					  return e1.first > e2.first;
				   });
		for (size_t i = 0; i < ordered.size(); ++i) {
			int v = ordered[i].second;
			if (stats_.var(v).freq() == 0) break;
			buf<<lang_.var_tostring(v)<<"\n";
			buf<<"freq: "<<stats_.var(v).freq()<<"/"<<stats_.var(v).n()<<"\n\n";
			buf<<lang_.drawn_var_tostring(xcvar, ycvar, lang_.lang().var(v))<<"\n\n";
		}
		return buf.str();
	}

}

#endif /* REEXP_PRINTER_I_H_ */
