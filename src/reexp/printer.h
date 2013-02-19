/*
 * print.h
 *
 *  Created on: Aug 25, 2011
 *      Author: arau
 */

#ifndef REEXP_PRINTER_H_
#define REEXP_PRINTER_H_

#include "reexp/learner.h"
#include "reexp/pred.h"
#include <iostream>
#include <sstream>
#include <string>

namespace reexp {

	inline char char_symbol_for(int i) {
		if (i < 10) {
			return '0' + i;
		}
		i -= 10;
		if (i <= ('z'-'a')) {
			return 'a' + i;
		}
		i -= ('z'-'a');
		if (i <= ('Z'-'A')) {
			return 'A' + i;
		}
		return '?';
	}



	struct bitmap {
		void ext_up() {
			bits_.resize(bits_.size() + w_);
			for (int i = bits_.size()-1; i >= w_; i--) {
				bits_[i] = bool(bits_[i - w_]);
			}
			for (int i = 0; i < w_; i++) {
				bits_[i] = false;
			}
			y_++;
			h_++;
		}
		void ext_down() {
			bits_.resize(bits_.size() + w_);
			h_++;
		}
		void ext_left() {
			bits_.resize(bits_.size() + h_);
			int at = (w_+1)*h_ -1;
			for (int i = h_-1; i >= 0 ; i--) {
				for (int j = w_-1; j >= 0 ; j--) {
					bits_[at] = bool(bits_[at-i+1]);
					at--;
				}
				bits_[at] = false;
			}
			w_++;
			x_++;
		}
		void ext_right() {
			bits_.resize(bits_.size() + h_);
			int at = (w_+1)*h_-1;
			for (int i = h_-1; i >= 0 ; i--) {
				bits_[at--] = false;
				for (int j = w_-1; j >= 0 ; j--) {
					bits_[at] = bool(bits_[at-i]);
					at--;
				}
			}
			w_++;
		}

		void set(int x, int y, bool v = true) {
			while (x < -x_) ext_left();
			while (y < -y_) ext_up();
			x += x_;
			y += y_;
			while (x >= w_) ext_right();
			while (y >= h_) ext_down();
			bits_[x+y*w_] = v;
		}
		bool get(int x, int y) const {
			return bits_[x+x_ + (y+y_) * w_];
		}
		bitmap() : w_(0), x_(), h_(0), y_(), bits_() {}
		int w_, x_;
		int h_, y_;
		reexp::bits bits_;
	};

	struct var_graphics {
		bitmap& b_;
		bitmap& b2_;
		int x_, y_;
		int xcvar_, ycvar_; // context variables that span x / y dimensions
		var_graphics(bitmap& b, bitmap& b2, int xcvar, int ycvar)
		: b_(b), b2_(b2), x_(), y_(), xcvar_(xcvar), ycvar_(ycvar) {}
		void set(int x, int y, bool b, bool b2) {
			b_.set(x+x_, y+y_, b);
			b2_.set(x+x_, y+y_, b2);
		}
		std::string to_string(char symbol) {
			std::ostringstream buf;
			for (int i = 0; i < b_.h_; ++i) {
				for (int j = 0; j < b_.w_; ++j) {
					if (b_.get(j, i)) {
						buf<<(b2_.get(j, i)?symbol:'-');
					} else {
						buf<<(".");
					}
				}
				buf<<std::endl;
			}
			return buf.str();
		}

	};

	struct translation_sentry {
		var_graphics& g_;
		int x_, y_;
		translation_sentry (var_graphics& g, int x, int y)
		:	g_(g), x_(x), y_(y) {
			g_.x_ += x;
			g_.y_ += y;
		}
		~translation_sentry () {
			g_.x_ -= x_;
			g_.y_ -= y_;
		}
	};

	/**
	 * information on the problem setting
	 */
	class pinfo {
		public:
			std::vector<std::string> vnames_; // bit variable names
			std::vector<std::string> cvnames_; // ctx variable names

			inline int vnameindex(const std::string& name) const {
				auto i = std::find(vnames_.begin(),
							  	   vnames_.end(),
							  	   name);
				if (i != vnames_.end()) {
					return i - vnames_.begin();
				} else {
					return -1;
				}
			}
	};

	template <typename P>
	class lang_info {
		public:
			lang_info(const pinfo& i, const reexp::lang<P>& l)
			: info_(i), lang_(l) {}

			std::string varid_tostring(size_t i) const {
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

			std::string shift_tostring(const cvec<P>& shift) const {
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

			std::string rel_tostring(const rel<P>& r, int state = 0xffff, const cvec<P>& shift = reexp::cvec<P>()) const {
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

			std::string rel_tostring(int rid, int state = 0xffff, const cvec<P>& shift = cvec<P>()) const {
				return rel_tostring(lang_.rel(rid), state, shift);
			}

			std::string var_tostring(const var<P>& v, const cvec<P>& shift = cvec<P>()) const {
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
						buf<<"v" + v.id();
						buf<<shift_tostring(shift);
					}
				}
				return buf.str();
			}
			std::string var_tostring(int i, const cvec<P>& shift = cvec<P>()) const {
				return var_tostring(lang_.var(i), shift);
			}
			std::string vars_tostring() const {
				std::ostringstream buf;
				for (int i = 0; i < lang_.var_count(); ++i) {
					buf<<var_tostring(i)<<std::endl;
				}
				return buf.str();
			}
			std::string rels_tostring() const {
				std::ostringstream buf;
				for (int i = 0; i < lang_.rel_count(); ++i) {
					buf<<rel_tostring(i)<<std::endl;
				}
				return buf.str();
			}

			static void draw_rel(var_graphics& g, const reexp::rel<P>& r, size_t state = 0xffff) {
				for (size_t i = 0; i < r.entries().size(); i++) {
					const reexp::rel_entry<P>& re = r.entries()[i];
					translation_sentry s(g, re.shift_[g.xcvar_], re.shift_[g.ycvar_]);
					draw_var(g, *re.var_, bool((state>>i)&1));
				}
			}

			static void draw_var(var_graphics& g, const reexp::var<P>& v, bool s) {
				const reexp::exp<P>* e = dynamic_cast<const reexp::exp<P>*>(&v);
				if (e) {
					draw_rel(g, e->rel(), s ? e->state() : 0);
				} else {
					g.set(0, 0, 1, s);
				}
			}

			std::string drawn_implmask_to_string(const implmask<P>& o) const {
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

			std::string drawn_var_rootmasks_to_string(const reexp::var<P>& var) const {
				std::ostringstream buf;
				for (const implmask<P>& o : var.rootmasks()) {
					buf<<drawn_implmask_to_string(o);
				}
				return buf.str();
			}

			std::string drawn_var_expmasks_to_string(const reexp::var<P>& var) const {
				std::ostringstream buf;
				for (const implmask<P>& o : var.expmasks()) {
					buf<<drawn_implmask_to_string(o);
				}
				return buf.str();
			}

			std::string drawn_var_implmasks_to_string(const reexp::var<P>& var) const {
				std::ostringstream buf;
				buf<<"rootmasks:\n";
				buf<<drawn_var_rootmasks_to_string(var);
				buf<<"\nexclmasks:\n";
				buf<<drawn_var_expmasks_to_string(var);
				return buf.str();
			}

			std::string drawn_var_tostring(int xcvar, int ycvar, const reexp::var<P>& v) const {
				bitmap bm;
				bitmap bm2;
				var_graphics g(bm, bm2, xcvar, ycvar);
				draw_var(g, v, true);
				return g.to_string(char_symbol_for(v.id()));
			}

			std::string drawn_rel_tostring(int xcvar, int ycvar, const reexp::rel<P>& r) const {
				bitmap bm;
				bitmap bm2;
				var_graphics g(bm, bm2, xcvar, ycvar);
				draw_rel(g, r);
				return g.to_string('X');
			}


			std::string drawn_vars_tostring(int xcvar, int ycvar) const {
				std::ostringstream buf;
				for (int i = 0; i < lang_.var_count(); ++i) {
					buf<<var_tostring(i)<<"\n\n";
					buf<<drawn_var_tostring(xcvar, ycvar, lang_.var(i))<<"\n\n";
				}
				return buf.str();
			}


			std::string drawn_rels_tostring(int xcvar, int ycvar) const {
				std::ostringstream buf;
				for (int i = 0; i < lang_.rel_count(); ++i) {
					buf<<rel_tostring(i)<<"\n\n";
					buf<<drawn_rel_tostring(xcvar, ycvar, lang_.rel(i))<<"\n\n";
				}
				return buf.str();
			}

			const reexp::lang<P>& lang() const {
				return lang_;
			}

		private:
			const pinfo& info_;
			const reexp::lang<P>& lang_;

	};

	template <typename P>
	class stats_info {
		public:
			stats_info(const pinfo& i, const reexp::stats<P>& s)
			: info_(i), lang_(i, s.data().lang()), stats_(s) {}

			const pinfo& info() const {
				return info_;
			}

			std::string drawn_data_var_tostring(int xvarid,
										  	    int yvarid,
										  	    const reexp::data_var<P>& pixels,
										  	   	const reexp::cvec<P>& at = reexp::cvec<P>()) {
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
							buf<<"?";
						}
					}
					buf<<"\n";
				}
				return buf.str();
			}

			std::string drawn_data_tostring(cvec<P> at, int xcvar, int ycvar) const {
				const reexp::data<P>& d = stats_.data();
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

			std::string drawn_data_tostring(cvec<P> at, int xcvar, int ycvar, int itercvar, int n) const {
				std::ostringstream buf;
				for (int i = 0; i < n; ++i) {
					buf<<i<<":\n";
					buf<<drawn_data_tostring(at, xcvar, ycvar)<<"\n";
					at[itercvar] = at[itercvar]+1;
				}
				return buf.str();
			}

			std::string var_tostring(const var_stats<P>& vs, const cvec<P>& shift = cvec<P>()) const {
				std::ostringstream buf;
				buf<<lang_.var_tostring(vs.data().var(), shift)<<"  "<<vs.freq()<<"/"<<vs.n()<<" ["<<vs.information()<<"]\n";
				return buf.str();
			}

			std::string var_tostring(int i, const cvec<P>& shift = cvec<P>()) const {
				return var_tostring(stats_.var(i), shift);
			}
			std::string vars_tostring() const {
				std::ostringstream buf;
				for (int i = 0; i < stats_.data().lang().var_count(); ++i) {
					buf<<var_tostring(i)<<std::endl;
				}
				return buf.str();
			}

			std::string detailed_vars_tostring(int xcvar, int ycvar) const {
				std::ostringstream buf;
				const lang<P>& lang = stats_.data().lang();
				for (int i = 0; i < lang.var_count(); ++i) {
					buf<<lang_.var_tostring(i)<<"\n\n";
					buf<<lang_.drawn_var_tostring(xcvar, ycvar, lang.var(i))<<"\n";
					buf<<var_deps_tostring(i)<<"\n\n";
				}
				return buf.str();
			}


			std::string rel_tostring(const rel_stats<P>& rs, int state = 0xffff, const cvec<P>& shift = cvec<P>()) const {
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
						double ne = rs.eStateP(i) * rs.eStateNaiveInfo(i);
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

			std::string rel_tostring(int i, int state = 0xffff, const cvec<P>& shift = cvec<P>()) const {
				return rel_tostring(stats_.rel(i), state, shift);
			}
			std::string rels_tostring() const {
				std::ostringstream buf;
				for (int i = 0; i < stats_.data().lang().rel_count(); ++i) {
					buf<<rel_tostring(i)<<std::endl;
				}
				return buf.str();
			}

			std::string scan_tostring(size_t max = 8, int details = 0) const {
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
						buf<<",naiveInfo(s)="<<rs.stateNaiveInfo(c.state_);
						buf<<",info(s)="<<rs.stateInfo(c.state_);
						buf<<",eNaiveInfo(s)="<<rs.eStateNaiveInfo(c.state_);
						buf<<",eInfo(s)="<<rs.eStateInfo(c.state_);
						buf<<"]";
					}
					buf<<"\n";
					cands.pop();
				}
				return buf.str();
			}

			void var_scan(const reexp::var<P>& var,
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

			double calc_influence(int from, const reexp::bits& to) {
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

			std::string top_influence_tostring(const reexp::bits& from,
											   const reexp::bits& to,
											   int cap) {
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

			std::string var_deps_tostring(const reexp::var<P>& var, int cap = 7) const {
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

			std::string var_deps_tostring(int var, int cap = 7) const {
				return var_deps_tostring(stats_.data().lang().var(var), cap);
			}
			std::string pred_tostring(const reexp::data<P>& d, int var, bool showbit = false) const {
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

			std::string pred_tostring(int var, bool showbit = false) const {
				return pred_tostring(stats_.data(), var, showbit);
			}

			std::string row_logdep_tostring(int var, int var2 = -1) const {
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
					pred.getInputStateLogDep(rs, riv, logdep, buf);
				}
				return buf.str();
			}

/*			std::string detailed_row_pred_tostring(int var) const {
				std::ostringstream buf;
				buf.setf(std::ios::fixed,std::ios::floatfield);
				buf.precision(2);
				explib::pred<P> pred(stats_);
				std::vector<double> p = pred.rowP(stats_.data(), var);

				const explib::data_var<P>& dv = stats_.data().var(var);
				explib::cvec<P> dim = dv.dim();
				for (size_t i = 0; i < p.size(); ++i) {
					buf<<"p("<<lang_.varid_tostring(var)<<"|"<<dim.at(i)<<") = "<<p[i]<<"%";
				}

				return buf.str();
			}*/

			std::string entropy_tostring(int var) const {
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

			std::string preds_tostring(int begin, int end, bool showbit = false) const {
				std::ostringstream buf;
				for (int i = begin; i < end; ++i) {
					buf<<pred_tostring(i, showbit);
				}
				return buf.str();
			}

			const reexp::lang_info<P>& lang_info() const {
				return lang_;
			}

			const reexp::stats<P>& stats() const {
				return stats_;
			}

			double var_naive_p(const std::vector<double>& origvarPs, int v) const {
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
			std::string drawn_vars_tostring_byexpbias(const std::vector<double>& origvarPs, int xcvar, int ycvar) const {
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



		private:

			const pinfo& info_;
			reexp::lang_info<P> lang_;
			const reexp::stats<P>& stats_;

	};

};


#endif /* REEXP_PRINTER_H_ */
