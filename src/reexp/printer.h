/*
 * print.h
 *
 *  Created on: Aug 25, 2011
 *      Author: arau
 */

#ifndef REEXP_PRINTER_H_
#define REEXP_PRINTER_H_

#include "reexp/reexp.h"
#include "reexp/io.h"

#include <iostream>
#include <sstream>
#include <string>

namespace reexp {

	struct bitmap {
		void ext_up();
		void ext_down();
		void ext_left();
		void ext_right();
		void set(int x, int y, bool v = true);
		bool get(int x, int y) const;
		bitmap();
		int w_, x_;
		int h_, y_;
		reexp::bits bits_;
	};

	struct var_graphics {
		bitmap& b_;
		bitmap& b2_;
		int x_, y_;
		int xcvar_, ycvar_; // context variables that span x / y dimensions
		var_graphics(bitmap& b, bitmap& b2, int xcvar, int ycvar);
		void set(int x, int y, bool b, bool b2);
		std::string to_string(char symbol);
	};

	struct translation_sentry {
		var_graphics& g_;
		int x_, y_;
		translation_sentry(var_graphics& g, int x, int y);
		~translation_sentry();
	};

	/**
	 * information on the problem setting
	 * TODO: should be renamed to names or labels
	 */
	class pinfo {
		public:
			std::vector<std::string> vnames_; // bit variable names
			std::vector<std::string> cvnames_; // ctx variable names

			int vnameindex(const std::string& name) const;
	};

	/**
	 * lang info can be used to print information about a lang<P> class
	 * and draw amazing ASCII visualizations of its data (e.g. various bit masks)
	 */
	template <typename P>
	class lang_info {
		public:
			lang_info(const pinfo& i, const reexp::lang<P>& l);
			std::string varid_tostring(size_t i) const;
			std::string shift_tostring(const cvec<P>& shift) const;
			std::string rel_tostring(const rel<P>& r, int state = 0xffff, const cvec<P>& shift = reexp::cvec<P>()) const;
			std::string rel_tostring(int rid, int state = 0xffff, const cvec<P>& shift = cvec<P>()) const;
			std::string var_tostring(const var<P>& v, const cvec<P>& shift = cvec<P>()) const;
			std::string var_tostring(int i, const cvec<P>& shift = cvec<P>()) const;
			std::string vars_tostring() const;
			std::string rels_tostring() const;
			static void draw_rel(var_graphics& g, const reexp::rel<P>& r, size_t state = 0xffff);
			static void draw_var(var_graphics& g, const reexp::var<P>& v, bool s);
			std::string drawn_implmask_to_string(const implmask<P>& o) const;
			std::string drawn_var_rootmasks_to_string(const reexp::var<P>& var) const;
			std::string drawn_var_expmasks_to_string(const reexp::var<P>& var) const;
			std::string drawn_var_implmasks_to_string(const reexp::var<P>& var) const;
			std::string drawn_var_tostring(int xcvar, int ycvar, const reexp::var<P>& v) const;
			std::string drawn_rel_tostring(int xcvar, int ycvar, const reexp::rel<P>& r) const;
			std::string drawn_vars_tostring(int xcvar, int ycvar) const;
			std::string drawn_rels_tostring(int xcvar, int ycvar) const;
			std::string invorder_diff_tostring(const reexp::data<P>& d1,
											   const reexp::data<P>& d2,
											   int maxprints = 100,
											   const index_over_var_bits<P>* d2KnownOffsets = 0,
											   const std::vector<int>* d2KnownAt = 0) const;
			const reexp::lang<P>& lang() const;
		private:
			const pinfo& info_;
			const reexp::lang<P>& lang_;

	};

	/**
	 * lang info can be used to print information about a stats<P> class and
	 * also some info of the classes it refers (data<P>, lang<P>)
	 */
	template <typename P>
	class stats_info {
		public:
			stats_info(const pinfo& i, const reexp::stats<P>& s);
			const pinfo& info() const;
			std::string drawn_data_var_tostring(int xvarid,
										  	    int yvarid,
										  	    const reexp::data_var<P>& pixels,
										  	   	const reexp::cvec<P>& at = reexp::cvec<P>()) const;
			std::string drawn_data_tostring(cvec<P> at, int xcvar, int ycvar, const data<P>& d) const;
			std::string drawn_data_tostring(cvec<P> at, int xcvar, int ycvar) const;
			std::string drawn_data_tostring(cvec<P> at, int xcvar, int ycvar, int itercvar, int n) const;
			std::string var_tostring(const var_stats<P>& vs, const cvec<P>& shift = cvec<P>()) const;
			std::string var_tostring(int i, const cvec<P>& shift = cvec<P>()) const;
			std::string vars_tostring() const;
			std::string detailed_vars_tostring(int xcvar, int ycvar) const;
			std::string rel_tostring(const rel_stats<P>& rs, int state = 0xffff, const cvec<P>& shift = cvec<P>()) const;
			std::string rel_tostring(int i, int state = 0xffff, const cvec<P>& shift = cvec<P>()) const;
			std::string rels_tostring() const;
			std::string scan_tostring(size_t max = 8, int details = 0) const;
			void var_scan(const reexp::var<P>& var,
						  std::priority_queue<candidate<P>>& cands) const;
			double calc_influence(int from, const reexp::bits& to) const;
			std::string top_influence_tostring(const reexp::bits& from,
											   const reexp::bits& to,
											   int cap) const;
			std::string var_deps_tostring(const reexp::var<P>& var, int cap = 7) const;
			std::string var_deps_tostring(int var, int cap = 7) const;
			std::string pred_tostring(const reexp::data<P>& d, int var, bool showbit = false) const;
			std::string pred_tostring(int var, bool showbit = false) const;
			std::string row_logdep_tostring(int var, int var2 = -1) const;
			std::string entropy_tostring(int var) const;
			std::string preds_tostring(int begin, int end, bool showbit = false) const;
			const reexp::lang_info<P>& lang_info() const;
			const reexp::stats<P>& stats() const;
			double var_naive_p(const std::vector<double>& origvarPs, int v) const;
			std::string drawn_vars_tostring_byexpbias(const std::vector<double>& origvarPs, int xcvar, int ycvar) const;
		private:
			const pinfo& info_;
			reexp::lang_info<P> lang_;
			const reexp::stats<P>& stats_;

	};


	extern template class lang_info<traits1d>;
	extern template class lang_info<traits2d>;
	extern template class lang_info<traits3d>;
	extern template class lang_info<traits4d>;

	extern template class stats_info<traits1d>;
	extern template class stats_info<traits2d>;
	extern template class stats_info<traits3d>;
	extern template class stats_info<traits4d>;

};


#endif /* REEXP_PRINTER_H_ */
