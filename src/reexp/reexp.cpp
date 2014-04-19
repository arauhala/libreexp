/*
 * reexp.cpp
 *
 *  Created on: Nov 23, 2013
 *      Author: arau
 */

#include "reexp/reexp.h"
#include "reexp/printer.h"
#include "reexp/pred.i.h"
#include "reexp/learner.i.h"
#include "reexp/io.i.h"
#include "reexp/printer.i.h"

namespace reexp {

	template struct bitmatrix<traits1d>;
	template struct bitmatrix<traits2d>;
	template struct bitmatrix<traits3d>;
	template struct bitmatrix<traits4d>;

	template struct ndim<traits1d>;
	template struct ndim<traits2d>;
	template struct ndim<traits3d>;
	template struct ndim<traits4d>;

	template class var<traits1d>;
	template class var<traits2d>;
	template class var<traits3d>;
	template class var<traits4d>;

	template class orig<traits1d>;
	template class orig<traits2d>;
	template class orig<traits3d>;
	template class orig<traits4d>;

	template class rel<traits1d>;
	template class rel<traits2d>;
	template class rel<traits3d>;
	template class rel<traits4d>;

	template class exp<traits1d>;
	template class exp<traits2d>;
	template class exp<traits3d>;
	template class exp<traits4d>;

	template bool operator < <traits1d>(const reexp::rel_entry<traits1d>& e1,
			 	 	 	 	 	 	 	const reexp::rel_entry<traits1d>& e2);
	template bool operator < <traits2d>(const reexp::rel_entry<traits2d>& e1,
			 	 	 	 	 	 	 	const reexp::rel_entry<traits2d>& e2);
	template bool operator < <traits3d>(const reexp::rel_entry<traits3d>& e1,
			 	 	 	 	 	 	 	const reexp::rel_entry<traits3d>& e2);
	template bool operator < <traits4d>(const reexp::rel_entry<traits4d>& e1,
			 	 	 	 	 	 	 	const reexp::rel_entry<traits4d>& e2);


	template int count_diff_in_datas<traits1d>(const reexp::data<traits1d>& d1,
											   const reexp::data<traits1d>& d2,
											   int& ddiff,
											   int& sdiff);
	template int count_diff_in_datas<traits2d>(const reexp::data<traits2d>& d1,
											   const reexp::data<traits2d>& d2,
											   int& ddiff,
											   int& sdiff);
	template int count_diff_in_datas<traits3d>(const reexp::data<traits3d>& d1,
											   const reexp::data<traits3d>& d2,
											   int& ddiff,
											   int& sdiff);
	template int count_diff_in_datas<traits4d>(const reexp::data<traits4d>& d1,
											   const reexp::data<traits4d>& d2,
											   int& ddiff,
											   int& sdiff);

	template class lang<traits1d>;
	template class lang<traits2d>;
	template class lang<traits3d>;
	template class lang<traits4d>;

	template class data<traits1d>;
	template class data<traits2d>;
	template class data<traits3d>;
	template class data<traits4d>;

	template struct rel_inputvars<traits1d>;
	template struct rel_inputvars<traits2d>;
	template struct rel_inputvars<traits3d>;
	template struct rel_inputvars<traits4d>;

	template struct rel_stats<traits1d>;
	template struct rel_stats<traits2d>;
	template struct rel_stats<traits3d>;
	template struct rel_stats<traits4d>;

	template struct exp_rel_stats<traits1d>;
	template struct exp_rel_stats<traits2d>;
	template struct exp_rel_stats<traits3d>;
	template struct exp_rel_stats<traits4d>;

	template class stats<traits1d>;
	template class stats<traits2d>;
	template class stats<traits3d>;
	template class stats<traits4d>;

	template void pred<traits1d>::getInputStateLogDep<std::ostream>(const rel_stats<traits1d>& rs,
											   	   	  	  	  	    const rel_inputvars<traits1d>& riv,
											   	   	  	  	  	    inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    std::ostream& explain) const;
	template void pred<traits2d>::getInputStateLogDep<std::ostream>(const rel_stats<traits2d>& rs,
											   	   	  	  	  	    const rel_inputvars<traits2d>& riv,
											   	   	  	  	  	    inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    std::ostream& explain) const;
	template void pred<traits3d>::getInputStateLogDep<std::ostream>(const rel_stats<traits3d>& rs,
											   	   	  	  	  	    const rel_inputvars<traits3d>& riv,
											   	   	  	  	  	    inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    std::ostream& explain) const;
	template void pred<traits4d>::getInputStateLogDep<std::ostream>(const rel_stats<traits4d>& rs,
											   	   	  	  	  	    const rel_inputvars<traits4d>& riv,
											   	   	  	  	  	    inputstate_logdep& inputStateLogDep,
											   	   	  	  	  	    std::ostream& explain) const;
	template class pred<traits1d>;
	template class pred<traits2d>;
	template class pred<traits3d>;
	template class pred<traits4d>;

	template class group_scaler<traits1d>;
	template class group_scaler<traits2d>;
	template class group_scaler<traits3d>;
	template class group_scaler<traits4d>;

	template class learner<traits1d>;
	template class learner<traits2d>;
	template class learner<traits3d>;
	template class learner<traits4d>;

	template class lang_info<traits1d>;
	template class lang_info<traits2d>;
	template class lang_info<traits3d>;
	template class lang_info<traits4d>;

	template class stats_info<traits1d>;
	template class stats_info<traits2d>;
	template class stats_info<traits3d>;
	template class stats_info<traits4d>;

	template class index_over_var_bits<traits1d>;
	template class index_over_var_bits<traits2d>;
	template class index_over_var_bits<traits3d>;
	template class index_over_var_bits<traits4d>;

	template class io<traits1d>;
	template class io<traits2d>;
	template class io<traits3d>;
	template class io<traits4d>;

}
