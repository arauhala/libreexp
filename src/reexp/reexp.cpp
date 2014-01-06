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
#include "reexp/printer.i.h"

namespace reexp {

	template struct bitmatrix<traits1d>;
	template struct bitmatrix<traits2d>;
	template struct bitmatrix<traits3d>;
	template struct bitmatrix<traits4d>;

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

}
