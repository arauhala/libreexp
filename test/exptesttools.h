/*
 * exptesttools.h
 *
 *  Created on: Aug 23, 2011
 *      Author: arau
 */

#ifndef EXPTESTTOOLS_H_
#define EXPTESTTOOLS_H_

#include "tester.h"
#include "reexp/all.h"

namespace exptest {

	template <typename P>
	void print_pixels(TestTool& t,
					  const explib::data_var<P>& pixels,
					  int xvarid,
					  int yvarid,
					  const explib::cvec<P>& at = explib::cvec<P>()) {
		explib::cvec<P> d = pixels.dim();
		explib::cvec<P> a( at );
		for (int i = at[yvarid]; i < d[yvarid]; i++) {
			a[yvarid] = i;
			for (int j = at[xvarid]; j < d[xvarid]; j++) {
				a[xvarid] = j;
				auto b = pixels[a];
				if (b) {
					if (*b) {
						t<<"X";
					} else {
						t<<".";
					}
				} else {
					t<<"?";
				}
			}
			t<<checkl;
		}
	}

	std::vector<std::string> read_split_line(std::istream& in);
}

#endif /* EXPTESTTOOLS_H_ */
