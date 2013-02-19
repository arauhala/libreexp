/*
 * exptest_tools.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: arau
 */

#include "exptesttools.h"

namespace exptest {

	std::vector<std::string> read_split_line(std::istream& in) {
		std::string line;
		std::vector<std::string> rv;
		if (std::getline(in, line)) {
			size_t pos = 0;
			while (pos < line.size()) {
				if (line[pos]=='"') {
					size_t begin = ++pos;
					while (line[pos] && line[pos] != '"') {
						pos++;
					}
					if (line[pos] != '"') throw std::runtime_error("non-closed citation");
					rv.push_back(line.substr(begin, pos-begin));
					pos++;
				} else {
					size_t begin = pos;
					while (line[pos] && line[pos] != ' ') pos++;
					rv.push_back(line.substr(begin, pos-begin));
				}
				if (line[pos] != ' ') break;
				pos++;
			}
		}
		return rv;
	}

}
