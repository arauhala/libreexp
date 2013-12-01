/*
 * bitmap.cpp
 *
 *  Created on: Dec 1, 2013
 *      Author: arau
 */

#include "reexp/printer.h"

namespace reexp {

	char char_symbol_for(int i) {
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

	void bitmap::ext_up() {
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
	void bitmap::ext_down() {
		bits_.resize(bits_.size() + w_);
		h_++;
	}
	void bitmap::ext_left() {
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
	void bitmap::ext_right() {
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

	void bitmap::set(int x, int y, bool v) {
		while (x < -x_) ext_left();
		while (y < -y_) ext_up();
		x += x_;
		y += y_;
		while (x >= w_) ext_right();
		while (y >= h_) ext_down();
		bits_[x+y*w_] = v;
	}
	bool bitmap::get(int x, int y) const {
		return bits_[x+x_ + (y+y_) * w_];
	}
	bitmap::bitmap() : w_(0), x_(), h_(0), y_(), bits_() {}

	var_graphics::var_graphics(bitmap& b, bitmap& b2, int xcvar, int ycvar)
	: b_(b), b2_(b2), x_(), y_(), xcvar_(xcvar), ycvar_(ycvar) {}
	void var_graphics::set(int x, int y, bool b, bool b2) {
		b_.set(x+x_, y+y_, b);
		b2_.set(x+x_, y+y_, b2);
	}
	std::string var_graphics::to_string(char symbol) {
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

	translation_sentry::translation_sentry(var_graphics& g, int x, int y)
	:	g_(g), x_(x), y_(y) {
		g_.x_ += x;
		g_.y_ += y;
	}
	translation_sentry::~translation_sentry () {
		g_.x_ -= x_;
		g_.y_ -= y_;
	}

	int pinfo::vnameindex(const std::string& name) const {
		auto i = std::find(vnames_.begin(),
						   vnames_.end(),
						   name);
		if (i != vnames_.end()) {
			return i - vnames_.begin();
		} else {
			return -1;
		}
	}
}
