#ifndef REEXP_BITS_H__
#define REEXP_BITS_H__

#include <stdlib.h>
#include <vector>
#include <stdint.h>
#include <string.h>
#include <sstream>
#include <stdexcept>
#include <assert.h>

namespace reexp {

	class bits;
	class cond_bits;

#if 0
	typedef unsigned long bchunk;
#else
	typedef unsigned long long bchunk;
#endif
	static const unsigned int bchunk_bsize = sizeof(bchunk)*8;
//	static const unsigned int bchunk_divmask = 0x1f;
	static const bchunk bchunk_mask = ~bchunk(0); // 64 bits

	struct sample {
		int class_;
		const char* text_;
	};


	template <class T>
	inline bool bit_of(int bit, T of) {
		return of&(1<<bit);
	}

	struct const_cond_bit {
		public:
			inline const_cond_bit(bool defined, bool state)
			: defined_ (defined),
			  state_(state) {}

			inline operator bool() {
				return defined_;
			}
			inline bool operator*() {
				return state_;
			}
		private:
			bool defined_ : 1;
			bool state_ : 1;
	};

	struct bit_ref {
		inline bit_ref(bits& b, int i) : bits_(b), i_(i) {}
		inline bit_ref& operator =(bool bit);
		inline bit_ref& operator |=(bool bit);
		inline bit_ref& operator &=(bool bit);
		inline operator bool () const;
		bits& bits_;
		int i_;
	};

	enum bit_tail_fill {
		false_tail_fill,
		true_tail_fill,
		no_tail_fill
	};

    struct bchunk_and_ostream {
		public:
			inline bchunk_and_ostream(bchunk* chunks, size_t i)
			: chunkp_(chunks + (i / bchunk_bsize)), offset_(i % bchunk_bsize) {}
			inline bchunk_and_ostream& operator<<(bchunk c) {
				c = ~c; // take this items' complement
				*(chunkp_++) &= ~(c << offset_);
				if (offset_) *chunkp_ &= ~(c >> (bchunk_bsize - offset_));
				//*chunkp_ &= ~(c >> ((bchunk_bsize-offset_)&bchunk_divmask));
				return *this;
			}
		private:
			bchunk* chunkp_;
			size_t offset_;
	};

    struct bchunk_or_ostream {
		public:
			inline bchunk_or_ostream(bchunk* chunks, int i)
			: chunkp_(chunks + (i / bchunk_bsize)), offset_(i % bchunk_bsize) {}
			inline bchunk_or_ostream& operator<<(bchunk c) {
				*(chunkp_++) |= (c << offset_);
				if (offset_) *chunkp_ |= (c >> (bchunk_bsize - offset_));
				//*chunkp_ |= (c >> ((bchunk_bsize-offset_)&bchunk_divmask));
				return *this;
			}
		private:
			bchunk* chunkp_;
			int offset_;
	};

    struct bchunk_istream {
		public:
			inline bchunk_istream(const bchunk* chunks, unsigned int i, unsigned int n, bit_tail_fill fill)
			: chunkp_(chunks + (i / bchunk_bsize)), offset_(i % bchunk_bsize), left_(n), fill_(fill) {}
			inline operator bool () const {
				return left_ > 0;
			}
			inline bchunk_istream& operator>>(bchunk& c) {
				c = *(chunkp_++) >> offset_;
				if (offset_) c |= *chunkp_ << (bchunk_bsize-offset_);
				if (left_ < int(bchunk_bsize)) {
					switch (fill_) { // let's hope fill_ get eliminated compile time, it should be possible
						case false_tail_fill:
							c &= ~(bchunk_mask<<left_); // fill non-existing tail with zeros
							break;
						case true_tail_fill:
							c |= (bchunk_mask<<left_); // fill non-existing tail with ones
							break;
						case no_tail_fill:  // nothing
							break;
					}
				}
				left_ -= bchunk_bsize;
				return *this;
			}
		private:
			const bchunk* chunkp_;
			unsigned int offset_;
			int left_;
			bit_tail_fill fill_;
	};

    template <int FILL>
    struct bchunk_istream2 {
		public:
			inline bchunk_istream2(const bchunk* chunks, unsigned int i, unsigned int n)
			: chunkp_(chunks + (i / bchunk_bsize)), offset_(i % bchunk_bsize), left_(n) {}
			inline operator bool () const {
				return left_ > 0;
			}
			inline bchunk_istream2& operator>>(bchunk& c) {
				c = *(chunkp_++) >> offset_;
				if (offset_) c |= *chunkp_ << (bchunk_bsize-offset_);
				if (left_ < int(bchunk_bsize)) {
					switch (FILL) {
						case false_tail_fill:
							c &= ~(bchunk_mask<<left_); // fill non-existing tail with zeros
							break;
						case true_tail_fill:
							c |= (bchunk_mask<<left_); // fill non-existing tail with ones
							break;
						case no_tail_fill:  // nothing
							break;
					}
				}
				left_ -= bchunk_bsize;
				return *this;
			}
		//private:
			const bchunk* chunkp_;
			unsigned int offset_;
			int left_;
	};

    struct bit_istream {
    private:
		const bchunk* chunkp_;
		size_t i_;
		size_t end_;
    public:
		bit_istream(const bchunk* chunkp, size_t offset, size_t length)
		: chunkp_(chunkp), i_(offset), end_(length) {}
		operator bool () const {
			return i_ != end_;
		}
		size_t pos() const {
			return i_;
		}
		inline bit_istream& operator>>(bool& b) {
			size_t chunk = i_ / bchunk_bsize;
			size_t bit = i_++ % bchunk_bsize;
			b = bool(chunkp_[chunk] & (bchunk(1) << bit));
			return *(this);
		}
    };
    struct bit_ostream {
    private:
		bchunk* chunkp_;
		size_t i_;
		size_t end_;
    public:
		bit_ostream(bchunk* chunkp, size_t bitoffset, size_t bitlength)
		: chunkp_(chunkp), i_(bitoffset), end_(bitlength) {}
		operator bool () const {
			return i_ != end_;
		}
		size_t pos() const {
			return i_;
		}
		inline bit_ostream& operator<<(bool b) {
			size_t chunk = i_ / bchunk_bsize;
			size_t bit = i_++ % bchunk_bsize;
			if (b) {
				chunkp_[chunk] |= (bchunk(1) << bit);
			} else {
				chunkp_[chunk] &= ~(bchunk(1) << bit);
			}
			return *(this);
		}
    };

	struct const_bits_ref {
	public:
		friend class bits;
		friend struct bits_ref;
		inline const_bits_ref(const bits& b, int i, int n)
		: bits_(b), i_(i), n_(n) {}
		inline bool operator[](int i) const;
		inline int size() const;
		inline bool true_overlap(const const_bits_ref& bits) const;
		inline bchunk_istream chunk_istream(bit_tail_fill fill) const;


	private:
		const bits& bits_;
		const int i_;
		const int n_;
	};

	struct bits_ref {
	public:
		friend class bits;
		inline bits_ref(bits& b, int i, int n) : bits_(b), i_(i), n_(n) {}
		inline bool operator[](int i) const;
		inline bit_ref operator[](int i);

		inline bits_ref& andNeg(const const_bits_ref& bits);
		inline bits_ref& operator&=(const const_bits_ref& bits);
		inline bits_ref& operator|=(const const_bits_ref& bits);
		inline bits_ref& fill(bool v);
/*		inline bits_ref& operator=(const const_bits_ref& bits);
		inline bits_ref& operator=(const bits_ref& bits);*/
		inline int size() const;
		inline operator const_bits_ref() const {
			return const_bits_ref(bits_, i_, n_);
		}
		inline bits_ref from(int from, int n) {
			return bits_ref(bits_, i_ + from, n_ - from - n);
		}
		inline bits_ref between(int from, int to) {
			return bits_ref(bits_, i_ + from, to - from);
		}
		inline const_bits_ref from(int from, int n) const {
			return bits_ref(bits_, i_ + from, n_ - from - n);
		}
		inline const_bits_ref between(int from, int to) const {
			return bits_ref(bits_, i_ + from, to - from);
		}
		inline bchunk_istream chunk_istream(bit_tail_fill fill) const;
		template <int FILL>
		inline bchunk_istream2<FILL> chunk_istream2() const;

	private:
		bits& bits_;
		const int i_;
		const int n_;
	};


	class bits {
		public:
			friend struct bits_ref;
			friend struct bit_ref;
			bits()
			:	size_(0), chunks_() {}
			bits(int size)
			:	size_(0), chunks_() {resize(size);}
			inline bool operator[](int i) const {
				assert(i >= 0 && i < size_);
				int chunk = i / bchunk_bsize;
				int bit = i % bchunk_bsize;
				bchunk b = chunks_[chunk];
				return b & (bchunk(1) << bit);
			}
			inline int popcount() const {
				int rv = 0;
				for (auto i = chunks_.begin(); i != chunks_.end(); ++i) {
					rv += __builtin_popcountl(*i);
				}
				return rv;
			}
			inline int sparse_popcount() const {
				int rv = 0;
				for (auto i = chunks_.begin(); i != chunks_.end(); ++i) {
					if (*i) rv += __builtin_popcountl(*i);
				}
				return rv;
			}
			inline bits& assignNeg(const bits& b) {
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] = ~b.chunks_[i];
				}
				return *this;
			}
			inline bits& andNeg(const bits& b) {
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] &= ~b.chunks_[i];
				}
				return *this;
			}
			inline bits& orNeg(const bits& b) {
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] |= ~b.chunks_[i];
				}
				return *this;
			}
			inline bits& operator^=(const bits& b) {
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] ^= b.chunks_[i];
				}
				return *this;
			}
			inline bits& operator&=(const bits& b) {
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] &= b.chunks_[i];
				}
				return *this;
			}
			inline bits& operator|=(const bits& b) {
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] |= b.chunks_[i];
				}
				return *this;
			}
			inline bits& copy(const bits& b) {
				assert(size() >= b.size());
				memcpy(chunks_.data(), b.chunks_.data(), sizeof(bchunk)*b.chunks_.size());
				return *this;
			}
			inline bits& operator=(const bits& b) {
				resize(b.size());
				return copy(b);
			}
			template <typename T>
			inline bits& copy_bits_ref(const T& b) {
				assert(size() >= b.size());
				// assert size() == b.size();
				if (b.i_ % bchunk_bsize) {
					const std::vector<bchunk>& data( b.bits_.chunks_ );
					size_t r = b.i_ / bchunk_bsize;
					size_t end = 1 + (b.i_ + b.n_ - 1) / bchunk_bsize;
					size_t boffset = b.i_ % bchunk_bsize;
					size_t w = 0;

					//
					// So,
					//   * for each write location, we have two read locations
					//   * one of the read locations may be empty
					//

					//
					bchunk chunk = data[r++];
					for (; w < chunks_.size(); r++) {
						chunks_[w] = chunk>>boffset;
						if (r == end) break;
						chunk = data[r];
						chunks_[w++] |= (chunk<<(bchunk_bsize-boffset));
					}
				} else {
					int offset = b.i_ / bchunk_bsize;
					int bytes = b.n_ / sizeof(bchunk);
					if (b.n_ % sizeof(bchunk)) bytes++;
					memcpy(chunks_.data(), b.bits_.chunks_.data() + offset, bytes);
				}
				clearTail();
				// this is needed so that popcount won't count extra bits
				return *this;
			}
			inline bits& copy(const const_bits_ref& b) {
				return copy_bits_ref(b);
			}
			inline bits& copy(const bits_ref& b) {
				return copy_bits_ref(b);
			}
			inline bits& operator=(const const_bits_ref& b) {
				resize(b.size());
				return copy(b);
			}
			inline bits& operator=(const bits_ref& b) {
				resize(b.size());
				return copy(b);
			}
			inline const_bits_ref from(int i, int n) const {
				return const_bits_ref(*this, i, n);
			}
			inline const_bits_ref from(int i) const {
				return const_bits_ref(*this, i, size()-i);
			}
			inline const_bits_ref between(int i, int j) const {
				return const_bits_ref(*this, i, j-i);
			}
			inline bits_ref from(int i, int n) {
				return bits_ref(*this, i, n);
			}
			inline bits_ref from(int i) {
				return bits_ref(*this, i, size()-i);
			}
			inline bits_ref between(int i, int j) {
				return bits_ref(*this, i, j-i);
			}
			inline bit_ref operator[](int i) {
				assert(i >= 0 && i < size_);
				return bit_ref(*this, i);
			}
			inline void fill(bool v) {
				bchunk f = v?bchunk_mask:0;
				for (size_t i = 0; i < chunks_.size(); ++i) {
					chunks_[i] = f;
				}
				clearTail();
			}
			inline void fill(bool v, int from) {
				bchunk f = v?bchunk_mask:0;
				int offset = from / bchunk_bsize;
				int boffset = from % bchunk_bsize;
				if (boffset) {
					if (v) {
						chunks_[offset] |= (bchunk_mask << boffset);
					} else {
						chunks_[offset] &= ~(bchunk_mask << boffset);
					}
					offset++;
				}
				for (size_t i = offset; i < chunks_.size(); ++i) {
					chunks_[i] = f;
				}
				clearTail();
			}
			inline void clearTail() {
				int lastchunkbit = size_%bchunk_bsize;
				if (lastchunkbit) {
					int lastchunk = size_/bchunk_bsize;
					chunks_[lastchunk] &= bchunk_mask >> (bchunk_bsize-lastchunkbit);
				}
			}
			inline void resize(int size) {
				if (size_ != size) {
					chunks_.resize((size/bchunk_bsize)+((size%bchunk_bsize)?1:0));
					this->size_ = size;
				}
			}

			inline void append(bool val) {
				chunks_.resize(size_+1);
				(*this)[size_++] = val;
			}
			inline int size() const {
				return size_;
			}
			inline std::vector<bchunk>& chunks() {
				return chunks_;
			}
			inline const std::vector<bchunk>& chunks() const {
				return chunks_;
			}
			inline bchunk chunk_from(int bitfrom) const {
				bchunk b;
				int firstidx = bitfrom / bchunk_bsize;
				int boffset = bitfrom % bchunk_bsize;
				b = chunks_[firstidx] >> boffset;
				if (boffset) b |= (chunks_[firstidx+1]
				                   << (bchunk_bsize - boffset));
				return b;
			}
			inline bit_istream istream(int offset = 0) {
				return bit_istream(chunks_.data(), offset, size_);
			}
			inline bit_ostream ostream(int offset = 0) {
				return bit_ostream(chunks_.data(), offset, size_);
			}
		private:
			int size_;
			std::vector<bchunk> chunks_;
	};

	bit_ref& bit_ref::operator=(bool value) {
		int chunk = i_ / bchunk_bsize;
		int bit = i_ % bchunk_bsize;
		if (value) {
			bits_.chunks_[chunk] |= bchunk(1) << bit;
		} else {
			bits_.chunks_[chunk] &= ~(bchunk(1) << bit);
		}
		return *this;
	}

	bit_ref& bit_ref::operator|=(bool value) {
		if (value) {
			int chunk = i_ / bchunk_bsize;
			int bit = i_ % bchunk_bsize;
			assert(chunk >= 0 && chunk < int(bits_.chunks_.size()));
			bits_.chunks_[chunk] |= bchunk(1) << bit;
		}
		return *this;
	}

	bit_ref& bit_ref::operator&=(bool value) {
		if (!value) {
			int chunk = i_ / bchunk_bsize;
			int bit = i_ % bchunk_bsize;
			assert(chunk >= 0 && chunk < int(bits_.chunks_.size()));
			bits_.chunks_[chunk] &= ~(bchunk(1) << bit);
		}
		return *this;
	}


	bit_ref::operator bool() const {
		int chunk = i_ / bchunk_bsize;
		int bit = i_ % bchunk_bsize;
		return bits_.chunks_[chunk] & (bchunk(1) << bit);
	}

	inline bool const_bits_ref::operator[](int i) const {
		return bits_[i_+i];
	}
	inline int const_bits_ref::size() const {
		return n_;
	}
	inline bool const_bits_ref::true_overlap(const const_bits_ref& bits) const {
		bchunk_istream tin( chunk_istream( false_tail_fill) );
		bchunk_istream oin( bits.chunk_istream( false_tail_fill) );
		while (oin) {
			bchunk tc, oc;
			tin>>tc;
			oin>>oc;
			if (tc & oc) return true;
		}
		return false;
	}

	inline bchunk_istream const_bits_ref::chunk_istream(bit_tail_fill fill) const {
		return bchunk_istream(bits_.chunks().data(), i_, n_, fill);
	}
	inline bool bits_ref::operator[](int i) const {
		return bits_[i_+i];
	}
	inline bit_ref bits_ref::operator[](int i) {
		return bits_[i_+i];
	}
	inline int bits_ref::size() const {
		return n_;
	}
	// todo: rename to setAndNeg
	inline bits_ref& bits_ref::andNeg(const const_bits_ref& bits) {
		bchunk_istream in( bits.chunk_istream( false_tail_fill ) );
		bchunk_and_ostream out( bits_.chunks().data(), i_ );
		while (in) {
			bchunk c;
			in>>c; c = ~c; out<<c;
		}
		return *this;
	}
	inline bits_ref& bits_ref::operator&=(const const_bits_ref& bits) {
		bchunk_istream in( bits.chunk_istream( true_tail_fill) );
		bchunk_and_ostream out( bits_.chunks().data(), i_ );
		while (in) {
			bchunk c;
			in>>c; out<<c;
		}
		return *this;
	}

	inline bits_ref& bits_ref::operator|=(const const_bits_ref& bits) {
		bchunk_istream in( bits.chunk_istream( false_tail_fill) );
		bchunk_or_ostream out( bits_.chunks().data(), i_ );
		while (in) {
			bchunk c;
			in>>c; out<<c;
		}
		return *this;
	}
	inline bits_ref& bits_ref::fill(bool v) {
		if (v) {
			bchunk_or_ostream out( bits_.chunks().data(), i_ );
			int left = n_;
			while (left >= int(bchunk_bsize)) { out<<bchunk_mask; left -= bchunk_bsize; }
			if (left) out<<(bchunk_mask>>(bchunk_bsize-left));
		} else {
			bchunk_and_ostream out( bits_.chunks().data(), i_ );
			int left = n_;
			while (left >= int(bchunk_bsize)) { out<<0; left -= bchunk_bsize; }
			out<<(bchunk_mask<<(left));
		}
		return *this;
	}


	inline bchunk_istream bits_ref::chunk_istream(bit_tail_fill fill) const {
		return bchunk_istream(bits_.chunks().data(), i_, n_, fill);
	}

	template <int FILL>
	inline bchunk_istream2<FILL> bits_ref::chunk_istream2() const {
		return bchunk_istream2<FILL>(bits_.chunks().data(), i_, n_);
	}

	class cond_bits;

	struct cond_bit_ref {
		inline cond_bit_ref(cond_bits& b, int i) : bits_(b), i_(i) {}
		inline operator bit_ref();
		inline operator bool() { return operator bit_ref(); }
		inline cond_bit_ref& operator =(bool b) { operator bit_ref().operator=(b); return *this; }
		inline bool operator * () const;
		inline bit_ref operator * ();
		cond_bits& bits_;
		int i_;
	};

	struct const_cond_bit_ref {
		inline const_cond_bit_ref(const cond_bits& b, int i) : bits_(b), i_(i) {}
		inline operator bool();
		inline bool operator * () const;
		const cond_bits& bits_;
		int i_;
	};

	struct cond_bits_ref {
		const cond_bits& bits_;
		int offset_;
		int n_;
		cond_bits_ref(const cond_bits& bits, int offset, int n)
		:   bits_(bits),
		    offset_(offset),
		    n_(n) {}
		inline operator const_bits_ref() const;
		inline const_bits_ref operator*() const;
	};

	class cond_bits {
			friend struct cond_bits_ref;
		public:
			friend struct cond_bit_ref;
			friend struct const_cond_bit_ref;
			inline cond_bit_ref operator[](int i) {
				return cond_bit_ref(*this, i);
			}
			inline const_cond_bit_ref operator[](int i) const {
				return const_cond_bit_ref(*this, i);
			}
			inline cond_bits_ref from(int offset, int n) const {
				return cond_bits_ref(*this, offset, n);
			}
			inline cond_bits& operator=(const cond_bits_ref& b) {
				defined_ = b.bits_.defined_.from(b.offset_, b.n_);
				states_ = b.bits_.states_.from(b.offset_, b.n_);
				return *this;
			}
			inline void resize(int size) {
				defined_.resize(size);
				states_.resize(size);
			}
			inline void append(bool val) {
				defined_.append(true);
				states_.append(val);
			}
			inline int size() const {
				return defined_.size();
			}
			inline bits& defined() {
				return defined_;
			}
			inline bits& states() {
				return states_;
			}
			inline const bits& defined() const {
				return defined_;
			}
			inline const bits& states() const {
				return states_;
			}
		private:
			bits defined_;
			bits states_;
	};

	cond_bit_ref::operator bit_ref () {
		return bits_.defined_[i_];
	}

	bit_ref cond_bit_ref::operator*() {
		return bits_.states_[i_];
	}

	bool cond_bit_ref::operator*() const {
		return bits_.states_[i_];
	}

	inline const_cond_bit_ref::operator bool() {
		return bits_.defined_[i_];
	}

	inline bool const_cond_bit_ref::operator * () const {
		return bits_.states_[i_];
	}

	cond_bits_ref::operator const_bits_ref() const {
		return bits_.defined_.from(offset_, n_);
	}
	const_bits_ref cond_bits_ref::operator*() const {
		return bits_.states_.from(offset_, n_);
	}


	template <typename T>
	std::string vector_todensestring(const T& a) {
		std::ostringstream buf;
		for (int i = 0; i < a.size(); ++i) {
			buf<<a[i];
		}
		return buf.str();
	}
	template <typename T>
	std::string vector_tolines(const T& a, int split) {
		std::ostringstream buf;
		for (int i = 0; i < a.size(); ++i) {
			if (i && !(i%split)) buf<<"\n";
			buf<<a[i];
		}
		return buf.str();
	}

}

#endif /*REEXP_BITS_H__*/
