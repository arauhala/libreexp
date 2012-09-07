/*
 * util.h
 *
 *  Created on: Nov 20, 2010
 *      Author: arauhala
 */

#ifndef UTIL_H_
#define UTIL_H_

#include "stdlib.h"
#include <memory>
#include <vector>

namespace explib {
	namespace util {

		int next_version();

		// As a difference to vector, the pointer to inserted item
		// is never changed. Also default constructor is not used internally.
		template <typename T>
		class arrays_list {
			public:
				arrays_list() : size_(), granularity_(100), arrays_() {}
				arrays_list(const arrays_list<T>& copy)
				: 	size_(),
				  	granularity_(copy.granularity_),
				  	arrays_() {
					for (int i = 0; i < copy.size(); i++) {
						push_back(copy[i]);
					}
				}
				~arrays_list() {
					for (size_t i = 0; i + 1 < arrays_.size(); i++) {
						for (size_t j = 0; j < granularity_; j++) {
							arrays_[i][j].~T();
						}
					}
					size_t n = size_ % granularity_;
					size_t k = arrays_.size()-1;
					for (size_t i = 0; i < n; i++) {
						arrays_[k][i].~T();
					}
					for (size_t i = 0; i < arrays_.size(); ++i) {
						free(arrays_[i]);
					}
				}
				void pop_back() {
					size_t a = array(--size_);
					size_t i = index(size_);
					arrays_[a][i].~T();
					if (i == 0 && arrays_.size()) { // last entry was deleted
						free(arrays_.back());
						arrays_.pop_back(); // free last array
					}
				}
				void push_back(const T& item) {
					size_t a = array(size_);
					while (a >= arrays_.size()) {
						arrays_.push_back( reinterpret_cast<T*>( malloc(sizeof(T) * granularity_) ) );
					}
					size_t i = index(size_++);
					new (&arrays_[a][i]) T(item);
				}
				T& operator[](int at) {
					return arrays_[array(at)][index(at)];
				}
				const T& operator[](int at) const {
					return arrays_[array(at)][index(at)];
				}
				T& back() {
					return operator[](size_-1);
				}
				const T& back() const {
					return operator[](size_-1);
				}
				inline int size() const {
					return size_;
				}
			protected:
				inline size_t array(size_t i) const {
					return i / granularity_;
				}
				inline size_t index(size_t i) const {
					return i % granularity_;
				}

			private:
				size_t size_;
				size_t granularity_;
				std::vector<T*> arrays_;
		};

		/**
		 * As a separation to traditional vector, this version never
		 * constructs anything with default constructor (with no parameters).
		 * Only copy constructors and destructors are used internally.
		 *
		 * NOTE: may not behave well on exceptions.
		 */
		template <typename T>
		class dense_vector {
			private:
				size_t size_;
				size_t allocSize_;
				T* array_;
			public:
				dense_vector() : size_(0), allocSize_(8), array_() {
					array_ = reinterpret_cast<T*>(malloc(allocSize_*sizeof(T)));
				}
				dense_vector(const dense_vector<T>& copy)
				: 	size_(0),
				  	allocSize_(copy.granularity_),
				  	array_() {
					array_ = malloc(allocSize_*sizeof(T));
					for (int i = 0; i < copy.size(); i++) {
						push_back(copy[i]);
					}
				}
				~dense_vector() {
					if (array_) {
						for (size_t i = 0; i < size_; i++) {
							array_[i].~T();
						}
						free(array_);
					}
				}
				void push_back(const T& item) {
					if (size_ == allocSize_) {
						allocSize_ *= 2;
						T* array = reinterpret_cast<T*>(malloc(allocSize_ * sizeof(T)));
						for (size_t i = 0; i < size_; i++) {
							new (&array[i]) T(array_[i]);
							array_[i].~T();
						}
						free(array_);
						array_ = array;
					}
					new (&array_[size_++]) T(item);
				}
				T& operator[](int at) {
					return array_[at];
				}
				const T& operator[](int at) const {
					return array_[at];
				}
				T& back() {
					return array_[size_-1];
				}
				const T& back() const {
					return array_[size_-1];
				}
				inline int size() const {
					return size_;
				}
		};

		struct hash_code {
			int hash_;
			hash_code() : hash_(5381) {}
			hash_code& operator << (int h) {
				hash_ = ((hash_ << 5) + hash_) + h;
				return *this;
			}
			int operator()() const {
				return hash_;
			}
		};

	}

}

#endif /* UTIL_H_ */
