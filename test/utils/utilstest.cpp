/*
 * utilstest.cpp
 *
 *  Created on: Sep 29, 2011
 *      Author: arau
 */

#include "reexp/all.h"

#include "tester.h"

void run_bit_ref_test(TestTool& t) {
	explib::bits b;
	b.resize(8);
	for (int i = 4; i < 8; ++i) {
		b[i] = true;
	}

	t<<"bits: "<<vector_todensestring(b)<<"\n";
	t<<"ref [0, 4]: "<<vector_todensestring(b.between(0, 4))<<"\n";
	t<<"ref [2, 6]: "<<vector_todensestring(b.between(2, 6))<<"\n";
	t<<"ref [4, 8]: "<<vector_todensestring(b.between(4, 8))<<"\n";
}

void copy_and_print(TestTool& t, explib::bits& dest, explib::bits& src, int offset) {
	dest = src.from(offset, dest.size());
	t<<"bits ["<<offset<<", "<<offset+dest.size()<<"]\n"<<vector_tolines(dest, 32)<<"\n";
}

void run_bit_copy_test(TestTool& t) {
	explib::bits b;
	b.resize(128);
	b.fill(false);
	for (int i = 32; i < b.size()-32; i+=2) {
		b[i] = true;
	}
	t<<"src:\n"<<vector_tolines(b, 32)<<"\n";
	explib::bits b2;
	b2.resize(64);
	copy_and_print(t, b2, b, 0);
	copy_and_print(t, b2, b, 3);
	copy_and_print(t, b2, b, 16);
	copy_and_print(t, b2, b, 32);
	copy_and_print(t, b2, b, 48);
	copy_and_print(t, b2, b, 61);
	copy_and_print(t, b2, b, 64);
}
void run_bits_bitop_test(TestTool& t) {
	explib::bits b1;
	explib::bits b2;
	explib::bits b3;
	b1.resize(256);
	b2.resize(256);
	b3.resize(256);

	for (int i = 1; i < b1.size(); i+=2) {
		b1[i] = true;
	}
	b2.copy(b1.from(1));

	t<<"b1, pop "<<b1.popcount()<<":\n"<<vector_tolines(b1, 64)<<"\n\n";
	t<<"b2, pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	b3 = b1;
	b3 &= b2;
	t<<"b1 & b2, pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3 |= b2;
	t<<"b1 | b2, pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.andNeg(b2);
	t<<"b1 & ~b2, pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.andNeg(b1);
	t<<"b1 & ~b2, pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";
}

void run_bits_ref_andneg_test(TestTool& t) {
	explib::bits b1;
	explib::bits b2;
	explib::bits b3;
	b1.resize(256);
	b2.resize(256);
	b3.resize(256);

	for (int i = 1; i < b1.size(); i+=2) {
		b1[i] = true;
	}
	b2.copy(b1.from(1));

	t<<"b1, pop "<<b1.popcount()<<":\n"<<vector_tolines(b1, 64)<<"\n\n";
	t<<"b2, pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	explib::bits b4;
	b4.resize(256);
	b4.fill(true);

	b3 = b1;
	b3.from(32, 128).andNeg(b4.from(96, 128));
	t<<"b1[32,160] & ~b4[92,224], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(20, 140).andNeg(b4.from(80, 140));
	t<<"b1[20, 160] & ~b4[80, 220], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";


	b3 = b1;
	b3.from(32, 64).andNeg(b4.from(92, 64));
	t<<"b1[32,98] & ~b4[92,156], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(10, 10).andNeg(b2.from(57, 10));
	t<<"b1[10,20] & ~b2[57,67], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(60, 10).andNeg(b2.from(57, 10));
	t<<"b1[60,10] & ~b2[57,67], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.between(0, 16).andNeg(b2.between(0, 16));
	t<<"b1[0,16] & ~b2[0,16], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.between(0, 16).andNeg(b2.between(1, 16));
	t<<"b1[0,16] & ~b2[1,16], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(10, 10).andNeg(b2.from(43, 10));
	t<<"b1[10,20] & ~b2[43,53], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(9, 10).andNeg(b2.from(43, 10));
	t<<"b1[9,19] & ~b2[43,53], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";
}

void run_bits_andneg_test(TestTool& t) {
	explib::bits b1;
	explib::bits b2;
	explib::bits b3;
	explib::bits b4;
	b1.resize(256);
	b2.resize(256);
	b3.resize(256);
	b4.resize(256);

	for (int i = 0; i < b1.size(); i+=2) {
		b1[i] = true;
	}
	for (int i = 1; i < b2.size(); i+=3) {
		b2[i] = true;
	}
	t<<"b1:\n"<<vector_tolines(b1, 64)<<"\n\n";
	t<<"b2:\n"<<vector_tolines(b2, 64)<<"\n\n";

	b3 = b1;
	b3.andNeg(b2);

	t<<"b1 & ~b2:\n"<<vector_tolines(b3, 64)<<"\n\n";

	b4 = b2;
	b4 &= b1;
	t<<"b2 & b1:\n"<<vector_tolines(b4, 64)<<"\n\n";

	b3 = b1;
	b3.andNeg(b4);
	t<<"b1 & ~(b2 & ~b1):\n"<<vector_tolines(b3, 64)<<"\n\n";

}


void run_bits_ref_and_test(TestTool& t) {
	explib::bits b1;
	explib::bits b2;
	explib::bits b3;
	b1.resize(256);
	b2.resize(256);
	b3.resize(256);

	for (int i = 1; i < b1.size(); i+=2) {
		b1[i] = true;
	}
	b2.copy(b1.from(1));

	t<<"b1, pop "<<b1.popcount()<<":\n"<<vector_tolines(b1, 64)<<"\n\n";
	t<<"b2, pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	explib::bits b4;
	b4.resize(256);
	b4.fill(0);

	b3 = b1;
	b3.from(32, 128) &= (b4.from(96, 128));
	t<<"b1[32,160] & b4[92,224], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(32, 64) &= (b4.from(92, 64));
	t<<"b1[32,98] & b4[92,156], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(10, 10) &= (b2.from(56, 10));
	t<<"b1[10,20] & b2[56,66], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(60, 10) &= (b2.from(56, 10));
	t<<"b1[60,10] & b2[56,66], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";


	b3 = b1;
	b3.between(0, 16) &= (b2.between(0, 16));
	t<<"b1[0,16] & b2[0,16], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.between(0, 16) &= (b2.between(1, 16));
	t<<"b1[0,16] & b2[1,16], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(10, 10) &= (b2.from(44, 10));
	t<<"b1[10,20] & b2[44,54], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(9, 10) &= (b2.from(44, 10));
	t<<"b1[9,19] & b2[44,54], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(32, 64) &= (b2.from(92, 64));
	t<<"b1[32,98] & b2[92,156], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";


}

void run_bits_ref_or_test(TestTool& t) {
	explib::bits b1;
	explib::bits b2;
	explib::bits b3;
	b1.resize(256);
	b2.resize(256);
	b3.resize(256);

	for (int i = 1; i < b1.size(); i+=2) {
		b1[i] = true;
	}
	b2.copy(b1.from(1));

	t<<"b1, pop "<<b1.popcount()<<":\n"<<vector_tolines(b1, 64)<<"\n\n";
	t<<"b2, pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	explib::bits b4;
	b4.resize(256);
	b4.fill(1);

	b3 = b1;
	b3.from(32, 128) |= (b4.from(96, 128));
	t<<"b1[32,160] | b4[92,224], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(32, 64) |= (b4.from(92, 64));
	t<<"b1[32,98] | b4[92,156], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(10, 10) |= (b2.from(56, 10));
	t<<"b1[10,20] | b2[56,66], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(60, 10) |= (b2.from(56, 10));
	t<<"b1[60,10] | b2[56,66], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";


	b3 = b1;
	b3.between(0, 16) |= (b2.between(0, 16));
	t<<"b1[0,16] | b2[0,16], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.between(0, 16) |= (b2.between(1, 16));
	t<<"b1[0,16] | b2[1,16], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(10, 10) |= (b2.from(44, 10));
	t<<"b1[10,20] | b2[44,54], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(9, 10) |= (b2.from(44, 10));
	t<<"b1[9,19] | b2[44,54], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";

	b3 = b1;
	b3.from(32, 64) |= (b2.from(92, 64));
	t<<"b1[32,98] | b2[92,156], pop "<<b3.popcount()<<":\n"<<vector_tolines(b3, 64)<<"\n\n";
}

void run_bits_ref_fill_test(TestTool& t) {
	explib::bits b1;
	explib::bits b;
	b1.resize(256);
	b.resize(256);

	for (int i = 1; i < b1.size(); i+=2) {
		b1[i] = true;
	}

	b = b1;
	b.from(32, 128).fill(true);
	t<<"b.from(32,128).fill(true), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(32, 128).fill(false);
	t<<"b.from(32,128).fill(false), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(32, 64).fill(true);
	t<<"b.from(32,64).fill(true), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(32, 64).fill(false);
	t<<"b.from(32,64).fill(false), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(10, 10).fill(true);
	t<<"b.from(10,10).fill(true), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(10, 10).fill(false);
	t<<"b.from(10,10).fill(false), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(0, 16).fill(true);
	t<<"b.from(0,16).fill(true), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(0,16).fill(false);
	t<<"b.from(0,16).fill(false), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(1, 250).fill(true);
	t<<"b.from(1, 250).fill(true), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

	b = b1;
	b.from(1, 250).fill(false);
	t<<"b.from(1, 250).fill(false), pop "<<b.popcount()<<":\n"<<vector_tolines(b, 64)<<"\n\n";

}

void run_bits_ref_perf_test(TestTool& t) {
	explib::bits b1;
	b1.resize(8*1024*1024);

	explib::bits b2;
	b2.resize(8*1024*1024);

	srand(0);
	for (int i = 0; i < b1.size(); i++) {
		b1[i] = rand()%2;
	}

	for (int i = 0; i < b2.size(); i++) {
		b2[i] = rand()%2;
	}
	int size = 1;
	char buf[24];
	for (int i = 0; i < 14; ++i) {
		int times = (1<<27)>>i;

		snprintf(buf, sizeof(buf),  "size:%d", size);

		int max1 = b1.size() - size;
		int max2 = b2.size() - size;

		{
			TimeSentry time;

			for (int j = 0; j < times; ++j) {
				int destpos = ((j<<3) + j)%max1;
				int srcpos = ((j<<3) + j)%max2;
				b1.from(destpos, size).andNeg(b2.from(srcpos, size));
			}
			double us = double(time.us())/times;
			t.record({buf, "property:&~_time(us)"}, double(us));
		}
		{
			TimeSentry time;
			for (int j = 0; j < times; ++j) {
				int destpos = ((j<<3) + j)%max1;
				int srcpos = ((j<<3) + j)%max2;
				b1.from(destpos, size) &= (b2.from(srcpos, size));
			}
			double us = double(time.us())/times;
			t.record({buf, "property:&_time(us)"}, double(us));
		}
		{
			TimeSentry time;
			for (int j = 0; j < times; ++j) {
				int destpos = ((j<<3) + j)%max1;
				int srcpos = ((j<<3) + j)%max2;
				b1.from(destpos, size) |= (b2.from(srcpos, size));
			}
			double us = double(time.us())/times;
			t.record({buf, "property:|time(us)"}, double(us));
		}
		size *= 2;
	}
	Table table(
		t.report(ToTable<Average>({}, "property:", "size:")));
	t.reported()<<table;

}

void run_popcount_test(TestTool& t) {
	explib::bits b1;
	explib::bits b2;
	b1.resize(256);
	b2.resize(7);

	for (int i = 1; i < b1.size(); i+=2) {
		b1[i] = true;
	}

	t<<"bits, pop "<<b1.popcount()<<":\n"<<vector_tolines(b1, 64)<<"\n\n";

	b2 = b1.from(0, 7);
	t<<"bits [0,7], pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	b2 = b1.from(30, 7);
	t<<"bits [30,37], pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	b2 = b1.from(60, 7);
	t<<"bits [60,67], pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	b2 = b1.from(64, 7);
	t<<"bits [64,71], pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	b2 = b1.from(122, 7);
	t<<"bits [122,129], pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";

	b2 = b1.from(129, 7);
	t<<"bits [129,136], pop "<<b2.popcount()<<":\n"<<vector_tolines(b2, 64)<<"\n\n";


}

void run_big_popcount_test(TestTool& t) {
	explib::bits b1;
	b1.resize(8*1024*1024);

	int ln = 13;
	for (int i = 0; i < b1.size(); i++) {
		ln = (ln * 1953) % 179;
		b1[i] = ln%2;
	}

	t<<"popcount: "<<b1.popcount()<<"\n\n";
}

void run_big_16b_popcount_test(TestTool& t) {
	static const int units = 16*1024;
	explib::bits b[units];

	int ln = 13;
	for (int i = 0; i < units; i++) {
		b[i].resize(16);
		for (int j = 0; j < 16; j++) {
			ln = (ln * 1953) % 179;
			b[i][j] = ln%2;
		}
	}

	int pop = 0;
	for (int i = 0; i < units; i++) {
		pop += b[i].popcount();
	}

	t<<"popcount: "<<pop<<"\n\n";
}

struct problem2d {
	static const int DIM = 2; // two context variables
	static const int MAX_REL_VARS = 2; // two max relation variables
};

struct problem3d {
	static const int DIM = 3; // two context variables
	static const int MAX_REL_VARS = 2; // two max relation variables
};

template <typename P>
std::string bitmatrix_tostring(const explib::bitmatrix<P>& m) {
	std::ostringstream buf;
	const explib::cvec<P>& dim = m.ndim_.dim_;

	for (explib::dim_row_iterator<P> i(dim, dim.rowdim()); i; ++i) {
		int offset = dim.offset(i.begin());
		for (int j = 0; j < i.length(); ++j) {
			if (m.bits_[offset+j]) {
				buf<<"X";
			} else {
				buf<<".";
			}
		}
		buf<<"\n";
	}
	return buf.str();
}

void run_ndim_test(TestTool& t) {
	typedef problem2d p;

	explib::ndim<p> v;

	auto fit= [&](const explib::cvec<p>& shift,
					const explib::cvec<p>& dim) {
		t<<"applying "<<dim<<" at "<<shift<<"\n";
		v.fit(shift, dim);
		t<<"shift: "<<v.shift_<<"  ndim: "<<v.dim_<<"  vol:"<<v.dim_.volume()<<"\n\n";
	};

	fit({0, 0}, {1, 1});
	fit({-1, -1}, {1, 1});
	fit({-1, -1}, {3, 1});
	fit({-1, -1}, {1, 3});
}

template <typename P>
void run_dim_bchunk_iteration_test_round(TestTool& t, int rowlen) {
	typedef P p;
	int size = 8*1024*1024;// 1 megabyte
	explib::cvec<p> dim;
	dim[0] = rowlen;
	if (P::DIM ==2 ) {
		int rows = size / rowlen;
		dim[1] = rows;
	} else {
		int rows = size / rowlen;
		rows /= 512;
		dim[1] = rows;
		dim[2] = 512;
	}
	explib::bits data;
	data.resize(dim.volume());

	for (int i = 0; i < data.size(); ++i) {
		if ((rand() % 200) == 0) { // fill around every third chunk
			data[i] = true;
		}
	}

	int true_chunks = 0;
	TimeSentry time;
	explib::dim_row_iterator<p> i(dim);
	while (i) {
		explib::cvec<p> beginv = i.begin();
		int begini = dim.offset(beginv);
		explib::bits_ref row = data.from(begini, i.length());
		auto j = row.chunk_istream2<explib::false_tail_fill>();
		while (j) {
			explib::bchunk c;
			j>>c;
			if (c) true_chunks++;
		}
		++i;
	}
	long us = time.us();

	t.record({sup()<<"rowlen:"<<rowlen, "prop:us"}, us);
	t.record({sup()<<"rowlen:"<<rowlen, "prop:true"}, true_chunks);
}

void run_bchunk_iteration_test(TestTool& t) {
	srand(0);
	int size = 8*1024*1024; // 1 megabyte
	explib::bits data;
	data.resize(size);
	for (int i = 0; i < data.size(); ++i) {
		if ((rand() % 200) == 0) { // fill around every third 64 bit chunk
			data[i] = true;
		}
	}

	{
		int true_chunks = 0;
		TimeSentry time;
		explib::bchunk_istream i = data.from(0).chunk_istream(explib::false_tail_fill);
		while (i) {
			explib::bchunk c;
			i>>c;
			if (c) true_chunks++;
		}
		long us = time.us();
		t.record({"iteration:bcistream", "prop:us"}, us);
		t.record({"iteration:bcistream", "prop:true"}, true_chunks);
	}

	{
		int true_chunks = 0;
		TimeSentry time;
		for (explib::bchunk c : data.chunks()) {
			if (c) true_chunks++;
		}
		long us = time.us();
		t.record({"iteration:bcvector", "prop:us"}, us);
		t.record({"iteration:bcvector", "prop:true"}, true_chunks);
	}
	Table table(
		t.report(ToTable<Average>({}, "prop:", "iteration:")));
	t.reported()<<table;
}

void run_dim_bchunk_iteration_test(TestTool& t) {
	srand(0);
	typedef problem2d p;
	run_dim_bchunk_iteration_test_round<p>(t, 1);
	run_dim_bchunk_iteration_test_round<p>(t, 4);
	run_dim_bchunk_iteration_test_round<p>(t, 16);
	run_dim_bchunk_iteration_test_round<p>(t, 64);
	run_dim_bchunk_iteration_test_round<p>(t, 256);
	run_dim_bchunk_iteration_test_round<p>(t, 1024);
	run_dim_bchunk_iteration_test_round<p>(t, 4096);

	Table table(
		t.report(ToTable<Average>({}, "prop:", "rowlen:")));
	t.reported()<<table;
}

void run_3dim_bchunk_iteration_test(TestTool& t) {
	srand(0);
	typedef problem3d p;
	run_dim_bchunk_iteration_test_round<p>(t, 1);
	run_dim_bchunk_iteration_test_round<p>(t, 4);
	run_dim_bchunk_iteration_test_round<p>(t, 16);
	run_dim_bchunk_iteration_test_round<p>(t, 64);
	run_dim_bchunk_iteration_test_round<p>(t, 256);
	run_dim_bchunk_iteration_test_round<p>(t, 1024);
	run_dim_bchunk_iteration_test_round<p>(t, 4096);

	Table table(
		t.report(ToTable<Average>({}, "prop:", "rowlen:")));
	t.reported()<<table;
}

void run_bitmatrix_test(TestTool& t) {
	typedef problem2d p;
	explib::bitmatrix<p> m;
	m.unite(explib::cvec<p>(), explib::cvec<p>(3, 3));
	m[explib::cvec<p>(1, 1)] = true;
	m[explib::cvec<p>(0, 0)] = true;

	t<<"\nmatrix:\n"<<bitmatrix_tostring(m);

	explib::bitmatrix<p> m2;
	m2.unite(explib::cvec<p>(), explib::cvec<p>(3, 3));
	m2[explib::cvec<p>(2, 2)] = true;
	m2[explib::cvec<p>(0, 2)] = true;
	m2[explib::cvec<p>(2, 0)] = true;

	t<<"\nmatrix 2:\n"<<bitmatrix_tostring(m2);
	m.blit(explib::cvec<p>(0, 0), m2);
	t<<"\nmatrix 1+2:\n"<<bitmatrix_tostring(m);
}

void run_bitmatrix_rzblit_test(TestTool& t) {
	typedef problem2d p;
	explib::bitmatrix<p> m;
	m.unite(explib::cvec<p>(), explib::cvec<p>(3, 3));
	m[explib::cvec<p>(1, 1)] = true;
	m[explib::cvec<p>(0, 0)] = true;

	explib::bitmatrix<p> m2;
	m2.unite(explib::cvec<p>(), explib::cvec<p>(3, 3));
	m2[explib::cvec<p>(2, 2)] = true;
	m2[explib::cvec<p>(0, 2)] = true;
	m2[explib::cvec<p>(2, 0)] = true;

	explib::bitmatrix<p> m3;

	t<<"\nmatrix m1:\n"<<bitmatrix_tostring(m);
	t<<"\nmatrix m2:\n"<<bitmatrix_tostring(m2);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(0, 0), m);
	t<<"\nmatrix m3+=m1:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(0, 0), m2);
	t<<"\nmatrix m3+=m2:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(-1, -1), m);
	t<<"\nmatrix m3+=m1@-1,-1:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(-3, 2), m2);
	t<<"\nmatrix m3+=m1@-1,-1:\n"<<bitmatrix_tostring(m3);
}

void run_bitmatrix_rzsetblit_test(TestTool& t) {
	typedef problem2d p;
	explib::bitmatrix<p> m;
	m.unite(explib::cvec<p>(), explib::cvec<p>(3, 3));
	m[explib::cvec<p>(1, 1)] = true;
	m[explib::cvec<p>(0, 0)] = true;

	explib::bitmatrix<p> m2;
	m2.unite(explib::cvec<p>(), explib::cvec<p>(3, 3));
	m2[explib::cvec<p>(2, 2)] = true;
	m2[explib::cvec<p>(0, 2)] = true;
	m2[explib::cvec<p>(2, 0)] = true;

	explib::bitmatrix<p> m3;
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_set({0, 0}, true);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_set({10, 0}, true);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_set({-10, 0}, true);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_set({0, 10}, true);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_set({0, -10}, true);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_set({10, 10}, true);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);

	t<<"\nmatrix m1:\n"<<bitmatrix_tostring(m);
	t<<"\nmatrix m2:\n"<<bitmatrix_tostring(m2);
	t<<"\nmatrix m3:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(0, 0), m);
	t<<"\nmatrix m3+=m1:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(0, 0), m2);
	t<<"\nmatrix m3+=m2:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(-1, -1), m);
	t<<"\nmatrix m3+=m1@-1,-1:\n"<<bitmatrix_tostring(m3);
	m3.resizing_blit(explib::cvec<p>(-3, 2), m2);
	t<<"\nmatrix m3+=m1@-1,-1:\n"<<bitmatrix_tostring(m3);

	m3.resizing_set({3, 3}, true);
	m3.resizing_set({-3, 3}, true);
	m3.resizing_set({0, 1000}, true);

	explib::bitmatrix<p> m4;
	t<<"\nmatrix m4:\n"<<bitmatrix_tostring(m4);
	m4.resizing_blit({1, 1}, m);
	t<<"\nmatrix m4:\n"<<bitmatrix_tostring(m4);

	explib::bitmatrix<p> m5;
	t<<"\nmatrix m5:\n"<<bitmatrix_tostring(m5);
	m5.resizing_blit({-1, -1}, m);
	t<<"\nmatrix m5:\n"<<bitmatrix_tostring(m5);

}

void addutilstest(TestRunner& runner) {
	runner.add("utils/bit_ref", 			 {"func"}, &run_bit_ref_test);
	runner.add("utils/bit_copy",			 {"func"}, &run_bit_copy_test);
	runner.add("utils/bits_bitop", 			 {"func"}, &run_bits_bitop_test);
	runner.add("utils/bits_ref_andneg", 	 {"func"}, &run_bits_ref_andneg_test);
	runner.add("utils/bits_ref_and", 	 	 {"func"}, &run_bits_ref_and_test);
	runner.add("utils/bits_ref_or", 	 	 {"func"}, &run_bits_ref_or_test);
	runner.add("utils/bits_ref_fill",    	 {"func"}, &run_bits_ref_fill_test);
	runner.add("utils/bits_ref_perf",    	 {"perf"}, &run_bits_ref_perf_test);
	runner.add("utils/bits_andneg", 	 	 {"func"}, &run_bits_andneg_test);
	runner.add("utils/popcount", 			 {"func"}, &run_popcount_test);
	runner.add("utils/big_popcount",		 {"func"}, &run_big_popcount_test);
	runner.add("utils/big_16b_popcount",	 {"func"}, &run_big_16b_popcount_test);
	runner.add("utils/ndim",			     {"func"}, &run_ndim_test);
	runner.add("utils/bchunk_iteration", 	 {"func"}, &run_bchunk_iteration_test);
	runner.add("utils/dim_bchunk_iteration", {"func"}, &run_dim_bchunk_iteration_test);
	runner.add("utils/3dim_bchunk_iteration",{"func"},&run_3dim_bchunk_iteration_test);
	runner.add("utils/bitmatrix",			 {"func"}, &run_bitmatrix_test);
	runner.add("utils/bitmatrix_rzblit",	 {"func"}, &run_bitmatrix_rzblit_test);
	runner.add("utils/bitmatrix_rzsetblit",  {"func"}, &run_bitmatrix_rzsetblit_test);
}
