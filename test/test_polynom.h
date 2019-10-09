/*
 * test_polynom.h
 *
 *  Created on: 5 Aug 2019
 *      Author: homac
 */

#ifndef TEST_POLYNOM_H_
#define TEST_POLYNOM_H_

#include <stdlib.h>

#include <random>

#include "perf_clock.h"

#include "math2d/polynom.h"

using namespace math2d;
namespace Test_polynom {


static inline bool validate_roots(double* roots, int count, double* expected_roots, int expected_count, double tolerance) {

	if (count != expected_count) {
		return false;
	}
	int found = 0;
	for (int u = 0; u < count; u++) {
		if (u>0 && ((roots[u-1] > roots[u]) || (expected_roots[u-1] > expected_roots[u]))) {
			// we expect roots to be sorted in ascending order
			break;
		}

		if (finite(roots[u]) && about_equal(roots[u], expected_roots[u], tolerance)) {
			found++;
		}
	}
	return (found == count);

}


static inline void assert_roots(double* roots, int count, double* expected_roots, int expected_count, double tolerance) {
	assert(validate_roots(roots, count, expected_roots, expected_count, tolerance));
}


static inline void test_polynom_N_root() {
	float tolerance = 0.00001;
	double params[6] = {+1,-2,-2,+3,+1,-1};
	double guess[5] = {-1,-0.8,0.5,1, 2.2};
	polynom_t<5> W(params);
	double t;
	for (int i = 0; i < 5; i++) {
		t = W.single_root(guess[i], tolerance);
		assert(finite(t));
	}
}
template <int N>
static inline void testsub_polynom_N_roots(double params[N+1], double guesses[N], int num_guesses) {
	float tolerance = 0.000000001f;
	polynom_t<N> W(params);
	double t;
	for (int i = 0; i < num_guesses; i++) {
		t = W.single_root(guesses[i], tolerance);
		assert(finite(t));
		guesses[i] = t;
	}

	double roots[N];
	POLYNOM_N_ROOTS_EVALUATION_RESET();
	int count = W.roots(roots, tolerance);
	POLYNOM_N_ROOTS_EVALUATION_REPORT();
	assert(num_guesses == count);
	assert_roots(roots, count, guesses, num_guesses, tolerance);

}


static inline void test_polynom4_roots_random() {
	// TODO: this will be obsolete once we have decided,
	// which of the two functions to keep.
	double tolerance = 0.0005f;
	double params[5];
	double roots_1[4];
	double roots_2[4];
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(-1,1);
	const int iterations = 10000000;

	for (int i = 0; i < iterations; i++) {
		printf("%05d: ", i);
		for (int p = 0; p < 5; p++) {
			params[p] = distribution(generator);
			printf("+(%g) x^%d ", params[p], 4-p);
		}
		printf(" = 0\n");
		int count_1 = polynom_N_roots<4>(params, roots_1, 0.0000001f);
		int count_2 = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots_2);
		if (!validate_roots(roots_1, count_1, roots_2, count_2, tolerance)) {

			fprintf(stderr,"roots_1:\n");
			fprintf(stderr,"double expected[] = {\n");
			for (int i = 0; i < count_1; i++) {
				fprintf(stderr,"\t%g,\n", roots_1[i]);
			}
			fprintf(stderr, "};\n");
			fprintf(stderr, "int num_expected = %d;\n", count_1);
			fprintf(stderr, "\n");
			fprintf(stderr,"roots_2:\n");
			fprintf(stderr,"double expected[] = {\n");
			for (int i = 0; i < count_2; i++) {
				fprintf(stderr,"\t%g,\n", roots_2[i]);
			}
			fprintf(stderr, "};\n");
			fprintf(stderr, "int num_expected = %d;\n", count_2);
			fprintf(stderr, "\n");

			// output test case
			fprintf(stderr, "double params[] = {\n");
			for (int i =0; i < 5; i++) {
				fprintf(stderr, "\t%g,\n", params[i]);
			}
			fprintf(stderr, "};\n");
			fprintf(stderr, "\n");
			assert(false);
		}
	}
}


static inline void test_polynom4_roots()
{
	{
		// 0 = t^4 - 2 t^3 - 2 t^2 + 3 t + 1
		double params[] = {+1,-2, -2, +3, +1};
		double guesses[] = {
				-1.193527085,
				-0.294962899,
				+1.294962899,
				+2.193527085,
		};
		double roots[4];
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, guesses, 4, 0.000001f);
	}

	{
		// 0 = t^4 - 2 t^3 - 2 t^2 + 3 t + 1
		double params[] = {+1, -6, -2, +3, +1};
		double guesses[] = {
				-0.63950,
				-0.33905,
				+0.73918,
				+6.2394,
		};
		double roots[4];
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, guesses, 4, 0.0001f);
	}

}


static inline void test_polynom4_roots_2()
{
	double roots[4];
	{
		double expected[] = {
			-1,
		};
		int num_expected = 1;

		double params[] = {
			1,
			1,
			0,
			1,
			1,
		};
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, 0.0001f);

	}
	{
		double expected[] = {
			-1,
			1,
		};
		int num_expected = 2;

		double params[] = {
			1,
			0,
			0,
			0,
			-1,
		};
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, 0.0001f);

	}

	{
		double expected[] = {
			-0.104452,
			3.10445,
		};
		int num_expected = 2;

		double params[] = {
			-8,
			48,
			-128,
			168,
			19,
		};
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, 0.0001f);

	}

	{
		double expected[] = {
			0.012914,
			0.64132,
		};
		int num_expected = 2;

		double params[] = {
			-206,
			136,
			-244,
			158,
			-2,
		};
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, 0.0001f);

	}


	{
		double params[] = {
			95.6276,
			220.383,
			13.7875,
			78.8065,
			103.01,
		};
		double expected[] = {
			-2.30923,
			-0.707767,
		};
		int num_expected = 2;

		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, 0.0001f);

	}

	{
		// 0 = t^4 - 2 t^3 - 2 t^2 + 3 t + 1
		double params[] = {+1,-2, -2, +3, +1};
		double guesses[] = {
				-1.193527085,
				-0.294962899,
				+1.294962899,
				+2.193527085,
		};
		int num_expected = 4;
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, guesses, num_expected, 0.0001f);
	}

	{
		// 0 = t^4 - 2 t^3 - 2 t^2 + 3 t + 1
		double params[] = {+1, -6, -2, +3, +1};
		double guesses[] = {
				-0.63950,
				-0.33905,
				+0.73918,
				+6.2394,
		};
		int num_expected = 4;
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, guesses, num_expected, 0.0001f);
	}

	{
		// 0 = t^4 - 2 t^3 - 2 t^2 + 3 t + 1
		double params[] = {+1, +1, +2, +3, +5};
		double expected[] = {};
		int num_expected = 0;
		int count = polynom4_roots_2(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, 0.0001f);
	}
}



static inline void test_polynom_N_roots() {
	{
		double params[] = {0,0,0,0,0,0};
		double guesses[] = {0};
		testsub_polynom_N_roots<5>(params, guesses, 0);
	}
	{
		double params[] = {+1,-2,-2,+3,+1,-1};
		double guesses[] = {-1,-0.8, 0.5, 1, 2.2};
		testsub_polynom_N_roots<5>(params, guesses, 5);
	}
	{
		double params[] = {0,-2,-2,+3,+1,-1};
		double guesses[] = {-1.6,-0.8,0.6,0.8};
		testsub_polynom_N_roots<5>(params, guesses, 4);
	}
	{
		double params[] = {1,0,0,0,0};
		double guesses[] = {0};
		testsub_polynom_N_roots<4>(params, guesses, 1);
	}
	{
		double params[] = {1,0,0,0,0,0};
		double guesses[] = {0};
		testsub_polynom_N_roots<5>(params, guesses, 1);
	}
	{
		double params[] = {1,0,0,0,0,1};
		double guesses[] = {-0.5};
		testsub_polynom_N_roots<5>(params, guesses, 1);
	}
}




static inline void all() {
	// TODO: test_polynom needs actual tests
	test_polynom4_roots();
	test_polynom_N_root();
	test_polynom_N_roots();

}



}



#endif /* TEST_POLYNOM_H_ */
