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





static inline void print_equation(FILE* out, double params[0], int degree) {
	for (int p = 0; p < degree+1; p++) {
		fprintf(out, "+(%.6f) x^%d ", params[p], degree-p);
	}
	fprintf(out, " = 0\n");
}
static inline void print_roots(FILE* out, double roots[0], int count) {
	fprintf(out,"double expected[] = {\n");
	for (int i = 0; i < count; i++) {
		fprintf(out,"\t%.6f,\n", roots[i]);
	}
	fprintf(out, "};\n");
	fprintf(out, "int num_expected = %d;\n", count);
}

static inline void print_params(FILE* out, double params[0], int degree) {
	fprintf(out, "double params[] = {\n");
	for (int i =0; i < 5; i++) {
		fprintf(out, "\t%.6f,\n", params[i]);
	}
	fprintf(out, "};\n");
}

static inline void print_polynom_N_roots_testcase(FILE* out, double params[0], int degree, double expected_roots[0], int num_expected, double roots[0], int count) {
	fprintf(out, "// ");
	print_equation(stderr, params, degree);
	fprintf(out, "\n");
	fprintf(out, "\n");
	fprintf(out,"expected roots:\n");
	print_roots(out, expected_roots, num_expected);
	fprintf(out, "\n");

	fprintf(out,"found roots:\n");
	print_roots(out, roots, count);
	fprintf(out, "\n");

	// output test case
	print_params(out, params, degree);
	fprintf(out, "\n");

}

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
template <int DEGREE>
static inline void testsub_polynom_N_roots(double params[DEGREE+1], double guesses[DEGREE], int num_guesses) {
	float tolerance = 0.00000001f;
	polynom_t<DEGREE> W(params);
	// improve guesses
	double t;
	for (int i = 0; i < num_guesses; i++) {
		t = guesses[i];
		if (W(t) != 0) t = W.single_root(guesses[i], tolerance/4);
		assert(finite(t));
		guesses[i] = t;
	}

	double roots[DEGREE];
	POLYNOM_N_ROOTS_EVALUATION_RESET();
	int count = W.roots(roots, tolerance);
	POLYNOM_N_ROOTS_EVALUATION_REPORT();
	assert(num_guesses == count);
	bool valid = validate_roots(roots, count, guesses, num_guesses, tolerance);
	if (!valid) {
		print_polynom_N_roots_testcase(stderr, params, DEGREE, guesses, num_guesses, roots, count);
		assert(false);
	}
}

static inline void test_polynom4_roots_random() {
	// TODO: this will become an evaluation function

	if (MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS) {
		fprintf(stderr, "This function can only work if MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS is set to 0");
		exit(-1);
	}


	double tolerance = 0.005;
	const int degree = 4;
	double params[degree+1];
	double roots_1[degree];
	double roots_2[degree];
	std::default_random_engine generator;
	// std::uniform_int_distribution<int> distribution(-1,1);
	std::uniform_int_distribution<int> distribution(-1000,1000);
	// std::uniform_real_distribution<double> distribution(-1000,1000);

#ifndef TEST_FEEDBACK
#define TEST_FEEDBACK 1
#endif
	POLYNOM_N_ROOTS_EVALUATION_RESET();
	const int iterations = 10000000;
	for (int i = 0; i < iterations; i++) {
		for (int p = 0; p < degree+1; p++) {
			params[p] = distribution(generator);
		}
#if TEST_FEEDBACK
		printf("%05d: ", i); print_equation(stdout, params, degree);
#endif

		int count_1 = polynom_N_roots<degree>(params, roots_1, 0.00000001);
		int count_2 = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots_2);
		if (!validate_roots(roots_1, count_1, roots_2, count_2, tolerance)) {
			print_polynom_N_roots_testcase(stderr, params, degree, roots_1, count_1, roots_2, count_2);
			assert(false);
		}
	}
	POLYNOM_N_ROOTS_EVALUATION_REPORT();
}



static inline void testsub_params_from_roots(double roots[0], double params[0], int num_roots) {
	params[0] = 1.0;
	for (int degree = 1; degree <= num_roots; degree++) {
		double x = roots[degree-1];
		params[degree] = params[degree-1]*(-x);
		for (int p = degree-1; p > 0; p--) {
			params[p] = params[p-1]*(-x) + params[p];
		}
	}
}

static inline void test_polynom4_roots_by_random_roots() {
	// TODO: this will become an evaluation function

	if (MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS) {
		fprintf(stderr, "This function can only work if MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS is set to 0");
		exit(-1);
	}


	double tolerance = 0.005;
	const int degree = 4;
	double params[degree+1];
	double roots_1[degree];
	double roots_2[degree];
	std::default_random_engine generator;
	//std::uniform_int_distribution<int> distribution(-1,1);
	std::uniform_int_distribution<int> distribution(-10,10);
	// std::uniform_real_distribution<double> distribution(-1000,1000);



#ifndef TEST_FEEDBACK
#define TEST_FEEDBACK 1
#endif
	POLYNOM_N_ROOTS_EVALUATION_RESET();
	const int iterations = 10000000;
	for (int i = 0; i < iterations; i++) {
		for (int x = 0; x < degree; x++) {
			roots_1[x] = distribution(generator);
		}

		testsub_params_from_roots(roots_1, params, degree);
		sort(roots_1, roots_1+degree);
		int count_1 = unique(roots_1, roots_1+degree) - roots_1;
#if TEST_FEEDBACK
		printf("%05d: ", i); print_equation(stdout, params, degree);
#endif
		int count_2 = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots_2);
		if (!validate_roots(roots_1, count_1, roots_2, count_2, tolerance)) {
			print_polynom_N_roots_testcase(stderr, params, degree, roots_1, count_1, roots_2, count_2);
			assert(false);
		}
	}
	POLYNOM_N_ROOTS_EVALUATION_REPORT();
}

static inline void test_polynom_N_roots_by_random_roots() {
	// TODO: this will become an evaluation function

	double tolerance = 0.000005;
	const int degree = 4;
	double params[degree+1];
	double roots_1[degree];
	double roots_2[degree];
	std::default_random_engine generator;
	//std::uniform_int_distribution<int> distribution(-1,1);
	std::uniform_int_distribution<int> distribution(-10,10);
	// std::uniform_real_distribution<double> distribution(-1000,1000);



#ifndef TEST_FEEDBACK
#define TEST_FEEDBACK 1
#endif
	POLYNOM_N_ROOTS_EVALUATION_RESET();
	const int iterations = 10000000;
	for (int i = 0; i < iterations; i++) {
		for (int x = 0; x < degree; x++) {
			roots_1[x] = distribution(generator);
		}

		testsub_params_from_roots(roots_1, params, degree);
		sort(roots_1, roots_1+degree);
		int count_1 = unique(roots_1, roots_1+degree) - roots_1;
#if TEST_FEEDBACK
		printf("%05d: ", i); print_equation(stdout, params, degree);
#endif
		int count_2 = polynom_N_roots<4>(params, roots_2, tolerance);
//		int count_2 = __internal__polynom_N_roots_highp(params, 4, roots_2, tolerance);
		if (!validate_roots(roots_1, count_1, roots_2, count_2, tolerance)) {
			print_polynom_N_roots_testcase(stderr, params, degree, roots_1, count_1, roots_2, count_2);
			assert(false);
		}
	}
	POLYNOM_N_ROOTS_EVALUATION_REPORT();
}



static inline void test_polynom4_roots()
{
	double tolerance = 0.0001f;
	double roots[4];




	{
		double expected[] = {
			-6.000000,
			8.000000,
			9.000000,
		};
		int num_expected = 3;

		double params[] = {
			1.000000,
			-19.000000,
			58.000000,
			672.000000,
			-3456.000000,
		};
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);

	}
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
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);

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
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);

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
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);

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
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);

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

		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);

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
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, guesses, num_expected, tolerance);
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
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, guesses, num_expected, tolerance);
	}

	{
		// 0 = t^4 - 2 t^3 - 2 t^2 + 3 t + 1
		double params[] = {+1, +1, +2, +3, +5};
		double expected[] = {};
		int num_expected = 0;
		int count = polynom4_roots(params[0], params[1], params[2], params[3], params[4], roots);
		assert_roots(roots, count, expected, num_expected, tolerance);
	}
}



static inline void test_polynom_N_roots() {
	{
		// double root at 0.0
		double params[] = {
			-149.000000,
			584.000000,
			-54.000000,
			0.000000,
			0.000000,
		};
		double expected[] = {
			0.000000,
			0.094757,
			3.824707,
		};
		int num_expected = 3;
		testsub_polynom_N_roots<4>(params, expected, num_expected);
	}
	{
		// double root at 0.0
		double params[] = {752,2172,0,0};
		double guesses[] = {-2.8883,0.0};
		testsub_polynom_N_roots<3>(params, guesses, 2);
	}
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
