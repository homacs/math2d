/*
 * test_polynom.h
 *
 *  Created on: 5 Aug 2019
 *      Author: homac
 */

#ifndef TEST_POLYNOM_H_
#define TEST_POLYNOM_H_

#include <math2d/polynom.h>
#include <stdlib.h>
#include "perf_clock.h"

using namespace math2d;
namespace Test_polynom {


static inline bool validate_roots(double* roots, int count, double* expected_roots, int expected_count, double tolerance) {

	assert(count == expected_count);
	int found = 0;
	for (int u = 0; u < count; u++) {
		for (int v = 0; v < count; v++) {
			if (about_equal(roots[u], expected_roots[v], tolerance)) {
				found++;
				expected_roots[v] = INFINITY;
				break;
			}
		}
	}
	return (found == count);

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
	float tolerance = 0.000001;
	polynom_t<N> W(params);
	double t;
	for (int i = 0; i < num_guesses; i++) {
		t = W.single_root(guesses[i], tolerance);
		assert(finite(t));
		guesses[i] = t;
	}

	double roots[N];
	int count = W.roots(roots, tolerance);
	validate_roots(roots, count, guesses, num_guesses, tolerance);
}

static inline void test_polynom_N_roots() {
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

	test_polynom_N_root();
	test_polynom_N_roots();

}



}



#endif /* TEST_POLYNOM_H_ */
