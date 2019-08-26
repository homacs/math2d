/*
 * test_polynom.h
 *
 *  Created on: 5 Aug 2019
 *      Author: homac
 */

#ifndef TEST_POLYNOM_H_
#define TEST_POLYNOM_H_

#include "math/polynom.h"
#include <stdlib.h>
#include "perf_clock.h"

using namespace math;
namespace Test_polynom {


static inline bool validate_roots(double roots[3], int count, double expected_roots[3], int expected_count, double tolerance) {

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

static inline void all() {
	// TODO: test_polynom needs actual tests
	// NOTE: we have tested against another implementation
	//       of polynom3_roots and found rare deviations in
	//       range 0.000001 in 1000000 tests. Accuracy of
	//       the reference implementation was unknown, though.
}



}



#endif /* TEST_POLYNOM_H_ */
