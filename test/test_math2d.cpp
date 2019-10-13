


#define MATH2D_EVALUATE   1



#include <iostream>

#include <string>
#include <stdio.h>


#include <assert.h>
#include <stdlib.h>

#include <vector>
#include <algorithm>
using namespace std;



#include <glm/gtc/constants.hpp>
using namespace glm;





#include "test_polynom.h"
#include "test_float-utils.h"
#include "test_line.h"
#include "test_bezier.h"
#include "test_Graph2DPlanar.h"
#include "Test_Triangulator.h"
#include "Test_MatrixMxM.h"
using namespace math2d;

#define PERF_INIT() \
	double perf_duration; \
	perf_time_t perf_start, perf_end; \
	perf_clock(perf_start)
#define PERF_REPORT(name) \
	perf_clock(perf_end); \
	perf_duration = perf_end - perf_start; \
	printf("TIME '%s': %f\n", #name, perf_duration/perf_time_t::NSEC_PER_SEC); \
	perf_start = perf_end

#define PERF_DONE()

void run_all_tests() {
	PERF_INIT();

	Test_MatrixMxM::all();
	PERF_REPORT(Test_MatrixMxM);

	Test_float_utils::all();
	PERF_REPORT(Test_float_utils);
	Test_polynom::all();
	PERF_REPORT(Test_polynom);
	test_line_all();
	PERF_REPORT(test_line_all);
	test_bezier_all();
	PERF_REPORT(test_bezier_all);

	test_Graph2DPlanar_all();
	PERF_REPORT(test_Graph2DPlanar_all);
	Test_Triangulator::all();
	PERF_REPORT(Test_Triangulator);

	PERF_DONE();
}



int main(int argc, char** argv) {
//	Test_polynom::test_polynom_N_roots();
	Test_polynom::test_polynom4_roots();
//	PERF_INIT();
//    Test_polynom::test_polynom_N_roots_by_random_roots();
//	PERF_REPORT(test_polynom4_roots_random);
//	PERF_DONE();

	// test_bezier_bezier_intersections();

//	run_all_tests();
	return 0;
}
