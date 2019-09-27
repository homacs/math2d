
#include <iostream>

#include <string>
#include <stdio.h>


#include <assert.h>
#include <stdlib.h>

#include <glm/gtc/constants.hpp>
#include <math2d/line.h>
using namespace glm;


#include "test_polynom.h"
#include "test_float-utils.h"
#include "test_line.h"
#include "test_bezier.h"
#include "test_Graph2DPlanar.h"
#include "Test_Triangulator.h"
#include "Test_MatrixMxM.h"
using namespace math2d;



#include <vector>
#include <algorithm>
using namespace std;



void run_all_tests() {
	Test_float_utils::all();
	Test_polynom::all();
	Test_MatrixMxM::all();
	test_line_all();
	test_Graph2DPlanar_all();
	Test_Triangulator::all();
	test_bezier_all();
}



int main(int argc, char** argv) {
	Test_polynom::test_polynom_N_roots();
	//  test_bezier_bezier_intersections();

	run_all_tests();
	return 0;
}
