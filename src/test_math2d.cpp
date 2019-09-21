
#include <iostream>

#include <string>
#include <stdio.h>


#include <assert.h>
#include "math/line.h"

#include <stdlib.h>

#include <glm/gtc/constants.hpp>
using namespace glm;


#include "test_polynom.h"
#include "test_float-utils.h"
#include "test_line.h"
#include "test_bezier.h"
#include "test_Graph2DPlanar.h"
#include "Test_Triangulator.h"
#include "Test_MatrixMxM.h"
using namespace math;



#include <vector>
#include <algorithm>
using namespace std;



void run_all_tests() {
	Test_float_utils::all();
	Test_MatrixMxM::all();
	test_line_all();
	test_Graph2DPlanar_all();
	Test_Triangulator::all();
	test_bezier_all();
}



int main(int argc, char** argv) {
	test_bezier_bezier_intersections();

	run_all_tests();
	return 0;
}
