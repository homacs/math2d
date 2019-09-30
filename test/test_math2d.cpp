


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
	test_bezier_point_closest_point();
	//  test_bezier_bezier_intersections();

	run_all_tests();
	return 0;
}
