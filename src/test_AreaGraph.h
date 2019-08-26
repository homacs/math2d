/*
 * test_MultiGraph.h
 *
 *  Created on: 21 Jun 2019
 *      Author: homac
 */

#ifndef MATH_TEST_AREAGRAPH_H_
#define MATH_TEST_AREAGRAPH_H_

#include <math/AreaGraph.h>
#include "perf_clock.h"

#include <glm/gtc/constants.hpp>

using namespace glm;
using namespace math;



struct Test_AreaGraph {

#define HOLE AreaGraph::AREA_HOLE
#define FACE AreaGraph::AREA_FACE

#define ASSERT_NEIGHBOURSHIP(_G_, ...) { \
	unsigned sizes[] = {__VA_ARGS__}; \
	assert_Neighbourships(_G_, sizeof(sizes)/sizeof(unsigned), sizes);\
}
static inline int assert_Neighbourships(AreaGraph& g, int n, unsigned sizes[0]) {
	int checks = n;

	AreaGraph::Areas* areas = &(g.areas);
	assert(areas->AREA_VOID->neighbours.size() == sizes[0]);
	sizes[0] = -1;
	checks--;

	for (AreaGraph::Areas::iterator it_a = areas->begin(); it_a != areas->end(); it_a++) {
		AreaGraph::Area* area = *it_a;
		unsigned neighboursize = area->neighbours.size();
		for (int i = 0; i < n; i++)	{
			if (neighboursize == sizes[i]) {
				sizes[i] = -1;
				checks--;
				break;
			}
		}
	}
	// is the number of areas as expected and have all expected sizes been found?
	assert(checks==0);
	return 0;
}

#define ASSERT_WEIGHTS_AND_SIZES(_G_, ...) { \
	int checks[] = {__VA_ARGS__}; \
	assert_weights_and_sizes(_G_, sizeof(checks)/sizeof(int), checks);\
}

static int assert_weights_and_sizes(AreaGraph& g, int n, int args[0]) {
	struct check_t {
		int weight;
		unsigned size;
	};

	assert(n%2 == 0);
	n = n / 2;

	check_t* values = (check_t*)args;
	int checks = n;

	AreaGraph::Areas& areas = g.areas;
	assert(areas.AREA_VOID->weight == values[0].weight);
	assert(areas.AREA_VOID->edges.size() == values[0].size);
	values[0] = {0};
	checks--;

	for (AreaGraph::Areas::iterator it = areas.begin(); it != areas.end(); it++) {
		AreaGraph::Area* area = (*it);
		for (int i = 0; i < n; i++) {
			if (values[i].size == area->edges.size() && values[i].weight == area->weight) {
				values[i] = {0};
				checks--;
				break;
			}
		}
	}
	assert(checks == 0);
	return 0;
}
#define ASSERT_TYPES_AND_SIZES(_G_, ...) { \
	int checks[] = {__VA_ARGS__}; \
	assert_types_and_sizes(_G_, sizeof(checks)/sizeof(int), checks);\
}

static int assert_types_and_sizes(AreaGraph& g, int n, int args[0]) {
	struct check_t {
		unsigned type;
		unsigned size;
	};

	assert(n%2 == 0);
	n = n / 2;

	check_t* values = (check_t*)args;
	int checks = n;

	AreaGraph::AreaRule rule = g.areaRule;
	AreaGraph::Areas& areas = g.areas;
	assert(areas.AREA_VOID->type(rule) == values[0].type);
	assert(areas.AREA_VOID->edges.size() == values[0].size);
	values[0] = {0xffffffff, 0xffffffff};
	checks--;

	for (AreaGraph::Areas::iterator it = areas.begin(); it != areas.end(); it++) {
		AreaGraph::Area* area = (*it);
		unsigned size = area->edges.size();
		AreaGraph::AreaType type = area->type(rule);
		for (int i = 0; i < n; i++) {
			if (values[i].size == size && values[i].type == type) {
				values[i] = {0xffffffff, 0xffffffff};
				checks--;
				break;
			}
		}
	}
	assert(checks == 0);
	return 0;
}

static int test_findAreas_diamond_shapes() {

	AreaGraph g;
	AreaGraph::Areas* areas;
	/////////////////////////////////////
	// clockwise diamond


	g.reset();
	vec2 a = vec2(0,0);
	vec2 b = vec2(-1,1);
	vec2 c = vec2(0,2);
	vec2 d = vec2(1,1);


	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,d);
	g.addLineSegment(d,a);

	areas = g.findAreas();

	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);



	/////////////////////////////////////
	// counter clockwise diamond
	g.reset();

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	g.addLineSegment(d,c);
	g.addLineSegment(c,b);
	g.addLineSegment(b,a);
	g.addLineSegment(a,d);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);



	return 0;
}

static int test_findAreas_square_shapes() {
	AreaGraph g;
	AreaGraph::Areas* areas;

	// clockwise square
	vec2 a(0,1);
	vec2 b(1,1);
	vec2 c(1,0);
	vec2 d(0,0);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,d);
	g.addLineSegment(d,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);

	g.mergeAreasOfSameType2();
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);


	/////////////////////////////////////
	// counter clockwise square
	g.reset();

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);

	g.addLineSegment(d,c);
	g.addLineSegment(c,b);
	g.addLineSegment(b,a);
	g.addLineSegment(a,d);

	areas = g.findAreas();

	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);


	return 0;
}

static int test_findAreas_rectangle() {
	AreaGraph g;
	AreaGraph::Areas* areas;

	// clockwise square
	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(1,1);
	vec2 d(2,1);
	vec2 e(2,0);
	vec2 f(1,0);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);

	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,d);
	g.addLineSegment(d,e);
	g.addLineSegment(e,f);
	g.addLineSegment(f,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 6);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 6);

	/////////////////////////////////////
	// counter clockwise square
	g.reset();

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);

	g.addLineSegment(a,f);
	g.addLineSegment(f,e);
	g.addLineSegment(e,d);
	g.addLineSegment(d,c);
	g.addLineSegment(c,b);
	g.addLineSegment(b,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, -1, 6);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, -1, 6);

	return 0;
}

static int test_findAreas_V_shape() {
	AreaGraph g;
	AreaGraph::Areas* areas;

	vec2 a(0,0);
	vec2 b(0,2);
	vec2 c(1,1); // center
	vec2 d(2,2);
	vec2 e(2,0);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);

	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,d);
	g.addLineSegment(d,e);
	g.addLineSegment(e,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 5, 1, 5);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 5, 1, 5);
	return 0;
}

static int test_findAreas_A_shape() {
	AreaGraph g;
	AreaGraph::Areas* areas;

	vec2 a(0,0);
	vec2 b(0,2);
	vec2 c(1,1); // center
	vec2 d(2,2);
	vec2 e(2,0);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);

	g.addLineSegment(a,b);
	g.addLineSegment(b,d);
	g.addLineSegment(d,e);
	g.addLineSegment(e,c);
	g.addLineSegment(c,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 5, 1, 5);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 5, 1, 5);
	return 0;
}

static int test_findAreas_square_hole() {

	AreaGraph g;
	AreaGraph::Areas* areas;

	vec2 a(0,0);
	vec2 b(0,3);
	vec2 c(3,3);
	vec2 d(3,0);

	vec2 e(1,1);
	vec2 f(2,1);
	vec2 h(2,2);
	vec2 i(1,2);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);
	g.addNode(h);
	g.addNode(i);

	// square CW
	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,d);
	g.addLineSegment(d,a);

	// inner square CCW
	g.addLineSegment(e,f);
	g.addLineSegment(f,h);
	g.addLineSegment(h,i);
	g.addLineSegment(i,e);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);

	return 0;
}
static int test_findAreas_diamond_hole() {

	AreaGraph g;
	AreaGraph::Areas* areas;

	vec2 a(0,0);
	vec2 b(0,3);
	vec2 c(3,3);
	vec2 d(3,0);

	vec2 e(1.5, 1.0);
	vec2 f(2.0, 1.5);
	vec2 h(1.5, 2.0);
	vec2 i(1.0, 1.5);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);
	g.addNode(h);
	g.addNode(i);

	// square CW
	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,d);
	g.addLineSegment(d,a);

	// inner square CCW
	g.addLineSegment(e,f);
	g.addLineSegment(f,h);
	g.addLineSegment(h,i);
	g.addLineSegment(i,e);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);


	return 0;
}



static int test_findAreas_square_figure_eight_lying() {
	AreaGraph g;
	AreaGraph::Areas* areas;



	// CW square aside CCW square, sharing middle edge
	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(1,1);
	vec2 d(2,1);
	vec2 e(2,0);
	vec2 f(1,0);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);

	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,f);
	g.addLineSegment(f,a);

	g.addLineSegment(f,e);
	g.addLineSegment(e,d);
	g.addLineSegment(d,c);
	g.addLineSegment(c,f);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 4, -1, 4);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);


	return 0;
}


static int test_findAreas_square_figure_eight_standing() {
	AreaGraph g;
	AreaGraph::Areas* areas;

	// CW square aside CCW square, sharing middle edge
	vec2 a(0,0);
	vec2 b(-1,0);
	vec2 c(-1,1);
	vec2 d(-1,2);
	vec2 e(2,0);
	vec2 f(1,0);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);

	g.addLineSegment(a,b);
	g.addLineSegment(b,c);
	g.addLineSegment(c,f);
	g.addLineSegment(f,a);

	g.addLineSegment(f,e);
	g.addLineSegment(e,d);
	g.addLineSegment(d,c);
	g.addLineSegment(c,f);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 4, -1, 4);

	g.mergeAreasOfSameType2();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	return 0;
}


static int test_findAreas_diamond_figure_eight_lying() {
	AreaGraph g;
	AreaGraph::Areas* areas;

	// CW square aside CCW square, sharing middle edge
	vec2 h(1,0);
	vec2 d(3,0);

	vec2 a(0,1);
	vec2 c(2,1); // center
	vec2 e(4,1);

	vec2 b(1,2);
	vec2 f(3,2);

	g.addNode(a);
	g.addNode(b);
	g.addNode(c);
	g.addNode(d);
	g.addNode(e);
	g.addNode(f);
	g.addNode(h);

	g.addLineSegment(c,h);
	g.addLineSegment(h,a);
	g.addLineSegment(a,b);
	g.addLineSegment(b,c);

	g.addLineSegment(c,d);
	g.addLineSegment(d,e);
	g.addLineSegment(e,f);
	g.addLineSegment(f,c);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 8, 1, 4, -1, 4);


	g.mergeAreasOfSameType2();
	assert(areas->size() == 2);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 8, FACE, 4, FACE, 4);

	return 0;
}


static int assert_memory() {
	assert(    AreaGraph::Node::total_nodes_in_system == 0
			&& AreaGraph::Edge::total_edges_in_system == 0
			&& AreaGraph::Area::total_areas_in_system == 0
			);
	return 0;
}


static int all ()
{


	test_findAreas_square_shapes();
	assert_memory();

	test_findAreas_diamond_shapes();
	assert_memory();

	test_findAreas_rectangle();
	assert_memory();

	test_findAreas_V_shape();
	assert_memory();

	test_findAreas_A_shape();
	assert_memory();

	test_findAreas_square_hole();
	assert_memory();

	test_findAreas_diamond_hole();
	assert_memory();

	test_findAreas_square_figure_eight_lying();
	assert_memory();

	test_findAreas_square_figure_eight_standing();
	assert_memory();

	test_findAreas_diamond_figure_eight_lying();
	assert_memory();

	return 0;
}

}; // struct Test_AreaGraph


#endif /* MATH_TEST_AREAGRAPH_H_ */
