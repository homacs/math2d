/*
 * test_MultiGraph.h
 *
 *  Created on: 21 Jun 2019
 *      Author: homac
 */

#ifndef MATH_TEST_AREAGRAPH_H_
#define MATH_TEST_AREAGRAPH_H_

#include "perf_clock.h"

#include <glm/gtc/constants.hpp>
#include <math2d/Triangulator.h>

using namespace glm;
using namespace math2d;

struct V_LT_X {
	bool operator () (const vec2& a, const vec2& b) const {
		return a.x < b.x || (a.x == b.x && a.y < b.y);
	}
};
static const V_LT_X vec2_lt_x;



struct Test_Triangulator {

	static const Triangulator::AreaType HOLE = Triangulator::AREA_HOLE;
	static const Triangulator::AreaType FACE = Triangulator::AREA_FACE;
	static const Triangulator::WindingOrder WINDING_CW  = Triangulator::WO_CW;
	static const Triangulator::WindingOrder WINDING_CCW = Triangulator::WO_CCW;

#define ASSERT_WEIGHTS_AND_SIZES(_G_, ...) { \
	int checks[] = {__VA_ARGS__}; \
	assert_weights_and_sizes(_G_, sizeof(checks)/sizeof(int), checks);\
}

static int assert_weights_and_sizes(Triangulator& g, int n, int args[0]) {
	struct check_t {
		int weight;
		unsigned size;
	};

	assert(n%2 == 0);
	n = n / 2;

	check_t* values = (check_t*)args;
	int checks = n;

	Triangulator::Areas& areas = g.areas;
	assert(areas.AREA_VOID->weight == values[0].weight);
	assert(areas.AREA_VOID->edges.size() == values[0].size);
	values[0] = {0};
	checks--;

	for (Triangulator::Areas::iterator it = areas.begin(); it != areas.end(); it++) {
		Triangulator::Area* area = (*it);
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

static int assert_types_and_sizes(Triangulator& g, int n, int args[0]) {
	struct check_t {
		unsigned type;
		unsigned size;
	};

	const unsigned INVALID = 0xffffffff;


	assert(n%2 == 0);
	n = n / 2;

	check_t* values = (check_t*)args;
	int checks = n;

	Triangulator::AreaRule rule = g.areaRule;
	Triangulator::Areas& areas = g.areas;
	Triangulator::AreaType type = areas.AREA_VOID->type(rule);
	unsigned size = areas.AREA_VOID->edges.size();
	assert(type == values[0].type);
	assert(size == values[0].size);
	values[0] = {INVALID, INVALID};
	checks--;

	for (Triangulator::Areas::iterator it = areas.begin(); it != areas.end(); it++) {
		Triangulator::Area* area = (*it);
		size = area->edges.size();
		type = area->type(rule);
		for (int i = 0; i < n; i++) {
			if (values[i].size == size && values[i].type == type) {
				values[i] = {INVALID, INVALID};
				checks--;
				break;
			}
		}
	}
	assert(checks == 0);
	return 0;
}


#define ASSERT_SPLITS(splits, ...) {\
	vec2 expected[] = { __VA_ARGS__}; \
	assert_splits(splits, sizeof(expected)/sizeof(vec2), expected);\
}

static void assert_splits(vector<Triangulator::AreaEdge*>& splits, int n, vec2 expected[0]) {
	typedef vector<Triangulator::AreaEdge*> AreaEdges;
	// two nodes per split edge!
	assert(n%2 == 0);

	const float INVALID = INFINITY;

	// correct number of expected splits?
	assert(unsigned(n)/2 == splits.size());

	unsigned found = 0;
	for (AreaEdges::iterator it_s = splits.begin(); it_s != splits.end(); it_s++) {
		Triangulator::AreaEdge* split = *it_s;
		for (int i = 0; i < n; i+=2) {
			vec2& a = expected[i];
			vec2& b = expected[i+1];
			// a has to be before b

			assert(a.y == INVALID || a.y < b.y || (a.y == b.y && a.x < b.x));
			if (split->a->pos == a && split->b->pos == b) {
				a.x = INVALID;
				a.y = INVALID;
				b.x = INVALID;
				b.y = INVALID;
				found++;
				break;
			}
		}
	}
	// all splits as expected?
	assert(found == splits.size());
}



#define POINTS(...) ({ \
	vec2 points[] = {__VA_ARGS__};\
	unsigned size = sizeof(points)/sizeof(vec2);\
	new vector<vec2>(points, points+size);\
})

#define ASSERT_POLYGONS(P, ...) {\
		vector<vec2>* expected[] = {__VA_ARGS__};\
		assert_polygons(P, sizeof(expected)/sizeof(void*), expected);\
}

static void assert_polygons(vector<Triangulator::Polygon*>& polygons, unsigned size, vector<vec2>* expected[0]) {
	assert(size == polygons.size());
	typedef Triangulator::Polygon Polygon;
	typedef vector<Polygon*> Polygons;

	for (Polygons::iterator it = polygons.begin(); it != polygons.end(); it++) {
		Polygon& poly = *(*it);

		Polygon::iterator it_v_begin = poly.begin();

		// search for a matching sequence in set of expected sequences
		for (unsigned i = 0; i < size; i++) {
			bool match = false;
			vector<vec2>* expected_points = expected[i];
			if (expected_points == NULL) continue;

			// search for match in this sequence of expected points
			for (vector<vec2>::iterator it_ep = expected_points->begin(); it_ep != expected_points->end(); it_ep++) {
				if (*it_v_begin == *it_ep) {
					Polygon::iterator it_v = it_v_begin;
					while (it_v != poly.end() && *it_v == *it_ep) {
//						printf("match: (%2.f,%2.f) == (%2.f,%2.f)\n", it_v->x, it_v->y, it_ep->x, it_ep->y);
						// advance in both lists
						it_v++;
						it_ep++;
						// cycle to start when hitting end
						if (it_ep == expected_points->end()) it_ep = expected_points->begin();
					}
					match = (it_v == poly.end()) && (poly.size() == expected_points->size());
					break;
				}
			}
			if (match) {
				// found match -> note and go to next polygon
				delete expected[i];
				expected[i] = NULL;
				break;
			}
		}



	}


	for (unsigned i = 0; i < size; i++) {
		assert(expected[i] == 0);
	}
}




#define CREATE_PATH(_G_, ...) {\
	vec2 nodes[] = {__VA_ARGS__};\
	create_path(_G_, sizeof(nodes)/sizeof(vec2), nodes);\
}


static void create_path(Triangulator& g, unsigned size, vec2 nodes[0]) {
	assert(size > 2);
	unsigned i = 0;
	for (; i < size-1; i++) {
		g.addLineSegment(nodes[i], nodes[i+1]);
	}
	g.addLineSegment(nodes[i], nodes[0]);
}


static void test_monotonize_no_splits(Triangulator& g) {
	Triangulator::AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits);
}


static int test_findAreas_diamond_shapes() {

	Triangulator g;
	Triangulator::Areas* areas;
	/////////////////////////////////////
	//    c   '
	//   / \  '
	//  b   d '
	//   \ /  '
	//    a   '
	g.reset();
	vec2 a = vec2(0,0);
	vec2 b = vec2(-1,1);
	vec2 c = vec2(0,2);
	vec2 d = vec2(1,1);

	CREATE_PATH(g, a,b,c,d);

	areas = g.findAreas();

	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);

	test_monotonize_no_splits(g);

	/////////////////////////////////////
	// counter clockwise diamond
	g.reset();
	CREATE_PATH(g, d,c,b,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);

	test_monotonize_no_splits(g);

	return 0;
}

static int test_findAreas_square_shapes() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a---b
	// |   |
	// d---c
	// clockwise square
	vec2 a(0,1);
	vec2 b(1,1);
	vec2 c(1,0);
	vec2 d(0,0);

	CREATE_PATH(g, a,b,c,d);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);

	g.mergeAreasOfSameType();
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 1, 4);

	test_monotonize_no_splits(g);



	/////////////////////////////////////
	// counter clockwise square
	g.reset();
	CREATE_PATH(g, d,c,b,a);

	areas = g.findAreas();

	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, -1, 4);

	test_monotonize_no_splits(g);

	return 0;
}

static int test_findAreas_rectangle() {
	Triangulator g;
	Triangulator::Areas* areas;
	// b---c---d
	// |       |
	// |       |
	// |       |
	// a---f---e
	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(1,1);
	vec2 d(2,1);
	vec2 e(2,0);
	vec2 f(1,0);

	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 6);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 6);

	test_monotonize_no_splits(g);

	/////////////////////////////////////
	// counter clockwise
	g.reset();
	CREATE_PATH(g, f,e,d,c,b,a);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, -1, 6);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, -1, 6);

	test_monotonize_no_splits(g);

	return 0;
}

static int test_findAreas_V_shape() {
	Triangulator g;
	Triangulator::Areas* areas;

	// b   d
	// |\ /|
	// | c |
	// |   |
	// a---e
	vec2 a(0,0);
	vec2 b(0,2);
	vec2 c(1,1); // center
	vec2 d(2,2);
	vec2 e(2,0);

	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 5, 1, 5);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	test_monotonize_no_splits(g);
	return 0;
}

static int test_findAreas_A_shape() {
	Triangulator g;
	Triangulator::Areas* areas;
	// b---d
	// |   |
	// | c |
	// |/ \|
	// a   e
	vec2 a(0,0);
	vec2 b(0,2);
	vec2 c(1,1); // center
	vec2 d(2,2);
	vec2 e(2,0);
	CREATE_PATH(g, a,b,d,e,c);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 5, 1, 5);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	test_monotonize_no_splits(g);
	return 0;
}

static int test_findAreas_square_hole() {

	Triangulator g;
	Triangulator::Areas* areas;
	// b-------c    '
	// | i---h |    '
	// | |   | |    '
	// | e---f |    '
	// a-------d    '

	vec2 a(0,0);
	vec2 b(0,3);
	vec2 c(3,3);
	vec2 d(3,0);

	vec2 e(1,1);
	vec2 f(2,1);
	vec2 h(2,2);
	vec2 i(1,2);

	// square CW
	CREATE_PATH(g, a,b,c,d);

	// inner square CCW
	CREATE_PATH(g, e,f,h,i);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);

	g.mergeAreasOfSameType();
	assert(areas->size() == 2);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 4, HOLE, 4, FACE, 8);


	Triangulator::AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, e,b, f,c);

	return 0;
}
static int test_findAreas_diamond_hole() {

	Triangulator g;
	Triangulator::Areas* areas;
	// b---------c    '
	// |    h    |    '
	// |   / \   |    '
	// |  i   f  |    '
	// |   \ /   |    '
	// |    e    |    '
	// a---------d    '
	vec2 a(0,0);
	vec2 b(0,3);
	vec2 c(3,3);
	vec2 d(3,0);

	vec2 e(1.5, 1.0);
	vec2 f(2.0, 1.5);
	vec2 h(1.5, 2.0);
	vec2 i(1.0, 1.5);

	// square CW
	CREATE_PATH(g, a,b,c,d);

	// inner square CCW
	CREATE_PATH(g, e,f,h,i);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);

	g.mergeAreasOfSameType();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 4, 0, 4, 1, 8);



	Triangulator::AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, i,b, f,c);



	return 0;
}



static int test_findAreas_square_figure_eight_lying() {
	Triangulator g;
	Triangulator::Areas* areas;
	// b---c---d  '
	// |   |   |  '
	// a---f---e  '
	// CW square aside CCW square, sharing middle edge
	vec2 a(0,0);
	vec2 b(0,1);
	vec2 c(1,1);
	vec2 d(2,1);
	vec2 e(2,0);
	vec2 f(1,0);

	CREATE_PATH(g, a,b,c,f,e,d,c,f);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 4, -1, 4);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	test_monotonize_no_splits(g);
	return 0;
}


static int test_findAreas_square_figure_eight_standing() {
	Triangulator g;
	Triangulator::Areas* areas;
	// d---e   '
	// |   |   '
	// c---f   '
	// |   |   '
	// b---a   '
	// CW square aside CCW square, sharing middle edge
	vec2 a(0,0);
	vec2 b(-1,0);
	vec2 c(-1,1);
	vec2 d(-1,2);
	vec2 e(0,2);
	vec2 f(0,1);


	CREATE_PATH(g, a,b,c,f,e,d,c,f)

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 6, 1, 4, -1, 4);

	g.mergeAreasOfSameType();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	test_monotonize_no_splits(g);

	return 0;
}


static int test_findAreas_diamond_figure_eight_lying() {
	Triangulator g;
	Triangulator::Areas* areas;
	//   b   f   '
	//  / \ / \  '
	// a   c   e '
	//  \ / \ /  '
	//   h   d   '
	vec2 b(1,2);
	vec2 f(3,2);

	vec2 a(0,1);
	vec2 c(2,1); // center
	vec2 e(4,1);

	vec2 h(1,0);
	vec2 d(3,0);
	CREATE_PATH(g, a,b,c,d,e,f,c,h);

	areas = g.findAreas();
	assert(areas->size() == 2);
	ASSERT_WEIGHTS_AND_SIZES(g, 0, 8, 1, 4, -1, 4);


	g.mergeAreasOfSameType();
	assert(areas->size() == 2);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 8, FACE, 4, FACE, 4);


	test_monotonize_no_splits(g);



	return 0;
}



static void test_monotonize_one_split_edge_neighbour() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a------b
	// |     /
	// |    c
	// |     \   '
	// e------d
	vec2 a(0,1);
	vec2 b(2,1);
	vec2 c(1,0); // split node
	vec2 d(2,-1);
	vec2 e(0,-1);
	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, c,a);

	vector<Triangulator::Polygon*> polygons;
	g.polygonizeFaces(polygons, Triangulator::WO_CW);
	ASSERT_POLYGONS(polygons, POINTS(a,b,c), POINTS(a,c,d,e));

	polygons.clear();
	g.polygonizeFaces(polygons, Triangulator::WO_CCW);
	ASSERT_POLYGONS(polygons, POINTS(c,b,a), POINTS(e,d,c,a));

}

static void test_monotonize_one_split_node_neighbour() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a------b
	// |     /
	// | f  c
	// |/ \  \   '
	// h   e--d
	vec2 a(-2,1);
	vec2 b(2,1);
	vec2 c(1,0); // split node
	vec2 d(2,-1);
	vec2 e(0,-1);
	vec2 f(-1,0); // closest connection node
	vec2 h(-2,-1);
	CREATE_PATH(g, a,b,c,d,e,f,h);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 7, FACE, 7);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, f,c);

}


static void test_monotonize_one_split_vertical() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a------b
	// |     /
	// |    c
	// |    |   '
	// e----d
	vec2 a(0,1);
	vec2 b(2,1);
	vec2 c(1,0); // split node
	vec2 d(1,-1);// vertical down
	vec2 e(0,-1);
	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits);


	g.reset();
	// a----b
	// |    |
	// |    c
	// |     \   '
	// e------d
	b = vec2(1,1); // vertical up
	d = vec2(2,-1); // set d to diagonal right
	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	splits.clear();
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, c,a);

}

static void test_monotonize_one_split_vertical_2() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a------b
	// |     /
	// |    c
	// |    |   '
	// |    d
	// |    |
	// f----e
	vec2 a(0,1);
	vec2 b(2,1);
	vec2 c(1,0); // split node
	vec2 d(1,-1);
	vec2 e(1,-2);
	vec2 f(0,-2);
	CREATE_PATH(g, a,b,c,d,e, f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits);


	g.reset();
	// a----b
	// |    |
	// |    c
	// |    |
	// |    d
	// |     \   '
	// f------e
	b = vec2(1,1); // vertical up
	e = vec2(2,-2); // set d to diagonal right
	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	splits.clear();
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, d,a);

}


static void test_monotonize_one_merge_edge_neighbour() {
	Triangulator g;
	Triangulator::Areas* areas;
	// b-----a
	//  \    |
	//   c   |
	//  /    |
	// d-----e
	vec2 a(2,1);
	vec2 b(0,1);
	vec2 c(1,0); // merge node, pointing right
	vec2 d(0,-1);
	vec2 e(2,-1);
	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, c, a);

}


static void test_monotonize_one_merge_node_neighbour() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a-----h
	//  \    |
	//   b e |
	//  / / \|
	// c-d   f
	vec2 a(-1,1);
	vec2 b(0,0); // merge node
	vec2 c(-1,-1);
	vec2 d(0,-1);
	vec2 e(1,0); // closest connection node
	vec2 f(2,-1);
	vec2 h(2,+1);

	// expected split edge (b,e)
	CREATE_PATH(g, a,b,c,d,e,f,h);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 7, FACE, 7);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, b, e);
}



static void test_monotonize_one_merge_vertical() {
	//
	// merge with vertical up (expect split (c,a))
	//

	Triangulator g;
	Triangulator::Areas* areas;
	//   b----a
	//   |    |
	//   c    |
	//  /     |
	// d------e
	vec2 a(2,1);
	vec2 b(1,1); // vertical up from merge node
	vec2 c(1,0); // merge node, pointing right
	vec2 d(0,-1);
	vec2 e(2,-1);
	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, c, a);

	//
	// merge with vertical down (expect no split edge)
	//

	g.reset();
	// b------a
	//  \     |
	//   c    |
	//   |    |
	//   d----e
	b = vec2(0,1); // back to diagonal
	d = vec2(1,-1); // vertical down
	CREATE_PATH(g, a,b,c,d,e);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 5, FACE, 5);

	splits.clear();
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits);
}


static void test_monotonize_one_merge_vertical_2() {
	//
	// merge with vertical up (expect split (c,a))
	//

	Triangulator g;
	Triangulator::Areas* areas;
	//   b----a
	//   |    |
	//   c    |
	//   |    |
	//   d    |
	//  /     |
	// e------f
	vec2 a(2,1);
	vec2 b(1,1);
	vec2 c(1,0);
	vec2 d(1,-1);
	vec2 e(0,-2);
	vec2 f(2,-2);
	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	vector<Triangulator::AreaEdge*> splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, d, a);

	//
	// merge with vertical down (expect no split edge)
	//

	g.reset();
	// b------a
	//  \     |
	//   c    |
	//   |    |
	//   d    |
	//   |    |
	//   e----f
	b = vec2(0,1); // diagonal
	e = vec2(1,-2); // vertical down
	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	splits.clear();
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits);
}



static void test_monotonize_merge_then_split_same_scanline() {
	Triangulator g;
	Triangulator::Areas* areas;
	// a------b
	//  \    /
	//   f  c
	//  /    \   '
	// e------d
	vec2 a(0,1);
	vec2 b(3,1);
	vec2 c(2,0); // split node
	vec2 d(3,-1);
	vec2 e(0,-1);
	vec2 f(1,0); // merge node

	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, f, c);
}

static void test_monotonize_split_then_merge_subsequent_scanline() {
	Triangulator g;
	Triangulator::Areas* areas;
	//  c------d '
	//   \    /  '
	//    b  /   '
	//   /  e    '
	//  /    \   '
	// a------f  '
	vec2 a(0,0);
	vec2 b(3,3); // merge node
	vec2 c(1,5);
	vec2 d(8,5);
	vec2 e(5,2); // split node
	vec2 f(7,0);

	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, e, b);
}

static void test_monotonize_merge_then_split_subsequent_scanline() {
	Triangulator g;
	Triangulator::Areas* areas;
	// c------d   '
	//  \    /    '
	//   \  e     '
	//    b  \    '
	//   /    \   '
	//  a------f  '
	vec2 a(1,0);
	vec2 b(3,2); // merge node
	vec2 c(0,5);
	vec2 d(7,5);
	vec2 e(5,3); // split node
	vec2 f(8,0);

	CREATE_PATH(g, a,b,c,d,e,f);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, b, e);
}

static void test_monotonize_split_with_conflict_above() {
	Triangulator g;
	Triangulator::Areas* areas;
	// 1   5-------4   '
	// \\ /       /    '
	//  \6       /     '
	//   \      3      '
	//    \      \     '
	//     \     \     '
	//      \    \     '
	//       \   \     '
	//        \   \    '
	//         \  \    '
	//          \ \    '
	//           \\    '
	//            2    '
	vec2 v1(0,0);
	vec2 v2(11,-12);
	vec2 v3(10,-2);
	vec2 v4(11,0);
	vec2 v5(2,0);
	vec2 v6(1,-1);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v3,v6);
}

static void test_monotonize_split_with_conflict_below() {
	Triangulator g;
	Triangulator::Areas* areas;
	//            2  '
	//           //    '
	//          / /    '
	//         /  /    '
	//        /   /     '
	//       /   /     '
	//      /    /      '
	//     /     /      '
	//    /      /       '
	//   /      3        '
	//  /6       \       '
	// // \       \      '
	// 1   5-------4     '
	vec2 v1(0,0);
	vec2 v2(11,12);
	vec2 v3(10,2);
	vec2 v4(11,0);
	vec2 v5(2,0);
	vec2 v6(1,1);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v6,v3);
}

static void test_monotonize_split_with_horizontal_conflict_above() {
	Triangulator g;
	Triangulator::Areas* areas;
	// 1           4   '
	// \\        //    '
	//  \6------5/     '
	//   \      3      '
	//    \      \     '
	//     \     \     '
	//      \    \     '
	//       \   \     '
	//        \   \    '
	//         \  \    '
	//          \ \    '
	//           \\    '
	//            2    '
	vec2 v1(0,0);
	vec2 v2(11,-12);
	vec2 v3(10,-2);
	vec2 v4(11,0);
	vec2 v5(10,-1);
	vec2 v6(1,-1);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v3,v6);


}

static void test_monotonize_split_with_horizontal_conflict_below() {
	Triangulator g;
	Triangulator::Areas* areas;
	//            2  '
	//           //    '
	//          / /    '
	//         /  /    '
	//        /   /     '
	//       /   /     '
	//      /    /      '
	//     /     /      '
	//    /      /       '
	//   /      3        '
	//  /6------5\       '
	// //        \\      '
	// 1          4     '
	vec2 v1(0,0);
	vec2 v2(11,12);
	vec2 v3(10,2);
	vec2 v4(11,0);
	vec2 v5(10,1);
	vec2 v6(1,1);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v6,v3);
}

static void test_monotonize_merge_with_conflict_above() {
	Triangulator g;
	Triangulator::Areas* areas;
	// 1-------2   4   '
	//  \       \ //   '
	//   \       3/    '
	//    6      /     '
	//   /      /      '
	//   /     /       '
	//   /    /        '
	//   /   /         '
	//  /   /          '
	//  /  /           '
	//  / /            '
	//  //             '
	//  5              '
	vec2 v1(0,0);
	vec2 v2(10,0);
	vec2 v3(11,-1);
	vec2 v4(12,0);
	vec2 v5(1,-12);
	vec2 v6(2,-2);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v6,v3);
}

static void test_monotonize_merge_with_conflict_below() {
	Triangulator g;
	Triangulator::Areas* areas;
	//  5              '
	//  \\             '
	//  \ \            '
	//  \  \           '
	//  \   \          '
	//   \   \         '
	//   \    \        '
	//   \     \       '
	//   \      \      '
	//    6      \     '
	//   /       3\    '
	//  /       / \\   '
	// 1-------2   4  '
	vec2 v1(0,0);
	vec2 v2(10,0);
	vec2 v3(11,1);
	vec2 v4(12,0);
	vec2 v5(1,12);
	vec2 v6(2,2);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v3,v6);
}

static void test_monotonize_merge_with_horizontal_conflict_above() {
	Triangulator g;
	Triangulator::Areas* areas;
	// 1           4   '
	//  \\        //   '
	//   \2------3/    '
	//    6      /     '
	//   /      /      '
	//   /     /       '
	//   /    /        '
	//   /   /         '
	//  /   /          '
	//  /  /           '
	//  / /            '
	//  //             '
	//  5              '
	vec2 v1(0,0);
	vec2 v2(2,-1);
	vec2 v3(11,-1);
	vec2 v4(12,0);
	vec2 v5(1,-12);
	vec2 v6(2,-2);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v6,v3);
}

static void test_monotonize_merge_with_horizontal_conflict_below() {
	Triangulator g;
	Triangulator::Areas* areas;
	//  5              '
	//  \\             '
	//  \ \            '
	//  \  \           '
	//  \   \          '
	//   \   \         '
	//   \    \        '
	//   \     \       '
	//   \      \      '
	//    6      \     '
	//   /2------3\    '
	//  //        \\   '
	// 1           4  '
	vec2 v1(0,0);
	vec2 v2(2,1);
	vec2 v3(11,1);
	vec2 v4(12,0);
	vec2 v5(1,12);
	vec2 v6(2,2);
	CREATE_PATH(g, v1,v2,v3,v4,v5,v6);

	areas = g.findAreas();
	assert(areas->size() == 1);
	ASSERT_TYPES_AND_SIZES(g, HOLE, 6, FACE, 6);

	typedef Triangulator::AreaEdges AreaEdges;
	AreaEdges splits;
	g.monotonize_X(splits);
	ASSERT_SPLITS(splits, v3,v6);
}


#define CREATE_POLYGON(P, WO, ...) {\
	vec2 v[] = { __VA_ARGS__};\
	vec2* v_end = v + (sizeof(v)/sizeof(vec2));\
	P.reset(WO);\
	P.insert(P.begin(), v,v_end);\
	P.setClosed(true);\
}


static void test_triangulate_init() {
	// a---b-c
	// |    /
	// e---d
	vec2 a(0,1);
	vec2 b(1,1);
	vec2 c(2,1);
	vec2 d(1,0);
	vec2 e(0,0);

	Triangulator::Polygon p;
	Triangulator::SimpleTriangulation T;

	// clockwise

	CREATE_POLYGON(p, WINDING_CW,  a,b,c,d,e);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == a && T.v[1] == b && T.v[2] == c); // upper
		assert(T.v[3] == d && T.v[4] == e); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CW, b,c,d,e,a);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == a && T.v[1] == b && T.v[2] == c); // upper
		assert(T.v[3] == d && T.v[4] == e); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CW, c,d,e,a,b);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == a && T.v[1] == b && T.v[2] == c); // upper
		assert(T.v[3] == d && T.v[4] == e); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CW, d,e,a,b,c);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == a && T.v[1] == b && T.v[2] == c); // upper
		assert(T.v[3] == d && T.v[4] == e); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CW, e,a,b,c,d);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == a && T.v[1] == b && T.v[2] == c); // upper
		assert(T.v[3] == d && T.v[4] == e); // lower
		assert(T.split == 3);
	}


	// counter clockwise

	CREATE_POLYGON(p, WINDING_CCW, e,d,c,b,a);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == c && T.v[1] == b && T.v[2] == a); // upper
		assert(T.v[3] == e && T.v[4] == d); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CCW, d,c,b,a,e);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == c && T.v[1] == b && T.v[2] == a); // upper
		assert(T.v[3] == e && T.v[4] == d); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CCW, c,b,a,e,d);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == c && T.v[1] == b && T.v[2] == a); // upper
		assert(T.v[3] == e && T.v[4] == d); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CCW, b,a,e,d,c);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == c && T.v[1] == b && T.v[2] == a); // upper
		assert(T.v[3] == e && T.v[4] == d); // lower
		assert(T.split == 3);
	}

	CREATE_POLYGON(p, WINDING_CCW, a,e,d,c,b);
	{
		vector<vec2> out;
		T.init(p);
		assert(T.v[0] == c && T.v[1] == b && T.v[2] == a); // upper
		assert(T.v[3] == e && T.v[4] == d); // lower
		assert(T.split == 3);
	}


}

/** tests the output of triangulations */
struct TriangleOutputTester : vector<vec2> {
	typedef vector<vec2> super;
	typedef pair<vec2,vec2> edge_t;

	struct E_LT {
		bool operator () (const edge_t& a, const edge_t& b) const {
			return vec2_lt_x(a.first, b.first)
					|| ((a.first == b.first) && vec2_lt_x(a.second, b.second));
		}
	};

	typedef set<edge_t, E_LT> edges_t;
	Triangulator::Polygon* polygon;
	bool clearOnNext;
	edges_t edges;
	const bool PRINT_INDEX = true;
	const bool DEBUG_OUTPUT = false;
	bool reverseIndex = false;
	bool printIndexOrderAware = false;

	TriangleOutputTester() {
		polygon = NULL;
		clearOnNext = true;
	}

	TriangleOutputTester(bool printIndexOrderAware) {
		polygon = NULL;
		clearOnNext = true;
		this->printIndexOrderAware = printIndexOrderAware;
	}

	void next(Triangulator::Polygon& p) {
		polygon = &p;
		if (clearOnNext) edges.clear();
		for (unsigned i = 1; i < p.size(); i++) {
			edge(p[i-1], p[i]);
		}
		edge(p[p.size()-1], p[0]);
		reverseIndex = (p.winding == Triangulator::WO_CCW);
	}

	void assert_non_intersecting(const vec2& a, const vec2& b) const {
		for (edges_t::iterator it = edges.begin(); it != edges.end(); it++) {
			const vec2& u = it->first;
			const vec2& v = it->second;
			int result = line_segments_relationship(a,b,u,v);
			if (result == 1) {
				print(a,b,u,v);
				assert(false);
			}
		}
	}

	void edge(vec2& a, vec2& b) {
		assert_non_intersecting(a,b);
		edges.emplace(edge_t(a,b));
	}


	int index(const vec2& v) const {
		unsigned i;
		for (i = 0; i < polygon->size(); i++) {
			if ((*polygon)[i] == v) break;
		}
		if (i == polygon->size()) {
			return -1;
		} else {
			if (reverseIndex) {
				return polygon->size()-1-i;
			} else {
				return i;
			}
		}
	}

	void print(const vec2& a, const vec2& b, const vec2& u, const vec2& v) const {
		if (!DEBUG_OUTPUT) return;
		// hit if they are truly intersecting and not touching
		if (PRINT_INDEX) {
			fprintf(stderr, "%d - %d  x  %d - %d\n",
					index(a), index(b),
					index(u), index(v));
		} else {
			fprintf(stderr, "(%.f,%.f) (%.f,%.f)  x  (%.f,%.f) (%.f,%.f)\n",
					a.x,a.y, b.x,b.y,
					u.x,u.y, v.x,v.y);
		}
	}
	void print(const vec2& a, const vec2& b, const vec2& c) const {
		if (!DEBUG_OUTPUT) return;
		if (PRINT_INDEX) {
			printf("%d - %d - %d\n", index(a), index(b), index(c));
		} else {
			printf("(%.f,%.f) (%.f,%.f) (%.f,%.f)\n", a.x, a.y , b.x, b.y, c.x, c.y);
		}
	}

	void push_back(const value_type& x) {
		// hit if edges of test polygon were not added
		assert(edges.size());
		super::push_back(x);
		if (size() % 3 == 0) {
			iterator it = end();
			vec2& c = *(--it);
			vec2& b = *(--it);
			vec2& a = *(--it);
			print(a,b,c);
			edge(a,b);
			edge(b,c);
			edge(c,a);
		}
	 }
};




static void test_triangulate(Triangulator::Polygon& p) {
	TriangleOutputTester out(true);
	Triangulator::SimpleTriangulation T;

	out.next(p);
	T.triangulate(p, out);
	assert(out.size() % 3 == 0);
	assert(out.size()/3 == p.size() - 2);


	p.reverseWinding();
	out.next(p);
	T.triangulate(p, out);

}



static void test_triangulate_simple() {
	// a---b-c
	// |    /
	// e---d
	vec2 a(0,1);
	vec2 b(1,1);
	vec2 c(2,1);
	vec2 d(1,0);
	vec2 e(0,0);
	Triangulator::Polygon p;
	CREATE_POLYGON(p, WINDING_CW, a,b,c,d,e);
	test_triangulate(p);
}


static void test_triangulate_case_1() {
	//     c       '
	//     |\      '
	//  ...b d...  '
    // a`````````e '
	vec2 a(0,0);
	vec2 b(4,1);
	vec2 c(5,9);
	vec2 d(6,1);
	vec2 e(9,0);
	Triangulator::Polygon p;
	CREATE_POLYGON(p, WINDING_CW, a,b,c,d,e);
	test_triangulate(p);
}

static void test_triangulate_case_2() {
	// 1---2   4---5
	// |    \ /    |
	// |     3     |
	// |           |
	// |     8     |
	// |    / \    |
	// 0---9   7---6
	vec2 v0(0,0);
	vec2 v1(0,3);
	vec2 v2(1,3);
	vec2 v3(2,2);
	vec2 v4(3,3);
	vec2 v5(4,3);
	vec2 v6(4,0);
	vec2 v7(3,0);
	vec2 v8(2,1);
	vec2 v9(1,0);
	Triangulator::Polygon p;
	CREATE_POLYGON(p, WINDING_CW, v0,v1,v2,v3,v4,v5,v6,v7,v8,v9);
	test_triangulate(p);
}

static void test_triangulate_case_3() {
	// 1---2
	// |   |
	// |   3
	// |   |
	// 0---4
	vec2 v0(0,0);
	vec2 v1(0,2);
	vec2 v2(1,2);
	vec2 v3(1,1);
	vec2 v4(1,0);
	Triangulator::Polygon p;
	CREATE_POLYGON(p, WINDING_CW, v0,v1,v2,v3,v4);
	test_triangulate(p);
}

static void test_triangulate_case_4() {
	// 2---3
	// |   |
	// 1   4
	// |   |
	// 0---5
	vec2 v0(0,0);
	vec2 v1(0,1);
	vec2 v2(0,2);
	vec2 v3(1,2);
	vec2 v4(1,1);
	vec2 v5(1,0);
	Triangulator::Polygon p;
	CREATE_POLYGON(p, WINDING_CW, v0,v1,v2,v3,v4,v5);
	test_triangulate(p);
}


static void test_triangulate_case_5() {

	// TODO: LÖSUNGEN: entweder strikt monoton
	//                 oder evtl. erst nur diagonalen und dann polygonize


	// 1---2 4
	// |   |/|
	// |   3 |
	// |     |
	// |   6 |
	// |   |\|
	// 0---7 5
	vec2 v0(0,0);
	vec2 v1(0,3);
	vec2 v2(1,3);
	vec2 v3(1,2);
	vec2 v4(2,3);
	vec2 v5(2,0);
	vec2 v6(1,1);
	vec2 v7(1,0);
	Triangulator::Polygon p;
	CREATE_POLYGON(p, WINDING_CW, v0,v1,v2,v3,v4,v5,v6,v7);
	test_triangulate(p);
}


static int assert_memory() {
	assert(    Triangulator::Node::total_nodes_in_system == 0
			&& Triangulator::Edge::total_edges_in_system == 0
			&& Triangulator::Area::total_areas_in_system == 0
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

	test_monotonize_one_split_edge_neighbour();
	assert_memory();

	test_monotonize_one_split_node_neighbour();
	assert_memory();

	test_monotonize_one_split_vertical();
	assert_memory();

	test_monotonize_one_split_vertical_2();
	assert_memory();

	test_monotonize_one_merge_edge_neighbour();
	assert_memory();

	test_monotonize_one_merge_node_neighbour();
	assert_memory();

	test_monotonize_one_merge_vertical();
	assert_memory();

	test_monotonize_one_merge_vertical_2();
	assert_memory();

	test_monotonize_merge_then_split_same_scanline();
	assert_memory();

	test_monotonize_split_then_merge_subsequent_scanline();
	assert_memory();

	test_monotonize_merge_then_split_subsequent_scanline();
	assert_memory();

	test_monotonize_split_with_conflict_above();
	assert_memory();

	test_monotonize_split_with_conflict_below();
	assert_memory();

	test_monotonize_split_with_horizontal_conflict_above();
	assert_memory();

	test_monotonize_split_with_horizontal_conflict_below();
	assert_memory();

	test_monotonize_merge_with_conflict_above();
	assert_memory();

	test_monotonize_merge_with_conflict_below();
	assert_memory();

	test_monotonize_merge_with_horizontal_conflict_above();
	assert_memory();

	test_monotonize_merge_with_horizontal_conflict_below();
	assert_memory();

	// TODO: test for polygonize only

	// TODO: test triangulate only

	test_triangulate_init();
	assert_memory();

	test_triangulate_simple();
	assert_memory();

	test_triangulate_case_1();
	assert_memory();

	test_triangulate_case_2();
	assert_memory();

	test_triangulate_case_3();
	assert_memory();

	test_triangulate_case_4();
	assert_memory();

	test_triangulate_case_5();
	assert_memory();
	return 0;
}

}; // struct Test_Triangulator


#endif /* MATH_TEST_AREAGRAPH_H_ */
