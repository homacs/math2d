/*
 * test_bezier.h
 *
 *  Created on: 26 Jul 2019
 *      Author: homac
 */

#ifndef TEST_BEZIER_H_
#define TEST_BEZIER_H_


#include <assert.h>


#include <math2d/bezier.h>


static inline void testsub_bezier_split_merge(double f, vec2 const & v0, vec2 const & v1, vec2 const & v2, vec2 const & v3) {
	vec2 v; // test variable
	vec2 p0 = v0;
	vec2 p1 = v1;
	vec2 p2 = v2;
	vec2 p3 = v3;

	vec2 q0;
	vec2 q1;
	vec2 q2;
	vec2 q3;


	// tolerance = 0.0001 % of its approximated total length
	float tolerance = glm::length(p1-p0) + glm::length(p2-p1) + glm::length(p3-p2);
	tolerance = 0.000001 * tolerance;

	bezier_split_highp(f, p0, p1, p2, p3, q0, q1, q2, q3);
	assert(p0 == v0);
	assert(p3 == q0);
	assert(q3 == v3);

	bezier_point_highp(f, v0,v1,v2,v3, v);
	assert(about_equal(p3, v, tolerance));


	bool merged = bezier_merge_highp(p0, p1, p2, p3, q0, q1, q2, q3, tolerance);
	assert(merged);

	assert(about_equal(p0, v0, tolerance));
	assert(about_equal(p1, v1, tolerance));
	assert(about_equal(p2, v2, tolerance));
	assert(about_equal(p3, v3, tolerance));
}

static inline void testsub_bezier_split_merge_bothway(vec2 const & v0, vec2 const & v1, vec2 const & v2, vec2 const & v3) {
	float f = 1.0/2;
	testsub_bezier_split_merge(f, v0,v1,v2,v3);
	testsub_bezier_split_merge(f, v3,v2,v1,v0);

	f = 1.0/4.0;
	testsub_bezier_split_merge(f, v0,v1,v2,v3);
	testsub_bezier_split_merge(f, v3,v2,v1,v0);

	f = 1.0/3.0;
	testsub_bezier_split_merge(f, v0,v1,v2,v3);
	testsub_bezier_split_merge(f, v3,v2,v1,v0);

	f = 0.1;
	testsub_bezier_split_merge(f, v0,v1,v2,v3);
	testsub_bezier_split_merge(f, v3,v2,v1,v0);

	f = 0.001;
	testsub_bezier_split_merge(f, v0,v1,v2,v3);
	testsub_bezier_split_merge(f, v3,v2,v1,v0);

	f = 0.0001;
	testsub_bezier_split_merge(f, v0,v1,v2,v3);
	testsub_bezier_split_merge(f, v3,v2,v1,v0);

}

static inline void test_bezier_split_merge() {

	// small

	// simply concave
	testsub_bezier_split_merge_bothway(vec2(0,0),vec2(0,1),vec2(1,1),vec2(1,0));
	// mit sattelpunkt
	testsub_bezier_split_merge_bothway(vec2(0,1),vec2(0,0),vec2(1,1),vec2(1,0));
	// closed loop
	testsub_bezier_split_merge_bothway(vec2(1,-1),vec2(0,0),vec2(0,0),vec2(1,1));

	// medium

	// simply concave
	testsub_bezier_split_merge_bothway(vec2(0,0),vec2(0,100),vec2(100,100),vec2(100,0));
	// mit sattelpunkt
	testsub_bezier_split_merge_bothway(vec2(0,100),vec2(0,0),vec2(100,100),vec2(100,0));
	// closed loop
	testsub_bezier_split_merge_bothway(vec2(100,-100),vec2(0,0),vec2(0,0),vec2(100,100));

	// large

	// simply concave
	testsub_bezier_split_merge_bothway(vec2(0,0),vec2(0,10000),vec2(10000,10000),vec2(10000,0));
	// mit sattelpunkt
	testsub_bezier_split_merge_bothway(vec2(0,10000),vec2(0,0),vec2(10000,10000),vec2(10000,0));
	// closed loop
	testsub_bezier_split_merge_bothway(vec2(10000,-10000),vec2(0,0),vec2(0,0),vec2(10000,10000));

}


static inline void test_bezier_inflection_point() {
	bool turning;
	float t;


	// with inflection point
	turning = bezier_inflection_point(vec2(0,7), vec2(5.5,6.5), vec2(-2,1.5), vec2(3.0,1.5), t);
	assert(turning);

	// with inflection point
	turning = bezier_inflection_point(vec2(1,0), vec2(1,7), vec2(3,-4), vec2(3,3),t);
	assert(turning);

	// with inflection point
	turning = bezier_inflection_point(vec2(1,0), vec2(1,100), vec2(3,-100), vec2(3,3), t);
	assert(turning);

	// straight line, equi-distant points
	turning = bezier_inflection_point(vec2(1,0), vec2(2,1), vec2(3,2), vec2(4,3), t);
	assert(!turning);

	// straight line, points in increasing distances (constant acceleration)
	turning = bezier_inflection_point(vec2(1,0), vec2(2,1), vec2(4,3), vec2(8,7), t);
	assert(!turning);

	// no inflection point
	turning = bezier_inflection_point(vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0), t);
	assert(!turning);


}


static inline void test_bezier_point_closest_point() {
	float tolerance = 0.000001;
	vec2 result_v;
	float t;
	// with inflection point
	t = bezier_point_closest_point_t(vec2(0,7), vec2(5.5,6.5), vec2(-2,1.5), vec2(3.0,1.5), vec2(0,7), tolerance, result_v);
	assert(about_equal((float)t, 0.f, tolerance));

	// with inflection point
	t = bezier_point_closest_point_t(vec2(0,7), vec2(5.5,6.5), vec2(-2,1.5), vec2(3.0,1.5), vec2(3.0,1.5), tolerance, result_v);
	assert(about_equal((float)t, 1.f, tolerance));

	// with inflection point
	t = bezier_point_closest_point_t(vec2(0,7), vec2(5.5,6.5), vec2(-2,1.5), vec2(3.0,1.5), vec2(1.0,4.5), tolerance, result_v);
	assert(about_equal((float)t, 0.484173596f, tolerance));

}


static inline void test_bezier_line_intersections_iterative() {

	float tolerance = 0.000001f;
	vec2 p[3];
	int count;

	// no turning point, two intersections
	count = bezier_line_segment_intersections_iterative(
			vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 2);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));

	// no turning point, one intersection
	count = bezier_line_segment_intersections_iterative(
			vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0),
			vec2(0,1), vec2(2,1), tolerance, p);
	assert(count == 1);
	assert(about_equal(p[0].y, 1, tolerance));


	// with turning point, three intersections
	count = bezier_line_segment_intersections_iterative(
			vec2(1,0), vec2(1,7), vec2(3,-4), vec2(3,3),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 3);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));
	assert(about_equal(p[2].y, 1, tolerance));

	// with turning point, three intersections
	count = bezier_line_segment_intersections_iterative(
			vec2(1,0), vec2(1,100), vec2(3,-100), vec2(3,3),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 3);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));
	assert(about_equal(p[2].y, 1, tolerance));



	count = bezier_line_segment_intersections_iterative(
			vec2(1,0), vec2(1,1000000), vec2(3,-1000000), vec2(3,3),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 3);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));
	assert(about_equal(p[2].y, 1, tolerance));
}



static inline void test_bezier_extrema() {
	dvec2 t_min, t_max;
	bool found;
	// no turning point
	found = bezier_extrema(vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0),	t_min, t_max);
	assert(found == true);
	assert(t_min.x == INFINITY);
	assert(t_max.x == INFINITY);
	assert(t_min.y == INFINITY);
	assert(t_max.y == 0.5);

	// with turning point
	found = bezier_extrema(vec2(1,0), vec2(1,7), vec2(3,-4), vec2(3,3), t_min, t_max);
	assert(found == true);
	assert(t_min.x == INFINITY);
	assert(t_max.x == INFINITY);
	assert(t_min.y != INFINITY);
	assert(t_max.y != INFINITY);

	// loop
	found = bezier_extrema(vec2(0,0), vec2(-3,1), vec2(3,1), vec2(0,0), t_min, t_max);
	assert(found == true);
	assert(t_min.x != INFINITY);
	assert(t_max.x != INFINITY);
	assert(t_min.y == INFINITY);
	assert(t_max.y != INFINITY);

	// diagonal loop
	found = bezier_extrema(vec2(0,0), vec2(-2,-1), vec2(2,1), vec2(0,0), t_min, t_max);
	assert(found == true);
	assert(t_min.x != INFINITY);
	assert(t_max.x != INFINITY);
	assert(t_min.y != INFINITY);
	assert(t_max.y != INFINITY);

	// saddle point
	found = bezier_extrema(vec2(0,0), vec2(2,0), vec2(0,2), vec2(2,2), t_min, t_max);
	assert(found == false);
	assert(t_min.x == INFINITY);
	assert(t_max.x == INFINITY);
	assert(t_min.y == INFINITY);
	assert(t_max.y == INFINITY);
}

static inline void test_bezier_line_intersections() {

	float tolerance = 0.000001f;
	vec2 p[3];
	int count;

	// no turning point, two intersections
	count = bezier_line_segment_intersections(
			vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 2);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));

	// no turning point, one intersection
	count = bezier_line_segment_intersections(
			vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0),
			vec2(0,1), vec2(2,1), tolerance, p);
	assert(count == 1);
	assert(about_equal(p[0].y, 1, tolerance));


	// with turning point, three intersections
	count = bezier_line_segment_intersections(
			vec2(1,0), vec2(1,7), vec2(3,-4), vec2(3,3),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 3);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));
	assert(about_equal(p[2].y, 1, tolerance));

	// with turning point, three intersections
	count = bezier_line_segment_intersections(
			vec2(1,0), vec2(1,100), vec2(3,-100), vec2(3,3),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 3);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));
	assert(about_equal(p[2].y, 1, tolerance));



	count = bezier_line_segment_intersections(
			vec2(1,0), vec2(1,1000000), vec2(3,-1000000), vec2(3,3),
			vec2(0,1), vec2(4,1), tolerance, p);
	assert(count == 3);
	assert(about_equal(p[0].y, 1, tolerance));
	assert(about_equal(p[1].y, 1, tolerance));
	assert(about_equal(p[2].y, 1, tolerance));
}




static inline void test_bezier_bezier_intersections() {

	float tolerance = 0.000001f;
	vec2 p[9];
	double t_p[9];
	double t_q[9];
	int count;

	// no turning point, same bezier mirrored
	// two intersections
	count = bezier_bezier_intersections_t(
			vec2(1,0), vec2(1,3), vec2(3,3), vec2(3,0),
			vec2(1,3), vec2(1,0), vec2(3,0), vec2(3,3),
			tolerance, t_p, t_q);
	assert(count == 2);
//	assert(about_equal(p[0].y, 1, tolerance));
//	assert(about_equal(p[1].y, 1, tolerance));

}


void test_bezier_all() {
	test_bezier_extrema();
	test_bezier_inflection_point();
	test_bezier_point_closest_point();
	test_bezier_split_merge();
	test_bezier_line_intersections_iterative();
	test_bezier_line_intersections();
	test_bezier_bezier_intersections();
}


#endif /* TEST_BEZIER_H_ */
