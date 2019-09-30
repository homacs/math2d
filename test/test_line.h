/*
 * test_math_line.h
 *
 *  Created on: 17 Jun 2019
 *      Author: homac
 */

#ifndef TEST_MATH_LINE_H_

#define TEST_MATH_LINE_H_



#include <stdio.h>
#include <stdlib.h>

#include <glm/gtx/vector_angle.hpp>
using namespace glm;


#include "math2d/line.h"
using namespace math2d;

#include "perf_clock.h"



void test_performance_pseudo_orientation() {
	float angle = 0;

	float rand_fac = ((float)1000.0)/RAND_MAX;
	vec2 v;
	vec2 e_y(1,0);
	perf_time_t start;
	perf_time_t end;
	int iterations = 100000000;

	// calibration: measuring loop time only
	perf_clock(start);
	for (int i = 0; i < iterations; i++) {
		v.x = rand()*rand_fac;
		v.y = rand()*rand_fac;
	}
	perf_clock(end);
	if (v.x == v.y) printf("%s", "");
	perf_time_t LOOP_TIME = end-start;
	printf("loop: %fs\n", LOOP_TIME.secs_nsecs());

	perf_clock(start);
	for (int i = 0; i < iterations; i++) {
		v.x = rand()*rand_fac;
		v.y = rand()*rand_fac;
//		angle = orientedAngle(v, e_y);
		angle = pseudo_orientation_e2(v);
	}

	perf_clock(end);
	perf_time_t total = end-start;
	printf("%03.2fs total, %03.2fns/call\n", total.secs_nsecs(), double(total-LOOP_TIME)/iterations );

	// use angle to prevent removable by optimisation
	if (angle == 0 && v.x == 0 && v.y == 0)  { printf("%s", ""); }

}


int test_sequence__pseudo_orientation() {

	// test if angle is continuous in range [-179.9999,180]
#define OUTPUT_PLOT 0
#if OUTPUT_PLOT
	FILE* fout = fopen("/tmp/angle.txt","w");
#endif

	float step = 0.000001f;
	float PI = glm::pi<double>();
	// NOTE: rotateZ(v, (-)180) actually turns v by an angle >(<) (-)180!


	vec3 e_y(0,1,0);
	vec3 v;
	vec3 prev_v;
	double angle = 0;
	double prev_angle = 0;
	float prev_alpha = -PI+step;
	v = glm::rotateZ(e_y, -prev_alpha);
	prev_angle = pseudo_orientation_e2(vec2(v.x,v.y));
	for (float alpha = prev_alpha+step; alpha < PI; alpha+=step) {
		// z axis interpreted in opposite direction (towards observer, as in OpenGL)
		v = glm::rotateZ(e_y, -alpha);
		angle = pseudo_orientation_e2(vec2(v.x,v.y));
		if (prev_angle >= angle) {
			fprintf(stderr, "%f->%f: %f -> %f ", prev_alpha, alpha, prev_angle, angle);
			fprintf(stderr, "(%f,%f) -> (%f,%f)\n", prev_v.x,prev_v.y, v.x,v.y);
			//assert(prev_angle < angle);
		}
		prev_alpha = alpha;
		prev_angle = angle;
		prev_v = v;
#if OUTPUT_PLOT
		fprintf(fout, "%2.8f %2.8f %2.8f %2.8f\n",alpha, angle, v.x, v.y);
		fflush(fout);
#endif
	}
#if OUTPUT_PLOT
	fclose(fout);
#endif

#undef OUTPUT_PLOT
	return 0;

}

void test_plot_pseudo_orientation() {




	FILE* fout = fopen("/tmp/angle.txt", "w");

	float angle = 0;
	float PI = glm::pi<float>();
	vec3 e_y(0,1,0);
	// NOTE: rotateZ(v, (-)180) actually turns v by an angle >(<) (-)180!
	for (int alpha = -179; alpha < 180; alpha++) {
		float alpha_radian = float(alpha) * 2.0 * PI/ 360.0 ;
		// z axis interpreted in opposite direction (towards observer, as in OpenGL)
		vec3 v = glm::rotateZ(e_y, -alpha_radian);
		angle = pseudo_orientation_e2(vec2(v.x,v.y));
		fprintf(fout, "%2.5f %2.5f %2.5f %2.5f\n", alpha_radian, angle, v.x, v.y);
	}

	fclose(fout);
}

/**
 * A segment is defined by a straight line between two points a and b.
 */
struct line_segment_t {
	glm::vec2 a;
	glm::vec2 b;
	line_segment_t(const line_segment_t& s) {
		this->a = s.a;
		this->b = s.b;
	}
	line_segment_t(const glm::vec2& a, const glm::vec2& b) {
		this->a = a;
		this->b = b;
	}
	line_segment_t reverse() {
		return line_segment_t(b,a);
	}
};



/** swaps two line segments with each other */
void testsub_swap(line_segment_t& s1, line_segment_t& s2) {
	line_segment_t tmp(s1);
	s1 = s2;
	s2 = tmp;
}




/**
 * helper function for test_sequence__segment_intersection()
 * tests all valid permutations of given segment vectors
 * for intersections.
 *
 */
template <typename T>
void testsub_segment_intersections(
		T m1, T m2,
		T n1, T n2,
		T u1, T u2,
		T v1, T v2,
		T p1, T p2,
		bool result)
{

	vec2 m(m1,m2);
	vec2 n(n1,n2);

	vec2 u(u1,u2);
	vec2 v(v1,v2);

	vec2 p(p1,p2);

	vec2 q;

	static struct Equal{
		const bool single;
		float_mantissa_mask_t precision;
		Equal() : single(true) {
			// remove 6 of 23 bits from float mantissa before comparison
			precision = float_mantissa_mask(6);
		}

		float trunc(float a) const {
			return float_mantissa_trunc(a, precision);
		}

		bool isEqual(const T a, const T b) const {
			if (single) {
				// single precision
				return trunc(a) == trunc(b);
			} else {
				// double precision
				return float(a) == float(b);
			}
		}

		bool operator() (const vec2& a, const vec2& b) const {
			return isEqual(a.x,b.x) && isEqual(a.y,b.y);
		}
		bool operator() (const T a, const T b) const {
			return isEqual(a,b);
		}
	} equal;

	line_segment_t s[] = {
			line_segment_t(m,n),
			line_segment_t(u,v),
	};

	float dist;
	// test all permutations of given test vectors
	for (int i = 0; i < 2; i++) {
		for (int k = 0; k < 2; k++) {
			bool intersects = line_segment_intersection(s[0].a,s[0].b,s[1].a,s[1].b,q, dist);
			assert (intersects == result);
			if (result) {
				assert (equal(q, p));
				assert (equal(dist, length(s[0].a - p)/length(s[0].b-s[0].a)));
			}

			s[k] = s[k].reverse();
		}

		bool intersects = line_segment_intersection(s[0].a,s[0].b,s[1].a,s[1].b,q, dist);
		assert (intersects == result);
		if (result) assert (equal(q, p));


		testsub_swap(s[0],s[1]);
	}

}

int test_sequence__segment_intersection() {
	// parallel and on top of each other
	testsub_segment_intersections(0.,0.,	1.,1.,		0.,0.,	1.,1.,		0.5,0.5, false);

	// obvious intersections
	testsub_segment_intersections(1.,0.,	0.,1.,		0.,0.,	1.,1.,		0.5,0.5, true);
	testsub_segment_intersections(1000.,0.,	0.,1000.,	0.,0.,	1000.,1000.,500.,500., true);
	testsub_segment_intersections(1.,0.,	0.,-1.,		0.,0.,	1.,-1.,		0.5,-0.5, true);
	testsub_segment_intersections(-1.,0.,	0.,-1.,		0.,0.,	-1.,-1.,	-0.5,-0.5, true);
	testsub_segment_intersections(-1.,0.,	0.,1.,		0.,0.,	-1.,1.,		-0.5,0.5, true);
	testsub_segment_intersections(-1.,0.,	+1.,0.,		0.,-1.,	0.,+1.,		0.,0., true);

	// lines intersect but segments do not
	testsub_segment_intersections(5.,0.,	0.,5.,		0.,0.,	2.,2.,		1.5,1.5, false);


	// segment a touches segment b in one point
	testsub_segment_intersections(4.,0.,	0.,4.,		0.,0.,	2.,2.,		2.0,2.0, true);



	testsub_segment_intersections(4461.66846,3018.90674,
			57557.0859,17640.0469,
			59054.2422,330.298187,
			20656.7598,20389.6426,
			36840.2812,11935.1602, true);

	return 0;
}




int test_sequence__segment_intersection_highp() {
	vec2 p;
	double s;
	line_segment_intersection_highp(
			vec2(4461.66846,3018.90674),
			vec2(57557.0859,17640.0469),
			vec2(59054.2422,330.298187),
			vec2(20656.7598,20389.6426),
			p, s);

	assert(p == vec2(36840.2812,11935.1611));


	line_segment_intersection_highp(
			vec2(4461.66846,3018.90674),
			vec2(57557.0859,17640.0469),
			vec2(59054.2422,330.298187),
			vec2(20656.7578,20389.6426),
			p, s);



	return 0;
}




int test_sequence__ray_left_intersects_segment_where_m_below() {
	//
	bool intersects;
	vec2 s,m;


	s = vec2(1,1); // ray start
	m = vec2(0,0); // below ray

	// above and crossing
	intersects = ray_left_intersects_segment_where_m_below(s, m, vec2(0,2));
	assert(intersects == true);

	// touching s
	intersects = ray_left_intersects_segment_where_m_below(s, m, vec2(2,2));
	assert(intersects == true);

	// above not crossing
	intersects = ray_left_intersects_segment_where_m_below(s, m, vec2(3,2));
	assert(intersects == false);



	s = vec2(1,1); // ray start
	m = vec2(0.5,0.5); // below ray

	// above and crossing
	intersects = ray_left_intersects_segment_where_m_below(s, m, vec2(0,2));
	assert(intersects == true);

	// touching s
	intersects = ray_left_intersects_segment_where_m_below(s, m, vec2(2,2));
	assert(intersects == true);

	// above not crossing
	intersects = ray_left_intersects_segment_where_m_below(s, m, vec2(3,2));
	assert(intersects == false);

	return 0;

}



int testsub_line_contains_point(float mx, float my, float nx, float ny, float px, float py, float expected_s, bool expected_result) {
	vec2 m(mx,my);
	vec2 n(nx,ny);
	vec2 p(px,py);
	float s;
	bool result;
	result = line_contains_point(
			m,
			n,
			p,
			s);
	assert(result == expected_result);
	if (result) assert(s == expected_s);

	result = line_contains_point(
			n,
			m,
			p,
			s);
	assert(result == expected_result);
	if (result) assert(s == (1-expected_s));


	return 0;
}


int testsub_line_contains_point_permute(float mx, float my, float nx, float ny, float px, float py, float expected_s, bool expected_result) {
	testsub_line_contains_point(
			mx,my,
			nx,ny,
			px,py,
			expected_s, expected_result);

	testsub_line_contains_point(
			mx,-my,
			nx,-ny,
			px,-py,
			expected_s, expected_result);

	testsub_line_contains_point(
			-mx,my,
			-nx,ny,
			-px,py,
			expected_s, expected_result);

	testsub_line_contains_point(
			-mx,-my,
			-nx,-ny,
			-px,-py,
			expected_s, expected_result);


	return 0;
}

int test_sequence__line_contains_point() {
	//
	//   ON THE LINE
	//

	// on line
	testsub_line_contains_point_permute(0,0, 4,2, 2,1, 0.5, true);
	testsub_line_contains_point_permute(4,2, 8,4, 6,3, 0.5, true);

	// on top of one of the points
	testsub_line_contains_point_permute(0,0, 4,2, 0,0, 0.0, true);
	testsub_line_contains_point_permute(0,0, 4,2, 4,2, 1.0, true);
	testsub_line_contains_point_permute(4,2, 8,4, 4,2, 0.0, true);
	testsub_line_contains_point_permute(4,2, 8,4, 8,4, 1.0, true);

	// horizontal lines
	testsub_line_contains_point_permute(0,0, 4,0, 2,0, 0.5, true);
	testsub_line_contains_point_permute(2,0, 4,0, 3,0, 0.5, true);
	testsub_line_contains_point_permute(0,0, 4,0, 0,0, 0.0, true);
	testsub_line_contains_point_permute(0,0, 4,0, 4,0, 1.0, true);
	testsub_line_contains_point_permute(2,0, 4,0, 2,0, 0.0, true);
	testsub_line_contains_point_permute(2,0, 4,0, 4,0, 1.0, true);

	// vertical lines
	testsub_line_contains_point_permute(0,0, 0,4, 0,2, 0.5, true);
	testsub_line_contains_point_permute(0,2, 0,4, 0,3, 0.5, true);
	testsub_line_contains_point_permute(0,0, 0,4, 0,0, 0.0, true);
	testsub_line_contains_point_permute(0,0, 0,4, 0,4, 1.0, true);
	testsub_line_contains_point_permute(0,2, 0,4, 0,2, 0.0, true);
	testsub_line_contains_point_permute(0,2, 0,4, 0,4, 1.0, true);


	//
	// OFF THAT LINE
	//

	// above
	testsub_line_contains_point_permute(0,0, 4,2, 2,2, 0, false);
	testsub_line_contains_point_permute(4,2, 8,4, 6,4, 0, false);
	// below
	testsub_line_contains_point_permute(0,0, 4,2, 2,0, 0, false);
	testsub_line_contains_point_permute(4,2, 8,4, 6,2, 0, false);

	// horizontal
	testsub_line_contains_point_permute(0,0, 4,0, 2,1, 0, false);
	testsub_line_contains_point_permute(0,0, 4,0, 2,-1, 0, false);
	testsub_line_contains_point_permute(2,0, 4,0, 3,1, 0, false);
	testsub_line_contains_point_permute(2,0, 4,0, 3,-1, 0, false);

	// vertical
	testsub_line_contains_point_permute(0,0, 0,4, 1,2, 0.5, false);
	testsub_line_contains_point_permute(0,0, 0,4, -1,2, 0.5, false);
	testsub_line_contains_point_permute(0,2, 0,4, 1,3, 0.5, false);
	testsub_line_contains_point_permute(0,2, 0,4, -1,3, 0.5, false);


	return 0;
}


int testsub_line_segment_superpose_permute(
		float mx, float my,
		float nx, float ny,
		float ux, float uy,
		float vx, float vy,
		float px, float py,
		float qx, float qy,
		bool expected_result)
{
	vec2 m(mx, my);
	vec2 n(nx,ny);
	vec2 u(ux,uy);
	vec2 v(vx,vy);

	vec2 expected_p(px,py);
	vec2 expected_q(qx,qy);

	vec2 p;
	vec2 q;
	bool result;

	//
	// this testsub tests all permutations by
	// swapping m,n and u,v
	//


	//
	// m,n , u,v
	//

	result = line_segment_superpose(m,n,u,v,p,q);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);

	// opposite direction -> swapped p and q
	result = line_segment_superpose(n,m,u,v,q,p);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);

	result = line_segment_superpose(m,n,v,u,p,q);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);

	result = line_segment_superpose(n,m,v,u,q,p);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);


	//
	//  u,v  ,  m,n
	//

	// swap(m,u)
	vec2 tmp = m;
	m = u;
	u = tmp;
	// swap (n,v)
	tmp = n;
	n = v;
	v = tmp;

	result = line_segment_superpose(m,n,u,v,p,q);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);

	// opposite direction -> swapped p and q
	result = line_segment_superpose(n,m,u,v,q,p);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);

	result = line_segment_superpose(m,n,v,u,p,q);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);

	result = line_segment_superpose(n,m,v,u,q,p);
	assert(result == expected_result);
	if (result) assert(expected_p == p && expected_q == q);



	return 0;
}

int test_sequence__line_segment_superpose() {

	testsub_line_segment_superpose_permute(0,0, 2,2, 4,4, 6,6, 0,0, 0,0, false);
	testsub_line_segment_superpose_permute(0,0, 4,4, 2,2, 6,6, 2,2, 4,4, true);
	testsub_line_segment_superpose_permute(0,0, 6,6, 2,2, 4,4, 2,2, 4,4, true);
	testsub_line_segment_superpose_permute(2,2, 4,4, 0,0, 6,6, 2,2, 4,4, true);
	testsub_line_segment_superpose_permute(2,2, 6,6, 0,0, 4,4, 2,2, 4,4, true);
	testsub_line_segment_superpose_permute(4,4, 6,6, 0,0, 2,2, 0,0, 0,0, false);

	return 0;
}


int test_line_all() {
#ifdef NDEBUG
	fprintf(stderr, "Compile *without* -dNDEBUG to run tests");
#endif
	test_sequence__ray_left_intersects_segment_where_m_below();
	test_sequence__segment_intersection();
	test_sequence__line_contains_point();
	test_sequence__line_segment_superpose();
	test_sequence__pseudo_orientation();
	return 0;
}


#endif /* TEST_MATH_LINE_H_ */
