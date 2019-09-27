/*
 * math_line.cpp
 *
 *  Created on: 17 Jun 2019
 *      Author: homac
 */

#include "float-utils.h"

#include <assert.h>
#include <math2d/line.h>
#include <math2d_config.h>




using namespace glm;


namespace math2d {


bool ray_left_intersects_segment_where_m_below(const vec2& s, const vec2& m, const vec2& n) {
	/*
		If m is below and n above the ray,
		then the angle from the direction
		vector (n-m) to the difference
		vector (s-m) must be greater or equal zero.
		Otherwise the segment crosses behind s.
	 */

	// m is supposed to be below the ray
	assert (m.y <= s.y);


	bool intersects = false;

	if (n.y >= s.y) {
		// m below and n above
		// angle from (n-m) to (s-m)
		float pseudo_angle = cross_z((s-m), (n-m));
		intersects = (pseudo_angle >= 0);
	}

	return intersects;
}

/**
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_x is set to intersection point p(p_x,s_y)
 */
int scanline_horizontal_intersection_m_below(const float s_y, vec2 m, vec2 n, float& p_x) {

	// m is supposed to be below the ray
	assert (m.y <= s_y);


	if (m.y > s_y || n.y < s_y) return -1;
	// ==>   m.y <= s_y   and    n.y >= s_y
	if (m.y == s_y && m.y == n.y) return 0;
	// ==>   m_y <= s_y and n.y >= s_y   but   m_y != n_y

	// calculate intersection point p(p_x, s_y)
	//
	// p_x = m_x + t*(n_x-m_x)
	// p_y = m_y + t*(n_y-m_y)
	// s_y = p_y
	// ==>
	// s_y = m_y + t*(n_y-m_y)
	// t = (s_y - m_y)/(n_y - m_y)
	// ==>
	// p_x = m_x + (s_y - m_y) * (n_x - m_x) / (n_y - m_y)

	p_x = m.x + (s_y - m.y) * (n.x - m.x) / (n.y - m.y);
	return 1;
}
/**
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_x is set to intersection point p(p_x,s_y)
 */
int scanline_horizontal_intersection(const float s_y, vec2 m, vec2 n, float& p_x) {
	if (m.y < s_y) return scanline_vertical_intersection_m_left(s_y, m, n, p_x);
	else return scanline_vertical_intersection_m_left(s_y, n, m, p_x);
}

/**
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_y is set to intersection point p(s_x,p_y)
 */
int scanline_vertical_intersection_m_left(const float s_x, vec2 m, vec2 n, float& p_y) {

	// m is supposed to be left of the vertical ray
	assert (m.x <= s_x);


	if (m.x > s_x || n.x < s_x) return -1;
	// ==>   m.x <= s_x   and    n.x >= s_x
	if (m.x == s_x && m.x == n.x) return 0;
	// ==>   m_x <= s_x and n.x >= s_x   but   m_x != n_x

	// calculate intersection point p(s_x, p_y)
	p_y = m.y + (s_x - m.x) * (n.y - m.y) / (n.x - m.x);
	return 1;
}

/**
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_y is set to intersection point p(s_x,p_y)
 */
int scanline_vertical_intersection(const float s_x, vec2 m, vec2 n, float& p_y) {
	if (m.x < s_x) return scanline_vertical_intersection_m_left(s_x, m, n, p_y);
	else return scanline_vertical_intersection_m_left(s_x, n, m, p_y);
}



bool line_intersection(const vec2& m, const vec2& n, const vec2& u, const vec2& v, vec2& p, float& s) {

	/*
		Given the equations
		s1 = m + (n-m) * s , 0 <= s <= 1
		s2 = u + (v-u) * t , 0 <= t <= 1
		We search for the intersection p,
		of both line equations, such that
			s1(s) = s2(t)   (=p)
		We solve the equation system to get
		s AND t, and test if both satisfy
		their boundary condition (meaning:
		p is on both segments).
	*/


	// pre calculate direction vectors
	vec2 v_u(v-u);
	vec2 n_m(n-m);
	vec2 u_m(u-m);

	// cross product of the direction vectors
	// of both lines. If 0 --> parallel.
	float denominator = cross_z(v_u,n_m);
	if (unlikely(denominator == 0)) {
		// lines are parallel to each other
		// or one of the segments is a point (length = 0)

		// Both cases are defined as non-intersecting.
		return false;
	} else {

		s = cross_z(v_u, u_m) / denominator;

		p = m + s * (n_m);

		return true;
	}
}



bool line_intersection_highp(const dvec2& m, const dvec2& n, const dvec2& u, const dvec2& v, dvec2& p, double& s) {
	// TODO: arithmetic precision guards
	// We do not test, whether the result is
	// incorrect due to the limited precision of float.
	// TODO: optimise
	// this can be further optimised,
	// - either by not using vec2 internally
	// - or not using const vec2& params

	/*
		Given the equations
		s1 = m + (n-m) * s , 0 <= s <= 1
		s2 = u + (v-u) * t , 0 <= t <= 1
		We search for the intersection p,
		of both line equations, such that
			s1(s) = s2(t)   (=p)
		We solve the equation system to get
		s AND t, and test if both satisfy
		their boundary condition (meaning:
		p is on both segments).
	*/


	// pre calculate direction vectors
#define vec_minus(a,b) a.x = a.x - b.x; a.y = a.y - b.y
	highp_dvec2 v_u(v);
	vec_minus(v_u,u);
	highp_dvec2 n_m(n);
	vec_minus(n_m,m);
	highp_dvec2 u_m(u);
	vec_minus(u_m,m);
#undef vec_minus
	// cross product of the direction vectors
	// of both lines. If 0 --> parallel.
	double denominator = cross_z(v_u,n_m);
	if (unlikely(denominator == 0)) {
		// lines are parallel to each other
		// or one of the segments is a point (length = 0)

		// Both cases are defined as non-intersecting.
		return false;
	} else {

		s = cross_z(v_u, u_m) / denominator;

		// g2(t):
		p.x = double(m.x) + s * n_m.x;
		p.y = double(m.y) + s * n_m.y;
		return true;
	}
}

/**
 *
 * Test if a given point p is on the line
 *
 *    l(s) = m + s (n - m) .
 *
 * Function returns true, if p is contained in line.
 * Function also returns the factor s, which represents the
 * distance from m to p, relative to the direction vector (n-m).
 * s gives information about order of these nodes in the following manner:
 *   - s < 0:          p -> m -> n
 *   - s = 0:          p = m
 *   - s = 1:          p = n
 *   - s > 1:          m -> n -> p
 *   - s > 0 && s < 1: m -> p -> n
 *
 *
 *
 * @param m start of line
 * @param n end of line
 * @param p point to be tested
 * @param s Output: distance from m to p relative to direction and length of line.
 * @return true, only if p on m->n
 */
bool line_contains_point(const vec2& m, const vec2& n, const vec2& p, float& s) {
	assert (m != n);


	vec2 r = n-m;
	// orthogonal vector
	vec2 o;
	o.x = -r.y;
	o.y = r.x;

	vec2 q;
	// intersection between m->n and p->o
	if (line_intersection(m,n,p, p+o, q, s)) {
		// q is the intersection point
		// if q == p then p on m -> n
		return p == q;
	}

	return false;
}

bool line_contains_point_highp(const dvec2& m, const dvec2& n, const dvec2& p, double& s, float tolerance) {
	assert(tolerance >= 0);
	// not a line?
	assert (m != n);

	/*
	 * We actually calculate the intersection point q of the given line
	 * and an orthogonal line through the point.
	 * Then we compare p and q under the given tolerance.
	 */

	dvec2 r = n-m;
	// orthogonal vector
	dvec2 o;
	o.x = -r.y;
	o.y = r.x;

	dvec2 q;
	// intersection between m->n and p->o
	if (line_intersection_highp(m,n,p, p+o, q, s)) {
		// q is the intersection point
		// if q == p then p on m -> n
		return about_equal(p,q, tolerance);
	}

	return false;
}


/**
 * For *parallel* *lines* only!
 *
 *
 * This function is meant to be called after
 * line_segment_intersection() returns false.
 *
 * Two lines superpose (overlap) each other
 * if they are parallel to each other
 * and share a common section [p,q].
 *
 * @return true if both lines have a common section. Section is stored in p and k
 */
bool line_segment_superpose(const vec2& m, const vec2& n, const vec2& u, const vec2& v, vec2& p, vec2& q) {
	// function does not accept points
	assert(m != n && u != v);
	// function does not accept non-parallel lines
	assert(0 == cross_z(v-u,n-m));


	bool superposing = false;

	p = m;
	q = n;


	float r_m = 0;
	float r_n = 1;
	float r_u;
	float r_v;


	if (unlikely(
			   line_contains_point(p, q, u, r_u)
			&& line_contains_point(p, q, v, r_v)
		))
	{
		// search the two middle points

		// We want to prevent inaccuracies due to floating point computation.
		// thus we cannot use r_i to calculate the points and have to test
		// the order of points manually.

		// TODO: optimise (that is a lot of branches)
		float r_p = r_m;
		float r_q = r_n;

		if (r_u < r_v) // m --> n and u --> v point in the same direction
		{
			// optimise: swap u and v instead of branching

			// now it comes down to range overlapping
			superposing = (r_m <= r_v && r_n >= r_u);
			// optimise: end here if !superposing?
			// optimise: move decision in branches instead?
			if (!superposing) return false;


			if (r_p < r_u) {
				p = u;
				if (r_q > r_v) {
					q = v;
				}
			}
			else /* r_p >= r_u */
			{
				if (r_q > r_v) {
					q = v;
				}
			}
		}
		else /*r_v < r_u */
		{
			// now it comes down to range overlapping
			superposing = (r_m <= r_u && r_n >= r_v);
			if (!superposing) return false;

			if (r_p < r_v) { // m -> v
				p = v;
				if (r_q > r_u) {
					q = u;
				}
			}
			else /* r_p >= r_v */
			{
				if (r_q > r_u) {
					q = u;
				}
			}
		}
	}
	return superposing;
}

LineSegmentRelationship line_segments_relationship(const vec2& m, const vec2& n, const vec2& u, const vec2& v) {
	// test result
	// -1: parallel to each other but not superposing
	// 0: unrelated (disjunct and not parallel)
	// 1: intersect each other
	// 2: superpose each other
	// 3: adjacent to each other
	// 4: one segment touches other segment in one point

	LineSegmentRelationship result;

	assert(m != n && u != v);


	bool adjacent = (m == u || m == v || n == u || n == v);

	// pre calculate direction vectors
	vec2 v_u(v-u);
	vec2 n_m(n-m);
	vec2 u_m(u-m);

	// cross product of the direction vectors
	// of both lines. If 0 --> parallel.
	float denominator = cross_z(v_u,n_m);
	if (unlikely(denominator == 0)) {
		// lines are parallel to each other
		float unused;
		bool u_on_mn = line_contains_point(m, n, u, unused);
		bool v_on_mn = line_contains_point(m, n, v, unused);
		if (u_on_mn && v_on_mn) {
			// segments superpose each other
			if (adjacent) {
				vec2 unused1;
				vec2 unused2;
				if (!line_segment_superpose(m,n,u,v, unused1, unused2)) {
					return LineSegmentRelationship::LINES_ADJACENT;
				}
			}
			return LineSegmentRelationship::LINES_SUPERPOSE;
		} else if (u_on_mn || v_on_mn) {
			return LineSegmentRelationship::LINES_TOUCHING;
		} else {
			// lines are truly parallel to each other.
			return LineSegmentRelationship::LINES_PARALLEL;
		}
	} else if (adjacent) {
		return LineSegmentRelationship::LINES_ADJACENT;
	} else {


		float t = cross_z(n_m, u_m) / denominator;
		float s = cross_z(v_u, u_m) / denominator;

		// truly intersect if result == 1
		result = LineSegmentRelationship(t >= 0 && t <= 1 && s >= 0 && s <= 1);

		if (result) {
			float unused;
			if (line_contains_point(m, n, u, unused)
				|| line_contains_point(m, n, v, unused)
				|| line_contains_point(u, v, m, unused)
				|| line_contains_point(u, v, n, unused))
			{
				return LineSegmentRelationship::LINES_TOUCHING;
			}
		}

		return result;
	}
}



bool line_segment_intersection(const vec2& m, const vec2& n, const vec2& u, const vec2& v, vec2& p, float& s) {
	// TODO: arithmetic precision guards
	// We do not test, whether the result is
	// incorrect due to the limited precision of float.
	// TODO: optimise
	// this can be further optimised,
	// - either by not using vec2 internally
	// - or not using const vec2& params

	/*
		Given the equations
		s1 = m + (n-m) * s , 0 <= s <= 1
		s2 = u + (v-u) * t , 0 <= t <= 1
		We search for the intersection p,
		of both line equations, such that
			s1(s) = s2(t)   (=p)
		We solve the equation system to get
		s AND t, and test if both satisfy
		their boundary condition (meaning:
		p is on both segments).
	*/

	// test result
	bool intersects;

	// pre calculate direction vectors
	vec2 v_u(v-u);
	vec2 n_m(n-m);
	vec2 u_m(u-m);

	// cross product of the direction vectors
	// of both lines. If 0 --> parallel.
	float denominator = cross_z(v_u,n_m);
	if (unlikely(denominator == 0)) {
		// lines are parallel to each other
		// or one of the segments is a point (length = 0)

		// Both cases are defined as non-intersecting.
		return false;
	} else {

		float t;
		t = cross_z(n_m, u_m) / denominator;
		s = cross_z(v_u, u_m) / denominator;

		intersects = (t >= 0 && t <= 1 && s >= 0 && s <= 1);
		if (unlikely(intersects)) {
			p = m + s * (n_m);

		}
		return intersects;
	}
}

bool line_segment_intersection_highp(const vec2& m, const vec2& n, const vec2& u, const vec2& v, vec2& p, double& s) {
	// TODO: arithmetic precision guards
	// We do not test, whether the result is
	// incorrect due to the limited precision of float.
	// TODO: optimise
	// this can be further optimised,
	// - either by not using vec2 internally
	// - or not using const vec2& params

	/*
		Given the equations
		s1 = m + (n-m) * s , 0 <= s <= 1
		s2 = u + (v-u) * t , 0 <= t <= 1
		We search for the intersection p,
		of both line equations, such that
			s1(s) = s2(t)   (=p)
		We solve the equation system to get
		s AND t, and test if both satisfy
		their boundary condition (meaning:
		p is on both segments).
	*/

	// test result
	bool intersects;

	// pre calculate direction vectors
#define vec_minus(a,b) a.x = a.x - b.x; a.y = a.y - b.y
	highp_dvec2 v_u(v);
	vec_minus(v_u,u);
	highp_dvec2 n_m(n);
	vec_minus(n_m,m);
	highp_dvec2 u_m(u);
	vec_minus(u_m,m);
#undef vec_minus
	// cross product of the direction vectors
	// of both lines. If 0 --> parallel.
	double denominator = cross_z(v_u,n_m);
	if (unlikely(denominator == 0)) {
		// lines are parallel to each other
		// or one of the segments is a point (length = 0)

		// Both cases are defined as non-intersecting.
		return false;
	} else {

		double t;
		t = cross_z(n_m, u_m) / denominator;
		s = cross_z(v_u, u_m) / denominator;

		intersects = (t >= 0 && t <= 1 && s >= 0 && s <= 1);
		if (unlikely(intersects)) {
			// g2(t):
			p.x = double(m.x) + s * n_m.x;
			p.y = double(m.y) + s * n_m.y;
		}

		return intersects;
	}
}
}; // namespace math
