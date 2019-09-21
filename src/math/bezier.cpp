/*
 * bezier.cpp
 *
 *  Created on: 26 Jul 2019
 *      Author: homac
 */


#include <queue>
#include <stack>
#include <math.h>

#include <stdio.h>


#include <glm/gtx/vector_angle.hpp>

#include "math/bezier.h"
#include "math/line.h"
#include "math/polynom.h"

using namespace glm;
using namespace std;



namespace math {



/**
 * checks for any kind of intersection (also superposition) and returns an
 * approximate intersection point.
 *
 * Accuracy is 1/4 of length(m, n), in case of a superposition.
 * Accuracy is way higher in case of actual intersections.
 *
 * Function returns false, if there is no such intersection or one of the
 * line segments is a point.
 *
 *
 */
static bool segment_intersection(vec2 const & m, vec2 const & n, vec2 const & u, vec2 const & v, vec2& p) {
	if (m == n || u == v) return false;

	float s;
	if (line_segment_intersection(m, n, u, v, p, s)) {
		return true;
	} else if (0 == cross_z(v-u,n-m)) {
		vec2 q;
		if (line_segment_superpose(m, n, u, v, p, q)) {
			p = mix(p,q,0.5);
			return true;
		}
	}
	return false;
}


/** internally used only */
struct __bezier_t {
	vec2 p0;
	vec2 p1;
	vec2 p2;
	vec2 p3;
	float approximateExtend;
	__bezier_t() : p0(),p1(),p2(),p3(), approximateExtend() {}
	__bezier_t(vec2 v0, vec2 v1, vec2 v2, vec2 v3) {
		p0 = v0;
		p1 = v1;
		p2 = v2;
		p3 = v3;
		approximateExtend = approx_ext();
	}

	/** Splits this bezier at f in two new bezier b1 (interval [0,f[) and b2 (interval [f,1]) */
	static void split(float f, __bezier_t& b1, __bezier_t& b2) {
		float previousExtend = b1.approximateExtend;
		bezier_split(f, b1.p0,b1.p1,b1.p2,b1.p3, b2.p0, b2.p1, b2.p2, b2.p3);
		b2.approximateExtend = previousExtend * (1-f);
		b1.approximateExtend = previousExtend * f;
	}

	/**
	 * Conservative calculates approximation for the maximum deviation
	 * of any given point inside the bounding box {p0,p1,p2,p3}
	 * and an actual intersection point on the curve.
	 */
	float approx_ext() {
		// grav. center c of bounding box
		vec2 c = (p0 + p1 + p2 + p3)/4.f;
		// maximum difference between grav. center and p0 times two
		return 2.0*std::max(std::max(length(p0-c), length(p1-c)), std::max(length(p2-c), length(p3-c)));

	}


	float maximumExtend() {
		return approximateExtend;
	}

	static bool intersects(const vec2& v0, const vec2& v1, const vec2& m, const vec2& n, vec2& v) {
		return (v0 != v1) &&  segment_intersection(v0, v1, m, n, v);
	}

	/**
	 * Checks whether the bounding box {p0,p1,p2,p3} intersects given
	 * line segment {m,n}.
	 * @param v vector to intersection point or approximate location (if one is found)
	 * @return
	 *   + 1 if p1->p2 intersects.
	 *   - 1 if primitive bounding box intersects
	 *     0 if no intersection is possible
	 */
	int bb_intersects(const vec2& m, const vec2& n, vec2& v) {
		if (intersects(p1, p2, m, n, v))
		{
			// We accept only intersections with
			// {p1,p2} as true intersections
			// This way, we can compare v to previous samples of v
			return 1;
		}
		else
		{
			// these intersections are *not* considered
			// to be true intersections.
			vec2 unused;
			return -(  intersects(p0, p1, m, n, unused)
					|| intersects(p2, p3, m, n, unused)
					|| intersects(p0, p3, m, n, unused)
					);
		}
	}

	/**
	 * Checks whether the bounding box {p0,p1,p2,p3} intersects the bounding
	 * box of the given bezier curve {q0,q1,q2,q3}.
	 *
	 * @param v vector to intersection point or approximate location (if one is found)
	 * @return
	 *   + 1 if p1->p2 intersects.
	 *   - 1 if primitive bounding box intersects
	 *     0 if no intersection is possible
	 */
	int bb_intersects(const vec2& q0, const vec2& q1, const vec2& q2, const vec2& q3, vec2& v) {
		vec2 results[3];
		float tolerance = 0;
		int count = bezier_line_segment_intersections(q0, q1, q2, q3, p1, p2, tolerance, results);
		if (1 == count)
		{
			// We accept only intersections between
			// {p1,p2} and {q1,q2} as true intersections
			// This way, we can compare v to previous samples of v
			v = results[0];
			return 1;
		}
		else
		{
			// these intersections are *not* considered
			// to be true intersections.
			return -(  0 != bezier_line_segment_intersections(q0, q1, q2, q3, p0, p1, tolerance, results)
					|| 0 != bezier_line_segment_intersections(q0, q1, q2, q3, p2, p3, tolerance, results)
					|| 0 != bezier_line_segment_intersections(q0, q1, q2, q3, p0, p3, tolerance, results)
					);
		}
	}

};






/**
 * Determines minima and maxima in x and y direction.
 * Returns true extrema only (no saddle points).
 *
 * @param t_min {t(min_x), t(min_y)}
 * @param t_max {t(max_x), t(max_y)}
 * @return true if at least one extrem was found
 */
bool bezier_extrema(vec2 const& p0, vec2 const& p1, vec2 const& p2, vec2 const& p3, dvec2& t_min, dvec2& t_max) {
	bool found = false;
	vec2 a = (p3 - 3.0f*p2 + 3.0f*p1 - p0);
	vec2 b = 3.f*(p2 - 2.0f*p1 + p0);
	vec2 c = 3.f*(p1 - p0);
	int analysis;

	auto bounds = CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE;

	analysis = polynom3_real_extrema(a.x,b.x,c.x, false, t_min.x, t_max.x, bounds);
	found = (analysis != 0);

	analysis = polynom3_real_extrema(a.y,b.y,c.y, false, t_min.y, t_max.y, bounds);
	found |= (analysis != 0);

	return found;

}

bool bezier_inflection_point(const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3, float& t) {
	// FIXME: missing proof that bezier can't have two inflection points in [0:1]
	dvec2 a = (p3 - 3.0f*p2 + 3.0f*p1 - p0);
	dvec2 b = (p2 - 2.0f*p1 + p0);
	dvec2 c = (p1 - p0);

	// solve   0 = 6 (b_a t^2 + c_a t + c_b)   for t
	// with    0 < t < 1    and    0 <> 6 (2 b_a t + c_a)
	double b_a = 6.f*cross_z(b,a);
	double c_a = 6.f*cross_z(c,a);
	double c_b = 6.f*cross_z(c,b);

	double t_i[2];
	int n = polynom2_real_roots(b_a, c_a, c_b, t_i, CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE);

	int count = 0;
	for (int i = 0; i < n; i++) {
		double derivate = polynom1_value(2.f*b_a, c_a, t_i[i]);
		if (derivate != 0) {
			t = t_i[i];
			count++;
		}
	}

	assert(count <= 1);
	return (count);
}


/* This function is used for comparison/debugging only.
 * Use faster (non-iterative) function bezier_line_segment_intersections instead. */
int bezier_line_segment_intersections_iterative(const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3, const vec2& m, const vec2& n, float tolerance, vec2 result[3]) {
	assert(tolerance >= 0.000001f);
	int count = 0;
	vec2 v;

#ifndef NDEBUG
	int iterations = 0;
#endif

	// work stack
	stack<__bezier_t> work;

	// temporary variables
	__bezier_t second_half;
	float f;

	vec2 INVALID(INFINITY,INFINITY);
	// Last sample of a valid intersection point
	// for the bezier part, currently worked on.
	// Valid means, the intersection was found in segment p1 and p2
	vec2 valid_sample = INVALID;

	if (bezier_inflection_point(p0, p1, p2, p3, f)) {
		// avoid special case: bezier with turning point
		__bezier_t b(p0, p1, p2, p3);
		// split in turning point
		__bezier_t::split(f, b, second_half);
		work.push(b);
		work.push(second_half);
	} else {
		work.push(__bezier_t(p0,p1,p2,p3));
	}
	// first split
	while (work.size()) {
#ifndef NDEBUG
		iterations++;
#endif
		__bezier_t& b = work.top();



		//
		// check first half of split
		//

		int intersecting = b.bb_intersects(m, n, v);
		if (intersecting) {
			if (b.maximumExtend() <= tolerance || length(valid_sample-v) <= (tolerance*2))
			{
				// accuracy ok -> store and proceed with next in queue
				result[count++] = v;
				valid_sample = INVALID;

				work.pop(); // done with this entry
			} else {
				// accuracy still too low -> split
				if (intersecting == 1) {
					valid_sample = v;
				}
				// TODO: performance: use reference to next item on stack instead of local variables?
				__bezier_t::split(0.5, b, second_half);
				work.push(second_half);
			}
		} else {
			// no intersection possible -> discard
			work.pop(); // loose end

			// NOTE: valid sample
			// Possible cases:
			// 1. valid_sample = INVALID (still)
			// 2. valid_sample of a predecessor of current bezier fraction -> reuse as valid sample for second_half of same predecessor
			// valid sample will always be kept until an intersection was written to p[]
		}


	}
#ifndef NDEBUG
	printf("%d\n", iterations);
#endif
	return count;
}

/**
 * This is an internal function used by bezier_line_segment_intersections.
 *
 *
 * Given bezier curve  B(t, p0, p1, p2, p3) and line segment L(s, m, n) = m + (n-m) s
 * calculate all t_i and corresponding s_i which satisfy
 *     B(t_i) - L(s_i) == 0
 *     and   0 <  t_i <  1  // TODO: why is that not less equal?!
 *     and   0 <= s_i <= 1
 *
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T>
int bezier_line_segment_intersections_t(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& m, const vec2& n,
		float tolerance, //<< @deprecated TODO: remove
		double t[3], // interpolation factor for bezier curve
		double s[3], // interpolation factor for line segment
		CO_DOMAIN_T interval = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE)
{
	int count = 0;
	// Known:
	//    B(t, p0, p1, p2, p3)       : bezier curve with control points {p0,p1,p2,p3}
	//    L(s, m, n) = m + (n-m) s   : line segment with start and end points {n,m}
	// Find all t_i for
	//      B(t_i) - L(s_i) == 0

	// Bezier coefficients for polynom form
	dvec2 a(p3 - 3.f*p2 + 3.f*p1 -     p0);
	dvec2 b(     3.f*p2 - 6.f*p1 + 3.f*p0);
	dvec2 c(              3.f*p1 - 3.f*p0);
	dvec2 d(                           p0);


	dvec2 n_m(n-m);
	dvec2 d_m(d-dvec2(m));

	// t_i for f(t_i) == 0

	double denominator;
	double A;
	double B;
	double C;
	double D;

	if (n_m.x != 0) {
		// use x coordinates to determine t_i
		denominator = n_m.x;
		A = (    cross_z(n_m,a))   / denominator;
		B = (    cross_z(n_m,b))   / denominator;
		C = (    cross_z(n_m,c))   / denominator;
		D = (    cross_z(n_m,d_m)) / denominator;

		count = polynom3_real_roots(A, B, C, D, t, interval);


		// prepare calculation of s_i(t_i) based on x coordinates
		A = a.x;
		B = b.x;
		C = c.x;
		D = d_m.x;

	} else if (n_m.y != 0) {
		// use y coordinates to determine t_i
		denominator = n_m.y;
		A = (    cross_z(a,  n_m)) / denominator;
		B = (    cross_z(b,  n_m)) / denominator;
		C = (    cross_z(c,  n_m)) / denominator;
		D = (    cross_z(d_m,n_m)) / denominator;
		count = polynom3_real_roots(A, B, C, D, t, interval);

		// prepare calculation of s_i(t_i) based on y coordinates
		A = a.y;
		B = b.y;
		C = c.y;
		D = d_m.y;

	}


	// Test if inside line segment bounds
	// -> calculate s_i(t_i)
	//    and check 0 <= s_i <= 1
	int total = count;
	for (int i = 0, j = 0; i < total; i++) {
		double s_i = polynom3_value(A, B, C, D, t[i])/denominator;
		if (CO_DOMAIN_REAL_IN_0_1_INCLUSIVE(s_i)) {
			s[j] = s_i;
			t[j] = t[i];
			j++;
		} else {
			count--;
		}
	}

	return count;
}



int bezier_line_segment_intersections(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& m, const vec2& n,
		float tolerance, ///< @deprecated TODO: remove
		vec2 result[3])
{
	int count = 0;

	// t_i for f(t_i) == 0
	double t[3];
	double s[3];

	count = bezier_line_segment_intersections_t(p0, p1, p2, p3, m, n, tolerance, t, s);

	for (int i = 0; i < count; i++) {
		bezier_point_highp(t[i], p0,p1,p2,p3, result[i]);
	}
	return count;
}


/* Computes all intersections of two bezier curves using
 * the recursive subdivision method.
 * Tolerance controls the maximum deviation of accepted
 * results (intersection points) of the true intersection
 * point. Thus, lower tolerance means longer processing time. */
int bezier_bezier_intersections(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& q0, const vec2& q1, const vec2& q2, const vec2& q3,
		float tolerance,
		vec2 result[9])
{
	assert(tolerance >= 0.000001f);
	int count = 0;
	vec2 v;



#ifndef NDEBUG
	int iterations = 0;
#endif

	// work stack
	stack<__bezier_t> work;

	// temporary variables
	__bezier_t second_half;
	float f;

	const vec2 INVALID(INFINITY,INFINITY);
	// Last sample of a valid intersection point
	// for the bezier part, currently worked on.
	// Valid means, the intersection was found in segment p1 and p2
	vec2 valid_sample = INVALID;

	if (bezier_inflection_point(p0, p1, p2, p3, f)) {
		// avoid special case: bezier with turning point
		__bezier_t b(p0, p1, p2, p3);
		// split in turning point
		__bezier_t::split(f, b, second_half);
		work.push(b);
		work.push(second_half);
	} else {
		work.push(__bezier_t(p0,p1,p2,p3));
	}
	// first split
	while (work.size()) {
#ifndef NDEBUG
		iterations++;
#endif
		__bezier_t& b = work.top();



		//
		// check first half of split
		//

		int intersecting = b.bb_intersects(q0, q1, q2, q3, v);
		if (intersecting) {
			if (b.maximumExtend() <= tolerance || length(valid_sample-v) <= (tolerance*2))
			{
				// accuracy ok -> store and proceed with next in queue
				result[count++] = v;
				valid_sample = INVALID;

				work.pop(); // done with this entry
			} else {
				// accuracy still too low -> split
				if (intersecting == 1) {
					valid_sample = v;
				}
				// TODO: performance: use reference to next item on stack instead of local variables?
				__bezier_t::split(0.5, b, second_half);
				work.push(second_half);
			}
		} else {
			// no intersection possible -> discard
			work.pop(); // loose end

			// NOTE: valid sample
			// Possible cases:
			// 1. valid_sample = INVALID (still)
			// 2. valid_sample of a predecessor of current bezier fraction -> reuse as valid sample for second_half of same predecessor
			// valid sample will always be kept until an intersection was written to p[]
		}


	}
#ifndef NDEBUG
	printf("%d\n", iterations);
#endif
	return count;
}












} // namespace math


