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

#include "math2d/bezier.h"
#include "math2d/line.h"
#include "math2d/polynom.h"

using namespace glm;
using namespace std;



namespace math2d {


int BEZIER_BEZIER_iterations = 0;

static const vec2 VEC2_INVALID(INFINITY,INFINITY);



/** internally used only */
struct __bezier_t__ {
	vec2 p0;
	vec2 p1;
	vec2 p2;
	vec2 p3;

	float approximateExtend;

	/*
	 * Split tracking:
	 *
	 * The original bezier O(t_o)={o0, o1, o2, o3} has scaling factor s = 1.0 .
	 * Splitting O(t_o) at  t_o = s1  gives two beziers P(t_p) and Q(t_q)
	 *         P(t_p)={p0,p1,p2,p3}    with s_p=s1, b_p = 0  and  p0 = O(f_p0) and p3 = O(f1)
	 *         Q(t_q)={q0,q1,q2,q3}    with s_q=1-s1, b_q = s1  and  q0 = O(f_q0) and q3 = O(1)
	 * Any P(t_p) can be mapped to O(t_o) by the following expression
	 *         t_o =  t_p/s_p + b_p    gives   O(t_o) == P(t_p)
	 * Any Q(t_q) can be mapped to O(t_o) by the following expression
	 *         t_o =  t_q/s_q + b_q    gives   O(t_o) == Q(t_q)
	 *
	 * Splitting Q again, now at t_q = f2 (= 0.5), gives P2 and Q2 with:
	 *         P2.b = Q.b + 0
	 *         Q2.b = P2.b + f2*s_p
	 *         f    = f2*f
	 * Using the same expression above to calculate t_o, will again satisfy
	 *         O(t_o) == P(t_p)
	 *       and
	 *         O(t_o) == Q(t_q)
	 */
	double s;
	double b;




	__bezier_t__() : p0(),p1(),p2(),p3(), approximateExtend(), s(1), b(0) {}

	__bezier_t__(vec2 v0, vec2 v1, vec2 v2, vec2 v3) {
		p0 = v0;
		p1 = v1;
		p2 = v2;
		p3 = v3;
		approximateExtend = approx_ext();
		s = 1.0;
		b = 0;
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

	/** Splits this bezier at f in two new bezier b1 (interval [0,f[) and b2 (interval [f,1]) */
	static void split(float t, __bezier_t__& b1, __bezier_t__& b2) {
		double f1 = b1.s * t;
		double f2 = b1.s * (1-t);
		double orig_b = b1.b;

		float previousExtend = b1.approximateExtend;

		bezier_split(t, b1.p0,b1.p1,b1.p2,b1.p3, b2.p0, b2.p1, b2.p2, b2.p3);

		b1.approximateExtend = previousExtend * t;
		b2.approximateExtend = previousExtend * (1-t);

		b1.s = f1;
		b1.b = orig_b + 0;
		b2.s = f2;
		b2.b = orig_b + f2;
	}

	double getOriginalFactor(double t_p) {
		return t_p*s + b;
	}

	float maximumExtend() {
		return approximateExtend;
	}

	/**
	 * checks for any kind of intersection (also superposition) and returns an
	 * approximate intersection point.
	 *
	 * Accuracy is 1/4 of length(m, n), in case of a superposition.
	 * Accuracy is way higher in case of actual intersections.
	 *
	 * Function returns false, if there is no such intersection or one of the
	 * line segments is a point.
	 * Note: not optimised
	 */
	static bool intersects(const vec2& m, const vec2& n, const vec2& u, const vec2& v, vec2& p) {
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

	/**
	 * Checks whether the bounding box {p0,p1,p2,p3} intersects given
	 * line segment {m,n}.
	 * @param v (out) vector to intersection point or approximate location (if one is found)
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
	 * Checks whether the bounding box {p0,p1,p2,p3} of this bezier curve intersects the bounding
	 * box of the given bezier curve Q(t)={q0,q1,q2,q3}.
	 *
	 * @param v (out) vector to intersection point or approximate location (if one is found)
	 * @param t (out) interpolation factor with Q(t)=v
	 * @return
	 *   + 1 if p1->p2 intersects.
	 *   - 1 if primitive bounding box intersects
	 *     0 if no intersection is possible
	 */
	int bb_intersects_t(const vec2& q0, const vec2& q1, const vec2& q2, const vec2& q3, vec2& v, double& t_q) {
		double result_s[3];
		double result_t[3];
		const float tolerance = 0;
		double t_q_samples[4] = {0};
		int t_count = 0;
		int sample_classification = 0;
		int total_count = 0;

		// base line intersection?
		int c03 = bezier_line_segment_intersections_t(q0, q1, q2, q3, p0, p3, tolerance, result_t, result_s);
		total_count += c03;
		if (c03 > 0) t_q_samples[t_count++] = result_t[0];

		int c01 = bezier_line_segment_intersections_t(q0, q1, q2, q3, p0, p1, tolerance, result_t, result_s);
		total_count += c01;
		if (c01 > 0) t_q_samples[t_count++] = result_t[0];

		int c12 = bezier_line_segment_intersections_t(q0, q1, q2, q3, p1, p2, tolerance, result_t, result_s);
		total_count += c12;
		if (c12 > 0) t_q_samples[t_count++] = result_t[0];

		int c23 = bezier_line_segment_intersections_t(q0, q1, q2, q3, p2, p3, tolerance, result_t, result_s);
		total_count += c23;
		if (c23 > 0) t_q_samples[t_count++] = result_t[0];

		if (c03 == 1 && total_count == 2) {
			// we accept two found intersections
			// if there are exactly two intersections in total
			// and one intersection is on the base line
			t_q = (t_q_samples[0] + t_q_samples[1])/2;
			bezier_point(t_q, q0,q1,q2,q3, v);
			// valid sample
			sample_classification = 1;
		} else {
			if (t_count) {
				// we have some samples for t_q, use them to get an average
				// value for t_q
				int i;
				t_q = 0;
				for (i = 0; i < t_count; i++) t_q += t_q_samples[i];
				t_q /= i;
				bezier_point(t_q, q0,q1,q2,q3, v);
			}
			// vague sample
			sample_classification = -(total_count > 0);
		}
		return sample_classification;
	}



};


/**
 * Determines minima and maxima in x and y direction.
 * Returns true extrema only (no saddle points).
 *
 * @param t_min {t(min_x), t(min_y)}
 * @param t_max {t(max_x), t(max_y)}
 * @return true if at least one extreme was found
 */
bool bezier_extrema(vec2 const& p0, vec2 const& p1, vec2 const& p2, vec2 const& p3, dvec2& t_min, dvec2& t_max) {
	bool found = false;
	vec2 a = (p3 - 3.0f*p2 + 3.0f*p1 - p0);
	vec2 b = 3.f*(p2 - 2.0f*p1 + p0);
	vec2 c = 3.f*(p1 - p0);
	int analysis;

	auto bounds = CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE;

	analysis = polynom3_extrema(a.x,b.x,c.x, false, t_min.x, t_max.x, bounds);
	found = (analysis != 0);

	analysis = polynom3_extrema(a.y,b.y,c.y, false, t_min.y, t_max.y, bounds);
	found |= (analysis != 0);

	return found;

}




bool bezier_inflection_point(const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3, float& t) {
	// FIXME: specification: missing proof that bezier can't have two inflection points in ]0:1[
	dvec2 a = (p3 - 3.0f*p2 + 3.0f*p1 - p0);
	dvec2 b = (p2 - 2.0f*p1 + p0);
	dvec2 c = (p1 - p0);

	// solve   0 = 6 (bxa t^2 + cxa t + cxb)   for t   --> same as  0 = bxa t^2 + cxa t + cxb
	// with    0 < t < 1    and    0 <> 6 (2 bxa t + cxa)
	double bxa = cross_z(b,a);
	double cxa = cross_z(c,a);
	double cxb = cross_z(c,b);

	double t_i[2];
	int n = polynom2_roots(bxa, cxa, cxb, t_i, CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE);

	int count = 0;
	for (int i = 0; i < n; i++) {
		double derivate = polynom1_value(2.f*bxa, cxa, t_i[i]);
		if (derivate != 0) {
			t = t_i[i];
			count++;
		}
	}

	assert(count <= 1);
	return (count);
}





template
<
	class CO_DOMAIN_T = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T,
	class CO_DOMAIN_S = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T
>
int bezier_line_segment_intersections_t(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& m, const vec2& n,
		float tolerance, //<< @deprecated TODO: remove
		double t[3], // interpolation factor for bezier curve
		double s[3], // interpolation factor for line segment
		CO_DOMAIN_T interval_t,
		CO_DOMAIN_S interval_s)
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

		count = polynom3_roots(A, B, C, D, t, interval_t);


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
		count = polynom3_roots(A, B, C, D, t, interval_t);

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
		if (interval_s(s_i)) {
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





/**
 * This function is used for comparison/debugging purposes only.
 * Use faster (non-iterative) function bezier_line_segment_intersections instead.
 */
int bezier_line_segment_intersections_iterative(const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3, const vec2& m, const vec2& n, float tolerance, vec2 result[3]) {
	assert(tolerance >= 0.000001f);
	int count = 0;
	vec2 v;


	// work stack
	stack<__bezier_t__> work;

	// temporary variables
	__bezier_t__ second_half;
	float f;

	// Last sample of a valid intersection point
	// for the bezier part, currently worked on.
	// Valid means, the intersection was found in segment p1 and p2
	vec2 valid_sample = VEC2_INVALID;

	if (bezier_inflection_point(p0, p1, p2, p3, f)) {
		// avoid special case: bezier with turning point
		__bezier_t__ b(p0, p1, p2, p3);
		// split in turning point
		__bezier_t__::split(f, b, second_half);
		work.push(b);
		work.push(second_half);
	} else {
		work.push(__bezier_t__(p0,p1,p2,p3));
	}
	// first split
	while (work.size()) {
		__bezier_t__& b = work.top();



		//
		// check first half of split
		//

		int intersecting = b.bb_intersects(m, n, v);
		if (intersecting) {
			if (b.maximumExtend() <= tolerance || length(valid_sample-v) <= (tolerance*2))
			{
				// accuracy ok -> store and proceed with next in queue
				result[count++] = v;
				valid_sample = VEC2_INVALID;

				work.pop(); // done with this entry
			} else {
				// accuracy still too low -> split
				if (intersecting == 1) {
					valid_sample = v;
				}
				__bezier_t__::split(0.5, b, second_half);
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
	return count;
}





int __internal__bezier_bezier_intersections_t(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& q0, const vec2& q1, const vec2& q2, const vec2& q3,
		float tolerance,
		double result_t_p[9],
		double result_t_q[9]
		);

int bezier_bezier_intersections_t(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& q0, const vec2& q1, const vec2& q2, const vec2& q3,
		float tolerance,
		double result_t_p[9],
		double result_t_q[9]
		)
{
	typedef __bezier_t__ bezier_t;

	// our internal bezier bezier intersections function expects the
	// second curve to be not self-intersecting. Thus, we have to deal
	// with this first.

	bool p_self_intersecting = bezier_self_intersecting(p0,p1,p2,p3);
	bool q_self_intersecting = bezier_self_intersecting(q0,q1,q2,q3);

	if (unlikely(p_self_intersecting && q_self_intersecting)) {
		// both self intersecting
		// split Q in half and test all combinations PxQ1 and PxQ2
		bezier_t Q1(q0,q1,q2,q3), Q2;
		bezier_t::split(0.5, Q1, Q2);

		// Since we split one curve, we may get duplicate entries at
		// the splitting location.
		// Memorise values of t_p at the boundary between Q1 and Q2
		// to check them in the second run.
		double t_p_boundary_1 = INFINITY;
		double t_p_boundary_2 = INFINITY;
		int total_count = 0;
		int count = __internal__bezier_bezier_intersections_t(p0, p1, p2, p3, Q1.p0, Q1.p1, Q1.p2, Q1.p3,
				tolerance,
				result_t_p, result_t_q);
		double *t_q_end = result_t_q + total_count + count;
		for (double *t_q = result_t_q + total_count; t_q < t_q_end; t_q++) {
			result_t_q[total_count] = Q1.getOriginalFactor(*t_q);
			if (about_equal(*t_q, 1.0, double(tolerance))) {
				// found a t_q at the splitting location
				if (t_p_boundary_1 == INFINITY) t_p_boundary_1 = result_t_p[total_count];
				else t_p_boundary_1 = result_t_p[total_count];
			}
			total_count++;
		}
		count = __internal__bezier_bezier_intersections_t(p0, p1, p2, p3, Q2.p0, Q2.p1, Q2.p2, Q2.p3,
				tolerance,
				result_t_p + total_count, result_t_q + total_count);
		t_q_end = result_t_q + total_count + count;
		for (double *t_q = result_t_q + total_count, *t_p = result_t_p + total_count;
				t_q < t_q_end; t_q++, t_p++) {
			result_t_q[total_count] = Q2.getOriginalFactor(*t_q);
			result_t_p[total_count] = *t_p;
			if (about_equal(*t_q, 0.0, double(tolerance)))
			{
				// found a t_q at the splitting location
				// -> check if we have an identical match from previous run
				if (about_equal(t_p_boundary_1, result_t_p[total_count], double(tolerance))) {
					t_p_boundary_1 = INFINITY;
					total_count--;
				} else if (about_equal(t_p_boundary_2, result_t_p[total_count], double(tolerance))) {
					t_p_boundary_2 = INFINITY;
					total_count--;
				}

			}
			total_count++;
		}
		return total_count;
	} else if (q_self_intersecting) {
		// P is not self intersecting
		// --> swap P and Q
		return __internal__bezier_bezier_intersections_t(q0,q1,q2,q3,p0,p1,p2,p3, tolerance, result_t_q, result_t_p);
	} else {
		// everything is fine -> go ahead
		return __internal__bezier_bezier_intersections_t(p0,p1,p2,p3,q0,q1,q2,q3, tolerance, result_t_p, result_t_q);
	}
	return 0;
}


int __internal__bezier_bezier_intersections_t(
		const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3,
		const vec2& q0, const vec2& q1, const vec2& q2, const vec2& q3,
		float tolerance,
		double result_t_p[9],
		double result_t_q[9]
		)
{

	typedef __bezier_t__ bezier_t;


	/*
	 * Procedure
	 * ---------
	 * This function implements the recursive subdivision method to
	 * determine intersection points of two bezier curves in an iterative
	 * fashion using a local working stack (internal stack management).
	 * One bezier curve cannot have intersection points with the other
	 * if there is no intersection of its bounding box {p0,p1,p2,p3}
	 * and the other bezier curve. Thus, an intersection of any of the
	 * line segments {p0,p1}, {p1,p2}, {p2,p3} or {p0,p3} with the other
	 * bezier curve {q0,q1,q2,q3} indicates a possible intersection.
	 * If a bounding box test fails, the curve (or sub-curve) gets removed from
	 * the work stack and will not be processed any further.
	 * A successful test may return a valid intersection point
	 * sample, that has to be further analysed to determine its accuracy
	 * (see below). If the accuracy of the sample does not satisfy the
	 * given tolerance requirement (see parameter) then the bezier curve
	 * is subdivided (split) at t=0.5 into two sub-curves, with its
	 * respective control points. Those two bezier sub-curves replace
	 * the previous curve on the working stack and will be tested in
	 * subsequent iterations of the main loop. Once, the accuracy
	 * of a sample has reached the given tolerance, the sample is
	 * stored as a result and the (sub-)curve is removed from the working stack.
	 *
	 * Intersection Point Samples
	 * --------------------------
	 * Intersection point samples may be classified as
	 * - valid:   Sample is an approximation for one specific intersection point,
	 *            which is known to exist.
	 * - vague:   Sample is an approximation for multiple possible intersection
	 *            points, which may also be non-existing.
	 * - invalid: Sample contains no information.
	 *
	 * Valid and vague samples will be further analysed to determine its accuracy.
	 *
	 * Accuracy is determined based on two methods:
	 * 1. Deviation (distance) between current and last valid sample.
	 * 2. Maximum extend of the bezier-curve estimated based on its bounding box.
	 *
	 * If one of these values is lower or equal to the given accuracy tolerance
	 * the sample is accepted as valid result.
	 *
	 * Initial sub-division
	 * --------------------
	 * The procedure considers only a simple bounding box based on segment
	 * {p0,p1}, {p1,p2}, {p2,p3} and {p0, p3}. This bounding box does not
	 * cover bezier curves with inflection point, where {p1,p2} and {p0,p3}
	 * cross each other. To prevent the algorithm from getting more complex
	 * we just sub-divide curves with inflection point at their inflection
	 * point, thereby turning them into two curves without inflection point.
	 *
	 * Main Loop and Working Order
	 * ---------------------------
	 * The main loop maintains (sub-)curves to be processed on a stack in right
	 * to left order. The right most sub-curve is always on top of the stack
	 * and will be processed first.
	 *
	 * The main loop keeps the last valid sample of a test to be able to
	 * estimate accuracy of the next sample. The sample is meaningful for
	 * the section of the curve, it was found on. Thus, every sub-curve of
	 * this section can be tested against this last sample. Every time a valid
	 * result was found, the last sample is set to INVALID, indicating that
	 * there is no last sample for the current sub-curve available. Thus, if
	 * there are two intersection points to be found on a given section, and
	 * one was already found, there will be one extra iteration of the loop
	 * to create a new sample for the other intersection point, before its
	 * accuracy can be determined. However, since there are multiple
	 * intersections on the section in question, the previous sample must
	 * have been quite fare away from the second intersection point anyway.
	 * Thus, that the additional test is wasted, is very unlikely.
	 *
	 */
	assert(tolerance >= 0.0000009f); // TODO: ?
	const double DOUBLE_INVALID = INFINITY;
	int count = 0;

	// work stack
	stack<bezier_t> work;

	// temporary variables
	bezier_t second_half;
	// interpolation factor for Q(t_q)
	double t_q;
	vec2 v_q;

	// Last sample of a valid intersection point
	// for the bezier part, currently worked on.
	// Valid means, the accuracy of the sample is at least higher than its
	// deviation to the next sample.
	vec2 valid_sample = VEC2_INVALID;

	float f;
	if (bezier_inflection_point(p0, p1, p2, p3, f)) {
		// avoid special case: bezier with turning point
		bezier_t b(p0, p1, p2, p3);
		// split in turning point
		bezier_t::split(f, b, second_half);
		work.push(b);
		work.push(second_half);
	} else {
		work.push(bezier_t(p0,p1,p2,p3));
	}

	double t_tolerance = tolerance/(1<<4);

	// first split
	while (work.size()) {
		BEZIER_BEZIER_COUNT();

		bezier_t& b = work.top();

		//
		// check for intersection
		//
		v_q = VEC2_INVALID;
		t_q = DOUBLE_INVALID;
		int intersecting = b.bb_intersects_t(q0, q1, q2, q3, v_q, t_q);
		if (intersecting) {
			// interpolation factor for P(t_p)
			double t_p = DOUBLE_INVALID;

			if (unlikely(b.maximumExtend() <= tolerance
					|| length(valid_sample-v_q) <= (tolerance*2)))
			{

				// accuracy ok
				// find intersection on original curve P
				// if successful -> store and proceed with next in queue

				// TODO: move all this to a separate function

				// interval of the split bezier projected onto the original bezier curve
				co_domain_t<> interval_on_P(b.getOriginalFactor(0), b.getOriginalFactor(1));


				// To find the factor t_p for the first bezier P(t)
				// we determine the intersection point between the tangent at Q(t_q)
				// and the bezier curve.
				vec2 o,m; // T(s) = o + m * s
				bezier_tangent(t_q, q0, q1, q2, q3, o, m);
				double s_i[3]; // <- receives the s_i for all intersections
				double t_i[3]; // <- receives the t_i for all intersections
				int num = bezier_line_segment_intersections_t(p0,p1,p2,p3,o,o+m, 0.0, t_i, s_i, interval_on_P, CO_DOMAIN_REAL);

				if (likely(num > 0)) {
					// we may get multiple intersection points
					// The one, which is closest to Q(t_q) has
					// an s_i[i] closest to 0.
					int i = 0;
					int best = i;
					s_i[best] *= s_i[best];
					for (i = 0; i < num; i++) {
						s_i[i] *= s_i[i];
						if (s_i[i] < s_i[best]) best = i;
					}
					// found closest

					// this will most certainly fail in rare but possible circumstances where the tangent varies a lot
					// at low deviations of t.
					vec2 v_p;
					bezier_point(t_i[best], p0, p1, p2, p3, v_p);
					if (length(v_q-v_p) <= tolerance) {
						t_p = t_i[best];
					}
				}

				if (unlikely(t_p == DOUBLE_INVALID)) {
					// didn't found a solution yet.
					// check if closest point on P satisfies (distance <= tolerance)
					vec2 v_p;
					double t = bezier_point_closest_point_t(p0,p1,p2,p3,v_q, tolerance, v_p, interval_on_P);
					if (length(v_q-v_p) <= tolerance) {
						t_p = t;
					}
				}
			}

			if (unlikely(t_p != DOUBLE_INVALID)) {
				// we have found a solution
				if (count
						&& about_equal(result_t_q[count-1], t_q, t_tolerance)
						&& about_equal(result_t_p[count-1], t_p, t_tolerance)
						&& bezier_equal_point(t_p, p0,p1,p2,p3, v_q, tolerance)
						)
				{
					// ignore duplicates

					// If an intersection point is exactly on a boundary between
					// two adjacent sub-curves, then it may be detected on both
					// sub-curves. Since the main loop is basically searching from
					// right to left (t=1 to t=0), the previously found identical
					// result is now the last element in the result sets.

					// TODO: performance: try to remove test on duplicates.
					// This branch here can be removed, if we can guarantee, that
					// an intersection point on a domain boundary is found exactly
					// once.

				}
				else
				{
					result_t_q[count] = t_q;
					result_t_p[count] = t_p;
					count++;
				}
				valid_sample = VEC2_INVALID;
				work.pop(); // done with this entry
			} else {
				// accuracy still too low -> split again
				if (intersecting == 1) {
					valid_sample = v_q;
				}
				// TODO: performance: use reference to next item on stack instead of local variables
				bezier_t::split(0.5, b, second_half);
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
	return count;
}










} // namespace math


