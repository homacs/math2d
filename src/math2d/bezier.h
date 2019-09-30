/*
 * Bezier.h
 *
 *  Created on: 9 Jul 2019
 *      Author: homac
 */

#ifndef MATH2D_BEZIER_H_
#define MATH2D_BEZIER_H_

#include "math2d-config.h"
#include <assert.h>
#include <glm/glm.hpp>
#include <glm/common.hpp>

#include "math2d/line.h"
#include "math2d/polynom.h"
#include "float-utils.h"

namespace math2d {


static const double_mantissa_mask_t BEZIER_DOUBLE_PRECISION = double_mantissa_mask(16);
static const float_mantissa_mask_t BEZIER_FLOAT_PRECISION = float_mantissa_mask(8);

static inline void bezier_point(float f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2& p);

static inline void bezier_point(float f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2& p);

static inline void bezier_point(float f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3, glm::vec2& p);

static inline void bezier_point_highp(double f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3, glm::vec2& p);


/**
 * Determines the parameters of a line L(s)=m*s + p, which is a tangent of the bezier curve P(f) at P(f_t).
 *
 * p = P(f_t)
 * m = P'(f_t)
 */
static inline void bezier_tangent(float f_t, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3, glm::vec2& p, glm::vec2& m);

bool bezier_inflection_point(const glm::vec2 p0, const glm::vec2 p1, const glm::vec2 p2, const glm::vec2 p3, float& t);


bool bezier_extrema(glm::vec2 const& p0, glm::vec2 const& p1, glm::vec2 const& p2, glm::vec2 const& p3, glm::dvec2& t_min, glm::dvec2& t_max);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T>
static inline double bezier_point_closest_point_t (
		const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3,
		const glm::vec2& v, float tolerance, glm::vec2& result_v,
		const CO_DOMAIN_T domain = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE);
/**
 * Intersections between bezier curve and a line segment from m to n.
 * L(s) = m + (n-m)* s
 */
int bezier_line_segment_intersections(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3, const glm::vec2& m, const glm::vec2& n, float tolerance, glm::vec2 results[3]);

/**
 * This is a function used by bezier_line_segment_intersections.
 * This function returns the interpolation factors for B(t) and L(s).
 *
 * Given bezier curve  B(t, p0, p1, p2, p3) and line segment L(s, m, n) = m + (n-m) s
 * calculate all t_i and corresponding s_i which satisfy
 *     B(t_i) - L(s_i) == 0
 *     and   0 <  t_i <  1  // TODO: why is that not less equal?!
 *     and   0 <= s_i <= 1
 *
 */
template
<
	class CO_DOMAIN_T = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T,
	class CO_DOMAIN_S = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T
>
int bezier_line_segment_intersections_t(
		const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3,
		const glm::vec2& m, const glm::vec2& n,
		float tolerance, //<< @deprecated TODO: remove
		double t[3], // interpolation factor for bezier curve
		double s[3], // interpolation factor for line segment
		CO_DOMAIN_T interval_t = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE,
		CO_DOMAIN_S interval_s = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE);

/**
 * Computes all intersections of two bezier curves using
 * the recursive subdivision method.
 * Parameter 'tolerance' controls the maximum deviation of accepted
 * results (intersection points) of the true intersection
 * point. Thus, lower tolerance effectively means higher
 * processing effort.
 */
int bezier_bezier_intersections_t(
		const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3,
		const glm::vec2& q0, const glm::vec2& q1, const glm::vec2& q2, const glm::vec2& q3,
		float tolerance, double t_p[9], double t_q[9]);



////////////////////////////////////////////////////////////////////////////////////
//          I M P L E M E N T A T I O N S     B E L O W
////////////////////////////////////////////////////////////////////////////////////

#ifdef MATH2D_EVALUATE
	extern int BEZIER_BEZIER_iterations;
#	define BEZIER_BEZIER_COUNT()             (BEZIER_BEZIER_iterations++)
#	define BEZIER_BEZIER_EVALUATION_RESET()  (BEZIER_BEZIER_iterations = 0)
#	define BEZIER_BEZIER_EVALUATION_REPORT() printf("bezier_bezier_intersections iterations: %d\n", BEZIER_BEZIER_iterations)

#else
#	define BEZIER_BEZIER_COUNT()
#	define BEZIER_BEZIER_EVALUATION_RESET()
#	define BEZIER_BEZIER_EVALUATION_REPORT()
#endif


/**
 * Intersections between bezier curve and a line segment from m to n.
 * L(s) = m + (n-m)* s
 */
int bezier_line_segment_intersections_iterative(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3, const glm::vec2& m, const glm::vec2& n, float tolerance, glm::vec2 results[3]);




static inline void bezier_point(float f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2& p) {
	assert(0 <= f && f <= 1.0);
	p = mix(p0, p1, f);
}
static inline void bezier_point_highp(double f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2& p) {
	assert(0 <= f && f <= 1.0);

	glm::dvec2 P0 = mix(p0, p1, f);
	glm::dvec2 P1 = mix(p1, p2, f);

	p = mix(P0, P1, f);
}
static inline void bezier_point(float f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2& p) {
	assert(0 <= f && f <= 1.0);

	glm::vec2 P0 = mix(p0, p1, f);
	glm::vec2 P1 = mix(p1, p2, f);

	p = mix(P0, P1, f);
}


static inline void bezier_point_highp(double f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3, glm::vec2& p) {
	assert(0 <= f && f <= 1.0);

	glm::dvec2 dp0 = p0;
	glm::dvec2 dp1 = p1;
	glm::dvec2 dp2 = p2;
	glm::dvec2 dp3 = p3;

	glm::dvec2 P1 = mix(dp0, dp1, f);
	glm::dvec2 P2 = mix(dp1, dp2, f);
	glm::dvec2 P3 = mix(dp2, dp3, f);

	glm::dvec2 P4 = mix(P1, P2, f);
	glm::dvec2 P5 = mix(P2, P3, f);

	p = mix(P4, P5, f);
}


static inline void bezier_point(float f, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3, glm::vec2& p) {
	assert(0 <= f && f <= 1.0);

	glm::vec2 P1 = mix(p0, p1, f);
	glm::vec2 P2 = mix(p1, p2, f);
	glm::vec2 P3 = mix(p2, p3, f);

	glm::vec2 P4 = mix(P1, P2, f);
	glm::vec2 P5 = mix(P2, P3, f);

	p = mix(P4, P5, f);
}


static inline bool bezier_equal_point(float t, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3,
			glm::vec2 const & v,float tolerance)
{

	glm::vec2 w;
	bezier_point(t, p0,p1,p2,p3,w);
	return about_equal(v,w,tolerance);
}

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T>
static inline double bezier_point_closest_point_t (const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3, const glm::vec2& v, float tolerance, glm::vec2& result_v, const CO_DOMAIN_T domain) {
	double result_t;
	double result_distance;
	glm::vec2 tmp;

	//
	// Initialise results by evaluating distance to the
	// edges of the domain. Closer edge wins.
	//

	result_t = domain.re_min;
	bezier_point(result_t, p0,p1,p2,p3, result_v);
	result_distance = length(v - result_v);

	bezier_point(domain.re_max, p0,p1,p2,p3, tmp);
	double distance = length(v - tmp);
	if (distance < result_distance) {
		result_t = domain.re_max;
		result_v = tmp;
		result_distance = distance;
	}

	//
	// Now calculate all closest points to v in the given
	// domain.
	//


	// A point p on B(t) closest to v defines line L(s)=p+s(v-p)
	// with an angle a(t) of 90° to the tangent in p.
	// Thus, to find p, we can solve
	//     a(t) = B'(t) * (B(t)-v) with a(t) = 0
	// This is what happens below.

	glm::vec2 a = p3 - 3.f*p2 + 3.f*p1 - p0;
	glm::vec2 b = 3.f * (p2 - 2.f*p1 + p0);
	glm::vec2 c = 3.f * (p1 - p0);

	glm::vec2 p0_v = p0-v;

	double param[6];

	param[0] = 3.f*scalar(a,a);                    // t^5
	param[1] = 5.f*scalar(a,b);                    // t^4
	param[2] = 4.f*scalar(a,c) + 2.f*scalar(b,b);  // t^3
	param[3] = 3.f*(scalar(b,c) + scalar(a,p0_v)); // t^2
	param[4] = scalar(c,c) + 2.f*scalar(b,p0_v);   // t^1
	param[5] = scalar(c,p0_v);                     // t^0

	polynom_t<5> W(param);


	// FIXME: calculate better tolerance and find method to improve if not met
	//
	// tolerance of W.roots() is expressed in respect to variable t
	// but we are looking for the tolerance in values of B(t).
	// Thus, tolerance range for t has to be scaled down by the
	// overall length of the curve.
	// t_tolerance = v_tolerance / length(B)
	// To shorten it up, we approximate the length by the
	// sum of lengths of the bounding box cage, which actually
	// gives us a higher accuracy than needed.
	float len = length(p1-p0)+length(p2-p1)+length(p3-p2);
	float t_tolerance = tolerance/len;


	double t_i[5];
	int count = W.roots(t_i, t_tolerance, domain);
	for (int i = 0; i < count; i++) {
		double t = t_i[i];
		t = domain.clip(t);
		bezier_point(t, p0,p1,p2,p3, tmp);
		double distance = length( v - tmp);
		if (distance < result_distance) {
			result_v = tmp;
			result_t = t;
			result_distance = distance;
		}
	}
	return result_t;
}


/**
 * Determines the parameters of a line L(s)=m*s + p, which is a tangent of the bezier curve P(f) at P(f_t).
 *
 * p = P(f_t)
 * m = P'(f_t)
 */
static inline void bezier_tangent(double f_t, glm::vec2 const & p0, glm::vec2 const & p1, glm::vec2 const & p2, glm::vec2 const & p3, glm::vec2& p, glm::vec2& m) {
	// p = P(f_t)
	bezier_point_highp(f_t, p0, p1, p2, p3, p);

	// m = P'(f_t)
	bezier_point_highp(f_t, p1-p0, p2-p1, p3-p2, m);
	// m *= 3.f; <-- doesn't change anything
}



static inline void bezier_split_highp(double f, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, glm::vec2& q0, glm::vec2& q1, glm::vec2& q2, glm::vec2& q3) {
	// using q3 as temporary variable
	assert(0 <= f && f <= 1.0);

	glm::dvec2 dp0 = p0;
	glm::dvec2 dp1 = p1;
	glm::dvec2 dp2 = p2;
	glm::dvec2 dp3 = p3;

	double g = 1.0-f;

	glm::dvec2 P1a = mix(dp0, dp1, f);
	glm::dvec2 P1b = mix(dp1, dp0, g);
	glm::dvec2 P1 = mix(P1a, P1b, 0.5);

	glm::dvec2 P2a = mix(dp1, dp2, f);
	glm::dvec2 P2b = mix(dp2, dp1, g);
	glm::dvec2 P2 = mix(P2a, P2b, 0.5);

	glm::dvec2 P3a = mix(dp2, dp3, f);
	glm::dvec2 P3b = mix(dp3, dp2, g);
	glm::dvec2 P3 = mix(P3a, P3b, 0.5);

	glm::dvec2 P4a = mix(P1, P2, f);
	glm::dvec2 P4b = mix(P2, P1, g);
	glm::dvec2 P4  = mix(P4a, P4b, 0.5);

	glm::dvec2 P5a = mix(P2, P3, f);
	glm::dvec2 P5b = mix(P3, P2, g);
	glm::dvec2 P5 = mix(P5a, P5b, 0.5);

	glm::dvec2 P6a = mix(P4, P5, f); // P6 -> splitting point
	glm::dvec2 P6b = mix(P5, P4, g); // P6 -> splitting point
	glm::dvec2 P6 = mix(P6a, P6b, 0.5); // P6 -> splitting point

	q0 = P6;
	q1 = P5;
	q2 = P3;
	q3 = p3; // q3 = p3
	p3 = P6; // p3 = P6
	p2 = P4;
	p1 = P1;
	// p0 unchanged

}





static inline void bezier_split(float f, glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, glm::vec2& q0, glm::vec2& q1, glm::vec2& q2, glm::vec2& q3) {
	// using q3 as temporary variable
	glm::vec2& px = q3;
	
	px = p1;
	p1 = mix(p0, px, f); // P1 = mix(p0, p1, f);
	px = mix(px, p2, f); // P2 = mix(p1, p2, f);
	q2 = mix(p2, p3, f); // P3 = mix(p2, p3, f);

	p2 = mix(p1, px, f); // P4 = mix(P1, P2, f);
	q1 = mix(px, q2, f); // P5 = mix(P2, P3, f);

	q0 = mix(p2, q1, f); // P6 -> splitting point

	// q0 = P6;
	// q1 = P5;
	// q2 = P3;
	q3 = p3; // q3 = p3
	p3 = q0; // p3 = P6
	// p2 = P4;
	// p1 = P1;
	// p0 unchanged

}



/**
 * Most accurate version of bezier_merge.
 *
 * Accuracy is about 10^4 times higher than bezier_merge().
 * Requires about 4 times more processing time than bezier_merge().
 *
 * @param tolerance Absolute tolerance of control point deviation (e.g. 0.5 if control points are ivec2 actually)
 */
static inline bool bezier_merge_highp(glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, glm::vec2 q0, glm::vec2 q1, glm::vec2 q2, glm::vec2 q3, float tolerance) {
	assert(tolerance >= 0);
	if (p3 != q0) {
		return false;
	}

	/*
	 * We calculate the control points of the merged bezier curve
	 * based on the given control points either in forward
	 * or reverse direction.
	 * We perform three tests to determine, whether the resulting
	 * curve actually represents a merge of the two given curves:
	 * 1. cannot merge, if p3 and q0 are not the same.
	 * 2. cannot merge, if p3 and q0 are not on a line between p2 and q1
	 * 3. cannot merge, if resulting control point p3/q0 is not equal to q3/p0.
	 *
	 * All these tests are quite stable against precision errors.
	 */

	glm::dvec2 dp0 = p0;
	glm::dvec2 dq3 = q3;
	glm::dvec2 P6 = p3; // == q0
	glm::dvec2 P5 = q1;
	glm::dvec2 P4 = p2;
	glm::dvec2 P3 = q2;
	glm::dvec2 P2; // = ?
	glm::dvec2 P1 = p1;

	double f;
	bool canMerge = line_contains_point_highp(P4,P5,P6,f, tolerance);
	canMerge = canMerge && (0.0 < f && f < 1.0);
	if (canMerge) {

		/*
		 * f can get very large if the split was performed with a
		 * very small interpolation factor. A very large factor
		 * causes very large deviations.
		 *
		 * Since we have the option to perform our calculations
		 * in reverse order, using the second set of control
		 * points,  we will perform calculations in that direction,
		 * which provides the most accuracy.
		 */
		if (f > 0.5) {
			f = 1.0/f;              // forward interpolation factor
			glm::dvec2 P2  = mix(P1,P4,f);
			glm::dvec2 dp1 = mix(dp0, P1, f); // p0 + (P1-p0) * f;
			glm::dvec2 dp2 = mix(dp1, P2, f); //
			glm::dvec2 dp3 = mix(dp2, P3, f);
			canMerge = about_equal(dp3,dq3, tolerance);
			p1 = dp1;
			p2 = dp2;
			p3 = q3;
		} else {
			f = 1.0/(1.0-f);        // reverse interpolation factor
			glm::dvec2 P2  = mix(P3,P5,f);
			glm::dvec2 dq2 = mix(dq3,P3,f); // p0 + (P1-p0) * f;
			glm::dvec2 dq1 = mix(dq2,P2,f); //
			glm::dvec2 dq0 = mix(dq1,P1,f);
			canMerge = about_equal(dq0,dp0, tolerance);
			p1 = dq1;
			p2 = dq2;
			p3 = q3;
		}
	}
	return canMerge;
}

static inline bool bezier_merge(glm::vec2& p0, glm::vec2& p1, glm::vec2& p2, glm::vec2& p3, glm::vec2 q0, glm::vec2 q1, glm::vec2 q2, glm::vec2 q3) {
	if (p3 != q0) return false;

	double f;
	glm::vec2 P6 = p3; // == q0
	glm::vec2 P5 = q1;
	glm::vec2 P4 = p2;
	glm::vec2 P3 = q2;
	glm::vec2 P2; // = ?
	glm::vec2 P1 = p1;
	bool canMerge = line_contains_point_highp(P4,P5,P6,f, 0.001);
	double g = 1.0/(1.0-f); // reverse direction factor
	f = 1.0/f;              // forward direction factor
	if (canMerge) {
		glm::vec2 P2a = mix(P1, P4, f);
		glm::vec2 P2b = mix(P3, P5, g);
		canMerge = about_equal(P2a, P2b, BEZIER_FLOAT_PRECISION);
		if (canMerge) {
			P2 = mix(P2a,P2b,0.5);
			p1 = mix(p0,P1,f); // p0 + (P1-p0) * f;
			p2 = mix(q3,P3,g); //
			p3 = q3;
		}
	}
	return canMerge;
}


}


#endif /* MATH2D_BEZIER_H_ */
