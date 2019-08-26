/*
 * Bezier.h
 *
 *  Created on: 9 Jul 2019
 *      Author: homac
 */

#ifndef MATH_BEZIER_H_
#define MATH_BEZIER_H_

#include "config.h"

#include <assert.h>
#include <glm/glm.hpp>
#include <glm/common.hpp>
using namespace glm;

#include "math/line.h"
#include "float-utils.h"

namespace math {


static const double_mantissa_mask_t BEZIER_DOUBLE_PRECISION = double_mantissa_mask(16);
static const float_mantissa_mask_t BEZIER_FLOAT_PRECISION = float_mantissa_mask(8);

static inline void bezier_point(float f, vec2 const & p0, vec2 const & p1, vec2& p) {
	assert(0 <= f && f <= 1.0);
	p = mix(p0, p1, f);
}

static inline void bezier_point(float f, vec2 const & p0, vec2 const & p1, vec2 const & p2, vec2& p) {
	assert(0 <= f && f <= 1.0);

	vec2 P0 = mix(p0, p1, f);
	vec2 P1 = mix(p1, p2, f);

	p = mix(P0, P1, f);
}

static inline void bezier_point_highp(double f, vec2 const & p0, vec2 const & p1, vec2 const & p2, vec2 const & p3, vec2& p) {
	assert(0 <= f && f <= 1.0);

	dvec2 dp0 = p0;
	dvec2 dp1 = p1;
	dvec2 dp2 = p2;
	dvec2 dp3 = p3;

	dvec2 P1 = mix(dp0, dp1, f);
	dvec2 P2 = mix(dp1, dp2, f);
	dvec2 P3 = mix(dp2, dp3, f);

	dvec2 P4 = mix(P1, P2, f);
	dvec2 P5 = mix(P2, P3, f);

	p = mix(P4, P5, f);
}


static inline void bezier_point(float f, vec2 const & p0, vec2 const & p1, vec2 const & p2, vec2 const & p3, vec2& p) {
	assert(0 <= f && f <= 1.0);

	vec2 P1 = mix(p0, p1, f);
	vec2 P2 = mix(p1, p2, f);
	vec2 P3 = mix(p2, p3, f);

	vec2 P4 = mix(P1, P2, f);
	vec2 P5 = mix(P2, P3, f);

	p = mix(P4, P5, f);
}

static inline void bezier_split_highp(double f, vec2& p0, vec2& p1, vec2& p2, vec2& p3, vec2& q0, vec2& q1, vec2& q2, vec2& q3) {
	// using q3 as temporary variable
	assert(0 <= f && f <= 1.0);

	dvec2 dp0 = p0;
	dvec2 dp1 = p1;
	dvec2 dp2 = p2;
	dvec2 dp3 = p3;

	double g = 1.0-f;

	dvec2 P1a = mix(dp0, dp1, f);
	dvec2 P1b = mix(dp1, dp0, g);
	dvec2 P1 = mix(P1a, P1b, 0.5);

	dvec2 P2a = mix(dp1, dp2, f);
	dvec2 P2b = mix(dp2, dp1, g);
	dvec2 P2 = mix(P2a, P2b, 0.5);

	dvec2 P3a = mix(dp2, dp3, f);
	dvec2 P3b = mix(dp3, dp2, g);
	dvec2 P3 = mix(P3a, P3b, 0.5);

	dvec2 P4a = mix(P1, P2, f);
	dvec2 P4b = mix(P2, P1, g);
	dvec2 P4  = mix(P4a, P4b, 0.5);

	dvec2 P5a = mix(P2, P3, f);
	dvec2 P5b = mix(P3, P2, g);
	dvec2 P5 = mix(P5a, P5b, 0.5);

	dvec2 P6a = mix(P4, P5, f); // P6 -> splitting point
	dvec2 P6b = mix(P5, P4, g); // P6 -> splitting point
	dvec2 P6 = mix(P6a, P6b, 0.5); // P6 -> splitting point

	q0 = P6;
	q1 = P5;
	q2 = P3;
	q3 = p3; // q3 = p3
	p3 = P6; // p3 = P6
	p2 = P4;
	p1 = P1;
	// p0 unchanged

}


static inline void bezier_split(float f, vec2& p0, vec2& p1, vec2& p2, vec2& p3, vec2& q0, vec2& q1, vec2& q2, vec2& q3) {
	// using q3 as temporary variable
	vec2& px = q3;
	
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
static inline bool bezier_merge_highp(vec2& p0, vec2& p1, vec2& p2, vec2& p3, vec2 q0, vec2 q1, vec2 q2, vec2 q3, float tolerance) {
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

	dvec2 dp0 = p0;
	dvec2 dq3 = q3;
	dvec2 P6 = p3; // == q0
	dvec2 P5 = q1;
	dvec2 P4 = p2;
	dvec2 P3 = q2;
	dvec2 P2; // = ?
	dvec2 P1 = p1;

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
			dvec2 P2  = mix(P1,P4,f);
			dvec2 dp1 = mix(dp0, P1, f); // p0 + (P1-p0) * f;
			dvec2 dp2 = mix(dp1, P2, f); //
			dvec2 dp3 = mix(dp2, P3, f);
			canMerge = about_equal(dp3,dq3, tolerance);
			p1 = dp1;
			p2 = dp2;
			p3 = q3;
		} else {
			f = 1.0/(1.0-f);        // reverse interpolation factor
			dvec2 P2  = mix(P3,P5,f);
			dvec2 dq2 = mix(dq3,P3,f); // p0 + (P1-p0) * f;
			dvec2 dq1 = mix(dq2,P2,f); //
			dvec2 dq0 = mix(dq1,P1,f);
			canMerge = about_equal(dq0,dp0, tolerance);
			p1 = dq1;
			p2 = dq2;
			p3 = q3;
		}
	}
	return canMerge;
}

static inline bool bezier_merge(vec2& p0, vec2& p1, vec2& p2, vec2& p3, vec2 q0, vec2 q1, vec2 q2, vec2 q3) {
	if (p3 != q0) return false;

	double f;
	vec2 P6 = p3; // == q0
	vec2 P5 = q1;
	vec2 P4 = p2;
	vec2 P3 = q2;
	vec2 P2; // = ?
	vec2 P1 = p1;
	bool canMerge = line_contains_point_highp(P4,P5,P6,f, 0.001);
	double g = 1.0/(1.0-f); // reverse direction factor
	f = 1.0/f;              // forward direction factor
	if (canMerge) {
		vec2 P2a = mix(P1, P4, f);
		vec2 P2b = mix(P3, P5, g);
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




int bezier_line_segment_intersections_iterative(const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3, const vec2& m, const vec2& n, float tolerance, vec2 p[3]);
int bezier_line_segment_intersections(const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3, const vec2& m, const vec2& n, float tolerance, vec2 p[3]);

bool bezier_inflection_point(const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3, float& t);
bool bezier_extrema(vec2 const& p0, vec2 const& p1, vec2 const& p2, vec2 const& p3, dvec2& t_min, dvec2& t_max);

}


#endif /* MATH_BEZIER_H_ */
