/*
 * math_line.h
 *
 *  Created on: 17 Jun 2019
 *      Author: homac
 */

#ifndef MATH_LINE_H_
#define MATH_LINE_H_


#include <glm/geometric.hpp>
#include <glm/glm.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <math.h>
#include <ieee754.h>

#include "math/general.h"
#include "math/math2d.h"
#include "float-utils.h"


namespace math {





#define PSEUDO_PI      (2.0)
#define PSEUDO_PI_HALF (1.0)

/**
 * This function calculates a fast pseudo ]-180°,180°] angle
 * between given dir and unit vector (0,1).
 *
 * The pseudo angle is a mapping of the 360 angle alpha to a
 * non-periodic continuously rising function f(alpha) in
 * the range ]-180°,180°]. This is useful to fast sort
 * vectors by angle, when the actual angle is
 * not needed.
 * <br/>
 * This function is 3 to 5 times faster than glm::orientedAngle()
 * and provides ~10^5 times more accurate differentiation
 * between angles of vec2<float>.
 */
static inline double pseudo_orientation_e2(glm::vec2 dir) {
	// reference vector (unit vector e2)
	// glm::vec2 e_y(0,1);

	// dir = normalize(dir);

	double len = sqrt(double(dir.x)*dir.x + double(dir.y)*dir.y);

	// sin of angle between dir and e_y
	// gives information if angle is in left (<0) or right half plane
	double _sin_angle = double(dir.x)/len;
	// cos of angle between dir and e_y (same as sin of cross_z(e_x, dir))
	// tells if dir points in upper or lower half plane
	double _cos_angle = double(dir.y)/len;


	//	 The mapping function f(alpha) is defined as
	//
	//	  					{-2 -sin(alpha), if alpha in ]-180, -90]
	//	  		f(alpha) = {	 sin(alpha), if alpha in ] -90, +90]
	//	 					{ 2 -sin(alpha), if alpha in ] +90,+180]
	//
	// Thus, f(alpha) can be written as
	//
	// 			f(alpha) = a + b * sin(alpha)
	//
	// where a and b have the following values based on sin(alpha) and cos(alpha).
	//
	//		a		b	range		sin			cos
	//		-2		-1	]-180,-90]	]-1, 0]		]0, -1]
	//		 0		 0	]-90,0,+90]	]0,1,0]		]-1,0,+1]
	//		+2		-1	] +90,180]	] 0,-1]		]+1, 0]
	//


	// To determine a and b *without branching*,
	// we use the sign of sin(alpha) and cos(alpha)
	// to calculate them.
	double sy = copysign_only(_sin_angle);
	double sx = copysign_only(_cos_angle);

	double a = sy-sx*sy;
	double b = sx;

	double angle = a + b *_sin_angle;

	return angle;
}
/**
 * reference vector: unit vector e1 = (1,0)
 */
static inline double pseudo_orientation_e1(glm::vec2 dir) {

	double len = sqrt(double(dir.x)*dir.x + double(dir.y)*dir.y);

	double _sin_angle = -double(dir.y)/len;
	double _cos_angle = double(dir.x)/len;

	double sy = copysign_only(_sin_angle);
	double sx = copysign_only(_cos_angle);

	double a = sy-sx*sy;
	double b = sx;

	double angle = a + b *_sin_angle;

	return angle;
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
 * @param s out: distance from m to p relative to direction and length of line (m->n).
 * @return true, only if p on m->n
 */
bool line_contains_point(const glm::vec2& m, const glm::vec2& n, const glm::vec2& p, float& s);
bool line_contains_point_highp(const glm::dvec2& m, const glm::dvec2& n, const glm::dvec2& p, double& s, float tolerance);



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
bool line_segment_superpose(const glm::vec2& m, const glm::vec2& n, const glm::vec2& u, const glm::vec2& v, glm::vec2& p, glm::vec2& q);



/**
 * For a given horizontal ray, starting in s and going to the left e(-1,0),
 * and a given segment, starting in m and going to n, where m is known to
 * be below s (m.y < s.y), this function tests, if the segment crosses the
 * ray.
 *
 * NOTE: This function does not test
 * - segments which lie upon the ray
 * - segments which cross s
 *
 */
bool ray_left_intersects_segment_where_m_below(const glm::vec2& s, const glm::vec2& m, const glm::vec2& n);

/**
 * test line segment m->n for intersection with horizontal scan line
 *
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_x is set to intersection point p(p_x,s_y)
 */
int scanline_horizontal_intersection_m_below(const float s_y, glm::vec2 m, glm::vec2 n, float& p_x);

/**
 * test line segment m->n for intersection with horizontal scan line
 *
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_x is set to intersection point p(p_x,s_y)
 */
int scanline_horizontal_intersection(const float s_y, glm::vec2 m, glm::vec2 n, float& p_x);

/**
 * test line segment m->n for intersection with horizontal scan line
 *
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_x is set to intersection point p(s_x,p_y)
 */
int scanline_vertical_intersection_m_left(const float s_x, glm::vec2 m, glm::vec2 n, float& p_y);

/**
 * test line segment m->n for intersection with horizontal scan line
 *
 * returns -1 if no intersection exists.
 * returns 0 if m->n is a horizontal line, on the scan line
 * returns 1 if intersection exists and p_x is set to intersection point p(s_x,p_y)
 */
int scanline_vertical_intersection(const float s_x, glm::vec2 m, glm::vec2 n, float& p_y);



enum LineSegmentRelationship {
	/** lines are parallel, but not otherwise related */
	LINES_PARALLEL  = -1,
	/** none of the relationships exists */
	LINES_UNRELATED = 0,
	/** lines intersect each other */
	LINES_INTERSECT = 1,
	/** line segments superpose each other */
	LINES_SUPERPOSE = 2,
	/** line segments are adjacent to each other (one point of both lines equal) */
	LINES_ADJACENT  = 3,
	/** one line segment touches other line segment in one point */
	LINES_TOUCHING  = 4,
};

/**
 * Determines realtionship between two line segments.
 * Result:
 * <ul>
 * <li>-1: parallel to each other but not superposing</li>
 * <li> 0: unrelated (disjunct)</li>
 * <li> 1: intersect each other</li>
 * <li> 2: superpose each other</li>
 * </ul>
 */
LineSegmentRelationship line_segments_relationship(const glm::vec2& m, const glm::vec2& n, const glm::vec2& u, const glm::vec2& v);

/**
 * Find intersection p of two segments (s1={m -> n} and s2={u -> v}) in R^2.
 *
 * ATTENTION ACCURACY: intersection(m,n, u,v, p1,s) and intersection(u,v, m,n, p2,s)
 *                      may give slightly different intersection points (p1 != p2)!
 *                      Use line_segment_intersection_highp for higher precision.
 *
 * A segment is defined by a straight line between two separate points.
 * Only if the intersection exists the function produces the following output:
 * Return value is true, intersection point is written to p, and s is set to s=length(p-m)/length(n-m).
 * s is a side product of internal calculation, so there is no extra effort to do so.
 *
 * There is no intersection between a point and a segment, and congruent (superpose)
 * segments do not create an intersection either (they superpose each other).
 * @return true if intersection exists, otherwise return false and also if segments are parallel to each other
 */
bool line_segment_intersection(const glm::vec2& m, const glm::vec2& n, const glm::vec2& u, const glm::vec2& v, glm::vec2& p, float& s);
/**
 * This function is a twin of line_segment_intersection.
 *
 * All internal calculations are performed on double precision floating
 * point registers, to achieve higher accuracy than line_segment_intersection().
 * Since inputs and intersection point p are given in single precision floating
 * point values, this function will give the same results for p, regardless of
 * the order of input parameters, unlike line_segments_intersection. Furthermore,
 * s will be the same, if casted to float afterwards.
 *
 * @see line_segment_intersection
 */
bool line_segment_intersection_highp(const glm::vec2& m, const glm::vec2& n, const glm::vec2& u, const glm::vec2& v, glm::vec2& p, double& s);
/**
 * Intersection point of two lines (not line segments!).
 */
bool line_intersection_highp(const glm::dvec2& m, const glm::dvec2& n, const glm::dvec2& u, const glm::dvec2& v, glm::dvec2& p, double& s);


}; // namespace math



#endif /* MATH_LINE_H_ */
