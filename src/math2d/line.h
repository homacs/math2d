/*
 * math_line.h
 *
 *  Created on: 17 Jun 2019
 *      Author: homac
 */

#ifndef MATH2D_LINE_H_
#define MATH2D_LINE_H_

#include "math2d-config.h"

#include <math.h>

#include <glm/geometric.hpp>
#include <glm/glm.hpp>

#include "float-utils.h"


#include "math2d/utils.h"



namespace math2d {





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



#endif /* MATH2D_LINE_H_ */
