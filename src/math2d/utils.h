/*
 * general.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH2D_UTILS_H_
#define MATH2D_UTILS_H_


#include "math2d_config.h"

#include <limits.h>
#include <float.h>
#include <math.h>

#include <glm/glm.hpp>

#include "math2d/codomain.h"


namespace math2d {


// higher accuracy than math.h : M_PI
static const long double CONSTANT_PI = 3.14159265358979323846264338327950288;



/**
 * Determines Z-component of the cross product of vectors in R^2
 * If a and b are normalised, then cross_z equals the cosine.
 */
template <typename T, glm::precision P>
static inline T cross_z(glm::tvec2<T, P> const & a, glm::tvec2<T, P> const & b)
{
	return a.x * b.y - a.y * b.x;
}

/**
 * scalar product (a*b).
 */
template <typename T, glm::precision P>
static inline T scalar(glm::tvec2<T, P> const & a, glm::tvec2<T, P> const & b)
{
	return a.x * b.x + a.y * b.y;
}

template <typename T, glm::precision P>
static inline bool equal(glm::tvec2<T, P> const & a, glm::tvec2<T, P> const & b)
{
	return a == b;
}

/**
 * signum function
 * signum(v) == 0 if v == 0
 * signum(v) == 1 if v > 0
 * signum(v) == -1 if v < 0
 * @return
 */
template <typename T>
T signum(T v) {
    return T((T(0) < v) - (v < T(0)));
}

}; // namespace math


#endif /* MATH2D_UTILS_H_ */
