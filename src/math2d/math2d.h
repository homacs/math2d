/*
 * math2d.h
 *
 *  Created on: 26 Jul 2019
 *      Author: homac
 */

#ifndef MATH2D_MATH2D_H_
#define MATH2D_MATH2D_H_
#include "float-utils.h"

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
 * Determines scalar product (a*b).
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


template <glm::precision P>
static inline glm::tvec2<float, P> precision_trunc(glm::tvec2<float, P> a, float_mantissa_mask_t precision)
{
	a.x = float_mantissa_trunc(a.x, precision);
	a.y = float_mantissa_trunc(a.y, precision);
	return a;
}

template <glm::precision P>
static inline glm::tvec2<float, P> precision_round(glm::tvec2<float, P> a, float_mantissa_mask_t precision)
{
	a.x = float_mantissa_round(a.x, precision);
	a.y = float_mantissa_round(a.y, precision);
	return a;
}

template <glm::precision P>
static inline glm::tvec2<double, P> precision_trunc(glm::tvec2<double, P> a, double_mantissa_mask_t precision)
{
	a.x = double_mantissa_trunc(a.x, precision);
	a.y = double_mantissa_trunc(a.y, precision);
	return a;
}

template <glm::precision P>
static inline glm::tvec2<double, P> precision_round(glm::tvec2<double, P> a, double_mantissa_mask_t precision)
{
	a.x = double_mantissa_round(a.x, precision);
	a.y = double_mantissa_round(a.y, precision);
	return a;
}


/**
 * Determines scalar product (a*b).
 */
template <glm::precision P>
static inline constexpr bool about_equal(glm::tvec2<float, P> const & a, glm::tvec2<float, P> const & b, float tolerance)
{
	return glm::length(a-b) <= tolerance;
}



/**
 * Determines scalar product (a*b).
 */
template <glm::precision P>
static inline bool about_equal(glm::tvec2<float, P> const & a, glm::tvec2<float, P> const & b, float_mantissa_mask_t precision)
{
	return float_mantissa_trunc(glm::length(a-b), precision) == 0;
}

/**
 * Determines scalar product (a*b).
 */
template <glm::precision P>
static inline bool about_equal(glm::tvec2<double, P> const & a, glm::tvec2<double, P> const & b, double_mantissa_mask_t precision)
{
	return double_mantissa_trunc(glm::length(a-b), precision) == 0;
}

template <glm::precision P>
static inline bool about_equal(glm::tvec2<double, P> const & a, glm::tvec2<double, P> const & b, double tolerance)
{
	return glm::length(a-b) <= tolerance;
}



#endif /* MATH2D_MATH2D_H_ */
