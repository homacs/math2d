/*
 * precision.h
 *
 *  Created on: 30 Sep 2019
 *      Author: homac
 */

#ifndef MATH2D_PRECISION_H_
#define MATH2D_PRECISION_H_


#include "math2d-config.h"

#include <glm/glm.hpp>


#include "float-utils.h"

namespace math2d {



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


template <glm::precision P>
static inline constexpr bool about_equal(glm::tvec2<float, P> const & a, glm::tvec2<float, P> const & b, float tolerance)
{
	return glm::length(a-b) <= tolerance;
}



template <glm::precision P>
static inline bool about_equal(glm::tvec2<float, P> const & a, glm::tvec2<float, P> const & b, float_mantissa_mask_t precision)
{
	return float_mantissa_trunc(glm::length(a-b), precision) == 0;
}

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


static inline bool about_equal(double a, double b, double const & tolerance) {
	assert(tolerance >= 0.0);
	double dist = fabs(a-b);
	return dist <= tolerance;
}

static inline bool about_equal(float a, float b, float const & tolerance) {
	assert(tolerance >= 0.0);
	float dist = fabs(a-b);
	return dist <= tolerance;
}


static inline bool about_equal(double a, double b, double_mantissa_mask_t const & precision) {
	return double_mantissa_round(a, precision) == double_mantissa_round(b, precision);
}

static inline bool about_equal(float a, float b, float_mantissa_mask_t const & precision) {
	return float_mantissa_round(a, precision) == float_mantissa_round(b, precision);
}


} // namespace math2d


#endif /* MATH2D_PRECISION_H_ */
