/*
 * general.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH2D_UTILS_H_
#define MATH2D_UTILS_H_

#include "math2d-config.h"

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
 * This function forces a reinterpret_cast, by
 * reinterpreting the address (pointer) of the given item.
 * Thus, it intentionally circumvents restrictions of C++ to
 * achieve compile-time casts.
 *
 * Useful for example, if you have a vector<A*> and a class B derived from A,
 * and you know, that every item in vector<A*> is actually a B*. Then you can cast
 * vector<B*>& vector_B = forced_cast<vector<B*>>(vector_A);
 */
template <class T>
static inline constexpr T& forced_cast(auto& item) {
	assert(sizeof(item) == sizeof(T));
	// cast reference to void pointer and then to T* and then returns T&
	// return (*((T*)( (void*)(&item) )));
	return (*reinterpret_cast<T*>( (void*)(&item) ));
}



}; // namespace math


#endif /* MATH2D_UTILS_H_ */
