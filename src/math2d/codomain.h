/*
 * codomain.h
 *
 *  Created on: 30 Sep 2019
 *      Author: homac
 */

#ifndef MATH2D_CODOMAIN_H_
#define MATH2D_CODOMAIN_H_

#include "math2d-config.h"


namespace math2d {





/**
 * Interval boundary declarator comparing given values against stored boundary.
 * operation: v >= bound
 */
template <typename T = double>
struct LOIN {
	T bound;
	LOIN(T v):bound(v){}
	bool operator ()(T v) const {
		return bound <= v;
	}
	operator T() const {return bound;}
	LOIN& operator = (T _value) {
		bound = _value;
		return *this;
	}
	LOIN& operator = (const LOIN& other) {
		return operator = (other.bound);
	}
};

/**
 * Interval boundary declarator comparing given values against stored boundary.
 * operation: v > bound
 */
template <typename T = double>
struct LOEX {
	T bound;
	LOEX(T v):bound(v){}
	bool operator ()(T v) const {
		return bound < v;
	}
	operator T() const {return bound;}
	LOEX& operator = (T _value) {
		bound = _value;
		return *this;
	}
	LOEX& operator = (const LOEX& other) {
		return operator = (other.bound);
	}
};

/**
 * Interval boundary declarator comparing given values against stored boundary.
 * operation: v <= bound
 */
template <typename T = double>
struct UPIN {
	T bound;
	UPIN(T v):bound(v){}
	bool operator ()(T v) const {
		return v <= bound;
	}
	operator T() const {return bound;}
	UPIN& operator = (T _value) {
		bound = _value;
		return *this;
	}
	UPIN& operator = (const UPIN& other) {
		return operator = (other.bound);
	}
};

/**
 * Interval boundary declarator comparing given values against stored boundary.
 * operation: v < value
 */
template <typename T = double>
struct UPEX {
	T bound;
	UPEX(T v):bound(v){}
	bool operator ()(T v) const {
		return v < bound;
	}
	operator T() const {return bound;}

	UPEX& operator = (T _value) {
		bound = _value;
		return *this;
	}
	UPEX& operator = (const UPEX& other) {
		return operator = (other.bound);
	}
};



/**
 *
 * Generally, a co-domain describes a set of valid values. We
 * use a simpler definition, where the domain is just an interval.
 *
 * In this implementation, we define boundaries by a value
 * and a comparison operator (one of GT, GE, LE, LT).
 * At initialisation, the comparator is given a value, which it
 * uses to check whether its boundary condition is met.
 *
 * A co-comain with default parameters accepts all finite double
 * values. Example:
 * <pre>
 *   co_domain_t<> finite_double;
 *   if (finite_double(0.0)) printf("yep");
 *   if (!finite_double(INFINITY)) printf("nope");
 * </pre>
 *
 * The type of accepted values defaults to double.
 * Thus, a simple [0,1] boundaries inclusive interval for
 * double values can be defined like this:
 * <pre>
 *    co_domain_t<> interval(LOIN<>(0),UPIN<>(1));
 * </pre>
 *
 * Since lower and upper boundaries default to LOIN and UPIN
 * the next declaration will does the same as the one above:
 * <pre>
 *    co_domain_t<> interval(0,1);
 * </pre>
 *
 * A boundary exclusive interval for const float values would
 * look like this:
 * <pre>
 *    const co_domain_t<const float> interval(LOEX<const float>(0), UPEX<const float>(1));
 * </pre>
 *
 * This can be made more readable when typedefs are used:
 * <pre>
 *    typedef const co_domain_t<const float, LOEX<const float>, UPEX<const float>> co_domain_const_float_excl_t;
 *    co_domain_const_float_excl_t interval(0,1);
 * </pre>
 */
template<typename T = double, class LOWER = LOIN<T>, class UPPER = UPIN<T>>
struct co_domain_t {
	LOWER lower;
	UPPER upper;

	co_domain_t(
			LOWER _min = LOIN<T>(-DBL_MAX),
			UPPER _max = UPIN<T>(DBL_MAX)
	)
	: lower(_min) , upper(_max)
	{}

	bool operator() (T re) const {
		return lower(re) && upper(re);
	}

	/** clip v to the boundary, which was exceeded. */
	T clip(T v) const {
		return !lower(v) ? T(lower) : !upper(v) ? T(upper) : v;
	}

};

/** co-domain for finite double values */
typedef co_domain_t<double> co_domain_std_t;
/** co-domain for finite double values including boundaries */
typedef co_domain_std_t co_domain_std_inclusive_t;
/** co-domain for finite double values excluding boundaries */
typedef co_domain_t<double, LOEX<>, UPEX<>> co_domain_std_exclusive_t;
/** const co-domain for finite double values */
typedef const co_domain_t<const double> co_domain_const_t;
/** const co-domain for finite double values including boundaries */
typedef const co_domain_t<const double, const LOIN<const double>, const UPIN<const double>> co_domain_const_inclusive_t;
/** const co-domain for finite double values excluding boundaries */
typedef const co_domain_t<const double, const LOEX<const double>, const UPEX<const double>> co_domain_const_exclusive_t;

/** const co-domain for finite double values */
typedef co_domain_const_t CO_DOMAIN_REAL_T;
const CO_DOMAIN_REAL_T CO_DOMAIN_REAL;

/** const co-domain for double values in interval [-1,1] inclusive */
typedef co_domain_const_inclusive_t CO_DOMAIN_REAL_IN_NEG_1_POS_1_T;
const CO_DOMAIN_REAL_IN_NEG_1_POS_1_T CO_DOMAIN_REAL_IN_NEG_1_POS_1(-1,1);

/** const co-domain for double values in interval [0,1] inclusive */
typedef co_domain_const_inclusive_t CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T;
const CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T CO_DOMAIN_REAL_IN_0_1_INCLUSIVE(0,1);

/** const co-domain for double values in interval ]0,1[ exclusive */
typedef co_domain_const_exclusive_t CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE_T;
const CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE_T CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE(0,1);



} // math2d


#endif /* MATH2D_CODOMAIN_H_ */
