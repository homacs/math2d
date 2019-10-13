/*
 * polynom.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH2D_POLYNOM_H_
#define MATH2D_POLYNOM_H_


#include "math2d-config.h"


#include <stdexcept>
#include <complex>
#include <string.h>
#include <algorithm>
#include <assert.h>
#include <math.h>


#include "float-utils.h"

#include "math2d/precision.h"
#include "math2d/utils.h"

namespace math2d {

template  <int ORDER, typename REAL_T = double>
struct polynom_t;

// TODO: beauty: change all f(t) to f(x), because it's more common
// TODO: beauty: change all polynom#_xxx to polynom_#_xxx
template<typename REAL_T = double>
static inline REAL_T polynom1_value(REAL_T a, REAL_T b, REAL_T t);

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom1_roots(REAL_T a, REAL_T b, REAL_T t[1], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template<typename REAL_T = double>
static inline REAL_T polynom2_value(REAL_T a, REAL_T b, REAL_T c, REAL_T t);

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_roots(REAL_T a, REAL_T b, REAL_T c, REAL_T roots[2], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

/**
 *
 * Searches extrema and saddle point
 * for a given polyom
 * 	  f(t) = a t^2 + b t + c
 * and writes corresponding values for t in t_min, t_max.
 *
 * Term c is omitted in the function parameters, because it is not used.
 *
 * Infinite values in t_min and t_max represent non-existing or invalid extrema.
 * @param a 1. term factor
 * @param b 2. term factor
 * @param t_min receives t_min, may be INFINITY to indicate invalid value.
 * @param t_max receives t_max, may be INFINITY to indicate invalid value.
 * @param domain valid range for t.
 * @return 0: no extrema, 1: minimum or maximum, -1 saddle point t_min==t_max
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_extrema(REAL_T a, REAL_T b, REAL_T& t_min, REAL_T& t_max, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template<typename REAL_T = double>
static inline REAL_T polynom3_value(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T t);

/**
 *
 * Searches extrema and saddle point
 * for a given polyom
 * 	  f(t) = a t^3 + b t^2 + c t + d
 * and writes corresponding values for t in t_min, t_max.
 *
 * Term d is omitted in the function parameters, because it is not used.
 *
 * Infinite values in t_min and t_max represent non-existing or invalid extrema.
 * @param a 1. term factor
 * @param b 2. term factor
 * @param c 3. term factor
 * @param accept_saddle_point whether saddle points have to be considered to be a double extrema or not.
 * @param t_min receives t_min, may be INFINITY to indicate invalid value.
 * @param t_max receives t_max, may be INFINITY to indicate invalid value.
 * @param domain valid range for t.
 * @return 0: no extrema, 1: minimum and/or maximum, -1 saddle point t_min==t_max
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_extrema(REAL_T a, REAL_T b, REAL_T c, bool accept_saddle_point, REAL_T& t_min, REAL_T& t_max, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

/**
 * General function to calculate real roots of a
 * cubic polynom of the form:
 *
 *  f(t) = a*t^3 + b*t^2 + c*t + d
 *
 * resulting in up to 3 t_i satisfying f(t_i) == 0.
 *
 * All results for t are written to t[] and number of
 * found roots is given in the return value of the function.
 *
 * @param t output parameter receiving the values for t where f(t) == 0
 * @return number of roots
 *
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_roots(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T roots[3], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

/**
 * Determine 1. derivative of F and store it in f .
 */
template <typename REAL_T = double>
static inline void polynom_N_derivative(const REAL_T F[0], int N, /* out */ REAL_T f[0]);

template <typename REAL_T = double>
static inline REAL_T polynom_N_value(const REAL_T f[0], int n, REAL_T t);

/**
 * Find up to N roots in a polynom F of degree N.
 *
 * This is the API level entry function to all root finding
 * functions, which tests inputs and then delegates to lower
 * level functions to solve the actual task.
 */
template <int ORDER, typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const REAL_T F[ORDER], REAL_T roots[ORDER], REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

/**
 * Find up to N roots in a polynom F of degree n,
 * having the first derivative function f of F.
 *
 * Preconditions:
 * 		F[0] != 0
 *
 * This is the main control function for root finding.
 * This function analysis the shape of F to make good
 * guesses before delegating to another function to
 * find a root for each guess.
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const REAL_T F[0], const REAL_T f[0], int n, REAL_T roots[0], REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

/**
 * Function determines all t_i for roots, extrema and inflection points of given function F(t).
 *
 * @param F polynom F(t) of degree n (F(t)=F[0]*t^n + F[1]*t^n-1 + ... + F[n]
 * @param f first derivative of F ( f(t) = F'(t) )
 * @param n degree of F
 * @param roots (out) t_i for F(t_i) = 0. Must have space for n values.
 * @param extrema (out) t_i for f(t_i) = 0. Must have space for n-1 values.
 * @param inflections (out) t_i for F''(t) = f'(t_i) = 0. Must have space for n-2 values.
 * @param sizes (out) number of roots, extrema and inflections found. must have space for 3 values.
 * @param tolerance maximum deviation of any t_i or F(t_i)
 * @param domain valid co-domain for any t_i.
 *
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline void polynom_N_features(const REAL_T F[0], const REAL_T f[0], int degree,
		REAL_T roots[0], REAL_T extrema[0], REAL_T inflections[0], int sizes[3],
		REAL_T tolerance, CO_DOMAIN_T domain);




/**
 * Find at most one root, based on the given guess
 * inside a given interval.
 *
 *
 * This function will only find the root in the given
 * interval, if the following conditions are met:
 *
 * - The interval must not contain an inflection point
 *   or extrema.
 * - The interval must contain a root (do boundary
 *   checks first).
 *
 * It is recommended to split the co-domain into intervals
 * between extrema, inflection points and given boundaries.
 * If this condition is met, then this function will always
 * converge on the root inside given interval.
 *
 * The algorithm will abort after MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS
 * iterations to avoid infinite loops. Return value in this case is NAN.
 * See also math2d-config.h.
 *
 * @param F function to be searched (polynom of degree n)
 * @param f derivative of F
 * @param n degree of F
 * @param guess guess for t_0 of a possible location of a root for F(t_0) = 0. guess must satisfy domain(guess).
 * @param tolerance deviation tolerance in respect to t_0 and F(t_0) for found root.
 * @param domain valid interval to be searched for the root
 * @return t for F(t) == 0 or NAN in case no root was found in MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS iterations.
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline REAL_T polynom_N_root(const REAL_T F[0], const REAL_T f[0], int n, REAL_T guess, const REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

////////////////////////////////////////////////////////////////////////////////////
//          I M P L E M E N T A T I O N S     B E L O W
////////////////////////////////////////////////////////////////////////////////////


// Stack local allocation is generally faster than malloc, but
// if there is something fishy going on, we can switch to
// malloc instead of alloca, to debug it.
#if 0
#define ALLOCA(__SIZE__) malloc(__SIZE__)
#define ALLOCA_INIT(__SYMBOL__) (__SYMBOL__ = NULL)
#define ALLOCA_FREE(__SYMBOL__) ((__SYMBOL__ == NULL) ? free(__SYMBOL__) : (void(0)))
#else
#define ALLOCA(__SIZE__) alloca(__SIZE__)
#define ALLOCA_INIT(__SYMBOL__) (void(0))
#define ALLOCA_FREE(__SYMBOL__) (void(0))
#endif


#ifdef MATH2D_EVALUATE
extern int POLYNOM_N_ROOTS_iterations;
extern int POLYNOM_N_ROOT_iterations_avg;
extern int POLYNOM_N_ROOT_calls_num;
#	define POLYNOM_N_ROOTS_COUNT()             (POLYNOM_N_ROOTS_iterations++)
#	define POLYNOM_N_ROOTS_EVALUATION_RESET()  {POLYNOM_N_ROOTS_iterations = 0;}
#	define POLYNOM_N_ROOT_EVALUATE_AVG_ITERATIONS(__COUNT__) {POLYNOM_N_ROOT_iterations_avg += __COUNT__; POLYNOM_N_ROOT_calls_num++;}
#	define POLYNOM_N_ROOTS_EVALUATION_REPORT() printf("polynom_N_root iterations eval result:\n" \
		"\ttot: %d\n"\
		"\tavg: %d\n"\
		"\tcnt: %d\n",\
	POLYNOM_N_ROOTS_iterations, POLYNOM_N_ROOT_calls_num?POLYNOM_N_ROOT_iterations_avg/POLYNOM_N_ROOT_calls_num:0, POLYNOM_N_ROOT_calls_num)
#else
#	define POLYNOM_N_ROOTS_COUNT()
#	define POLYNOM_N_ROOTS_EVALUATION_RESET()
#	define POLYNOM_N_ROOTS_EVALUATION_REPORT()
#	define POLYNOM_N_ROOT_EVALUATE_AVG_ITERATIONS(__COUNT__)
#endif


// TODO: decide whether we should generally allow double precision tolerance in internal functions
/**
 * WARNING: intervals cannot contain extrema!
 * checks if there is a root between interval boundaries (boundaries excluded).
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline REAL_T __internal__polynom_N_interval_contains_root(const REAL_T F[0], int n, REAL_T tolerance, CO_DOMAIN_T two_extrema_interval);

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots(const REAL_T F[0], int n, REAL_T roots[0], REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots_highp(const REAL_T F[0], int degree, REAL_T roots[0], REAL_T tolerance, CO_DOMAIN_T domain  = CO_DOMAIN_REAL);

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots_highp(const REAL_T F[0], const REAL_T f[0], int degree, REAL_T roots[0], const REAL_T tolerance, CO_DOMAIN_T domain);

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline REAL_T __internal_polynom_N_root_fine_tune(const REAL_T F[0], const REAL_T f[0], int degree, REAL_T guess, REAL_T& tolerance, int& iterations, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template <typename REAL_T = double>
static inline int __internal__sorted_unique(REAL_T list[0], int size, REAL_T tolerance) {
	std::sort(list, list + size);
	REAL_T* end = list+size;
	REAL_T* prev = list;
	for (REAL_T* item = prev+1; item != end; item++) {
		if (about_equal(*prev, *item, tolerance)) {
			*prev = (*prev + *item)/2;
			size--;
		} else {
			prev++;
			*prev = *item;
		}
	}
	return size;
}


/**
 * Defines a polynom f(t) of degree DEGREE.
 * Function f(t) := p[0] * t^N + p[1] * t^(N-1) + ... + p[N]
 *
 * Required: DEGREE > 0
 */
template  <int ORDER, typename REAL_T>
struct polynom_t {
	REAL_T p[ORDER+1];
	polynom_t() {assert(ORDER > 0);}
	polynom_t(const REAL_T _p[ORDER+1])
	{
		assert(ORDER > 0);
		memcpy(p, _p, sizeof(REAL_T)*(ORDER+1));
	}

	int degree() {
		return ORDER;
	}

	/** returns value of the polynom at t */
	REAL_T operator () (REAL_T t) const {
		return polynom_N_value(p, ORDER, t);
	}
	/**
	 * Find one single root close to 'guess' using the Newton-Raphson method.
	 *
	 * Only use this method, if you are looking for just one single root.
	 * Don't call this method repeatedly to find multiple roots.
	 *
	 * Either returns t of a root or +-INFINITY or NAN to indicate, that
	 * there is no root to be found. Use ::finite(t) to check whether it's
	 * a valid result.
	 */
	template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
	REAL_T single_root(REAL_T guess, REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) const {
		// TODO: remove method single_root
		polynom_t<ORDER-1> dfdt;
		derivative(dfdt);
		return polynom_N_root(p, dfdt.p, ORDER, guess, tolerance, domain);
	}

	template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
	int roots(REAL_T roots[ORDER], REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) const {
		return polynom_N_roots<ORDER>(p, roots, tolerance, domain);
	}

	void derivative(polynom_t<ORDER-1>& f) const {
		polynom_N_derivative(p, ORDER, f.p);
	}

};


template<typename REAL_T>
static inline REAL_T polynom1_value(REAL_T a, REAL_T b, REAL_T t) {
	return a*t + b;
}

/**
 * Formally, a polynomial of N-th degree has N distinguishable
 * roots. According to that, this function will always return 0
 * even if it is not the truth.
 *
 * A polynomial of degree 0 (i.e. f(x): x -> c) is a special case
 * in regards of root finding.
 *
 * A polynomial of degree 0 is actually a constant. Formally,
 * in this case, there should be N = 0 roots. But practically,
 * that is not true for the case c == 0, because every x will
 * yield f(x) = 0. A constant is either c = 0 and every x is a
 * root of f(x) or c != 0 and no root exists.
 *
 * Thus, special care has to be taken for this case.
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom0_roots(REAL_T a) {
	return 0;
}



template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom1_roots(REAL_T a, REAL_T b, REAL_T t[0], CO_DOMAIN_T domain) {
	int count = 0;
	// 0 = a*t + b
	// t = -b/a

	if (a == 0) {
		return polynom0_roots(b);
	} else {
		REAL_T t_i = -b/a;
		if (domain(t_i)) t[count++] = t_i;
	}
	return count;
}

template <typename REAL_T = double>
static inline REAL_T polynom2_value(REAL_T a, REAL_T b, REAL_T c, REAL_T t) {
	return t * polynom1_value(a, b, t) + c;
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_roots(REAL_T a, REAL_T b, REAL_T c, REAL_T roots[2], CO_DOMAIN_T domain)
{
	int count = 0;

	if (a == 0) {
		count = polynom1_roots(b, c, roots, domain);
	} else if (c == 0) {
		// one root at 0.0 + remaining in g(x) = f(x)/(x-0)
		if (domain(0.0)) roots[count++] = 0.0;
		count += polynom1_roots(a, b, roots + count, domain);
	} else {

		// 0 = at^2 + bt + c
		// t_1,2 = -b/(2a) +- sqrt ( (b/(2a))^2 - c/a )

		REAL_T t_i;

		REAL_T b_a_2 = b/(a*2.0f);
		REAL_T inner = b_a_2*b_a_2 - c/a;
		if (inner >= 0) {
			REAL_T sqrt_inner = std::sqrt(inner);
			t_i = -b_a_2 + sqrt_inner;
			if (domain(t_i)) roots[count++] = t_i;
			// prevent double root
			if (sqrt_inner)	{
				t_i = -b_a_2 - sqrt_inner;
				if (domain(t_i)) roots[count++] = t_i;
			}
		}
	}
	count = __internal__sorted_unique(roots, count, REAL_T(double_smallest_increment(1.0)));
	return count;
}

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_extrema(REAL_T a, REAL_T b, REAL_T& t_min, REAL_T& t_max, CO_DOMAIN_T domain) {
	REAL_T roots[1];
	int n;
	int result = 0;
	t_min = INFINITY;
	t_max = INFINITY;


	if (a != 0) { 		// polynome of degree 1 has no extrema.
		// roots of 1. derivate
		// f(t) = at^2 + bt + c
		// f'(t) = 2at + b
		// f''(t) = 2a

		n = polynom1_roots(2.f*a,b, roots, domain);

		for (int i = 0; i < n; i++) {

			// check value of 2. derivate
			REAL_T f_t = 2.f*a;
			if (f_t < 0) {
				// maxima
				t_max = roots[i];
				result = 1;
			} else /* (f_t > 0) */ {
				// minima
				t_min = roots[i];
				result = 1;
			}
		}
	}
	return result;
}

template <typename REAL_T = double>
static inline REAL_T polynom3_value(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T t) {
	return t * polynom2_value(a, b, c, t) + d;
}



template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_roots(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T roots[3], CO_DOMAIN_T domain)
{
	static const REAL_T PI =  math2d::CONSTANT_PI;
	int count = 0;

	if (a == 0) {
		return polynom2_roots(b,c,d, roots, domain);
	} else if (d == 0) {
		// one root at 0.0 + remaining in g(x) = f(x)/(x-0)
		if (domain(0.0)) roots[count++] = 0.0;
		count += polynom2_roots(a,b,c, roots + count, domain);
	} else {
		// a != 0
		REAL_T A=b/a;
		REAL_T B=c/a;
		REAL_T C=d/a;

		REAL_T p = B - A*A/REAL_T(3.0);
		REAL_T q = REAL_T(2.0)*A*A*A/27.0 - A*B/3.0 + C;

		REAL_T D = q*q/4.0 + p*p*p/27.0; // discriminant

		REAL_T t_i;

		// temp
		REAL_T A_3 = A/3.0;


		if (D > 0) {
			REAL_T sqrt_D = std::sqrt(D);
			REAL_T neg_q_2 = -q/2.0;
			REAL_T u = std::cbrt(neg_q_2 + sqrt_D);
			REAL_T v = std::cbrt(neg_q_2 - sqrt_D);

			t_i = u + v - A_3;
			if (domain(t_i)) roots[count++] = t_i;
			/* 2 complex roots
				t_i = -(u+v)/2.0 - A_3 + Imaginary((u-v)/2.0 * sqrt(3.0));
				t_i = -(u+v)/2.0 - A_3 - Imaginary((u-v)/2.0 * sqrt(3.0));
			*/
		} else if (D == 0) {
			if (p == 0) {
				t_i = -A_3;  // triple
				if (domain(t_i)) roots[count++] = t_i;
			} else {
				REAL_T q_3_p = 3.0*q/p;
				t_i =   q_3_p       - A_3;
				if (domain(t_i)) roots[count++] = t_i;

				t_i = - q_3_p/(2.0) - A_3; // double
				if (domain(t_i)) roots[count++] = t_i;
			}
		} else /* (D < 0) */ {
			REAL_T sqrt_4_3_p = std::sqrt(REAL_T(-4.0)/3.0*p);
			REAL_T onethird_arccos = REAL_T(1.0)/3.0 * std::acos(-q/2.0 * std::sqrt(REAL_T(-27.0)/(p*p*p)));
			REAL_T PI_3 = PI/3.0;
			t_i = - sqrt_4_3_p * std::cos(onethird_arccos + PI_3) - A_3;
			if (domain(t_i)) roots[count++] = t_i;

			t_i =   sqrt_4_3_p * std::cos(onethird_arccos       ) - A_3;
			if (domain(t_i)) roots[count++] = t_i;

			t_i = - sqrt_4_3_p * std::cos(onethird_arccos - PI_3) - A_3;
			if (domain(t_i)) roots[count++] = t_i;
		}
	}

	count = __internal__sorted_unique(roots, count, REAL_T(double_smallest_increment(1.0)));
    return count;
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_extrema(REAL_T a, REAL_T b, REAL_T c, bool accept_saddle_point, REAL_T& t_min, REAL_T& t_max, CO_DOMAIN_T domain) {
	REAL_T roots[2];
	int n;
	int result = 0;
	t_min = INFINITY;
	t_max = INFINITY;


	if (a == 0) {
		return polynom2_extrema(b,c, t_min, t_max, domain);
	} else {

		// f(t)=a t^3 + b t^2 + c t + d
		// f'(t)= 3at^2 + 2bt + c
		// f''(t)= 6at + 2b

		REAL_T A = 3.f * a;
		REAL_T B = 2.f * b;
		n = polynom2_roots( A, B, c, roots, domain);

		for (int i = 0; i < n; i++) {

			// check value of 2. derivate
			REAL_T f_t_i = polynom1_value(2.f*A, B, roots[i]);
			if (f_t_i < 0) {
				// maxima
				t_max = roots[i];
				result = 1;
			} else if (f_t_i > 0) {
				// minima
				t_min = roots[i];
				result = 1;
			} else if (accept_saddle_point) {
				// saddle point
				t_min = t_max = roots[i];
				result = -1;
			}
		}
		return result;
	}
}

template <typename REAL_T = double>
static inline REAL_T polynom4_value(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T e, REAL_T t) {
	return t * polynom3_value(a, b, c, d, t) + e;
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom4_roots_old_and_buggy(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T e, REAL_T roots[4], CO_DOMAIN_T domain = CO_DOMAIN_REAL)
{

	int count = 0;
	if (a == 0.0) {
		// actually 3rd order polynom
		return polynom3_roots(b,c,d,e, roots, domain);
	} else if (b == 0.0 && d == 0.0) {
		// substitute: x = t^2
		// to: a*x^2 + c*x + e = a*t^4 + c*t^2 + e
		count = polynom2_roots(a, c, e, roots, domain);

		// substitute back: t = sqrt(x)
		for (int i = 0; i < count; i++) {
			roots[i] = sqrt(roots[i]);
		}
		return count;
	} else if (e == 0.0) {
		// f(0) is a root and  g(x) := f(x)/(x-0) contains the remaining roots
		roots[count++] = 0;
		count += polynom3_roots(a,b,c,d, roots + count, domain);
	} else {

		// This implementation is based on the ferrari method
		// as shown in https://en.wikipedia.org/wiki/Quartic_function

		REAL_T DELTA_0 = c*c - 3.0*b*d + 12.0*a*e;
		REAL_T DELTA_1 = 2.0*c*c*c -9.0*b*c*d + 27.0*b*b*e + 27.0*a*d*d - 72.0*a*c*e;
		REAL_T DELTA   = - 1.0/27.0 * (DELTA_1*DELTA_1 - 4.0*DELTA_0*DELTA_0*DELTA_0);

		// FIXME: check if DELTA has correct sign by comparing it to its original form

		REAL_T P = 8.0*a*c - 3.0*b*b;
		REAL_T D = 64.0*a*a*a*e - 16.0*a*a*c*c + 16.0*a*b*b*c - 16.0*a*a*b*d - 3.0*b*b*b*b;

		REAL_T p = P/(8.0*a*a);                                 // always finite since a != 0
		REAL_T q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d)/(8.0*a*a*a); // always finite since a != 0

		REAL_T S;

		if (DELTA > 0.0) {
			// Q is complex if DELTA > 0 and all roots are either complex or real
			// if all roots are real, then there is a way to circumvent complex radicals.
			// if either P > 0 or D > 0, then all roots are complex
			if (P > 0 || D > 0) {
				// complex roots only
				return 0;
			} else {
				// real roots only
				// use trigonometric functions to circumvent complex radicals
				REAL_T PHI;
				if (DELTA_0 != 0) {
					PHI = acos(DELTA_1/(2.0*sqrt(DELTA_0*DELTA_0*DELTA_0)));
					// this looks like a cubic root, because of cos(PHI/3)
					S = 1.0/2.0 * sqrt(-2.0/3.0*p + 2.0/(3.0*a) * sqrt(DELTA_0) * cos(PHI/3.0));
				} else {
					// sqrt(DELTA_0) := 0 -> PHI not necessary
					S = 1.0/2.0 * sqrt(-2.0/3.0*p);
				}
				// since we are supposed to find only real roots in this branch,
				// S always has to be a real number
			}
		} else {
			// Q and S cannot be a complex number
			//
			// If S was a complex number, all roots would be complex numbers too
			// because the second root of any complex number is a complex number too.
			//
			// Thus, S cannot be a complex number here, because not all roots
			// cannot be complex in this branch, due to DELTA <= 0.

			REAL_T Q;
			REAL_T sqrt_inner;

			if (DELTA != 0 && DELTA_0 == 0) {
				// --> (DELTA_1 != 0)

				// making sure that Q != 0
				// DELTA_1 :=  +sqrt(DELTA_1*DELTA_1)
				// thus:
				//    0.5*(DELTA_1 + sqrt(DELTA_1*DELTA_1)) := DELTA_1
				// thus:
				Q = cbrt(DELTA_1);
				// equation is reduced, considering (DELTA_0 == 0) and (Q != 0)
				sqrt_inner = -2.0/3.0*p + 1.0/(3.0*a) * Q;
				// should not fire (see comment at parent branch entry)
				assert(sqrt_inner != 0);
				S = 1.0/2.0 * sqrt(sqrt_inner);
				// S may be zero, which is fine here
			} else {
				// preconditions according to branches above:
				// a != 0
				// ( b != 0  ||  d != 0 )  // either one but not both zero

				for (int i = 0; i < 2; i++) {
					REAL_T sign = 1 - (i%2)*2;

					// Possible cases:
					// 1) DELTA != 0 && DELTA_0 != 0 --> DELTA_1 ?    --> Q != 0
					// 2) DELTA == 0 && DELTA_0 != 0 --> DELTA_1 != 0 --> Q != 0
					// 3) DELTA == 0 && DELTA_0 == 0 --> DELTA_1 == 0 --> Q == 0
					// --> Q != 0 in all cases
					// Generally: if and only if DELTA_0 == 0, then Q == 0 possible
					Q = cbrt((DELTA_1 + sign * sqrt(-27.0*DELTA)) / 2.0);
					if (DELTA_0 == 0) {
						// since DELTA_0 == 0 we can omit the last quotient
						sqrt_inner = -2.0/3.0*p + 1.0/(3.0*a) * (Q);
					} else {
						sqrt_inner = -2.0/3.0*p + 1.0/(3.0*a) * (Q+DELTA_0/Q);
					}
					// should not fire (see comment at parent branch entry)
					assert (sqrt_inner >= 0);
					S = 1.0/2.0 * sqrt(sqrt_inner);
					// If S equals 0, we check another result to be safe.
					if (S != 0) break;
				}
			}
		}

		REAL_T t_i;
		REAL_T sqrt_inner;
		REAL_T sqrt_inner_last_term;
		REAL_T sqrt_inner_first_terms = -4.0*S*S - 2.0*p;

		if (S != 0)
		{
			sqrt_inner_last_term = q/S;
		}
		else
		{
			// FIXME: proof that q == 0 if S == 0
			assert(q == 0);
			sqrt_inner_last_term = 0;
		}
		sqrt_inner = sqrt_inner_first_terms + sqrt_inner_last_term;
		if (sqrt_inner >= 0) {
			// 2 real roots
			t_i = - b/(4.0*a) - S + 1.0/2.0 * sqrt(sqrt_inner);
			if (domain(t_i))roots[count++] = t_i;
			t_i = - b/(4.0*a) - S - 1.0/2.0 * sqrt(sqrt_inner);
			if (domain(t_i))roots[count++] = t_i;
		}

		sqrt_inner = sqrt_inner_first_terms - sqrt_inner_last_term;
		if (sqrt_inner >= 0) {
			// 2 real roots
			t_i = - b/(4.0*a) + S + 1.0/2.0 * sqrt(sqrt_inner);
			if (domain(t_i))roots[count++] = t_i;
			t_i = - b/(4.0*a) + S - 1.0/2.0 * sqrt(sqrt_inner);
			if (domain(t_i))roots[count++] = t_i;
		}

	}
	std::sort(roots, roots + count);

	return count;
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom4_roots(REAL_T a, REAL_T b, REAL_T c, REAL_T d, REAL_T e, REAL_T roots[4], CO_DOMAIN_T domain = CO_DOMAIN_REAL)
{
	// This version is based on the ferrari method.
	// It differs from other implementations in that it
	// delegates to other functions, rather than making
	// the solution and decision making unreadable.


	int count = 0;
	if (a == 0.0) {
		// actually 3rd order polynom
		return polynom3_roots(b,c,d,e, roots, domain);
	} else if (b == 0.0 && d == 0.0) {
		// substitute: x = t^2
		// to: a*x^2 + c*x + e = a*t^4 + c*t^2 + e
		REAL_T x_i[2];
		int num_x = polynom2_roots(a, c, e, x_i, domain);

		// substitute back: t = sqrt(x)
		double x;
		for (int i = 0; i < num_x; i++) {
			if (x_i[i] == 0) {
				x = 0.0;
				if (domain(x)) roots[count++] = x;
			}
			else if (x_i[i] > 0)
			{
				x = sqrt(x_i[i]);
				if (domain(+x)) roots[count++] = +x;
				if (domain(-x)) roots[count++] = -x;
			}
		}
	} else if (e == 0.0) {
		// f(0) is a root and  g(x) := f(x)/(x-0) contains the remaining roots
		if (domain(0.0)) roots[count++] = 0.0;
		count += polynom3_roots(a,b,c,d, roots + count, domain);
	} else {
		typedef long double real_t;
		// Substitute according to Tschirnhaus
		real_t alpha = - real_t(3.0)*b*b/(real_t(8.0)*a*a)   + real_t(c)/a;
		real_t beta  = + real_t(b)*b*b/(real_t(8.0)*a*a*a) - (b*c)/(real_t(2.0)*a*a) + real_t(d)/a;
		real_t gamma = - real_t(3.0)*b*b*b*b/(real_t(256.0)*a*a*a*a) + real_t(b)*b*c/(real_t(16.0)*a*a*a) - real_t(b)*d/(real_t(4.0)*a*a) + real_t(e)/a;
		// resolves in:
		//     0 = u^4 + alpha u^2 + beta u + gamma
		// with
		//     x_i = u_i - b/(4 a)
		real_t x_u_offset = - b/(real_t(4.0) * a);



		if (beta == 0.0) {
			// special case:
			//     0 = u^4 + alpha u^2 + gamma
			// substitute:
			//     z = u^2
			// and solve
			//     0 = z^2 + alpha z + gamma
			// with
			//     x_i = +/- sqrt(z_1,2) - b/(4 a)
			real_t z_i[2];
			// we do not care about complex solutions in z,
			// because sqrt(z) gives another complex number
			// and we are not interested in complex roots in x.
			int n = polynom2_roots(real_t(1.0), alpha, gamma, z_i, domain);
			real_t x;
			for (int i = 0; i < n; i++) {
				if (z_i[i] >= 0) {
					x = std::sqrt(z_i[i]) + x_u_offset;
					if (domain(x)) roots[count++] = x;

					x = -std::sqrt(z_i[i]) + x_u_offset;
					if (domain(x)) roots[count++] = x;
				}
			}
		} else {
			// beta != 0

			// Common case
			// Solving with Ferrari's method

			// substitutions result in:
			//     0 = y^3 + p1 y^2 + p2 y + p3

			real_t y_i[3];
			real_t p0 = 1;
			real_t p1 = real_t(5.0)/2.0*alpha;
			real_t p2 = (real_t(2.0)*alpha*alpha - gamma);
			real_t p3 = (real_t(4.0)*alpha*(alpha*alpha - gamma)-beta*beta)/8.0;
			int n = polynom3_roots(p0, p1, p2, p3, y_i, domain);

			// NOTE: A third degree polynomial always has at least one real root (not complex).
			//       To determine the roots in x, we can use any root in y
			//       and will always get the same results.
			//       Since there is at least one real y, we will use that.

			// this will fire only, if there is an error in polynom_3_root(1, ..) (see reason above)
			assert(n);
			real_t y = y_i[0];
			// now determine the u_1,2,3,4
			real_t sqrt_inner_1 = alpha + real_t(2.0)*y;
			if (sqrt_inner_1 >= 0) {

				real_t w = std::sqrt(sqrt_inner_1);

				// According to en wiki on ferrari's solution:
				// w is never gonna be zero, because if (w == 0)
				// then also (beta == 0) and we would have been
				// in the branch above (see beta == 0).
				// But keep that assert here!
				assert(w != 0);
				real_t z = beta/(real_t(2.0)*w);

				for (int i = 0; i < 4; i++) {
					real_t sign_1 = +1 - 2 * (i/2);
					real_t sign_2 = +1 - 2 * (i%2);

					real_t sqrt_inner_2 = w*w - real_t(4.0) * (alpha + y + sign_1 * -z);
					if (sqrt_inner_2 < 0) continue; // -> 2 complex roots

					// real root
					real_t u_i = real_t(0.5)*(sign_1 * -w + sign_2 * std::sqrt(sqrt_inner_2));
					real_t x_i = u_i + x_u_offset;
					if (domain(x_i)) roots[count++] = x_i;
				}
			} else {
				// the first term is a complex number, which will result in real roots, if
				// (and only if) the second term is a complex number too.
				// Thus, we have to deal with complex numbers here.
				typedef std::complex<real_t> complex_t;

				// we need a little bit of tolerance, because the
				// imaginary value can be slightly off due to
				// trigonometric functions, resulting in seemingly
				// complex roots.
				int tolerance_exp = int(-real_t(sizeof(real_t))*8 * 4.0/5.0);
				const real_t tolerance = std::pow(real_t(2.0),tolerance_exp);

				complex_t w = sqrt(complex_t(sqrt_inner_1));

				// Keep that assert here! See reasoning in branch above.
				assert(w != complex_t(0.0));
				complex_t z = beta/(w*real_t(2.0));

				for (int i = 0; i < 4; i++) {
					real_t sign_1 = +1 - 2 * (i/2);
					real_t sign_2 = +1 - 2 * (i%2);

					complex_t sqrt_inner_2 = w*w - real_t(4.0) * (alpha + y + sign_1 * -z);

					// real root
					complex_t u_i = real_t(0.5)*(sign_1 * -w + sign_2 * sqrt(sqrt_inner_2));
					complex_t x_i = u_i  + x_u_offset;
					// accept all complex with (|Im - 0| < 2^-16) as real numbers
					if (domain(x_i.real()) && about_equal(real_t(x_i.imag()), real_t(0), real_t(tolerance))) roots[count++] = x_i.real();
				}

			}
		}


	}

	count = __internal__sorted_unique(roots, count, double_smallest_increment(1.0));
	return count;
}

template <typename REAL_T = double>
static inline void polynom_N_derivative(const REAL_T F[0], int N, REAL_T f[0]) {
	for (int i = 0; i < N; i++) {
		f[i] = (N-i) * F[i];
	}
}


template <typename REAL_T = double>
static inline REAL_T polynom_N_value(const REAL_T f[0], int n, REAL_T t)
{
	long double result = f[0];
	for (int i = 1; i < n+1; i++) {
		result *= t;
		result += f[i];
	}
	return result;
}

/**
 * Determine maximum interval in which roots of given function can appear.
 *
 * Note: For real roots only.
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline bool polynom_N_root_bounds(const REAL_T F[0], int n, REAL_T& t_lower, REAL_T& t_upper, REAL_T tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) {


	// We initialise bounds with the domain bounds
	// and then try to improve them, using known methods.
	t_lower = domain.lower;
	t_upper = domain.upper;

	if (n > 1) {
		if (F[0] == 0) {
			// function requires (F[0] != 0)
			return polynom_N_root_bounds(F+1, n-1, t_lower, t_upper, tolerance, domain);
		}

		// TODO: optimise: remove ineffective methods and merge remaining methods.

		// Currently we are using all well known methods and
		// take the best bounds. Each method has its own block
		// to keep them distinguishable.

		{
			// Cauchy's Bound
			REAL_T B = ::fabs(F[1]/F[0]);
			for (int i = 2; i <= n; i++) {
				B = std::max(B, ::fabs(F[i]/F[0]));
			}
			B += 1;

			assert(finite(B));
			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);

		}

		{
			// Lagrange's Bound
			REAL_T B = 0;
			for (int i = 1; i <= n; i++) {
				B += ::fabs(F[i]/F[0]);
			}
			B = std::max(1.0, B);

			assert(finite(B));
			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);
		}

		{
			// Zassenhausen's Bound
			REAL_T B = ::fabs(F[1]/F[0]);
			for (int i = 2; i <= n; i++) {
				B = std::max(B, ::pow(::fabs(F[i]/F[0]), 1.0/i));
			}
			B = 2 * B;

			assert(finite(B));

			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);

		}

		{
			// Zassenhausen's Bound improved by Lagrange
			REAL_T first_largest  = -INFINITY;
			REAL_T second_largest = -INFINITY;
			REAL_T B              = -INFINITY;
			for (int i = 1; i <= n; i++) {
				B = std::max(B, ::pow(::fabs(F[i]/F[0]), 1.0/i));

				if (B > first_largest) {
					second_largest = first_largest;
					first_largest = B;
				} else if (B > second_largest) {
					second_largest = B;
				}
			}
			B = first_largest + second_largest;

			assert(finite(B));

			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);
		}

		{
			// Other methods:

			// "Samuelson's inequality"
			// see https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Polynomials_whose_roots_are_all_real
			// required: polynomial without any complex roots (no appropriate test available)

		}

	}

	bool improved = false;
	if (t_lower > domain.lower || t_upper < domain.upper) {
		improved = true;

		// in some cases the interval is really narrow
		// so we open it up by a little notch
		t_lower = std::max(REAL_T(domain.lower), t_lower-tolerance);
		t_upper = std::min(REAL_T(domain.upper), t_upper+tolerance);

#ifndef NDEBUG
		// weak test if there is actually no root between root boundary and domain boundary
		bool contains_root;
		if (t_lower > domain.lower) {
			co_domain_const_t off_boundary(REAL_T(domain.lower), t_lower);
			contains_root = __internal__polynom_N_interval_contains_root(F, n, tolerance, off_boundary);
			assert(!contains_root);
		}
		if (t_upper < domain.upper) {
			co_domain_const_t off_boundary(t_upper, REAL_T(domain.upper));
			contains_root = __internal__polynom_N_interval_contains_root(F, n, tolerance, off_boundary);
			assert(!contains_root);
		}
#endif
	}
	return improved;
}

template <int ORDER, typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const REAL_T F[ORDER], REAL_T roots[ORDER], REAL_T tolerance, CO_DOMAIN_T domain) {
	int i;
	for (i = 0; unlikely(i <= ORDER && F[i] == 0); i++);
	int n = ORDER-i;
	if (n <= 0) return 0;
	const REAL_T* W = F + i;
	return __internal__polynom_N_roots(W, n, roots, tolerance, domain);
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const REAL_T F[0], const REAL_T f[0], int n, REAL_T roots[0], REAL_T tolerance, CO_DOMAIN_T domain)
{
	int sizes[3];
	REAL_T extrema[n>0?n-1:0];
	REAL_T inflections[n>1?n-2:0];
	int& num_roots = sizes[0];
	REAL_T lower_bound = domain.lower;
	REAL_T upper_bound = domain.upper;
	polynom_N_root_bounds(F, n, lower_bound, upper_bound, tolerance, domain);
	co_domain_const_t root_domain(lower_bound, upper_bound);
	polynom_N_features(F, f, n, roots, extrema, inflections, sizes, tolerance, root_domain);
	return num_roots;
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots(const REAL_T F[0], int n, REAL_T roots[0], REAL_T tolerance, CO_DOMAIN_T domain) {
	assert(F[0] != 0);
	int count = 0;
	switch(n) {
	case 0:
		count = polynom0_roots(F[0]);
		break;
	case 1:
		count = polynom1_roots(F[0], F[1], roots, domain);
		break;
	case 2:
		count = polynom2_roots(F[0], F[1], F[2], roots, domain);
		break;
#if MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS
	case 3:
		count = polynom3_roots(F[0], F[1], F[2], F[3], roots, domain);
		break;
	case 4:
		count = polynom4_roots(F[0], F[1], F[2], F[3], F[4], roots, domain);
		break;
#endif
	default:
		{
			REAL_T* f;
			ALLOCA_INIT(f);
			f = (REAL_T*)ALLOCA(sizeof(REAL_T)*n);
			polynom_N_derivative(F, n, f);
			int count = polynom_N_roots(F, f, n, roots, tolerance, domain);
			ALLOCA_FREE(f);
			return count;
		}
	}

	return count;
}

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots_highp(const REAL_T F[0], const REAL_T f[0], int degree, REAL_T roots[0], const REAL_T tolerance, CO_DOMAIN_T domain)
{
	int count = 0;

	if (degree <= 4) {
		count = __internal__polynom_N_roots(F, degree, roots, tolerance, domain);
		if (3 <= degree) {
			for (int i = 0; i < count; i++) {
				REAL_T achieved_tolerance = tolerance;
				int iterations = MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS;
				REAL_T x = __internal_polynom_N_root_fine_tune(F, f, degree, roots[i], achieved_tolerance, iterations, domain);
				assert(finite(x));
				roots[i] = x;
			}
		}
	} else {
		count = polynom_N_roots(F, f, degree, roots, tolerance, domain);
	}
	return count;
}


template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots_highp(const REAL_T F[0], int degree, REAL_T roots[0], REAL_T tolerance, CO_DOMAIN_T domain)
{
	REAL_T* f;
	ALLOCA_INIT(f);
	f = (REAL_T*)ALLOCA(sizeof(REAL_T)*degree);
	polynom_N_derivative(F, degree, f);

	int count = __internal__polynom_N_roots_highp(F, f, degree, roots, tolerance, domain);
	ALLOCA_FREE(f);
	return count;
}


/**
 * This function checks if there is a root between
 * interval boundaries (boundaries excluded).
 *
 * WARNING: intervals cannot contain true extrema!
 */
template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline REAL_T __internal__polynom_N_interval_contains_root(const REAL_T F[0], int n, REAL_T tolerance, CO_DOMAIN_T two_extrema_interval) {
	// Interval contains root if F(t) for lower and upper
	// boundary have different sign and none of it equals zero.
	// The latter case means, that there is a root on a boundary
	// which includes, that there cannot be a root inside the interval.

	// F(t) on lower boundary
	REAL_T v_min = polynom_N_value(F, n, REAL_T(two_extrema_interval.lower));
	// F(t) on upper boundary
	REAL_T v_max = polynom_N_value(F, n, REAL_T(two_extrema_interval.upper));

	if (v_min*v_max < 0) {
		return true;
	}
	return false;
}

template <typename REAL_T = double>
static inline REAL_T polynom_N_root_binary_search(const REAL_T F[0], int n, REAL_T guess, const REAL_T tolerance, co_domain_t<REAL_T>& interval, int& iterations)
{
	// t
	REAL_T t = guess;

	// we cannot accept 0 tolerance, because
	// at some point the result cannot be improved,
	// due to floating point inaccuracies in evaluation
	// of F(t).
	assert(tolerance > 0);
	assert(finite(guess));
	assert(interval(guess));
	assert(finite(interval.lower));
	assert(finite(interval.upper));

	// binary search
	do {
		POLYNOM_N_ROOTS_COUNT();
		t = REAL_T(interval.lower)/2 + REAL_T(interval.upper)/2;

		REAL_T F_lower = polynom_N_value(F, n, REAL_T(interval.lower));
		REAL_T F_upper = polynom_N_value(F, n, REAL_T(interval.upper));

		REAL_T F_t     = polynom_N_value(F, n, t);

		if (F_t*F_lower < 0) {
			interval.upper = t;
		} else if (F_t*F_upper < 0){
			interval.lower = t;
		} else {
			interval.lower = (REAL_T(interval.lower) + t)/2;
			interval.upper = (REAL_T(interval.upper) + t)/2;
			return t;
		}

	} while (--iterations && interval.size() > tolerance);
	return REAL_T(interval.lower)/2 + REAL_T(interval.upper)/2;
}

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline REAL_T polynom_N_root(const REAL_T F[0], const REAL_T f[0], int n, REAL_T guess, const REAL_T tolerance, CO_DOMAIN_T domain)
{
	// F(t)
	REAL_T F_t;
	// f(t)
	REAL_T f_t;

	co_domain_t<REAL_T> interval(REAL_T(domain.lower), REAL_T(domain.upper));

	// guess = polynom_N_root_binary_search(F, n, guess, tolerance*32, interval);


	// t
	REAL_T t = guess;

	// TODO: maximum number of iterations as optional parameter
	int iterations = MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS;

	// we cannot accept 0 tolerance, because
	// at some point the result cannot be improved,
	// due to floating point inaccuracies in evaluation
	// of F(t).
	assert(tolerance > 0);
	assert(finite(guess));
	assert(interval(guess));

	REAL_T previous_accuracy = INFINITY;
	// Newton-Raphson method
	do {
		POLYNOM_N_ROOTS_COUNT();
		guess = t;

		F_t = polynom_N_value(F, n, guess);

		f_t = polynom_N_value(f, n-1, guess);
		// this will only be null, if the interval contains an extrema
		// which is forbidden
		// FIXME: this can also happen, if accuracy is too low to determine the slope
		assert (f_t != 0);

		t = guess - F_t / f_t;
		if (t >= interval.upper) {
			t = (interval.upper + guess)/2;
		} else if (t <= interval.lower) {
			t = (guess + interval.lower)/2;
		}
		assert(finite(t) && interval(t));
		REAL_T accuracy = fabs(guess-t);

		if ((accuracy <= tolerance) && about_equal(F_t, 0.0, REAL_T(tolerance))) {
			guess = t;
			// -> will end the loop and proceed with improvement on t
		}
		if (accuracy == previous_accuracy) {
			// does not converge: algorithm bounces between two values of t
			// -> switch to different method
			// TODO: do we give up too early here?
			guess = t = polynom_N_root_binary_search(F, n, guess, tolerance, interval, iterations);
			// -> will end the loop and proceed with improvement on t
		}

		previous_accuracy = accuracy;
		// abort if either we run out of iterations
		// or we can't improve
	} while (--iterations && t != guess);



	double accuracy = tolerance;
	if (t == guess) {
		// found root, going to improve it
		t = __internal_polynom_N_root_fine_tune(F, f, n, guess, accuracy, iterations, interval);
		if (accuracy > tolerance) {
			t = NAN;
		}
	} else {
		// no root found
		t = NAN;
	}

	// aborted to avoid infinite loop or we can't reach requested accuracy
	POLYNOM_N_ROOT_EVALUATE_AVG_ITERATIONS(MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS - iterations);
	return t;
}

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline REAL_T __internal_polynom_N_root_fine_tune(const REAL_T F[0], const REAL_T f[0], int degree, REAL_T guess, REAL_T& tolerance, int& iterations, CO_DOMAIN_T domain) {
	// F(t)
	REAL_T F_t;
	// f(t)
	REAL_T f_t;
	// t
	REAL_T t = guess;

	assert(finite(guess));
	assert(domain(guess));

	// This function improves the given guess for t_0 by determining the narrowest
	// interval [t_lower, t_upper] which satisfies
	//       F(t_lower) == 0 and F(t_upper) == 0
	// meaning, that the values at the boundaries have different signs.
	// This function assumes, that the given guess is very close to
	// one of these boundaries. Therefore it uses a single step hill-climbing
	// method to determine the boundaries. Once the boundaries were found
	// t_0 is computed as the linear interpolation between them:
	//     t_0 = (t_lower + t_upper)/2


	REAL_T t_lower = guess;
	REAL_T t_upper = guess;

	int boundaries_to_be_checked = 2;

	// compute a minimal step-width
	// If the step is too small for the current t so that
	// (t + step == t), then we use the smallest possible
	// increment instead.

	REAL_T min_step = double_smallest_increment(guess);

	if (tolerance == 0) {
		tolerance = min_step;
	}


	REAL_T step = tolerance;

	if (guess + step == guess) {
		// we can't improve beyond given tolerance
		return guess;
	}

	REAL_T F_t_previous = polynom_N_value(F, degree, t - min_step);
	F_t = polynom_N_value(F, degree, t);

	// determine in which direction we have to search.
	f_t = polynom_N_value(f, degree-1, t);
	// f_t = F_t - F_t_previous;

	if (F_t < 0 && f_t < 0) step = -step;
	if (F_t > 0 && f_t > 0) step = -step;

	do {
		POLYNOM_N_ROOTS_COUNT();

		F_t_previous = F_t;
		guess = t;

		// new value for t
		t = guess + step;

		// check guess and boundaries
		if (t == guess) {
			// cannot improve -> step is too small for current t
			break;
		}
		assert(finite(t));

		// check F(t)
		F_t = polynom_N_value(F, degree, t);

		REAL_T F_eval = F_t*F_t_previous;
		if (F_eval < 0) {
			// We have crossed the X-axis in one single step.
			// t_lower and t_upper are guess and t (or vice versa)
			if (step > 0) {
				// Forward direction
				t_lower = guess;
				t_upper = t;
			} else {
				// Backward direction
				t_lower = t;
				t_upper = guess;
			}
			// assuming that we are done
			boundaries_to_be_checked = 0;

			// prepare to resuming searching in the opposite
			// direction with higher accuracy
			step = -step;

			// check the approximate accuracy of F((t+guess)/2)
			if (about_equal((F_t+F_t_previous)/2, 0, tolerance)
					|| (guess+step/2 == guess))
			{
				// Either the accuracy is ok
				// or we cannot improve accuracy of F_t any further.
				// For example, if f(t) very large (almost infinite)
				// it will be very hard to find a t with F(t) = 0.
				// This is actually not too important, because
				// tolerance is in respect to t and not F(t).

				REAL_T achieved_accuracy = fabs(step);
				tolerance = achieved_accuracy;
				t = (guess + t)/2.0;
				if (!domain(t)) {
					t = NAN;
				}
				return t;
			} else {
				// keep searching in opposite direction (is prepared)
				// from current position
				step = step/2.0;
				boundaries_to_be_checked = 1;
			}
		}
		else if (F_eval == 0)
		{
			if (F_t_previous != 0) {
				// --> F_t == 0
				// entering interval from left or right
				// t is the first t with F(t) == 0
				if (step > 0) t_lower = t;
				else t_upper = t;
				boundaries_to_be_checked--;
			} else if (F_t != 0) {
				// F_t_previous == 0
				// leaving interval to the left or right
				// guess is the last t with F(t) == 0
				if (step > 0) {
					t_upper = guess;
					t = t_upper; // prepare to improve lower boundary next
				}
				else
				{
					t_lower = guess;
					t = t_lower; // prepare to improve upper boundary next
				}
				// reverse search direction
				step = -step;
				F_t = F_t_previous;
				boundaries_to_be_checked--;
			}

			if (!boundaries_to_be_checked) {
				// assume t_0 in the middle
				guess = (t_upper + t_lower)/2;
				// check accuracy we have achieved
				REAL_T achieved_accuracy = t_upper-t_lower;
				if (likely(achieved_accuracy <= tolerance)) {
					tolerance = achieved_accuracy;
					if (!domain(guess)) {
						guess = NAN;
					}
					return guess;
				} else {
					// improve opposite boundary with increased accuracy
					boundaries_to_be_checked = 1;
					step = step/2.0;
				}
			}

		}
		else if (F_t*F_t >= F_t_previous*F_t_previous)
		{
			// distance to zero is getting larger
			// --> we have passed the root
			// Turn around and try again with higher accuracy
			if (step > 0) {
				t_upper = t;
			}
			else
			{
				t_lower = t;
			}
			// reverse search direction
			step = -step/2;
			// provide proper value F(t)
			F_t = polynom_N_value(F, degree, t);
		}

	} while (--iterations);

	// aborted to avoid infinite loop or we can't reach requested accuracy
	return NAN;
}

template <typename REAL_T = double, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline void polynom_N_features(const REAL_T F[0], const REAL_T f[0], int degree, REAL_T roots[0], REAL_T extrema[0], REAL_T inflections[0], int sizes[3], REAL_T tolerance, CO_DOMAIN_T domain) {
	//
	// Searching roots using Newton-Raphson method
	// -------------------------------------------
	// Newton-raphson has some caveats to be avoided.
	//     a. Algorithm can run into infinite loop, jumping around one particular root.
	//     b. Algorithm might fail to find certain root and returns another one instead.
	//
	// This method first splits the domain into search intervals, where each interval is
	// guaranteed to contain at most one root, based on known features of the function.
	// Then it checks the boundaries of each interval and uses the Newton_raphson
	// algorithm only, if an interval is guaranteed to contain a root.
	//
	// Search intervals are determined based on the following conclusions:
	//     1. there is at most one root between two extrema.
	//     2. there is at most one root between domain boundary and the extreme next to it inside the domain.
	//     3. there is at most one root between domain boundaries, if no extreme was found
	//
	// For each search interval, the following statements are true:
	//     4. there is no root in a search interval to be found, if one of the boundaries is a root.
	//     5. there is exactly one root to be found in a search interval, if one boundary is below and the other above zero
	//
	// And finally, regarding the newton-raphson method:
	//     6. The search algorithm cannot run into an infinite loop,
	//        if there is no inflection point inside a given search interval
	//
	// So basically, we need to determine the extrema and inflection points of the function,
	// before we can start searching roots of it.
	assert(tolerance >= 0);

	int& num_roots = sizes[0];
	int& num_extrema = sizes[1];
	int& num_inflections = sizes[2];

	double* dfdt;
	ALLOCA_INIT(dfdt);
	if (2 < degree) {
		// F''(t)
		dfdt = (REAL_T*)ALLOCA(sizeof(REAL_T)*(degree-1));
		polynom_N_derivative(f, degree-1, dfdt);
	}

	if (4 < degree)
	{
		// determine extrema and inflections by recursively calling this same function again
		REAL_T* unused = roots; // using 'roots' as temp
		polynom_N_features(f, dfdt, degree-1, extrema, inflections, unused, sizes, tolerance, domain);
		num_inflections = num_extrema;
		num_extrema = num_roots;
	}
	else /* degree <= 4 */
	{

		if (2 < degree) {
			num_inflections = __internal__polynom_N_roots_highp(dfdt, degree-2, inflections, tolerance, domain);
			//assert(!num_inflections || domain(inflections[0]));
		} else {
			num_inflections = 0;
		}
		if (1 < degree) {
			if (2 < degree) num_extrema = __internal__polynom_N_roots_highp(f, dfdt, degree-1, extrema, tolerance, domain);
			else num_extrema = __internal__polynom_N_roots_highp(f, degree-1, extrema, tolerance, domain);
			//assert(!num_extrema || domain(extrema[0]));
		} else {
			num_extrema = 0;
		}
	}
	ALLOCA_FREE(dfdt);


	if (degree <= 4 && MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS) {
		// lucky -> delegate to simpler root function
		num_roots = __internal__polynom_N_roots_highp(F, f, degree, roots, tolerance, domain);
		assert(!num_roots || domain(roots[0]));
	} else {
		// determine roots using newton-raphson method as explained in introduction
		num_roots = 0;


		// TODO: guesses for first and last roots can still be improved
		// since boundaries are actually guesses for the first and last root.
		// But we have to know, whether the boundary is computed or a
		// given boundary of the application. Which means, the boundary
		// computation has to move here, maybe.
		// See polynom_N_root_bounds().



		struct merger {
			REAL_T* e;
			REAL_T* e_end;
			REAL_T* i;
			REAL_T* i_end;

			merger(REAL_T* extrema, int num_extrema, REAL_T* inflections, int num_inflections)
			: e(extrema),
			  e_end(extrema+num_extrema),
			  i(inflections),
			  i_end(inflections+num_inflections)
			{}

			REAL_T next() {
				if (e == e_end) {
					assert(i != i_end);
					return *i++;
				}
				else if (i == i_end){
					assert(e != e_end);
					return *e++;
				}
				else {
					assert(e != e_end && i != i_end);
					return *e < *i ? *e++ : *i++;
				}
			}

			bool hasMore() {
				return (e != e_end) || (i != i_end);
			}

		} features(extrema, num_extrema, inflections, num_inflections);


		// FIXME: we have to fine tune roots on boundaries
		//        to check whether its actually on the boundary
		//        and not inside the interval (extrema and inflections)
		// TODO:  ACTUALLY, THE WHOLE THING SHOULD BE REDESIGNED
		co_domain_std_t interval;

		interval.lower = domain.lower.bound;
		REAL_T F_lower = polynom_N_value(F, degree, interval.lower.bound);
		REAL_T F_upper = -INFINITY;
		if (F_lower == 0) {
			// root on lower boundary
			roots[num_roots++] = interval.lower.bound;
		}

		// consider all features and upper bound
		do {
			interval.upper = features.hasMore() ? features.next() : domain.upper.bound;

			if (interval.lower.bound < interval.upper.bound) {

				F_upper = polynom_N_value(F, degree, interval.upper.bound);

				if(F_lower*F_upper < 0) {
					// interval contains root

					REAL_T guess = (interval.upper + interval.lower)/2;
					// This should not happen because we have estimated upper
					// and lower bounds for roots and the boundaries should be
					// out of computational range.
					assert(finite(guess*guess));
					guess = polynom_N_root(F, f, degree, guess, tolerance, interval);
					assert(finite(guess));
					assert (!num_roots || roots[num_roots-1] != guess);
					// not a duplicate of previous root
					assert(num_roots < degree);
					roots[num_roots++] = guess;
				} else if (F_upper == 0) {
					// should not happen, unless boundaries are equal
					assert (!num_roots || roots[num_roots-1] != interval.upper.bound);
					assert(num_roots < degree);
					roots[num_roots++] = interval.upper;
				} else if (F_lower != 0 && about_equal(F_lower, 0.0, tolerance)) {
					// According to the given tolerance, the lower end is about zero.
					// Since the interval did not contain a root, the uncertainty
					// about an inaccurate duplicate root on the boundary and
					// a found root in the interval, has been ruled out here.
					// This case can only happen, if the lower boundary is an extreme.
					assert(num_roots < degree);
					roots[num_roots++] = interval.lower.bound;
				}
			}
			// prepare next interval
			interval.lower = interval.upper;
			F_lower = F_upper;
		} while (interval.lower.bound < domain.upper.bound); // FIXME: have to check last boundary too, because of application given boundaries
	}
}


}; // namespace math

#endif /* MATH2D_POLYNOM_H_ */
