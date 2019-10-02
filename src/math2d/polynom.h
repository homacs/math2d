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
#include <string.h>
#include <algorithm>
#include <assert.h>



#include "float-utils.h"

#include "math2d/precision.h"
#include "math2d/utils.h"

namespace math2d {

template  <int DEGREE>
struct polynom_t;


static inline double polynom1_value(double a, double b, double t);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom1_roots(double a, double b, double t[1], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

static inline double polynom2_value(double a, double b, double c, double t);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_roots(double a, double b, double c, double roots[2], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_extrema(double a, double b, double& t_min, double& t_max, CO_DOMAIN_T domain = CO_DOMAIN_REAL);


static inline double polynom3_value(double a, double b, double c, double d, double t);

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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_extrema(double a, double b, double c, bool accept_saddle_point, double& t_min, double& t_max, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_roots(double a, double b, double c, double d, double roots[3], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

/**
 * Determine 1. derivative of F and store it in f .
 */
static inline void polynom_N_derivative(const double F[0], int N, /* out */ double f[0]);

static inline double polynom_N_value(const double f[0], int n, double t);

/**
 * Find up to N roots in a polynom F of degree N.
 *
 * This is the API level entry function to all root finding
 * functions, which tests inputs and then delegates to lower
 * level functions to solve the actual task.
 */
template <int DEGREE, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[DEGREE], double roots[DEGREE], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[0], double f[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline void polynom_N_features(const double F[0], const double f[0], int n,
		double roots[0], double extrema[0], double inflections[0], int sizes[3],
		float tolerance, CO_DOMAIN_T domain);




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
 * converge on the root inside given interval and never abort.
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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline double polynom_N_root(const double F[0], const double f[0], int n, double guess, float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

////////////////////////////////////////////////////////////////////////////////////
//          I M P L E M E N T A T I O N S     B E L O W
////////////////////////////////////////////////////////////////////////////////////

#ifdef MATH2D_EVALUATE
	extern int POLYNOM_N_ROOTS_iterations;
#	define POLYNOM_N_ROOTS_COUNT()             (POLYNOM_N_ROOTS_iterations++)
#	define POLYNOM_N_ROOTS_EVALUATION_RESET()  (POLYNOM_N_ROOTS_iterations = 0)
#	define POLYNOM_N_ROOTS_EVALUATION_REPORT() printf("polynom_N_root iterations: %d\n", POLYNOM_N_ROOTS_iterations)
#else
#	define POLYNOM_N_ROOTS_COUNT()
#	define POLYNOM_N_ROOTS_EVALUATION_RESET()
#	define POLYNOM_N_ROOTS_EVALUATION_REPORT()
#endif




template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots(const double F[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);


/**
 * Defines a polynom f(t) of degree DEGREE.
 * Function f(t) := p[0] * t^N + p[1] * t^(N-1) + ... + p[N]
 *
 * Required: DEGREE > 0
 */
template  <int DEGREE>
struct polynom_t {
	double p[DEGREE+1];
	polynom_t() {assert(DEGREE > 0);}
	polynom_t(const double _p[DEGREE+1])
	{
		assert(DEGREE > 0);
		memcpy(p, _p, sizeof(double)*(DEGREE+1));
	}

	int degree() {
		return DEGREE;
	}

	/** returns value of the polynom at t */
	double operator () (double t) const {
		return polynom_N_value(p, DEGREE, t);
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
	double single_root(double guess, float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) const {
		// TODO: remove method single_root
		polynom_t<DEGREE-1> dfdt;
		derivative(dfdt);
		return polynom_N_root(p, dfdt.p, DEGREE, guess, tolerance, domain);
	}

	template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
	int roots(double roots[DEGREE], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) const {
		return polynom_N_roots<DEGREE>(p, roots, tolerance, domain);
	}

	void derivative(polynom_t<DEGREE-1>& f) const {
		polynom_N_derivative(p, DEGREE, f.p);
	}

};



static inline double polynom1_value(double a, double b, double t) {
	return a*t + b;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom1_roots(double a, double b, double t[1], CO_DOMAIN_T domain) {
	int count = 0;
	// 0 = a*t + b
	// t = -b/a

	if (a != 0) {
		double t_i = -b/a;
		if (domain(t_i)) t[count++] = t_i;
	}
	return count;
}


static inline double polynom2_value(double a, double b, double c, double t) {
	return t * polynom1_value(a, b, t) + c;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_roots(double a, double b, double c, double roots[2], CO_DOMAIN_T domain) {
	int count = 0;

	if (a == 0) {
		count = polynom1_roots(b, c, roots, domain);
	} else {

		// 0 = at^2 + bt + c
		// t_1,2 = -b/(2a) +- sqrt ( (b/(2a))^2 - c/a )

		double t_i;

		double b_a_2 = b/(a*2.0f);
		double inner = b_a_2*b_a_2 - c/a;
		if (inner >= 0) {
			double sqrt_inner = sqrt(inner);
			t_i = -b_a_2 + sqrt_inner;
			if (domain(t_i)) roots[count++] = t_i;
			// prevent double root
			if (sqrt_inner)	{
				t_i = -b_a_2 - sqrt_inner;
				if (domain(t_i)) roots[count++] = t_i;
			}
		}
	}
	std::sort(roots, roots + count);
	return count;
}

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom2_extrema(double a, double b, double& t_min, double& t_max, CO_DOMAIN_T domain) {
	double roots[1];
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
			double f_t = 2.f*a;
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


static inline double polynom3_value(double a, double b, double c, double d, double t) {
	return t * polynom2_value(a, b, c, t) + d;
}



template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_roots(double a, double b, double c, double d, double roots[3], CO_DOMAIN_T domain)
{
	static const double PI =  math2d::CONSTANT_PI;
	int count = 0;

	if (a == 0) {
		return polynom2_roots(b,c,d, roots, domain);
	} else {
		// a != 0
		double A=b/a;
		double B=c/a;
		double C=d/a;

		double p = B - A*A/3.0;
		double q = 2.0*A*A*A/27.0 - A*B/3.0 + C;

		double D = q*q/4.0 + p*p*p/27.0; // discriminant

		double t_i;

		// temp
		double A_3 = A/3.0;


		if (D > 0) {
			double sqrt_D = ::sqrt(D);
			double neg_q_2 = -q/2.0;
			double u = ::cbrt(neg_q_2 + sqrt_D);
			double v = ::cbrt(neg_q_2 - sqrt_D);

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
				double q_3_p = 3.0*q/p;
				t_i =   q_3_p       - A_3;
				if (domain(t_i)) roots[count++] = t_i;

				t_i = - q_3_p/(2.0) - A_3; // double
				if (domain(t_i)) roots[count++] = t_i;
			}
		} else /* (D < 0) */ {
			double sqrt_4_3_p = ::sqrt(-4.0/3.0*p);
			double onethird_arccos = 1.0/3.0 * ::acos(-q/2.0 * sqrt(-27.0/(p*p*p)));
			double PI_3 = PI/3.0;
			t_i = - sqrt_4_3_p * ::cos(onethird_arccos + PI_3) - A_3;
			if (domain(t_i)) roots[count++] = t_i;

			t_i =   sqrt_4_3_p * ::cos(onethird_arccos       ) - A_3;
			if (domain(t_i)) roots[count++] = t_i;

			t_i = - sqrt_4_3_p * ::cos(onethird_arccos - PI_3) - A_3;
			if (domain(t_i)) roots[count++] = t_i;
		}
	}

	std::sort(roots, roots + count);

    return count;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom3_extrema(double a, double b, double c, bool accept_saddle_point, double& t_min, double& t_max, CO_DOMAIN_T domain) {
	double roots[2];
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

		double A = 3.f * a;
		double B = 2.f * b;
		n = polynom2_roots( A, B, c, roots, domain);

		for (int i = 0; i < n; i++) {

			// check value of 2. derivate
			double f_t_i = polynom1_value(2.f*A, B, roots[i]);
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



static inline void polynom_N_derivative(const double F[0], int N, double f[0]) {
	for (int i = 0; i < N; i++) {
		f[i] = (N-i) * F[i];
	}
}



static inline double polynom_N_value(const double f[0], int n, double t) {
	long double result = f[0];
	for (int i = 1; i < n+1; i++) {
		result *= t;
		result += f[i];
	}
	return result;
}


template <int DEGREE, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[DEGREE], double roots[DEGREE], float tolerance, CO_DOMAIN_T domain) {
	int i;
	for (i = 0; unlikely(i <= DEGREE && F[i] == 0); i++);
	int n = DEGREE-i;
	const double* W = F + i;
	return __internal__polynom_N_roots(W, n, roots, tolerance, domain);
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[0], double f[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	int sizes[3];
	double extrema[n>0?n-1:0];
	double inflections[n>1?n-2:0];
	int& num_roots = sizes[0];
	polynom_N_features(F, f, n, roots, extrema, inflections, sizes, tolerance, domain);
	return num_roots;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots(const double F[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	assert(F[0] != 0);

	switch(n) {
	case 0:
		return 0;
	case 1:
		return polynom1_roots(F[0], F[1], roots, domain);
	case 2:
		return polynom2_roots(F[0], F[1], F[2], roots, domain);
	case 3:
		return polynom3_roots(F[0], F[1], F[2], F[3], roots, domain);
	default:
		{
			double* f = (double*)alloca(sizeof(double)*n);
			polynom_N_derivative(F, n, f);
			return polynom_N_roots(F, f, n, roots, tolerance, domain);
		}
	}
}


/**
 * WARNING: intervals cannot contain extrema!
 * checks if there is a root between interval boundaries (boundaries excluded).
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline double __internal__polynom_N_interval_contains_root(const double F[0], int n, double tolerance, CO_DOMAIN_T two_extrema_interval) {
	double v_min = polynom_N_value(F, n, two_extrema_interval.re_min);
	double v_max = polynom_N_value(F, n, two_extrema_interval.re_max);

	// if at least one of the boundaries is a root,
	// then there is no more root to be found.
	if (about_equal(v_min, 0.0, tolerance) || about_equal(v_max, 0.0, tolerance)) {
		return false;
	}

	// if both are above zero or both are below zero,
	// then there is no more root to be found
	if (v_min*v_max > 0) {
		return false;
	}
	return true;
}

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline double polynom_N_root(const double F[0], const double f[0], int n, double guess, float tolerance, CO_DOMAIN_T domain) {
	// F(t)
	double F_t;
	// f(t)
	double f_t;
	// t
	double t = guess;

	// maximum number of iterations
	int iterations = MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS;

	assert(tolerance > 0);
	assert(domain(guess));

	while (--iterations) {
		POLYNOM_N_ROOTS_COUNT();
		guess = t;

		F_t = polynom_N_value(F, n, guess);
		if (F_t == 0) return t;

		f_t = polynom_N_value(f, n-1, guess);
		// this will only be null, if the interval contains an extrema
		assert (f_t != 0);

		t = guess - F_t / f_t;
		if (t >= domain.re_max) {
			t = (domain.re_max + guess)/2;
		} else if (t <= domain.re_min) {
			t = (guess + domain.re_min)/2;
		}
		assert(domain(t));


		if (!finite(t) || (fabs(guess-t) <= tolerance && about_equal(F_t, 0.0, double(tolerance)))) {
			return t;
		}

	}

	// aborted to avoid infinite loop
	return NAN;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline void polynom_N_features(const double F[0], const double f[0], int n, double roots[0], double extrema[0], double inflections[0], int sizes[3], float tolerance, CO_DOMAIN_T domain) {
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
	assert(tolerance > 0);

	int& num_roots = sizes[0];
	int& num_extrema = sizes[1];
	int& num_inflections = sizes[2];

	double* dfdt;
	if (2 < n) {
		// F''(t)
		dfdt = (double*)alloca(sizeof(double)*(n-1));
		polynom_N_derivative(f, n-1, dfdt);
	}
	if (3 < n) {
		// determine extrema and inflections by recursively calling this same function again
		double* unused = roots; // using 'roots' as temp
		polynom_N_features(f, dfdt, n-1, extrema, inflections, unused, sizes, tolerance, domain);
		num_inflections = num_extrema;
		num_extrema = num_roots;
	} else {
		// determine inflection points and extrema using the simpler root finding functions
		if (2 < n) {
			num_inflections = __internal__polynom_N_roots(dfdt, n-2, inflections, tolerance, domain);
			assert(!num_inflections || domain(inflections[0]));
		} else {
			num_inflections = 0;
		}
		if (1 < n) {
			num_extrema = __internal__polynom_N_roots(f, n-1, extrema, tolerance, domain);
			assert(!num_extrema || domain(extrema[0]));
		} else {
			num_extrema = 0;
		}
	}



	if (n < 4) {
		// lucky -> delegate to simpler root function
		num_roots = __internal__polynom_N_roots(F, n, roots, tolerance, domain);
		assert(!num_roots || domain(roots[0]));
	} else {
		// determine roots using newton-raphson method as explained in introduction
		num_roots = 0;
		double t_tolerance = tolerance / 4;

		// merge features
		int num_feats = num_inflections + num_extrema;

		co_domain_dynamic_t interval;

		interval.re_min = domain.re_min;
		double F_feature = polynom_N_value(F, n, interval.re_min);
		if (about_equal(F_feature, 0, (double)tolerance)) {
			assert(num_roots < n);
			roots[num_roots++] = interval.re_min;
		}

		for (int i = 0, i_e = 0, i_i = 0; i < num_feats+1; i++) {
			if (i < num_feats) {
				if (i_e < num_extrema && (i_i == num_inflections || extrema[i_e] < inflections[i_i])) {
					interval.re_max = extrema[i_e++];
				} else {
					assert(i_i < num_inflections);
					interval.re_max = inflections[i_i++];
				}
			} else {
				interval.re_max = domain.re_max;
			}

			assert(interval.re_min <= interval.re_max);

			F_feature = polynom_N_value(F, n, interval.re_max);

			if (about_equal(F_feature, 0, (double)tolerance)) {
				if (!num_roots || !about_equal(roots[num_roots-1], interval.re_max, t_tolerance)) {
					assert(num_roots < n);
					roots[num_roots++] = interval.re_max;
				}
			} else {
				if(__internal__polynom_N_interval_contains_root(F, n, tolerance, interval)) {
					double guess = (interval.re_max + interval.re_min)/2;
					if (!finite(guess*guess)) {
						// TODO: we can get better guesses in case of infinite boundaries, I guess :P
						// for now, we just make sure to search in proper direction
						if (finite(interval.re_min*interval.re_min)) {
							guess = interval.re_min + 100*tolerance;
						} else if (finite(interval.re_max*interval.re_max)) {
							guess = interval.re_max - 100*tolerance;
						} else {
							guess = 0;
						}
					}
					guess = polynom_N_root(F, f, n, guess, tolerance, interval);
					assert (finite(guess));
					if (!num_roots || !about_equal(roots[num_roots-1], guess, t_tolerance)) {
						assert(num_roots < n);
						roots[num_roots++] = guess;
					}
				}
			}

			interval.re_min = interval.re_max;
		}
	}
}


}; // namespace math

#endif /* MATH2D_POLYNOM_H_ */
