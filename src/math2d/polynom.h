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

template  <int ORDER>
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
template <int ORDER, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[ORDER], double roots[ORDER], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

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
static inline void polynom_N_features(const double F[0], const double f[0], int degree,
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



/**
 * WARNING: intervals cannot contain extrema!
 * checks if there is a root between interval boundaries (boundaries excluded).
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline double __internal__polynom_N_interval_contains_root(const double F[0], int n, double tolerance, CO_DOMAIN_T two_extrema_interval);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots(const double F[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);


/**
 * Defines a polynom f(t) of degree DEGREE.
 * Function f(t) := p[0] * t^N + p[1] * t^(N-1) + ... + p[N]
 *
 * Required: DEGREE > 0
 */
template  <int ORDER>
struct polynom_t {
	double p[ORDER+1];
	polynom_t() {assert(ORDER > 0);}
	polynom_t(const double _p[ORDER+1])
	{
		assert(ORDER > 0);
		memcpy(p, _p, sizeof(double)*(ORDER+1));
	}

	int degree() {
		return ORDER;
	}

	/** returns value of the polynom at t */
	double operator () (double t) const {
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
	double single_root(double guess, float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) const {
		// TODO: remove method single_root
		polynom_t<ORDER-1> dfdt;
		derivative(dfdt);
		return polynom_N_root(p, dfdt.p, ORDER, guess, tolerance, domain);
	}

	template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
	int roots(double roots[ORDER], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) const {
		return polynom_N_roots<ORDER>(p, roots, tolerance, domain);
	}

	void derivative(polynom_t<ORDER-1>& f) const {
		polynom_N_derivative(p, ORDER, f.p);
	}

};



static inline double polynom1_value(double a, double b, double t) {
	return a*t + b;
}

/**
 * A polynomial of degree 0 is a constant.
 * A constant cannot have a distinguishable x for a root.
 * A constant is either 0.0 and every x is a root of f(x)
 * or not 0.0 and no root exists.
 *
 * This function just exists to provide reasoning for this
 * case.
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom0_roots(double a) {
	return 0;
}



template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom1_roots(double a, double b, double t[1], CO_DOMAIN_T domain) {
	int count = 0;
	// 0 = a*t + b
	// t = -b/a

	if (a == 0) {
		return polynom0_roots(b);
	} else {
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
	} else if (c == 0) {
		// one root at 0.0 + remaining in g(x) = f(x)/(x-0)
		if (domain(0.0)) roots[count++] = 0.0;
		count += polynom1_roots(a, b, roots + count, domain);
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
	double* end = std::unique(roots, roots + count);
	count = end - roots;

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
	} else if (d == 0) {
		// one root at 0.0 + remaining in g(x) = f(x)/(x-0)
		if (domain(0.0)) roots[count++] = 0.0;
		count += polynom2_roots(a,b,c, roots + count, domain);
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
	double* end = std::unique(roots, roots + count);
	count = end - roots;

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


static inline double polynom4_value(double a, double b, double c, double d, double e, double t) {
	return t * polynom3_value(a, b, c, d, t) + e;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom4_roots(double a, double b, double c, double d, double e, double roots[4], CO_DOMAIN_T domain = CO_DOMAIN_REAL)
{
	typedef long double real_t;

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

		real_t DELTA_0 = c*c - 3.0*b*d + 12.0*a*e;
		real_t DELTA_1 = 2.0*c*c*c -9.0*b*c*d + 27.0*b*b*e + 27.0*a*d*d - 72.0*a*c*e;
		real_t DELTA   = - 1.0/27.0 * (DELTA_1*DELTA_1 - 4.0*DELTA_0*DELTA_0*DELTA_0);

		// FIXME: check if DELTA has correct sign by comparing it to its original form

		real_t P = 8.0*a*c - 3.0*b*b;
		real_t D = 64.0*a*a*a*e - 16.0*a*a*c*c + 16.0*a*b*b*c - 16.0*a*a*b*d - 3.0*b*b*b*b;

		real_t p = P/(8.0*a*a);                                 // always finite since a != 0
		real_t q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d)/(8.0*a*a*a); // always finite since a != 0

		real_t S;

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
				real_t PHI;
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

			real_t Q;
			real_t sqrt_inner;

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
					real_t sign = 1 - (i%2)*2;

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

		double t_i;
		real_t sqrt_inner;
		real_t sqrt_inner_last_term;
		real_t sqrt_inner_first_terms = -4.0*S*S - 2.0*p;

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


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom4_roots_2(double a, double b, double c, double d, double e, double roots[4], CO_DOMAIN_T domain = CO_DOMAIN_REAL)
{
	// This version is based on the ferrari method.
	// It differs from other implementations in that it
	// delegates to other functions, rather than making
	// the solution and decision making unreadable.

	typedef long double real_t;

	int count = 0;
	if (a == 0.0) {
		// actually 3rd order polynom
		return polynom3_roots(b,c,d,e, roots, domain);
	} else if (b == 0.0 && d == 0.0) {
		// substitute: x = t^2
		// to: a*x^2 + c*x + e = a*t^4 + c*t^2 + e
		double x_i[2];
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
		// Substitute according to Tschirnhaus
		real_t alpha = - 3.0*b*b/(8.0*a*a)   + c/a;
		real_t beta  = + (b*b*b)/(8.0*a*a*a) - (b*c)/(2.0*a*a) + d/a;
		real_t gamma = - (3.0*b*b*b*b)/(256.0*a*a*a*a) + (b*b*c)/(16.0*a*a*a) - (b*d)/(4.0*a*a) + e/a;
		// resolves in:
		//     0 = u^4 + alpha u^2 + beta u + gamma
		// with
		//     x_i = u_i - b/(4 a)
		real_t x_u_offset = - b/(4.0 * a);



		if (beta == 0.0) {
			// special case:
			//     0 = u^4 + alpha u^2 + gamma
			// substitute:
			//     z = u^2
			// and solve
			//     0 = z^2 + alpha z + gamma
			// with
			//     x_i = +/- sqrt(z_1,2) - b/(4 a)
			double z_i[2];
			// we do not care about complex solutions in z,
			// because sqrt(z) gives another complex number
			// and we are not interested in complex roots in x.
			int n = polynom2_roots(1, alpha, gamma, z_i, domain);
			double x;
			for (int i = 0; i < n; i++) {
				if (z_i[i] >= 0) {
					x = sqrt(z_i[i]) + x_u_offset;
					if (domain(x)) roots[count++] = x;

					x = -sqrt(z_i[i]) + x_u_offset;
					if (domain(x)) roots[count++] = x;
				}
			}
		} else {
			// beta != 0

			// Common case
			// Solving with Ferrari's method

			// substitutions result in:
			//     0 = y^3 + p1 y^2 + p2 y + p3

			double y_i[3];
			double p0 = 1;
			double p1 = 5.0/2.0*alpha;
			double p2 = (2.0*alpha*alpha - gamma);
			double p3 = (4.0*alpha*(alpha*alpha - gamma)-beta*beta)/8.0;
			int n = polynom3_roots(p0, p1, p2, p3, y_i, domain);

			// NOTE: A third degree polynomial always has at least one real root (not complex).
			//       To determine the roots in x, we can use any root in y
			//       and will always get the same results.
			//       Since there is at least one real y, we will use that.

			// this will fire only, if there is an error in polynom_3_root(1, ..) (see reason above)
			assert(n);
			real_t y = y_i[0];
			// now determine the u_1,2,3,4
			real_t sqrt_inner_1 = alpha + 2.0*y;
			if (sqrt_inner_1 >= 0) {

				real_t w = sqrt(sqrt_inner_1);

				// According to en wiki on ferrari's solution:
				// w is never gonna be zero, because if (w == 0)
				// then also (beta == 0) and we would have been
				// in the branch above (see beta == 0).
				// But keep that assert here!
				assert(w != 0);
				real_t z = beta/(2.0*w);

				for (int i = 0; i < 4; i++) {
					real_t sign_1 = +1 - 2 * (i/2);
					real_t sign_2 = +1 - 2 * (i%2);

					real_t sqrt_inner_2 = w*w - 4.0 * (alpha + y + sign_1 * -z);
					if (sqrt_inner_2 < 0) continue; // -> 2 complex roots

					// real root
					real_t u_i = 0.5*(sign_1 * -w + sign_2 * sqrt(sqrt_inner_2));
					real_t x_i = u_i + x_u_offset;
					if (domain(x_i)) roots[count++] = x_i;
				}
			} else {
				// the first term is a complex number, which will result in real roots, if
				// (and only if) the second term is a complex number too.
				// Thus, we have to deal with complex numbers here.
				typedef std::complex<real_t> complex_t;

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
					// TODO: when introducing different real types, adjust tolerance for (Im == 0)
					// accept all complex with (|Im - 0| < 2^-16) as real numbers
					if (domain(x_i.real()) && about_equal(x_i.imag(),0, 1e-16)) roots[count++] = x_i.real();
				}

			}
		}


	}
	std::sort(roots, roots + count);
	double* end = std::unique(roots, roots + count);
	count = end - roots;
	return count;
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

/**
 * Determine maximum interval in which roots of given function can appear.
 *
 * Note: For real roots only.
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline bool polynom_N_root_bounds(const double F[0], int n, double& t_lower, double& t_upper, float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) {


	// We initialise bounds with the domain bounds
	// and then try to improve them, using known methods.
	t_lower = domain.lower;
	t_upper = domain.upper;

	if (n > 1) {
		if (F[0] == 0) return polynom_N_root_bounds(F+1, n-1, t_lower, t_upper, tolerance, domain);

		// TODO: optimise: evaluate, remove ineffective methods and merge remaining methods.

		// Currently we are using all well known methods and
		// take the best bounds. Each method has its own block
		// to keep them distinguishable.

		{
			// Cauchy's Bound
			double B = ::fabs(F[1]/F[0]);
			for (int i = 2; i <= n; i++) {
				B = std::max(B, ::fabs(F[i]/F[0]));
			}
			B += 1;

			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);

		}

		{
			// Lagrange's Bound
			double B = 0;
			for (int i = 1; i <= n; i++) {
				B += ::fabs(F[i]/F[0]);
			}
			B = std::max(1.0, B);

			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);
		}

		{
			// Zassenhausen's Bound
			double B = ::fabs(F[1]/F[0]);
			for (int i = 2; i <= n; i++) {
				B = std::max(B, ::pow(::fabs(F[i]/F[0]), 1.0/i));
			}
			B = 2 * B;

			t_lower = std::max(t_lower, -B);
			t_upper = std::min(t_upper, +B);

		}

		{
			// Zassenhausen's Bound improved by Lagrange
			double first_largest = - INFINITY;
			double second_largest = -INFINITY;
			double B;
			for (int i = 1; i <= n; i++) {
				B = std::max(B, ::pow(::fabs(F[i]/F[0]), 1.0/i));
				if (B > first_largest) {
					first_largest = B;
				} else if (B > second_largest) {
					second_largest = B;
				}
			}
			B = first_largest + second_largest;

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

#ifndef NDEBUG
		// test if there is no root between root boundary and domain boundary
		bool contains_root;
		if (t_lower > domain.lower) {
			co_domain_const_t off_boundary(double(domain.lower), t_lower);
			contains_root = __internal__polynom_N_interval_contains_root(F, n, tolerance, off_boundary);
			assert(!contains_root);
		}
		if (t_upper < domain.upper) {
			co_domain_const_t off_boundary(t_upper, double(domain.upper));
			contains_root = __internal__polynom_N_interval_contains_root(F, n, tolerance, off_boundary);
			assert(!contains_root);
		}
#endif
	}
	return improved;
}

template <int ORDER, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[ORDER], double roots[ORDER], float tolerance, CO_DOMAIN_T domain) {
	int i;
	for (i = 0; unlikely(i <= ORDER && F[i] == 0); i++);
	int n = ORDER-i;
	if (n <= 0) return 0;
	const double* W = F + i;
	return __internal__polynom_N_roots(W, n, roots, tolerance, domain);
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[0], double f[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	int sizes[3];
	double extrema[n>0?n-1:0];
	double inflections[n>1?n-2:0];
	int& num_roots = sizes[0];
	double lower_bound = domain.lower;
	double upper_bound = domain.upper;
	polynom_N_root_bounds(F, n, lower_bound, upper_bound, tolerance, domain);
	co_domain_const_t root_domain(lower_bound, upper_bound);
	polynom_N_features(F, f, n, roots, extrema, inflections, sizes, tolerance, root_domain);
	return num_roots;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int __internal__polynom_N_roots(const double F[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	assert(F[0] != 0);

	switch(n) {
	case 0:
		return polynom0_roots(F[0]);
	case 1:
		return polynom1_roots(F[0], F[1], roots, domain);
	case 2:
		return polynom2_roots(F[0], F[1], F[2], roots, domain);
#if MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS
	case 3:
		return polynom3_roots(F[0], F[1], F[2], F[3], roots, domain);
	case 4:
		return polynom4_roots_2(F[0], F[1], F[2], F[3], F[4], roots, domain);
#endif
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
	double v_min = polynom_N_value(F, n, two_extrema_interval.lower);
	double v_max = polynom_N_value(F, n, two_extrema_interval.upper);

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

	// we cannot accept 0 tolerance, because
	// at some point the result cannot be improved,
	// due to floating point inaccuracies in evaluation
	// of F(t).
	assert(tolerance > 0);
	assert(finite(guess));
	assert(domain(guess));

	// Newton-Raphson method
	do {
		POLYNOM_N_ROOTS_COUNT();
		guess = t;

		F_t = polynom_N_value(F, n, guess);
		if (F_t == 0) return t;

		f_t = polynom_N_value(f, n-1, guess);
		// this will only be null, if the interval contains an extrema
		// which is forbidden
		assert (f_t != 0);

		t = guess - F_t / f_t;
		if (t >= domain.upper) {
			t = (domain.upper + guess)/2;
		} else if (t <= domain.lower) {
			t = (guess + domain.lower)/2;
		}
		assert(finite(t) && domain(t));

		if ((about_equal(guess, t, double(tolerance)) && about_equal(F_t, 0.0, double(tolerance)))) {
			guess = t;
			// -> will end the loop and proceed with improvement on t
		}
		// abort if either we run out of iterations
		// or we can't improve
	} while (--iterations && t != guess);



	if (t == guess) {
		// There are three reasons to improve on the results we have got from
		// the newton-raphson method:
		// 1. t didn't change in last iteration of newton-raphson method.
		// 2. f(t) is close to 0 which results in a large interval for t in which F(t) is about equal 0.
		// 3. newton-raphson method closes in on a root, but never crosses it.
		//
		// The first case means, that the accuracy of the function to
		// compute F(t) is not high enough to reach the required accuracy of F_t.
		// The second case means, that even if the last step between guess and t
		// is smaller than tolerance, the actual distance from t_0 may be much
		// larger. The third issue means, that t will move closer to the root
		// but never reaches it. Crossing the root once will improve the accuracy
		// significantly.
		//
		// To solve this, we are now going to improve the given t by stepping
		// towards F(t) = 0 until we cross it. If we have found the two values
		// t_1 and t_2 which are left and right of the actual root t_0, we
		// then interpolate between them.

		// compute a minimal step-width
		// If the step is too small for the current t so that
		// (t + step == t), then we use the smallest possible
		// increment instead.
		double step = tolerance*2;

		double F_t_previous;

		// check the actual values for F(t) and f(t)
		// to determine in which direction we have to search.
		F_t = polynom_N_value(F, n, t);
		f_t = polynom_N_value(f, n-1, t);
		if (F_t == 0) return t;
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
			if (t >= domain.upper) {
				t = (domain.upper + guess)/2;
			} else if (t <= domain.lower) {
				t = (guess + domain.lower)/2;
			}
			assert(finite(t) && domain(t));

			// check F(t)
			F_t = polynom_N_value(F, n, guess);
			if (F_t*F_t_previous < 0) {
				// We have crossed the X-axis.
				// The actual t_0 is somewhere between guess and t
				return (t+guess)/2.0;
			}

		} while (--iterations);
	}

	// aborted to avoid infinite loop or we can't reach requested accuracy
	return NAN;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline void polynom_N_features(const double F[0], const double f[0], int degree, double roots[0], double extrema[0], double inflections[0], int sizes[3], float tolerance, CO_DOMAIN_T domain) {
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
	if (2 < degree) {
		// F''(t)
		dfdt = (double*)alloca(sizeof(double)*(degree-1));
		polynom_N_derivative(f, degree-1, dfdt);
	}
	if (4 < degree) {
		// determine extrema and inflections by recursively calling this same function again
		double* unused = roots; // using 'roots' as temp
		polynom_N_features(f, dfdt, degree-1, extrema, inflections, unused, sizes, tolerance, domain);
		num_inflections = num_extrema;
		num_extrema = num_roots;
	} else {
		// determine inflection points and extrema using the simpler root finding functions
		if (2 < degree) {
			num_inflections = __internal__polynom_N_roots(dfdt, degree-2, inflections, tolerance, domain);
			assert(!num_inflections || domain(inflections[0]));
		} else {
			num_inflections = 0;
		}
		if (1 < degree) {
			num_extrema = __internal__polynom_N_roots(f, degree-1, extrema, tolerance, domain);
			assert(!num_extrema || domain(extrema[0]));
		} else {
			num_extrema = 0;
		}
	}



	if (degree < 4) {
		// lucky -> delegate to simpler root function
		num_roots = __internal__polynom_N_roots(F, degree, roots, tolerance, domain);
		assert(!num_roots || domain(roots[0]));
	} else {
		// determine roots using newton-raphson method as explained in introduction
		num_roots = 0;
		double t_tolerance = tolerance / 4;


		// TODO: guesses for first and last roots can still be improved
		// since boundaries are actually guesses for the first and last root.
		// But we have to know, whether the boundary is computed or a
		// given boundary of the application. Which means, the boundary
		// computation has to move here, maybe.
		// See polynom_N_root_bounds().
		double lower_bound = domain.lower;
		double upper_bound = domain.upper;


		co_domain_std_t interval;

		interval.lower = lower_bound;
		double F_feature = polynom_N_value(F, degree, interval.lower);
		if (about_equal(F_feature, 0, (double)tolerance)) {
			assert(num_roots < degree);
			roots[num_roots++] = interval.lower;
		}

		// consider all features and upper bound
		int num_feats = num_inflections + num_extrema;
		for (int i = 0, i_e = 0, i_i = 0; i < num_feats+1; i++) {
			if (i < num_feats) {
				if (i_e < num_extrema && (i_i == num_inflections || extrema[i_e] < inflections[i_i])) {
					interval.upper = extrema[i_e++];
				} else {
					assert(i_i < num_inflections);
					interval.upper = inflections[i_i++];
				}
			} else {
				interval.upper = upper_bound;
			}

			assert(interval.lower <= interval.upper);

			F_feature = polynom_N_value(F, degree, interval.upper);

			if (about_equal(F_feature, 0, (double)tolerance)) {
				if (!num_roots || !about_equal(roots[num_roots-1], interval.upper, t_tolerance)) {
					assert(num_roots < degree);
					roots[num_roots++] = interval.upper;
				}
			} else {
				if(__internal__polynom_N_interval_contains_root(F, degree, tolerance, interval)) {
					double guess = (interval.upper + interval.lower)/2;
					// This should not happen because we have estimated upper
					// and lower bounds for roots and the boundaries should be
					// out of computational range.
					assert(finite(guess*guess));
					guess = polynom_N_root(F, f, degree, guess, tolerance, interval);
					assert(finite(guess));
					if (!num_roots || !about_equal(roots[num_roots-1], guess, t_tolerance)) {
						assert(num_roots < degree);
						roots[num_roots++] = guess;
					}
				}
			}

			interval.lower = interval.upper;
		}
	}
}


}; // namespace math

#endif /* MATH2D_POLYNOM_H_ */
