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

template  <int N>
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
template <int N, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[N], double roots[N], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

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
 * This function finds roots based on given guesses.
 *
 * Precondition:
 * 		F[0] != 0
 *
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[0], const double f[0], int n, double guesses[0], int num_guesses, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);


/**
 * Find at most one root, based on the given guess
 * inside a given interval. Roots on interval borders
 * will be ignored.
 *
 *
 * Runtime of the function depends on guess and the
 * interval. Two tips are important in this regard:
 *
 * - You should check first, if your interval actually
 *   contains a root (boundary checks).
 *
 * - Use extrema for intervals, because there is at most
 *   one root between two extrema.
 *
 *
 * It is recommended to choose intervals between
 * extrema or one extreme and +/- DBL_MAX. If this
 * condition is met, then this function will always
 * converge on at most one root and never abort to
 * avoid an infinite loop.
 *
 * The function will return NAN in two cases:
 * (1) There is no root to be found in the given interval.
 * (2) The function did not converge and would have ended
 *     up in an infinite loop. Avoid this by using proper
 *     intervals (see above).
 *
 * @param F function to be searched (polynom of degree n)
 * @param f derivative of F
 * @param n degree of F
 * @param guess guess for t_0 of a possible location of a root for F(t_0) = 0
 * @param tolerance deviation tolerance for found root.
 * @param domain valid interval to be search for the root
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
 * Defines a polynom f(t) of degree N.
 * Function f(t) := p[0] * t^N + p[1] * t^(N-1) + ... + p[N]
 */
template  <int N>
struct polynom_t {
	double p[N+1];
	polynom_t() {assert(N > 0);}
	polynom_t(const double _p[N+1])
	{
		assert(N > 0);
		memcpy(p, _p, sizeof(double)*(N+1));
	}

	int degree() {
		return N;
	}

	/** returns value of the polynom at t */
	double operator () (double t) const {
		return polynom_N_value(p, N, t);
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
	double single_root(double guess, float tolerance) {
		polynom_t<N-1> dfdt;
		derivative(dfdt);
		return polynom_N_root(p, dfdt.p, N, guess, tolerance);
	}

	template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
	int roots(double roots[N], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) {
		return polynom_N_roots<N>(p,roots, tolerance, domain);
	}

	void derivative(polynom_t<N-1>& f) const {
		polynom_N_derivative(p, N, f.p);
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


template <int N, class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline int polynom_N_roots(const double F[N], double roots[N], float tolerance, CO_DOMAIN_T domain) {
	int i;
	for (i = 0; unlikely(i <= N && F[i] == 0); i++);
	int n = N-i;
	const double* W = F + i;
	return __internal__polynom_N_roots(W, n, roots, tolerance, domain);
}

/**
 * WARNING: for intervals between extremas only!
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
static inline int polynom_N_roots(const double F[0], double f[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	// WICHTIG: Newton-Verfahren hat Schwächen: Konvergiert gegen falsche Nullstelle oder gar nicht (Ausnahme)
	// Konvergenz gegen falsche Nullstelle kann durch Bestimmung
	// von Intervallen zwischen Extrema verhindert werden.
	//
	// Zwischen zwei Extrema befindet sich maximal eine Nullstelle!
	// Das gilt auch für doppelte Extrema. Auf einem k-fachen Extrema befinden sich k-1 Nullstellen.
	//
	// Extrema müssen manuell auf Nullstellen getestet werden,
	// weil das Newton-Verfahren in der Nähe von Extream sehr schlecht konvergiert (zu geringer Steigung).
	// Schätzung sollte deshalb immer zwischen Extrema ansetzen.
	//
	// Sollte Verfahren gar nicht konvergieren, dann muss zu anderem Verfahren gewechselt werden.


	// TODO: Kurvendiskussion zur Bestimmung der Intervalle in denen Nullstellen sein können.
	//       1. Bestimme Extrema über Nullstellen 1. Ableitung
	//       2. Stelle fest, ob vor dem 1. Extrema und hinter dem letzten Extrema eine Nullstelle sein kann
	//          Wenn (F(t_e_1,2) == 0) dann nicht
	//          Wenn (F(t_e_1) > 0 && f(t_e_1 - delta) > 0) dann beispielsweise ja. Es gibt 4 Fälle.
	//       3. Bestimme linke und rechte Grenze für Nullstellen anhand Ergebnissen aus 2
	//          Finde die Nullstellen der Ränder durch Wahl eines t != t_e mit einem f(t) von ca. 1.0.
	//          Lässt sich aus den Extrema eine ungefähre Breite abschätzen, kann sie als dt verwendet
	//          werden. Existiert nur ein einziges Extrema, dann kann stattdessen ein sehr großes dt
	//          verwendet werden, beispielsweise dt = 10^5.
	//       4. Aus der geordneten Reihe der Extrema ergibt sich eine Liste von Intervallen, in
	//          denen weitere Nullstellen gesucht werden können. Der Schätzwert ist dabei
	//          jeweils die Mitte des Intervalls.

	// polynom is expected to be already normalised
	assert(F[0] != 0);

	if (n < 4 || n > 5) {
		// FIXME: make root guesses for any degree of a function
		throw std::range_error("cannot deal with n < 4 or n > 5");
	}


	// create the second derivative
	double* dfdt = (double*)alloca(sizeof(double)*(n-1));
	polynom_N_derivative(f, n-1, dfdt);


	int num_extrema;
	double* extrema = (double*)alloca(sizeof(double)*(n-1));
	if (n-1 < 4) {
		num_extrema = __internal__polynom_N_roots(f, n-1, extrema, tolerance);
	} else {
		num_extrema = polynom_N_roots(f, dfdt, n-1, extrema, tolerance, domain);
	}
	// note: extrema are sorted because roots get sorted


	int count = 0;
	double F_extrema;

	if (num_extrema == 0) {
		// function has no extrema

		// find single root

		co_domain_t interval (domain.re_min, domain.re_max, 0, 0);
		// TODO: we can get better guesses I guess :P
		// for now, we just make sure to search in the proper direction
		double guess = extrema[0] - 100*tolerance;
		if(__internal__polynom_N_interval_contains_root(F, n, tolerance, interval)) {
			guess = polynom_N_root(F, f, n, guess, tolerance, interval);
			assert (finite(guess));
			roots[count++] = guess;
		}

	} else {
		// we have at least one extrema

		// determine root between -INFINITY and extrema[0]

		// check if extrema is a root
		F_extrema = polynom_N_value(F, n, extrema[0]);

		if (about_equal(F_extrema, 0, (double)tolerance)) {
			roots[count++] = extrema[0];
		} else if (domain.re_min < extrema[0]) {

			co_domain_t interval (domain.re_min, extrema[0], 0, 0);
			// TODO: we can get better guesses I guess :P
			// for now, we just make sure to search in the proper direction
			double guess = extrema[0] - 100*tolerance;
			if(__internal__polynom_N_interval_contains_root(F, n, tolerance, interval)) {
				guess = polynom_N_root(F, f, n, guess, tolerance, interval);
				assert (finite(guess));
				roots[count++] = guess;
			}
		}

		// generate guesses from extrema
		for (int i = 1; i < num_extrema; i++) {
			F_extrema = polynom_N_value(F, n, extrema[i]);
			if (about_equal(F_extrema, 0, (double)tolerance)) {
				roots[count++] = extrema[i];
			} else if (domain(extrema[i-1]) || domain(extrema[i])) {
				co_domain_t interval (extrema[i-1], extrema[i], 0, 0);
				double guess = (extrema[i-1] + extrema[i])/2;

				if(__internal__polynom_N_interval_contains_root(F, n, tolerance, interval)) {
					guess = polynom_N_root(F, f, n, guess, tolerance, interval);
					assert (finite(guess));
					roots[count++] = guess;
				}
			}
		}

		if (domain(extrema[num_extrema-1])) {
			co_domain_t interval (extrema[num_extrema-1], domain.re_max, 0, 0);
			// TODO: we can get better guesses I guess :P
			// for now, we just make sure to search in the proper direction
			double guess = extrema[num_extrema-1] + 100*tolerance;
			if(__internal__polynom_N_interval_contains_root(F, n, tolerance, interval)) {
				guess = polynom_N_root(F, f, n, guess, tolerance, interval);
				assert (finite(guess));
				roots[count++] = guess;
			}
		}
	}
	return count;

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


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_T>
static inline double polynom_N_root(const double F[0], const double f[0], int n, double guess, float tolerance, CO_DOMAIN_T domain) {


	double F_t;
	double f_t;
	double t = guess;
	double t_avg = t + tolerance;
	double F_t_previous = INFINITY;
	do {
		POLYNOM_N_ROOTS_COUNT();
		guess = t;
		F_t = polynom_N_value(F, n, guess);
		f_t = polynom_N_value(f, n-1, guess);
		if (f_t != 0) {
			t = guess - F_t / f_t;
			if (t >= domain.re_max) {
				t = guess + (domain.re_max-guess)/2;
			} else if (t <= domain.re_min) {
				t = guess - (guess-domain.re_min)/2;
			}
		} else if (F_t == 0) {
			return t;
		}

		if (!finite(t) || (fabs(guess-t) <= tolerance && fabs(F_t-F_t_previous) <= tolerance)) {
			return t;
		} else if (fabs(t-t_avg) < tolerance/2) {
			// will not converge
			return NAN;
		}
		F_t_previous = F_t;
		t_avg = (t_avg + t)/2;
	} while (true);

	return t;
}


}; // namespace math

#endif /* MATH2D_POLYNOM_H_ */
