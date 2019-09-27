/*
 * polynom.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH2D_POLYNOM_H_
#define MATH2D_POLYNOM_H_

#include <stdexcept>
#include <math2d/general.h>
#include <string.h>

namespace math2d {

template  <int N>
struct polynom_t;


static inline double polynom1_value(double a, double b, double t);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom1_real_roots(double a, double b, double t[1], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

static inline double polynom2_value(double a, double b, double c, double t);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom2_real_roots(double a, double b, double c, double t[2], CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom2_real_extrema(double a, double b, bool saddle, double& t_min, double& t_max, CO_DOMAIN_T domain = CO_DOMAIN_REAL);


static inline double polynom3_value(double a, double b, double c, double d, double t);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom3_real_extrema(double a, double b, double c, bool saddle, double& t_min, double& t_max, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom3_real_roots(double a, double b, double c, double d, double t[3], CO_DOMAIN_T domain = CO_DOMAIN_REAL);


static inline void polynom_N_derivative(const double F[0], double f[0], int N);

static inline double polynom_N_value(const double f[0], int n, double t);

/**
 * Find one root close to 'guess' using the Newton-Raphson method.
 * Either returns t of a root or +-INFINITY or NAN to indicate, that
 * there is no root to be found. Use ::finite(t) to check whether it's
 * a valid result.
 */
static inline double polynom_N_root(const double F[0], const double f[0], int n, double guess, float tolerance);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom_N_roots(const double F[0], const double f[0], int n, double guesses[0], int num_guesses, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom_N_roots(const double F[0], double f[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

template <int N, class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom_N_roots(const double F[N], double roots[N], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);

////////////////////////////////////////////////////////////////////////////////////
//          I M P L E M E N T A T I O N S     B E L O W
////////////////////////////////////////////////////////////////////////////////////



template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int __internal__polynom_N_roots(const double W[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL);



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

	template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
	int roots(double roots[N], float tolerance, CO_DOMAIN_T domain = CO_DOMAIN_REAL) {
		return polynom_N_roots<N>(p,roots, tolerance, domain);
	}

	void derivative(polynom_t<N-1>& f) const {
		polynom_N_derivative(p, f.p, N);
	}

};



static inline double polynom1_value(double a, double b, double t) {
	return a*t + b;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom1_real_roots(double a, double b, double t[1], CO_DOMAIN_T domain) {
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


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom2_real_roots(double a, double b, double c, double t[2], CO_DOMAIN_T domain) {
	int count = 0;

	if (a == 0) {
		count = polynom1_real_roots(b, c, t, domain);
	} else {

		// 0 = at^2 + bt + c
		// t_1,2 = -b/(2a) +- sqrt ( (b/(2a))^2 - c/a )

		double t_i;

		double b_a_2 = b/(a*2.0f);
		double inner = b_a_2*b_a_2 - c/a;
		if (inner >= 0) {
			double sqrt_inner = sqrt(inner);
			t_i = -b_a_2 + sqrt_inner;
			if (domain(t_i)) t[count++] = t_i;
			// prevent double root
			if (sqrt_inner)	{
				t_i = -b_a_2 - sqrt_inner;
				if (domain(t_i)) t[count++] = t_i;
			}
		}
	}
	return count;
}

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
 *
 * @return 0: no extrema, 1: minimum and/or maximum, -1 saddle point t_min==t_max
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom2_real_extrema(double a, double b, bool saddle, double& t_min, double& t_max, CO_DOMAIN_T domain) {
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

		n = polynom1_real_roots(2.f*a,b, roots, domain);

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
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom3_real_roots(double a, double b, double c, double d, double t[3], CO_DOMAIN_T domain)
{
	static const double PI =  math2d::CONSTANT_PI;
	int count = 0;

	if (a == 0) {
		return polynom2_real_roots(b,c,d, t, domain);
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
			if (domain(t_i)) t[count++] = t_i;
			/* 2 complex roots
				t_i = -(u+v)/2.0 - A_3 + Imaginary((u-v)/2.0 * sqrt(3.0));
				t_i = -(u+v)/2.0 - A_3 - Imaginary((u-v)/2.0 * sqrt(3.0));
			*/
		} else if (D == 0) {
			if (p == 0) {
				t_i = -A_3;  // triple
				if (domain(t_i)) t[count++] = t_i;
			} else {
				double q_3_p = 3.0*q/p;
				t_i =   q_3_p       - A_3;
				if (domain(t_i)) t[count++] = t_i;

				t_i = - q_3_p/(2.0) - A_3; // double
				if (domain(t_i)) t[count++] = t_i;
			}
		} else /* (D < 0) */ {
			double sqrt_4_3_p = ::sqrt(-4.0/3.0*p);
			double onethird_arccos = 1.0/3.0 * ::acos(-q/2.0 * sqrt(-27.0/(p*p*p)));
			double PI_3 = PI/3.0;
			t_i = - sqrt_4_3_p * ::cos(onethird_arccos + PI_3) - A_3;
			if (domain(t_i)) t[count++] = t_i;

			t_i =   sqrt_4_3_p * ::cos(onethird_arccos       ) - A_3;
			if (domain(t_i)) t[count++] = t_i;

			t_i = - sqrt_4_3_p * ::cos(onethird_arccos - PI_3) - A_3;
			if (domain(t_i)) t[count++] = t_i;
		}
	}

    return count;
}


/**
 *
 * Searches extrema and saddle point
 * for a given polynom
 * 	  f(t) = a t^3 + b t^2 + c t + d
 * and writes corresponding values for t in t_min, t_max.
 *
 * Term d is omitted in the function parameters, because it is not used.
 *
 * Infinite values in t_min and t_max represent non-existing or invalid extrema.
 *
 * @return 0: no extrema, 1: minimum and/or maximum, -1 saddle point t_min==t_max
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom3_real_extrema(double a, double b, double c, bool saddle, double& t_min, double& t_max, CO_DOMAIN_T domain) {
	double roots[2];
	int n;
	int result = 0;
	t_min = INFINITY;
	t_max = INFINITY;


	if (a == 0) {
		return polynom2_real_extrema(b,c, saddle, t_min, t_max, domain);
	} else {

		// f(t)=a t^3 + b t^2 + c t + d
		// f'(t)= 3at^2 + 2bt + c
		// f''(t)= 6at + 2b

		double A = 3.f * a;
		double B = 2.f * b;
		n = polynom2_real_roots( A, B, c, roots, domain);

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
			} else if (saddle) {
				// saddle point
				t_min = t_max = roots[i];
				result = -1;
			}
		}
		return result;
	}
}



static inline void polynom_N_derivative(const double F[0], double f[0], int N) {
	for (int i = 0; i < N; i++) {
		f[i] = (N-i) * F[i];
	}
}



static inline double polynom_N_value(const double f[0], int n, double t) {
	double result = f[0];
	for (int i = 1; i < n+1; i++) {
		result *= t;
		result += f[i];
	}
	return result;
}

static inline double polynom_N_root(const double F[0], const double f[0], int n, double guess, float tolerance) {
	double t = guess;
	do {
		guess = t;
		t = guess - polynom_N_value(F, n, guess) / polynom_N_value(f, n-1, guess);
	} while (::fabs(guess-t) > tolerance);
	return t;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom_N_roots(const double F[0], const double f[0], int n, double guesses[0], int num_guesses, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	int count = 0;
	double t;
	double min = +INFINITY;
	double max = -INFINITY;
	for (int i = 0; i < num_guesses; i++) {
		t = polynom_N_root(F, f, n, guesses[i], tolerance);
		if (domain(t)) {
			roots[count++] = t;
			if (t < min) min = t;
			if (t > max) max = t;
		}
	}


	if (num_guesses > 1) {
		// the shape suggests, that there may be other roots at interval edges

		// check lower boundary
		t = polynom_N_root(F, f, n, domain.re_min, tolerance);
		if (domain(t)) {
			roots[count++] = t;
		} else {
			// re_min was a very bad guess
			t = min - (max-min);
			t = polynom_N_root(F, f, n, t, tolerance);
			if (domain(t)) roots[count++] = t;

		}

		// check upper boundary
		t = polynom_N_root(F, f, n, domain.re_max, tolerance);
		if (domain(t)) {
			roots[count++] = t;
		} else {
			// re_max was a very bad guess
			t = max + (max-min);
			t = polynom_N_root(F, f, n, t, tolerance);
			if (domain(t)) roots[count++] = t;
		}
	}
	return count;
}


template <int N, class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom_N_roots(const double F[N], double roots[N], float tolerance, CO_DOMAIN_T domain) {
	int i;
	for (i = 0; i <= N && F[i] == 0; i++);
	int n = N-i;
	const double* W = F + i;
	return __internal__polynom_N_roots(W, n, roots, tolerance, domain);
}




template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom_N_roots(const double F[0], double f[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	if (n < 3 || n > 5) {
		// FIXME: make root guesses for any degree of a function
		throw std::range_error("cannot deal with n < 3 or n > 5");
	}

	// generate guesses using the second derivate
	double* dfdt = (double*)alloca(sizeof(double)*(n-1));
	polynom_N_derivative(f, dfdt, n-1);

	double* guesses = (double*)alloca(sizeof(double)*(n-2));
	int num_guesses = __internal__polynom_N_roots( dfdt, n-2, guesses, tolerance);

	// calculate roots
	return polynom_N_roots(F, f, n, guesses, num_guesses, roots, tolerance, domain);

}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int __internal__polynom_N_roots(const double W[0], int n, double roots[0], float tolerance, CO_DOMAIN_T domain) {
	assert(W[0] != 0);

	switch(n) {
	case 0:
		return 0;
	case 1:
		return polynom1_real_roots(W[0], W[1], roots, domain);
	case 2:
		return polynom2_real_roots(W[0], W[1], W[2], roots, domain);
	case 3:
		return polynom3_real_roots(W[0], W[1], W[2], W[3], roots, domain);
	default:
		{
			double* w = (double*)alloca(sizeof(double)*n);
			polynom_N_derivative(W, w, n);
			return polynom_N_roots(W, w, n, roots, tolerance, domain);
		}
	}
}



}; // namespace math

#endif /* MATH2D_POLYNOM_H_ */