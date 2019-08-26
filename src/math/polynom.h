/*
 * polynom.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH_POLYNOM_H_
#define MATH_POLYNOM_H_

#include <math/general.h>

namespace math {



static inline double polynom1_value(double a, double b, double t) {
	return a*t + b;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom1_real_roots(double a, double b, double t[1], CO_DOMAIN_T interval = CO_DOMAIN_REAL) {
	int count = 0;
	// 0 = a*t + b
	// t = -b/a

	if (a != 0) {
		double t_i = -b/a;
		if (interval(t_i)) t[count++] = t_i;
	}
	return count;
}


static inline double polynom2_value(double a, double b, double c, double t) {
	return t * polynom1_value(a, b, t) + c;
}


template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom2_real_roots(double a, double b, double c, double t[2], CO_DOMAIN_T interval = CO_DOMAIN_REAL) {
	int count = 0;

	if (a == 0) {
		count = polynom1_real_roots(b, c, t, interval);
	} else {

		// 0 = at^2 + bt + c
		// t_1,2 = -b/(2a) +- sqrt ( (b/(2a))^2 - c/a )

		double t_i;

		double b_a_2 = b/(a*2.0f);
		double inner = b_a_2*b_a_2 - c/a;
		if (inner >= 0) {
			double sqrt_inner = sqrt(inner);
			t_i = -b_a_2 + sqrt_inner;
			if (interval(t_i)) t[count++] = t_i;
			// prevent double root
			if (sqrt_inner)	{
				t_i = -b_a_2 - sqrt_inner;
				if (interval(t_i)) t[count++] = t_i;
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
static inline int polynom2_real_extrema(double a, double b, bool saddle, double& t_min, double& t_max, CO_DOMAIN_T interval = CO_DOMAIN_REAL) {
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

		n = polynom1_real_roots(2.f*a,b, roots, interval);

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
 *
 * @return 0: no extrema, 1: minimum and/or maximum, -1 saddle point t_min==t_max
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom3_real_extrema(double a, double b, double c, bool saddle, double& t_min, double& t_max, CO_DOMAIN_T interval = CO_DOMAIN_REAL) {
	double roots[2];
	int n;
	int result = 0;
	t_min = INFINITY;
	t_max = INFINITY;


	if (a == 0) {
		return polynom2_real_extrema(b,c, saddle, t_min, t_max, interval);
	} else {

		// f(t)=a t^3 + b t^2 + c t + d
		// f'(t)= 3at^2 + 2bt + c
		// f''(t)= 6at + 2b

		double A = 3.f * a;
		double B = 2.f * b;
		n = polynom2_real_roots( A, B, c, roots, interval);

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
 *
 */
template <class CO_DOMAIN_T = CO_DOMAIN_REAL_NUMBERS_T>
static inline int polynom3_real_roots(double a, double b, double c, double d, double t[3], CO_DOMAIN_T domain = CO_DOMAIN_REAL)
{
	// NEW VERSION
	static const double PI =  math::CONSTANT_PI;
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

static inline double polynom3_value(double a, double b, double c, double d, double t) {
	return t * polynom2_value(a, b, c, t) + d;
}


}; // namespace math

#endif /* MATH_POLYNOM_H_ */
