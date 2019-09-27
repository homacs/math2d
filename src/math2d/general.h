/*
 * general.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH2D_GENERAL_H_
#define MATH2D_GENERAL_H_
#include <limits.h>


namespace math2d {


// higher accuracy than math.h : M_PI
static const long double CONSTANT_PI = 3.14159265358979323846264338327950288;



/**
 * IMPORTANT: this is currently used only to specify valid ranges for results!
 *
 * Template to describe a co-domain for mathematical functions.
 *
 * Generally the co-domain describes the range of valid numbers.
 * Bounds of value ranges are inclusive, meaning values equal to
 * min or max are accepted as well.
 * To limit the domain to real numbers, im_min and im_max have to
 * be 0 . This will accept all real numbers and complex numbers
 * with im == 0.
 *
 * A co-domain can be described by other means too. For example
 * a function like
 *
 * <pre>
 * 		bool mydomain(double re, double im = 0.0) {
 * 			return (im==0) && 0 < re <= 256;
 * 		}
 * </pre>
 *
 * will accept all real numbers in range ]0:256].
 *
 * For most cases where bounds are unlimited, there predefined
 * structures, which have better runtime performance.
 *
 *		CO_DOMAIN_REAL       all real numbers.
 *		CO_DOMAIN_COMPLEX    all complex numbers.
 *
 *
 * @param T scalar number type (float, double, long, int, ...)
 * @param re_min lower bound for real numbers
 * @param re_max upper bound for real numbers
 * @param im_min lower bound for imaginary part of complex numbers
 * @param im_max upper bound for imaginary part of complex numbers
 */
template <typename T = double>
struct co_domain_t {
	const T re_min = -DBL_MAX;
	const T re_max = DBL_MAX;
	const T im_min = -DBL_MAX;
	const T im_max = DBL_MAX;

	co_domain_t(const double _re_min = -DBL_MAX, const double _re_max = DBL_MAX, double _im_min = -DBL_MAX, double _im_max = DBL_MAX)
	:
		re_min(_re_min),
		re_max(_re_max),
		im_min(_im_min),
		im_max(_im_max)
	{}

	bool operator() (T re) const {
		return re_min <= re && re <= re_max;
	}
	bool operator() (T re, T im) const {
		return (im_min <= im && im <= im_max) && *this(re);
	}
};

const struct CO_DOMAIN_REAL_NUMBERS_T {
	const double re_min = -DBL_MAX;
	const double re_max = DBL_MAX;
	const double im_min = 0;
	const double im_max = 0;
	bool operator() (double re) const {
		return finite(re);
	}
	bool operator() (double re, double im) const {
		return finite(re) && (im == 0);
	}
} CO_DOMAIN_REAL;

const struct CO_DOMAIN_COMPLEX_NUMBERS_T {
	const double re_min = -DBL_MAX;
	const double re_max = DBL_MAX;
	const double im_min = -DBL_MAX;
	const double im_max = DBL_MAX;
	bool operator() (double re) const {
		return finite(re);
	}
	bool operator() (double re, double im) const {
		return finite(re) && finite(im);
	}
} CO_DOMAIN_COMPLEX;

const struct CO_DOMAIN_REAL_IN_NEG_1_POS_1_T {
	const double re_min = -1;
	const double re_max = 1;
	const double im_min = 0;
	const double im_max = 0;
	bool operator() (double value, double im = 0.0) const {
		return re_min <= value && value <= re_max && (im == 0.0);
	}
} CO_DOMAIN_REAL_IN_NEG_1_POS_1; /**< All real numbers in range [-1:1] */

const struct CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T {
	const double re_min = 0;
	const double re_max = 1;
	const double im_min = 0;
	const double im_max = 0;
	bool operator() (double value) const {
		return re_min <= value && value <= re_max;
	}
	bool operator() (double value, double im) const {
		return re_min <= value && value <= re_max && (im == 0.0);
	}
} CO_DOMAIN_REAL_IN_0_1_INCLUSIVE; /**< All real numbers in range [0:1] */
static const struct CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE_T {
	const double re_min = 0;
	const double re_max = 1;
	const double im_min = 0;
	const double im_max = 0;
	bool operator() (double value) const {
		return 0 < value && value < 1;
	}
} CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE; /**< All real numbers in range ]0:1[ */


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

}; // namespace math


#endif /* MATH2D_GENERAL_H_ */
