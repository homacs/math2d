/*
 * general.h
 *
 *  Created on: 3 Aug 2019
 *      Author: homac
 */

#ifndef MATH_GENERAL_H_
#define MATH_GENERAL_H_



namespace math {


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
template <typename T, T re_min = -INFINITY, T re_max = INFINITY, T im_min = 0, T im_max = 0>
struct co_domain_t {
	bool operator() (T re) const {
		return re_min <= re && re <= re_max;
	}
	bool operator() (T re, T im) const {
		return (im_min <= im && im <= im_max) && *this(re);
	}
};

const struct CO_DOMAIN_REAL_NUMBERS_T {
	bool operator() (double re) const {
		return true;
	}
	bool operator() (double re, double im) const {
		return (im == 0);
	}
} CO_DOMAIN_REAL;

const struct CO_DOMAIN_COMPLEX_NUMBERS_T {
	bool operator() (double re) const {
		return true;
	}
	bool operator() (double re, double im) const {
		return true;
	}
} CO_DOMAIN_COMPLEX;


const struct CO_DOMAIN_REAL_IN_0_1_INCLUSIVE_T {
	bool operator() (double value, double im = 0.0) const {
		return 0 <= value && value <= 1 && (im == 0.0);
	}
} CO_DOMAIN_REAL_IN_0_1_INCLUSIVE; /**< All real numbers in range [0:1] */
static const struct CO_DOMAIN_REAL_IN_0_1_EXCLUSIVE_T {
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


#endif /* MATH_GENERAL_H_ */
