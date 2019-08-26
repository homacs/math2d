/*
 * precision.h
 *
 *  Created on: 24 Jul 2019
 *      Author: homac
 */

#ifndef MATH_FLOAT_UTILS_H_
#define MATH_FLOAT_UTILS_H_

#include "config.h"
#include <ieee754.h>
#include <stdint.h>
#include <math.h>



const int32_t FLOAT_MANTISSA_SIZE = 23;
const int32_t FLOAT_EXPONENT_SIZE = 8;



/** mantissa most significant bit */
const uint32_t FLOAT_MANTISSA_MSB = (1u<<22);
/** mantissa maximum value */
const uint32_t FLOAT_MANTISSA_MAX = (FLOAT_MANTISSA_MSB<<1)-1;


static inline int float_exponent_biased(float f) {
	ieee754_float result = {.f = f};
	return result.ieee.exponent - IEEE754_FLOAT_BIAS;
}





const int32_t DOUBLE_MANTISSA0_SIZE = 20;
const int32_t DOUBLE_MANTISSA1_SIZE = 32;
const int32_t DOUBLE_EXPONENT_SIZE  = 11;

/** high-word mantissa most significant bit */
const uint32_t DOUBLE_MANTISSA0_MSB = (1u<<19);
/** high-word mantissa maximum value */
const uint32_t DOUBLE_MANTISSA0_MAX = (DOUBLE_MANTISSA0_MSB<<1)-1;
/** low-word mantissa maximum value */
const uint32_t DOUBLE_MANTISSA1_MAX = 0xFFFFFFFFu;

/** combined double mantissa maximum value */
const uint64_t DOUBLE_MANTISSA_FULL_MSB = uint64_t(DOUBLE_MANTISSA0_MSB)<<32;
/** combined double mantissa most significant bit */
const uint64_t DOUBLE_MANTISSA_FULL_MAX = uint64_t(DOUBLE_MANTISSA0_MAX)<<32 | DOUBLE_MANTISSA1_MAX;

const int32_t DOUBLE_MANTISSA_FULL_SIZE = DOUBLE_MANTISSA0_SIZE + DOUBLE_MANTISSA1_SIZE;



static inline int double_exponent_biased(double d) {
	ieee754_double result = {.d = d};
	return result.ieee.exponent - IEEE754_DOUBLE_BIAS;
}


static inline float copysign_only(float f) {
	ieee754_float result = {.f = 1};
	ieee754_float& bits = *((ieee754_float*)&f);
	result.ieee.negative = bits.ieee.negative;
	return result.f;
}

static inline double copysign_only(double f) {
	ieee754_double result = {.d = 1};
	ieee754_double& bits = *((ieee754_double*)&f);
	result.ieee.negative = bits.ieee.negative;
	return result.d;
}

/**
 * Data structure to deal with 52 bit unsigned integer mantissa of
 * ieee753 double precision floating point number.
 *
 * Supports type conversions to/from uint64_t in constexpr.
 */
struct double_mantissa_full_t {
#if	__BYTE_ORDER == __BIG_ENDIAN
	uint32_t mantissa0  :20;
	uint32_t mantissa1  :32;
#endif
#if	__BYTE_ORDER == __LITTLE_ENDIAN
	uint32_t mantissa1  :32;
	uint32_t mantissa0  :20;
#endif



	constexpr double_mantissa_full_t() :
#if	__BYTE_ORDER == __BIG_ENDIAN
	mantissa0(),
	mantissa1()
#else
	mantissa1(),
	mantissa0()
#endif
	{}

	constexpr double_mantissa_full_t(ieee754_double f)
	: double_mantissa_full_t(f.ieee.mantissa0, f.ieee.mantissa1)
	{}


	constexpr double_mantissa_full_t(uint32_t _mantissa0, uint32_t _mantissa1)
	:
#if	__BYTE_ORDER == __BIG_ENDIAN
		mantissa0(_mantissa0),
		mantissa1(_mantissa1)
#endif
#if	__BYTE_ORDER == __LITTLE_ENDIAN
		mantissa1(_mantissa1),
		mantissa0(_mantissa0)
#endif
	{}


	constexpr double_mantissa_full_t(uint64_t mantissa_full) : double_mantissa_full_t() {
		*((uint64_t*)this) = mantissa_full;
	}

	constexpr operator uint64_t() const {
		return *((uint64_t*)this);
	}

	constexpr operator uint64_t&() {
		return *((uint64_t*)this);
	}

};

static inline void ieee754_double_mantissa_add(ieee754_double& f, const uint64_t mantissa_add) {
	double_mantissa_full_t mantissa = double_mantissa_full_t(f);
	uint64_t& full_mantissa = mantissa;
	full_mantissa = full_mantissa + mantissa_add;

	const uint64_t OVERFLOW_BIT = DOUBLE_MANTISSA_FULL_MSB<<1; // one bit beyond upper limit
	const bool overflow = (full_mantissa & OVERFLOW_BIT) != 0;

	f.ieee.exponent += overflow;
	f.ieee.mantissa0 = mantissa.mantissa0;
	f.ieee.mantissa1 = mantissa.mantissa1;

}




typedef double_mantissa_full_t double_mantissa_mask_t;

struct float_mantissa_mask_t {
	uint32_t mantissa;
};




/**
 * Tests if the given value does not exceed the
 * precision identified by precision_key.
 */
static inline bool precision_guard(float a, float_mantissa_mask_t precision_key) {
	ieee754_float f = {.f = a};
	return 0 == (f.ieee.mantissa & ~(precision_key.mantissa));
}

/**
 * Tests if the given value does not exceed the
 * precision identified by precision_key.
 */
static inline bool precision_guard(double a, double_mantissa_mask_t precision_key) {
	ieee754_double f = {.d = a};
	return 0 == ((f.ieee.mantissa0 & ~precision_key.mantissa0) | (f.ieee.mantissa1 & ~precision_key.mantissa1));
}


/**
 * Creates a precision key for the given mantissa reduction.
 * A precision key is required by function precision_trunc().
 *
 * @param mantissa_reduction Number of bits to be removed
 * 	from the least significant bits of a float mantissa
 * 	when used with precision_trunc().
 */
constexpr float_mantissa_mask_t float_mantissa_mask(unsigned char mantissa_reduction) {
	const uint32_t mantissa = FLOAT_MANTISSA_MAX & (FLOAT_MANTISSA_MAX << mantissa_reduction);
	return float_mantissa_mask_t ({mantissa});
}

constexpr double_mantissa_mask_t double_mantissa_mask(unsigned char mantissa_reduction) {
	const int shift = mantissa_reduction;
	const uint64_t mask = DOUBLE_MANTISSA_FULL_MAX << shift;
	return double_mantissa_mask_t(mask);
}

static constexpr float_mantissa_mask_t STD_HALF_PRECISION = float_mantissa_mask(FLOAT_MANTISSA_SIZE-10);


/**
 * Truncates mantissa of given float value to comply with given precision key.
 */
static inline float float_mantissa_trunc(const float& a, const float_mantissa_mask_t& precision_key) {
	ieee754_float f = {.f = a};
	f.ieee.mantissa &= precision_key.mantissa;
	return f.f;
}

static inline double double_mantissa_trunc(const double& a, const double_mantissa_mask_t& precision_key) {
	ieee754_double f = {.d = a};
	f.ieee.mantissa0 &= precision_key.mantissa0;
	f.ieee.mantissa1 &= precision_key.mantissa1;
	return f.d;
}

static inline void float_mantissa_add(ieee754_float& f, const int32_t mantissa_add) {
	uint32_t f_mantissa = f.ieee.mantissa + mantissa_add;

	const uint32_t OVERFLOW_BIT = FLOAT_MANTISSA_MSB<<1; // one bit beyond upper limit
	const bool overflow = (f_mantissa & OVERFLOW_BIT) != 0;

	f.ieee.exponent += overflow;
	f.ieee.mantissa = f_mantissa;
}


/**
 * rounds mantissa of given float value to comply with given precision key.
 */
static inline float float_mantissa_round(const float& a, const float_mantissa_mask_t& precision_key) {
	ieee754_float f = {.f = a};

	// 1/2 * LSB of the given mantissa mask
	uint32_t precision_mantissa_half = ((~precision_key.mantissa) & (precision_key.mantissa>>1)) | (0==precision_key.mantissa);
	float_mantissa_add(f, precision_mantissa_half);
	return float_mantissa_trunc(f.f, precision_key);
}



/**
 * rounds mantissa of given float value to comply with given precision key.
 */
static inline float float_mantissa_floor(const float& a, const float_mantissa_mask_t& precision_key) {
	ieee754_float f = {.f = a};

	if (f.ieee.negative) {
		uint32_t precision_mantissa_lsb = (((~precision_key.mantissa)<<1) & (precision_key.mantissa)) | (0==precision_key.mantissa);
		float_mantissa_add(f, precision_mantissa_lsb);
		return float_mantissa_trunc(f.f, precision_key);
	} else {
		return float_mantissa_trunc(a, precision_key);
	}
}

/**
 * rounds mantissa of given float value to comply with given precision key.
 */
static inline float float_mantissa_ceil(const float& a, const float_mantissa_mask_t& precision_key) {
	return   -1.0 * float_mantissa_floor(-1.0 * a, precision_key);
}

static inline double double_mantissa_round(const double& a, const double_mantissa_mask_t& precision_key) {
	const uint64_t& precision_key_mantissa = (const uint64_t&)(precision_key);
	// 1/2 * LSB of the given mantissa mask
	uint64_t precision_key_half = precision_key_mantissa;
	precision_key_half = ((~precision_key_half) & (precision_key_half>>1)) | (0==precision_key_half);

	ieee754_double f = {.d = a};
	ieee754_double_mantissa_add(f, precision_key_half);
	return double_mantissa_trunc(f.d, precision_key);
}


/**
 * rounds mantissa of given float value to comply with given precision key.
 */
static inline double double_mantissa_floor(const double& a, const double_mantissa_mask_t& precision_key) {
	ieee754_double f = {.d = a};

	if (f.ieee.negative) {
		const uint64_t& precision_key_mantissa = (const uint64_t&)(precision_key);
		uint64_t precision_mantissa_lsb = (((~precision_key_mantissa)<<1) & (precision_key_mantissa)) | (0==precision_key_mantissa);
		ieee754_double_mantissa_add(f, precision_mantissa_lsb);
		return double_mantissa_trunc(f.d, precision_key);
	} else {
		return double_mantissa_trunc(a, precision_key);
	}
}

/**
 * rounds mantissa of given float value to comply with given precision key.
 */
static inline double double_mantissa_ceil(const double& a, const double_mantissa_mask_t& precision_key) {
	return   -1.0 * double_mantissa_floor(-1.0 * a, precision_key);
}

/**
 * Truncates value to a given precision.
 *
 * Calculates the nearest multiple of 'precision'
 * from 'value' towards 0.
 *
 * Unlike ceil() and floor(), trunc() simply cuts off
 * any fractional part below precision.
 *
 * Let 'integer' be an infinite integer type, then
 * the pseudo code for this function would be:
 *
 *    trunc_frac := double(integer(value/precision)) * precision
 *
 * @param v (value)
 * @param p (precision) Smallest accepted positive value, not null.
 * @return
 */
static inline float trunc_frac(float value, float precision) {
	value = value / precision;
	modff(value, &value); // trunc()
	value = value * precision;
	return value;
}


static inline double trunc_frac(double value, double precision) {
	value = value / precision;
	modf(value, &value); // trunc()
	value = value * precision;
	return value;
}


/** round value towards multiples of precision */
static inline float round_frac(float value, float precision) {
	float neg = 1.0f - ((value < 0)*2.0f);
	value = value + neg*precision/2.0;
	return trunc_frac(value, precision);
}

/** round value towards multiples of precision */
static inline double round_frac(double value, double precision) {
	double neg = 1.0f - ((value < 0)*2.0f);
	value = value + neg*precision/2.0;
	return trunc_frac(value, precision);
}




/**
 * truncates to fit in so-called half precision floating point
 */
static inline float half(float a) {
	return float_mantissa_trunc(a,STD_HALF_PRECISION);
}

static inline bool about_equal(double a, double b, double const & tolerance) {
	assert(tolerance >= 0.0);
	double dist = fabs(a-b);
	return dist <= tolerance;
}

static inline bool about_equal(float a, float b, float const & tolerance) {
	assert(tolerance >= 0.0);
	float dist = fabs(a-b);
	return dist <= tolerance;
}


static inline bool about_equal(double a, double b, double_mantissa_mask_t const & precision) {
	return double_mantissa_round(a, precision) == double_mantissa_round(b, precision);
}

static inline bool about_equal(float a, float b, float_mantissa_mask_t const & precision) {
	return float_mantissa_round(a, precision) == float_mantissa_round(b, precision);
}



#endif /* MATH_FLOAT_UTILS_H_ */
