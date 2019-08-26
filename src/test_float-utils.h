/*
 * test_float-utils.h
 *
 *  Created on: 27 Jul 2019
 *      Author: homac
 */

#ifndef TEST_FLOAT_UTILS_H_
#define TEST_FLOAT_UTILS_H_

#include "float-utils.h"
#include "perf_clock.h"


namespace Test_float_utils {

static inline ieee754_float ieee754_float_init(bool negative, uint32_t mantissa, int8_t exponent) {
	ieee754_float f;
	f.ieee.negative = negative;
	f.ieee.mantissa = mantissa;
	f.ieee.exponent = exponent;
	return f;
}

static inline ieee754_double ieee754_double_init(bool negative, uint32_t mantissa0, int32_t mantissa1, uint16_t exponent) {
	ieee754_double d;
	d.ieee.negative = negative;
	d.ieee.mantissa0 = mantissa0;
	d.ieee.mantissa1 = mantissa1;
	d.ieee.exponent = exponent;
	return d;
}


static inline void relative_float_precision() {
	float_mantissa_mask_t relative_precision;

	relative_precision = float_mantissa_mask(1);
	assert(relative_precision.mantissa == FLOAT_MANTISSA_MAX-1);
	relative_precision = float_mantissa_mask(FLOAT_MANTISSA_SIZE-1);
	assert(relative_precision.mantissa == FLOAT_MANTISSA_MSB);

	ieee754_float f;
	ieee754_float result;
	ieee754_float expected;

	f = ieee754_float_init(0, FLOAT_MANTISSA_MSB | (FLOAT_MANTISSA_MSB>>1), IEEE754_FLOAT_BIAS-1);
	result.f = float_mantissa_trunc(f.f, relative_precision);
	expected = ieee754_float_init(f.ieee.negative, FLOAT_MANTISSA_MSB, f.ieee.exponent);
	assert(result.f == expected.f);

	result.f = float_mantissa_round(f.f, relative_precision);
	expected = ieee754_float_init(f.ieee.negative, 0, f.ieee.exponent+1);
	assert(result.f == expected.f);
}

static inline void relative_double_precision() {
	double_mantissa_mask_t relative_precision;

	relative_precision = double_mantissa_mask(1);
	assert(relative_precision.mantissa0 == DOUBLE_MANTISSA0_MAX);
	assert(relative_precision.mantissa1 == DOUBLE_MANTISSA1_MAX-1);
	relative_precision = double_mantissa_mask(DOUBLE_MANTISSA_FULL_SIZE-1);
	assert(relative_precision.mantissa0 == DOUBLE_MANTISSA0_MSB);
	assert(relative_precision.mantissa1 == 0x0u);

	ieee754_double f;
	ieee754_double result;
	ieee754_double expected;




	f = ieee754_double_init(0,DOUBLE_MANTISSA0_MSB | (DOUBLE_MANTISSA0_MSB>>1),	0, uint16_t(IEEE754_DOUBLE_BIAS-1));
	result.d = double_mantissa_trunc(f.d, relative_precision);
	expected = ieee754_double_init(f.ieee.negative, DOUBLE_MANTISSA0_MSB, 0, f.ieee.exponent);
	assert(result.d == expected.d);

	result.d = double_mantissa_round(f.d, relative_precision);
	expected = ieee754_double_init(f.ieee.negative, 0, 0, f.ieee.exponent+1);
	assert(result.d == expected.d);
}


static inline void test_round_frac_double() {

	double value;
	double result;
	double expected;
	double precision = 0.025;

	result = round_frac(value = 0.024,precision);
	expected = precision;
	assert(result == expected);

	// unfortunately value = p/2 does not round to 0
	// due to lost accuracy in the least significant bits
	// This value is the upper limit for p = 0.025
	ieee754_double v = {.d = precision/2};
	v.ieee.mantissa1 -= 2;
	result = round_frac(value = v.d,precision);
	expected = 0.0;
	assert(result == expected);

	result = round_frac(0.006,precision);
	expected = 0.0;
	assert(result == expected);


	// same with negative values

	result = round_frac(value = -0.024,precision);
	expected = -precision;
	assert(result == expected);

	v = ieee754_double({.d = -precision/2});
	v.ieee.mantissa1 -= 2;
	result = round_frac(value = v.d,precision);
	expected = -0.0;
	assert(result == expected);

	result = round_frac(-0.006,precision);
	expected = -0.0;
	assert(result == expected);





	// some variations

	result = round_frac(value = 10000.06, precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	result = round_frac(value = 10.06,precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	precision = 1.0;
	result = round_frac(value = 10.05,precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	precision = 1.0;
	result = round_frac(value = 101.0,precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	result = round_frac(value = 1.87123e+27, precision = 0.025);
	expected = value;
	assert(result == expected);
}



static inline void test_round_frac_float() {

	float value;
	float result;
	float expected;
	float precision = 0.025;

	result = round_frac(value = 0.024,precision);
	expected = precision;
	assert(result == expected);

	ieee754_float v = {.f = precision/2};
	v.ieee.mantissa -= 1;
	result = round_frac(value = v.f,precision);
	expected = 0.0;
	assert(result == expected);

	result = round_frac(value = 0.006,precision);
	expected = 0.0;
	assert(result == expected);


	result = round_frac(value = 10000.06, precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	result = round_frac(value = 10.06,precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	precision = 1.0;
	result = round_frac(value = 10.05,precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	precision = 1.0;
	result = round_frac(value = 101.0,precision);
	expected = double(int(value/precision))*precision;
	assert(result == expected);

	result = round_frac(value = 1.87123e+27, precision = 0.025);
	expected = value;
	assert(result == expected);




	result = round_frac(value = -0.024,precision);
	expected = -precision;
	assert(result == expected);

	v = ieee754_float({.f = -precision/2});
	v.ieee.mantissa -= 1;
	result = round_frac(value = v.f,precision);
	expected = 0.0;
	assert(result == expected);

	result = round_frac(value = -0.006,precision);
	expected = 0.0;
	assert(result == expected);

}


static inline float float_precision_round(float v, float p) {
	ieee754_float f = {.f = v};
	int e = f.ieee.exponent;
	f.ieee.exponent = IEEE754_FLOAT_BIAS;
	f.f = round_frac(f.f, p);
	f.ieee.exponent = e + (f.ieee.exponent - IEEE754_FLOAT_BIAS);
	return f.f;
}


static inline void test_relative_float() {
	float v = 0.0001117234;
	float lsb = 1.0f/8;
	v = float_precision_round(v, lsb);

	v = float_precision_round(0, lsb);
}


static inline void all() {

	test_relative_float();


	relative_float_precision();
	relative_double_precision();
	test_round_frac_float();
	test_round_frac_double();
}


};




#endif /* TEST_FLOAT_UTILS_H_ */
