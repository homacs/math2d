/*
 * config.h
 *
 *  Created on: 28 Jun 2019
 *      Author: homac
 */

#ifndef MATH2D_CONFIG_H_
#define MATH2D_CONFIG_H_



#ifndef MATH2D_EVALUATE
#define MATH2D_EVALUATE 1
#endif


#ifdef __GNUC__
#define unlikely(BOOLEAN_EXPR) __builtin_expect(BOOLEAN_EXPR, false)
#define likely(BOOLEAN_EXPR)  __builtin_expect(BOOLEAN_EXPR, true)
#else
#define unlikely(BOOLEAN_EXPR) (BOOLEAN_EXPR)
#define likely(BOOLEAN_EXPR)  (BOOLEAN_EXPR)
#endif


#define MATH2D_POLYNOM_N_ROOT_MAX_ITERATIONS 1000000

#ifndef MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS
	/**
	 * This flag determines, whether the root finding
	 * function polynom_N_roots uses arithmetics for all
	 * polynomials with degree <= 4.
	 *
	 * It can be turned off for debugging purposes.
	 */
#	define MATH2D_POLYNOM_N_ROOT_USE_ARITHMETICS   1
#endif




#endif /* MATH2D_CONFIG_H_ */
