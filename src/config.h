/*
 * config.h
 *
 *  Created on: 28 Jun 2019
 *      Author: homac
 */

#ifndef CONFIG_H_
#define CONFIG_H_



// TODO: move to config.h
#ifdef __GNUC__
#define unlikely(BOOLEAN_EXPR) __builtin_expect(BOOLEAN_EXPR, false)
#define likely(BOOLEAN_EXPR)  __builtin_expect(BOOLEAN_EXPR, true)
#else
#define unlikely(BOOLEAN_EXPR) (BOOLEAN_EXPR)
#define likely(BOOLEAN_EXPR)  (BOOLEAN_EXPR)
#endif




/**
 * This function forces a reinterpret_cast, by
 * reinterpreting the address (pointer) of the given item.
 * Thus, it intentionally circumvents restrictions of C++ to
 * achieve compile-time casts.
 *
 * Useful for example, if you have a vector<A*> and a class B derived from A,
 * and you know, that every item in vector<A*> is actually a B*. Then you can cast
 * vector<B*>& vector_B = forced_cast<vector<B*>>(vector_A);
 */
template <class T>
static inline constexpr T& forced_cast(auto& item) {
	assert(sizeof(item) == sizeof(T));
	// cast reference to void pointer and then to T* and then returns T&
	// return (*((T*)( (void*)(&item) )));
	return (*reinterpret_cast<T*>( (void*)(&item) ));
}



#endif /* CONFIG_H_ */
