/*
 * perf_clock.h
 *
 *  Created on: 23 Jun 2019
 *      Author: homac
 */

#ifndef __PERF_CLOCK_H_
#define __PERF_CLOCK_H_

#include <stdlib.h>
#include <math.h>

struct perf_time_t : timespec {
	static const long NSEC_PER_SEC = 1000000000;


	perf_time_t() {}

	perf_time_t(const long long& nsec) {
		this->operator=(nsec);
	}

	perf_time_t(const double& _secs_nsecs) {
		secs_nsecs(_secs_nsecs);
	}

	perf_time_t& operator=(const long long& nsec) {
		tv_sec = nsec/NSEC_PER_SEC;
		tv_nsec = abs(nsec%NSEC_PER_SEC);
		return *this;
	}

	operator long long() const {
		return ((long long)tv_sec)*NSEC_PER_SEC + tv_nsec;
	}

	perf_time_t& secs_nsecs(const double& _secs_nsecs) {
		double secs;
		double nsecs = modf(_secs_nsecs, &secs);
		tv_sec = secs;
		tv_nsec = nsecs<0?-nsecs:nsecs;
		return *this;
	}

	double secs_nsecs() {
		return double(tv_sec) + double(tv_nsec)/NSEC_PER_SEC;
	}

};



perf_time_t operator-(const perf_time_t& a, const perf_time_t& b) {
	const long NSEC_PER_SEC = 1000000000;
	long long la = a;
	long long lb = b;
	la = la - lb;
	perf_time_t result = la;
	result.tv_sec  = la/NSEC_PER_SEC;
	result.tv_nsec = la%NSEC_PER_SEC;
	return result;
}


#define PERF_CLOCK CLOCK_REALTIME


void perf_clock(perf_time_t& t) {
	clock_gettime(PERF_CLOCK, &t);
}



#endif /* __PERF_CLOCK_H_ */
