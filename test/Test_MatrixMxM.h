/*
 * TestMatrixMxM.h
 *
 *  Created on: 9 Jul 2019
 *      Author: homac
 */

#ifndef MATH_TEST_MATRIXMXM_H_
#define MATH_TEST_MATRIXMXM_H_


#include <math2d/MatrixMxM.h>



namespace math2d {

struct Test_MatrixMxM {



	static bool all() {
		MatrixMxM<int> m;

		m.setSize(10);


		for (unsigned i = 0; i < m.getSize(); i++) {
			int v = i;
			m.set(i,i,v);
			int j = m.get(i,i);
			assert(v == j);
		}

		for (unsigned i = 0; i < m.getSize(); i++) {
			int v = i;
			int j = m.get(i,i);
			assert(v == j);
		}

		for (unsigned col = 0; col < m.getSize(); col++) {
			MatrixMxM<int>::RowEnumerator e = m.getRow(col);
			for (unsigned i = 0; e.hasMore(); i++, e++) {
				int j = *e;
				if (i==col) assert(j==(int)col);
				else assert(j==0);
			}
		}

		for (unsigned row = 0; row < m.getSize(); row++) {
			MatrixMxM<int>::ColumnEnumerator e = m.getColumn(row);
			for (unsigned i = 0; e.hasMore(); i++, e++) {
				int j = *e;
				if (i==row) assert(j==(int)row);
				else assert(j==0);
			}
		}

		return true;
	}

};






} /* namespace math */

#endif /* MATH_TEST_MATRIXMXM_H_ */
