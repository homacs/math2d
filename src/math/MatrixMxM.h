/*
 * AdjacencyMatrix.h
 *
 *  Created on: 9 Jul 2019
 *      Author: homac
 */

#ifndef MATH_MATRIXMXM_H_
#define MATH_MATRIXMXM_H_

#include <string.h>
#include <stdlib.h>
#include <assert.h>
namespace math {

template <class T>
class MatrixMxM {
	unsigned capacity;
	bool autoshrink;
	unsigned growth;
	unsigned size;
	T* data;
public:
	/**
	 * Create m x m matrix.
	 * @param _capacity: Expected number of rows or columns.
	 */
	MatrixMxM(unsigned _capacity, bool _autoshrink) {
		growth = 8;
		size = 0;
		capacity = 0;
		data = NULL;
		setAutoshrink(_autoshrink);
		setCapacity(_capacity);
	}

	MatrixMxM() : MatrixMxM(8,false) {}


	virtual ~MatrixMxM() {
		delete data;
	}

	unsigned getCapacity() {
		return capacity;
	}

	void setCapacity(unsigned _capacity) {
		if (_capacity > capacity) {
			resize(_capacity);
		} else if (_capacity < capacity) {
			assert(size <= _capacity);
			resize(_capacity);
		}
	}

	unsigned getSize() {
		return size;
	}


	void setSize(unsigned _size) {
		if (_size > capacity) {
			resize(_size);
			size = _size; // set size afterwards, to copy only those we had
		} else if(autoshrink && _size < capacity - growth) {
			size = _size; // set size before, to copy only those we keep
			resize(_size);
		} else if (_size < size) {
			// reset all cells, which are now invalid

			unsigned r = 0;
			// reset remainder of rows
			for (; r < size; r++) {
				memset(cell(data, capacity, r, size), 0, sizeof(T) * (size - _size));
			}
			// reset remaining old rows
			memset(row(data, capacity, r), 0, sizeof(T) * (size * (size-r)));
			size = _size;
		}
	}


	void resetColumn(int c) {
		for (unsigned r = 0; r < size; r++) {
			cell(data, capacity, r, c) = 0;
		}
	}

	void resetRow(int r) {
		memset(row(data, capacity, r), 0, size);
	}


	void setAutoshrink(bool _autoshrink) {
		this->autoshrink = _autoshrink;
	}


	T& get(unsigned row, unsigned column) {
		assert (row <= size);
		assert (column <= size);
		return *cell(data, capacity, row, column);
	}


	T set(unsigned row, unsigned column, const T& v) {
		assert (row <= size);
		assert (column <= size);
		T& c = get(row, column);
		T result = c;
		c = v;
		return result;
	}

	class RowEnumerator {
		T* e;
		T* end;
	public:
		RowEnumerator(const MatrixMxM* _M, unsigned r) {
			e = row(_M->data, _M->capacity, r);
			end = (e + _M->size);
		}
		RowEnumerator(const RowEnumerator& it) {
			e = it.e;
			end = it.end;
		}

		bool hasMore() {
			return e != end;
		}

		bool operator != (RowEnumerator& it) {
			return e != it->e;
		}

		RowEnumerator& operator ++ (int) {
			e++;
			return *this;
		}

		T& operator ->() {
			return *e;
		}

		T& operator *() {
			return *e;
		}
	};

	class ColumnEnumerator {
		unsigned step;
		T* e;
		T* end;

	public:

		ColumnEnumerator(const MatrixMxM* _M, unsigned _c) {
			step = _M->capacity;
			e = cell(_M->data, _M->capacity, 0, _c);
			end = e+(step*_M->size);
		}
		ColumnEnumerator(const ColumnEnumerator& it) {
			step = it.step;
			e = it.e;
			end = it.end;
		}

		bool hasMore() {
			return e != end;
		}

		bool operator != (ColumnEnumerator& it) {
			return e != it->e;
		}

		ColumnEnumerator& operator ++ (int) {
			e += step;
			return *this;
		}

		T& operator ->() {
			return *e;
		}

		T& operator *() {
			return *e;
		}

	};


	ColumnEnumerator getColumn(unsigned _c) {
		return ColumnEnumerator(this, _c);
	}

	RowEnumerator getRow(unsigned _r) {
		return RowEnumerator(this, _r);
	}


private:


	static T* row(T* M, unsigned capacity, unsigned r) {
		return &(M[r*capacity]);
	}

	static T* cell(T* M, unsigned capacity, unsigned _row, unsigned _column) {
		return row(M, capacity, _row) + _column;
	}

	T* copy(T* dest, T* src) {
		assert(dest != NULL);
		for (unsigned r = 0; r < size; r++) {
			T* row_old = row(src, capacity, r);
			T* row_new = row(dest, capacity, r);
			memcpy(row_new, row_old, size * sizeof(T*));
		}
		return dest;
	}

	int align(unsigned _capacity, unsigned alignment) {
		unsigned alignedCapacity = _capacity / alignment;
		alignedCapacity += ((_capacity % alignment) != 0);
		return alignedCapacity * alignment;
	}

	void resize(unsigned _capacity) {
		unsigned newCapacity = align(_capacity, growth);
		T* newData = (T*) calloc(newCapacity*newCapacity, sizeof(T));
		data = copy(newData, data);
		capacity = newCapacity;
	}

};

} /* namespace math */

#endif /* MATH_MATRIXMXM_H_ */
