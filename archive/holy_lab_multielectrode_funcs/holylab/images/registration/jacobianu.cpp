/*
 * The Jacobian matrix and the determinant of Jacobian matrix is calculated
 * basing on
 * 		"Topology preservation and regularity in estimated deformation fields"
 * 		by Bilge Karacah and Christos Davatzikos, 2003.
 *
 * The gradient of the sum of Jacobian determinant with respect to each u_i_new
 * component is simply implemented basing on the partial derivative formula and
 * the successive derivative formula.
 */

#include <cstddef>
#include <cmath> 		/* for log */
#include <limits> 	/* for infinity */
#include <stdlib.h>

namespace Jacobianu {

/* struct pow2: compile-time computation of 2^n */

template <int n>
struct pow2 {
	enum { value = 2*pow2<n-1>::value };
};

template <>
struct pow2<0> {
	enum { value = 1 };
};

/* class jacobianu: computes the Jacobian matrix of the composed deformation
 * u(x) and the determinant of the Jacobian matrix */

template <class T, int n_dims>
class jacobianu {

private:

	T LDetJ; 							/* log(det(J)) */
	T cellLDetJ; 					/* the sum of log(det(J)) in each cell */

	/* element offset along 1D array for each vortex of the cell */
	size_t memoffset[pow2<n_dims>::value];

	/* three independent components of u_composed */
	T uf[pow2<n_dims>::value]; /* uf */
	T ug[pow2<n_dims>::value]; /* ug */
	T uh[pow2<n_dims>::value]; /* uh */

	/* Variable for mapping uf/ug/uh in the Jacobian matrix calculation.
	 * Each row separately counts for the index of cell vertex chosen for
	 * computing each Jacobian matrix. The vertex index is the same as that of
	 * pIu index.
	 * For 2D, each row index corresponds to the 1, 2, 3 vertexes.
	 * For 3D, each row index corresponds to the 1, 2, 3, 4 vertexes. */
	int uMap[pow2<n_dims>::value][n_dims+1];
	/* relative x/y/z coordinate difference map for each vertex of cell relative
	 * to the 0 index vertex of cell */
	T xyzMap[pow2<n_dims>::value][n_dims];

public:

	jacobianu(void) {
		;
	}
	~jacobianu(void) {
		;
	}

	void init(const int szu[]);

	/* snip out a single cell of ucomp from "A" and then store the components in
	 * uf, ug & uh for further usage */
	void stuff(size_t offset, int nArrays, const T* const A[]);

	void calcLDetJ(size_t offset, int dim, T c, T* grad[], T* ucompgrad[]);

	T getCellLDetJ() const {
		return this->cellLDetJ;
	}

};

/* g(x) = x + u_composed(x)
 * J is in column-major order
 * u_composed is stored in grid-then-component order (and grid is column-major)
 * */

template <class T,int n_dims>
void jacobianu<T,n_dims>::init(const int szu[])
{
	this->LDetJ = 0.0;
	this->cellLDetJ = 0.0;

	/* initialize memoffset */
	int dimIdx, nbrIdx, bitmask;
	/* accumulative number of elements from the lower order of dimension */
	size_t cumsz[n_dims];
	cumsz[0] = 1;
	for (dimIdx = 1; dimIdx < n_dims; ++dimIdx) {
		cumsz[dimIdx] = cumsz[dimIdx-1] * szu[dimIdx-1];
	}
	for (nbrIdx = 0; nbrIdx < pow2<n_dims>::value; ++nbrIdx) {
		size_t offset = 0;
		for (dimIdx = 0, bitmask = 1; dimIdx < n_dims;
				dimIdx++, bitmask = bitmask << 1) {
			if (nbrIdx & bitmask) {
				offset += cumsz[dimIdx];
			}
		}
		this->memoffset[nbrIdx] = offset;
	}

	/* initialize uMap */
	int dims = n_dims;
	if (dims == 2) {
		int dummy_uMap[pow2<2>::value][3] = {
				{2,0,1},
				{0,1,3},
				{3,2,0},
				{1,3,2}};
		for (int i = 0; i < pow2<2>::value; ++i) {
			for (int j = 0; j < 3; ++j) {
				this->uMap[i][j] = dummy_uMap[i][j];
			}
		}
	}
	else if (dims == 3) {
		int dummy_uMap[pow2<3>::value][4] = {
				{1,4,2,0},
				{0,3,5,1},
				{0,6,3,2},
				{1,2,7,3},
				{0,5,6,4},
				{1,7,4,5},
				{2,4,7,6},
				{3,6,5,7}};
		for (int i = 0; i < pow2<3>::value; ++i) {
				for (int j = 0; j < 4; ++j) {
					this->uMap[i][j] = dummy_uMap[i][j];
				}
		}
	}

	/* initialize xyzMap */
	if (dims == 2) {
		T dummy_xyzMap[pow2<2>::value][2] = {
				{0.0, 0.0},
				{1.0, 0.0},
				{0.0, 1.0},
				{1.0, 1.0}};
		for (int i = 0; i < pow2<2>::value; ++i) {
			for (int j = 0; j < 2; ++j) {
				this->xyzMap[i][j] = dummy_xyzMap[i][j];
			}
		}
	}
	else if (dims == 3) {
		T dummy_xyzMap[pow2<3>::value][3] = {
				{0.0, 0.0, 0.0},
				{1.0, 0.0, 0.0},
				{0.0, 1.0, 0.0},
				{1.0, 1.0, 0.0},
				{0.0, 0.0, 1.0},
				{1.0, 0.0, 1.0},
				{0.0, 1.0, 1.0},
				{1.0, 1.0, 1.0}};
		for (int i = 0; i < pow2<3>::value; ++i) {
			for (int j = 0; j < 3; ++j) {
				this->xyzMap[i][j] = dummy_xyzMap[i][j];
			}
		}
	}
}

template <class T,int n_dims>
void jacobianu<T,n_dims>::stuff(size_t offset, int nArrays, const T * const A[])
{
	for (int aIdx = 0; aIdx < nArrays; ++aIdx) {
		for (int vIdx = 0; vIdx < pow2<n_dims>::value; ++vIdx) {
			if (aIdx == 0) {
				this->uf[vIdx] = A[aIdx][offset + this->memoffset[vIdx]];
			}
			else if (aIdx == 1) {
				this->ug[vIdx] = A[aIdx][offset + this->memoffset[vIdx]];
			}
			else if (aIdx == 2) {
				this->uh[vIdx] = A[aIdx][offset + this->memoffset[vIdx]];
			}
		}
	}
}

template <class T,int n_dims>
void jacobianu<T,n_dims>::calcLDetJ(size_t offset, int dim, T c, T* grad[],
		T* ucompgrad[])
{
	this->LDetJ = 0.0;
	this->cellLDetJ = 0.0;

	for (int vIdx = 0; vIdx < pow2<n_dims>::value; ++vIdx) {

		T detJ = 0.0;
		int uIdx = 0;
		int uOffset = 0;
		if (dim == 2) {
			uIdx = this->uMap[vIdx][0];
			T f1 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g1 = this->ug[uIdx] + this->xyzMap[uIdx][1];

			uIdx = this->uMap[vIdx][1];
			T f2 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g2 = this->ug[uIdx] + this->xyzMap[uIdx][1];

			uIdx = this->uMap[vIdx][2];
			T f3 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g3 = this->ug[uIdx] + this->xyzMap[uIdx][1];

			detJ = f1*g2 - f2*g1 - f1*g3 + f3*g1 + f2*g3 - f3*g2;

			if (detJ < 0.0) { /* avoid calculating gradient */
					this->LDetJ = std::numeric_limits<T>::infinity();
					this->cellLDetJ += this->LDetJ;
			}
			else {
				this->LDetJ = log(detJ / c);
				this->cellLDetJ += this->LDetJ * this->LDetJ;

				T df1 = g2 - g3;
				T dg1 = f3 - f2;
				T df2 = g3 - g1;
				T dg2 = f1 - f3;
				T df3 = g1 - g2;
				T dg3 = f2 - f1;

				if (grad[0] != NULL) {
					uOffset = offset + this->memoffset[this->uMap[vIdx][0]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df1 * ucompgrad[0][uOffset*2] + dg1 * ucompgrad[0][uOffset*2+1]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df1 * ucompgrad[1][uOffset*2] + dg1 * ucompgrad[1][uOffset*2+1]);

					uOffset = offset + this->memoffset[this->uMap[vIdx][1]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df2 * ucompgrad[0][uOffset*2] + dg2 * ucompgrad[0][uOffset*2+1]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df2 * ucompgrad[1][uOffset*2] + dg2 * ucompgrad[1][uOffset*2+1]);

					uOffset = offset + this->memoffset[this->uMap[vIdx][2]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df3 * ucompgrad[0][uOffset*2] + dg3 * ucompgrad[0][uOffset*2+1]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df3 * ucompgrad[1][uOffset*2] + dg3 * ucompgrad[1][uOffset*2+1]);
				}
			}
		}
		else if (dim == 3) {
			uIdx = this->uMap[vIdx][0];
			T f1 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g1 = this->ug[uIdx] + this->xyzMap[uIdx][1];
//			T h1 = -1.0*(this->uh[uIdx] + this->xyzMap[uIdx][2]);/*right-handed coord*/
			T h1 = this->uh[uIdx] + this->xyzMap[uIdx][2];

			uIdx = this->uMap[vIdx][1];
			T f2 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g2 = this->ug[uIdx] + this->xyzMap[uIdx][1];
//			T h2 = -1.0*(this->uh[uIdx] + this->xyzMap[uIdx][2]);
			T h2 = this->uh[uIdx] + this->xyzMap[uIdx][2];

			uIdx = this->uMap[vIdx][2];
			T f3 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g3 = this->ug[uIdx] + this->xyzMap[uIdx][1];
//			T h3 = -1.0*(this->uh[uIdx] + this->xyzMap[uIdx][2]);
			T h3 = this->uh[uIdx] + this->xyzMap[uIdx][2];

			uIdx = this->uMap[vIdx][3];
			T f4 = this->uf[uIdx] + this->xyzMap[uIdx][0];
			T g4 = this->ug[uIdx] + this->xyzMap[uIdx][1];
//			T h4 = -1.0*(this->uh[uIdx] + this->xyzMap[uIdx][2]);
			T h4 = this->uh[uIdx] + this->xyzMap[uIdx][2];

			detJ = f1*g3*h2 - f1*g2*h3 + f2*g1*h3 - f2*g3*h1 - f3*g1*h2 + f3*g2*h1 +
					   f1*g2*h4 - f1*g4*h2 - f2*g1*h4 + f2*g4*h1 + f4*g1*h2 - f4*g2*h1 -
					   f1*g3*h4 + f1*g4*h3 + f3*g1*h4 - f3*g4*h1 - f4*g1*h3 + f4*g3*h1 +
					   f2*g3*h4 - f2*g4*h3 - f3*g2*h4 + f3*g4*h2 + f4*g2*h3 - f4*g3*h2;

			if (detJ < 0.0) {
				this->LDetJ = std::numeric_limits<T>::infinity();
				this->cellLDetJ += this->LDetJ;
			}
			else {
				this->LDetJ = log(detJ / c);
				this->cellLDetJ += this->LDetJ * this->LDetJ;

				T df1 = g3*h2 - g2*h3 + g2*h4 - g4*h2 - g3*h4 + g4*h3;
				T dg1 = f2*h3 - f3*h2 - f2*h4 + f4*h2 + f3*h4 - f4*h3;
				T dh1 = f3*g2 - f2*g3 + f2*g4 - f4*g2 - f3*g4 + f4*g3;

				T df2 = g1*h3 - g3*h1 - g1*h4 + g4*h1 + g3*h4 - g4*h3;
				T dg2 = f3*h1 - f1*h3 + f1*h4 - f4*h1 - f3*h4 + f4*h3;
				T dh2 = f1*g3 - f3*g1 - f1*g4 + f4*g1 + f3*g4 - f4*g3;

				T df3 = g2*h1 - g1*h2 + g1*h4 - g4*h1 - g2*h4 + g4*h2;
				T dg3 = f1*h2 - f2*h1 - f1*h4 + f4*h1 + f2*h4 - f4*h2;
				T dh3 = f2*g1 - f1*g2 + f1*g4 - f4*g1 - f2*g4 + f4*g2;

				T df4 = g1*h2 - g2*h1 - g1*h3 + g3*h1 + g2*h3 - g3*h2;
				T dg4 = f2*h1 - f1*h2 + f1*h3 - f3*h1 - f2*h3 + f3*h2;
				T dh4 = f1*g2 - f2*g1 - f1*g3 + f3*g1 + f2*g3 - f3*g2;

				if (grad[0] != NULL) {
					uOffset = offset + this->memoffset[this->uMap[vIdx][0]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df1 * ucompgrad[0][uOffset*3] + dg1 * ucompgrad[0][uOffset*3+1] +
							 dh1 * ucompgrad[0][uOffset*3+2]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df1 * ucompgrad[1][uOffset*3] + dg1 * ucompgrad[1][uOffset*3+1] +
							 dh1 * ucompgrad[1][uOffset*3+2]);
					grad[2][uOffset] += 2.0 * this->LDetJ / detJ *
							(df1 * ucompgrad[2][uOffset*3] + dg1 * ucompgrad[2][uOffset*3+1] +
							 dh1 * ucompgrad[2][uOffset*3+2]);

					uOffset = offset + this->memoffset[this->uMap[vIdx][1]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df2 * ucompgrad[0][uOffset*3] + dg2 * ucompgrad[0][uOffset*3+1] +
							 dh2 * ucompgrad[0][uOffset*3+2]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df2 * ucompgrad[1][uOffset*3] + dg2 * ucompgrad[1][uOffset*3+1] +
							 dh2 * ucompgrad[1][uOffset*3+2]);
					grad[2][uOffset] += 2.0 * this->LDetJ / detJ *
							(df2 * ucompgrad[2][uOffset*3] + dg2 * ucompgrad[2][uOffset*3+1] +
							 dh2 * ucompgrad[2][uOffset*3+2]);

					uOffset = offset + this->memoffset[this->uMap[vIdx][2]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df3 * ucompgrad[0][uOffset*3] + dg3 * ucompgrad[0][uOffset*3+1] +
							 dh3 * ucompgrad[0][uOffset*3+2]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df3 * ucompgrad[1][uOffset*3] + dg3 * ucompgrad[1][uOffset*3+1] +
							 dh3 * ucompgrad[1][uOffset*3+2]);
					grad[2][uOffset] += 2.0 * this->LDetJ / detJ *
							(df3 * ucompgrad[2][uOffset*3] + dg3 * ucompgrad[2][uOffset*3+1] +
							 dh3 * ucompgrad[2][uOffset*3+2]);

					uOffset = offset + this->memoffset[this->uMap[vIdx][3]];
					grad[0][uOffset] += 2.0 * this->LDetJ / detJ *
							(df4 * ucompgrad[0][uOffset*3] + dg4 * ucompgrad[0][uOffset*3+1] +
							 dh4 * ucompgrad[0][uOffset*3+2]);
					grad[1][uOffset] += 2.0 * this->LDetJ / detJ *
							(df4 * ucompgrad[1][uOffset*3] + dg4 * ucompgrad[1][uOffset*3+1] +
							 dh4 * ucompgrad[1][uOffset*3+2]);
					grad[2][uOffset] += 2.0 * this->LDetJ / detJ *
							(df4 * ucompgrad[2][uOffset*3] + dg4 * ucompgrad[2][uOffset*3+1] +
							 dh4 * ucompgrad[2][uOffset*3+2]);
				}
			}
		}
	}
}

} /* end namespace Jacobianu */


