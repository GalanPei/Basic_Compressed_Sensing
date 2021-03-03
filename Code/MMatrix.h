#ifndef MMATRIX_H // the 'include guard'
#define MMATRIX_H

#include <vector>
#include "MVector.h"
#include <ostream>
#include <iomanip>

// Class that represents a mathematical matrix
class MMatrix{
public:
    // constructors
    MMatrix() : nRows(0), nCols(0) {}
    MMatrix(int n, int m, double x = 0) : nRows(n), nCols(m), A(n * m, x) {}

    // set all matrix entries equal to a double
    MMatrix &operator=(double x){
        for (unsigned i = 0; i < nRows * nCols; i++) A[i] = x;
        return *this;
    }

    MMatrix &operator=(MVector v){
        for (int i = 0; i < nRows*nCols; ++i) {
            A[i] = v[i];
        }
        return *this;
    }

    // access element, indexed by (row, column) [rvalue]
    double operator()(int i, int j) const{
        return A[j + i*nCols];
    }

    // access element, indexed by (row, column) [lvalue]
    double &operator()(int i, int j){
        return A[j + i*nCols];
    }

    MMatrix &operator()(int i, char c){
        assert(c == ':');
        MMatrix B(1, nCols, 0.0);
        for (int j = 0; j < nCols; ++j) {
            B(0, j) = A[j + i*nCols];
        }
        return B;
    }

    MVector &operator()(char c, int j){
        assert(c == ':');
        MVector B(nRows, 0.0);
        for (int i = 0; i < nRows; ++i) {
            B[i] = A[j + i*nCols];
        }
        return B;
    }

    // Transpose of MMatrix
    MMatrix transpose() const{
        MMatrix B(nCols, nRows, 0.0);
        for (int i = 0; i < nCols; ++i) {
            for (int j = 0; j < nRows; ++j) {
                B(i, j) = A[i + j*nCols];
            }
        }
        return B;
    }

    // return the last column of matrix as a vector
    MVector last_column() const{
        MVector v(nRows, 0.0);
        for (int i = 0; i < nRows; ++i) {
            v[i] = A[nCols - 1 + i*nCols];
        }
        return v;
    }

    MMatrix alter_size(int new_row, int new_col) const{
        MMatrix B(new_row, new_col, 0.0);
        for (int i = 0; i < new_row; ++i) {
            for (int j = 0; j < new_col; ++j) {
                B(i, j) = A[j + i*nCols];
            }
        }
        return B;
    }

    void vector_evalue(int col, MVector& v){
        for(int i = 0; i < nRows; i++){
            A[col + i*nCols] = v[i];
        }
    }

    // size of matrix
    int Rows() const { return nRows; }
    int Cols() const { return nCols; }

private:
    unsigned int nRows, nCols;
    std::vector<double> A;
};

// Operator overload for "matrix(m*n) * matrix(n*p)"
inline MMatrix operator*(const MMatrix& lm, const MMatrix& rm){
    MMatrix temp(lm.Rows(), rm.Cols(), 0.0);
    if(lm.Cols() == rm.Rows()){ // Continue the calculation if the number lm's col = rm's row
        for(int i = 0; i < lm.Rows(); ++i){
            for(int j = 0; j < rm.Cols(); ++j){
                for(int k = 0; k < lm.Cols(); ++k){
                    temp(i, j) += lm(i, k)*rm(k, j);
                }
            }
        }
    }
    return temp;
}

// MMatrix * MVector
inline MVector operator*(const MMatrix &m, const MVector &v){
    assert(m.Cols() == v.size());
    MVector r(m.Rows());
    for (int i = 0; i < m.Rows(); i++){
        for (int j = 0; j < m.Cols(); j++){
            r[i] +=m (i, j)*v[j];
        }
    }
    return r;
}

// transpose(MMatrix) * MVector
MVector TransposeTimes(const MMatrix &m, const MVector &v){
    assert(m.Rows() == v.size());
    MVector r(m.Cols());
    for (int i = 0; i < m.Cols(); i++){
        for (int j=0; j<m.Rows(); j++){
            r[i] += m(j, i)*v[j];
        }
    }
    return r;
}

// MMatrix = MVector <outer product> MVector
// M = a <outer product> b
MMatrix OuterProduct(const MVector &a, const MVector &b){
    MMatrix m(a.size(), b.size());
    for (int i = 0; i < a.size(); i++){
        for (int j=0; j<b.size(); j++){
            m(i, j) = a[i]*b[j];
        }
    }
    return m;
}

// double * MMatrix
inline MMatrix operator*(double d, const MMatrix &m){
    MMatrix r(m);
    for (int i = 0; i < m.Rows(); i++){
        for (int j = 0; j < m.Cols(); j++){
            r(i,j) *= d;
        }
    }
    return r;
}

// MMatrix -= MMatrix
inline MMatrix operator-=(MMatrix &m1, const MMatrix &m2){
    assert (m1.Rows() == m2.Rows() && m1.Cols() == m2.Cols());
    for (int i = 0; i < m1.Rows(); i++) {
        for (int j = 0; j < m1.Cols(); j++) {
            m1(i, j) -= m2(i, j);
        }
    }
    return m1;
}

// Output function for MMatrix
inline std::ostream &operator<<(std::ostream &os, const MMatrix &A){
    int m = A.Rows(), n = A.Cols();
    os << "[";
    for(int i = 0; i < m; ++i) {
        if (i > 0){
            os << " ";
        }
        for(int j = 0; j < n; ++j) {
            os << std::setiosflags(std::ios::left) << std::setw(9) << A(i, j);
        }
        if (i < m - 1){
            os << std::endl;
        }
    }
    os << "]" << std::endl;
    return os;
}

MMatrix Transpose(const MVector& A){
    MMatrix AT(1, A.size(), 0.0);
    for (int i = 0; i < A.size(); ++i) {
        AT(0, i) = A[i];
    }
    return AT;
}

// When size of MMatrix A is 1*1, we can treat it as a double variable in division
inline double operator/(const double& x, const MMatrix& A){
    double temp = 0.0;
    if(A.Rows() == 1 && A.Cols() == 1){
        temp = x/A(0, 0);
    }
    return temp;
}

// Transformation of a n*1 MVector to a n*1 MMatrix
inline MMatrix vector_to_matrix(MVector& A){
    MMatrix temp(A.size(), 1, 0.0);
    for(int i = 0; i < A.size(); ++i){
        temp(i, 0) = A[i];
    }
    return temp;
}

// initialize_normal_matrix function sets each entry to a random value,
// normally distributed, with zero mean and variance 1/m.
// n, m are the row and column number of the MMatrix respectively.
MMatrix initialize_normal_matrix(int n, int m){
    // Create normal distributed random number
    static std::default_random_engine e;
    static std::normal_distribution<double> norm(0.0, 1.0/m);
    // Assign the number to an n*m MMatrix
    MMatrix A(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A(i, j) = norm(e);
        }
    }
    return A;
}

#endif