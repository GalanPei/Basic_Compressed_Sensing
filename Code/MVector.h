
#ifndef MVECTOR_H // the 'include guard'
#define MVECTOR_H // see C++ Primer Sec. 2.9.2

#include <vector>
#include <ostream>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <random>

// Class that represents a mathematical vector
class MVector
{
public:
    // constructors
    MVector() {}
    explicit MVector(int n) : v(n) {}
    MVector(int n, double x) : v(n, x) {}
    MVector(std::initializer_list<double> l) : v(l) {}

//    MVector &operator=(std::vector<double> vec){
//        for (int i = 0; i < vec.size(); ++i) {
//            v[i] = vec[i];
//        }
//        return *this;
//    }

    // access element (lvalue) (see example sheet 5, q5.6)
    double &operator[](int index){
        return v[index];
    }

    // access element (rvalue) (see example sheet 5, q5.7)
    double operator[](int index) const {
        return v[index];
    }

    bool operator==(MVector& vec){
        if(v.size() != vec.size()) return false;
        for (int i = 0; i < v.size(); ++i) {
            if(v[i] != vec[i]) return false;
        }
        return true;
    }

    double LTwoNorm() const{ // 2-norm of vector v
        double result = 0.0;
        std::size_t s = size();
        for (int i = 0; i < s; i++){
            result +=  v[i]*v[i];
        }
        return sqrt(result);
    }

    double mean() const{
        double mean_value = 0;
        for (double i : v) {
            mean_value += i;
        }
        mean_value = mean_value/v.size();
        return mean_value;
    }

    int size() const {
        return v.size();
    } // number of elements

    void resize(int n){
        v.resize(n);
    }

    // This is the threshold function of MVector
    // Input positive integer k and it sets all but the largest k elements
    // (by modulus) to zero
    void threshold(int k){
        assert(k > 0); // k must be positive integer
        // Create vector temp to store each element's modulus in v
        std::vector<double> temp = v;
        for(double & i : temp){
            i = fabs(i);
        }
        // Sort temp's elements in descending order
        std::sort(temp.begin(), temp.end(), std::greater<>());
        // Introduce count_zero to make sure that when there exists several elements
        // equal to the kst largest modulus, the function can still function well.
        int count_zero = 0;
        for(double & j : v){ // Set element less than kst largest to 0
            if(fabs(j) < temp[k - 1]){
                j = 0;
                count_zero++; // Count the number of element we set to 0
            }
        }
        // If count_zero is incorrect, we set several elements equal to kst largest to 0
        // from the end of vector
        if(count_zero < v.size() - k){
            for (unsigned long l = v.size() - 1; l >= 0; --l) {
                if(fabs(v[l]) == temp[k - 1]){
                    v[l] = 0;
                    count_zero++;
                }
                if(count_zero == v.size() - k) break; // When count_zero is correct, break.
            }
        }
    }

private:
    std::vector<double> v;
};

// MVector + MVector
inline MVector operator+(const MVector &lhs, const MVector &rhs){
    assert(lhs.size() == rhs.size());
    MVector r(lhs);
    for (int i=0; i<lhs.size(); i++) {
        r[i] += rhs[i];
    }
    return r;
}

// MVector - MVector
inline MVector operator-(const MVector &lhs, const MVector &rhs){
    assert(lhs.size() == rhs.size());
    MVector r(lhs);
    for (int i=0; i<lhs.size(); i++) {
        r[i] -= rhs[i];
    }
    return r;
}

// Hadamard product
inline MVector operator*(const MVector &a, const MVector &b){
    assert(a.size() == b.size());
    MVector r(a.size());
    for (int i=0; i<a.size(); i++) {
        r[i] = a[i] * b[i];
    }
    return r;
}

// double * MVector
inline MVector operator*(double d, const MVector &v){
    MVector r(v);
    for (int i=0; i<v.size(); i++) {
        r[i] *= d;
    }
    return r;
}


// MVector -= MVector
inline MVector operator+=(MVector &v1, const MVector &v){
    assert(v1.size()==v.size());
    for (int i=0; i<v1.size(); i++) {
        v1[i] += v[i];
    }
    return v1;
}

// MVector -= MVector
inline MVector operator-=(MVector &v1, const MVector &v){
    assert(v1.size()==v.size());
    for (int i=0; i<v1.size(); i++) {
        v1[i] -= v[i];
    }
    return v1;
}

// Output function for MVector
inline std::ostream &operator<<(std::ostream &os, const MVector &rhs){
    std::size_t n = rhs.size();
    os << "(";
    for (std::size_t i=0; i<n; i++){
        os << rhs[i];
        if (i!=(n-1)) os << ", ";
    }
    os << ")";
    return os;
}

// initialize_normal function initialises the vector to a random k-sparse vector
// where each non-zero component has zero mean and unit variance.
MVector initialize_normal(int n, int k){
    assert(k < n); // k must be less than n
    // Create norm distribution random numbers
    static std::default_random_engine e(time(0));
    static std::normal_distribution<double> norm(0.0, 1.0);
    std::vector<double> v(n, 0.0);
    for (int i = 0; i < k; ++i) {
        v[i] = norm(e);
    }
    // Shuffle the vector to make the order of 0 element random
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle (v.begin(), v.end(), std::default_random_engine (seed));
    MVector v_new(n, 0.0); // Create a MVector variable to store the result
    for (int j = 0; j < n; ++j) {
        v_new[j] = v[j];
    }
    return v_new;
}

#endif