
#ifndef MVECTOR_H  // the 'include guard'
#define MVECTOR_H  // see C++ Primer Sec. 2.9.2

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include <vector>

namespace mmath {

// Class that represents a mathematical vector
template <typename dtype>
class MVector {
   public:
    using size_t = std::size_t;
    // constructors
    MVector() {}
    explicit MVector(const int &n) : vector_(n) {}
    MVector(const int &n, const dtype &x) : vector_(n, x) {}
    MVector(const std::initializer_list<dtype> &list) : vector_(list) {}

    //    MVector &operator=(std::vector<double> vec){
    //        for (int i = 0; i < vec.size(); ++i) {
    //            v[i] = vec[i];
    //        }
    //        return *this;
    //    }

    // access element (lvalue) (see example sheet 5, q5.6)
    dtype &operator[](unsigned index) { return vector_[index]; }

    // access element (rvalue) (see example sheet 5, q5.7)
    dtype operator[](unsigned index) const { return vector_[index]; }

    bool operator==(const MVector &vec) {
        if (vector_.size() != vec.size()) {
            return false;
        }

        for (unsigned idx{0}; idx < vector_.size(); ++idx) {
            if (vector_[idx] != vec[idx]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief
     *
     * @return double
     */
    double LTwoNorm() const {  // 2-norm of vector v
        double result{0.};
        for (unsigned idx{0}; idx < size(); ++idx) {
            result += vector_[idx] * vector_[idx];
        }
        return sqrt(result);
    }

    /**
     * @brief
     *
     * @return double
     */
    double mean() const {
        double mean_value{0.};
        for (auto &ele : vector_) {
            mean_value += ele;
        }
        mean_value = mean_value / size();
        return mean_value;
    }

    size_t size() const { return vector_.size(); }  // number of elements

    void resize(const int &n) { vector_.resize(n); }

    // This is the threshold function of MVector
    // Input positive integer k and it sets all but the largest k elements
    // (by modulus) to zero
    void threshold(const int &k) {
        assert(k > 0);  // k must be positive integer
        // Create vector temp to store each element's modulus in v
        std::vector<dtype> temp(vector_);
        for (auto &i : temp) {
            i = fabs(i);
        }
        // Sort temp's elements in descending order
        std::sort(temp.begin(), temp.end(), std::greater<>());
        // Introduce count_zero to make sure that when there exists several
        // elements equal to the kst largest modulus, the function can still
        // function well.
        int count_zero = 0;
        for (auto &j : vector_) {  // Set element less than kst largest to 0
            if (fabs(j) >= temp[k - 1]) {
                continue;
            }
            j = 0;
            count_zero++;  // Count the number of element we set to 0
        }
        // If count_zero is incorrect, we set several elements equal to kst
        // largest to 0 from the end of vector
        if (count_zero >= size() - k) {
            return;
        }

        for (unsigned l = size() - 1; l >= 0; --l) {
            if (fabs(vector_[l]) == temp[k - 1]) {
                vector_[l] = 0;
                count_zero++;
            }
            if (count_zero == vector_.size() - k)
                break;  // When count_zero is correct, break.
        }

        return;
    }

   private:
    std::vector<dtype> vector_;
    size_t size_{0};
};

template <typename dtype>
MVector<dtype> &operator+(const MVector<dtype> &lhs,
                          const MVector<dtype> &rhs) {
    if (lhs.size() != rhs.size() && lhs.size() != 1 && lhs.size() != 1) {
        throw std::invalid_argument("Invalid size of two vectors!");
    }

    std::vector<dtype> ret;
    size_t size{std::max(lhs.size(), rhs.size())};

    bool l_flag{1 == lhs.size()};
    bool r_flag{1 == rhs.size()};

    switch (l_flag | r_flag << 1) {
        case 0:
            for (unsigned idx{0}; idx < size; ++idx) {
                ret.emplace_back(lhs[idx] + rhs[idx]);
            }
            break;
        case 1:
            for (unsigned idx{0}; idx < size; ++idx) {
                ret.emplace_back(lhs[0] + rhs[idx]);
            }
            break;

        default:
            break;
    }

    for (unsigned idx{0}; idx < size; ++idx) {
    }
    return MVector<dtype>(ret);
}

// MVector - MVector
template <typename dtype>
inline MVector<dtype> operator-(const MVector<dtype> &lhs,
                                const MVector<dtype> &rhs) {
    assert(lhs.size() == rhs.size());
    MVector r(lhs);
    for (int i = 0; i < lhs.size(); i++) {
        r[i] -= rhs[i];
    }
    return r;
}

// Hadamard product
template <typename dtype>
inline MVector<dtype> operator*(const MVector<dtype> &a,
                                const MVector<dtype> &b) {
    assert(a.size() == b.size());
    MVector r(a.size());
    for (int i = 0; i < a.size(); i++) {
        r[i] = a[i] * b[i];
    }
    return r;
}

// double * MVector
template <typename dtype>
inline MVector<dtype> operator*(double d, const MVector<dtype> &v) {
    MVector r(v);
    for (int i = 0; i < v.size(); i++) {
        r[i] *= d;
    }
    return r;
}

// MVector -= MVector
template <typename dtype>
inline MVector<dtype> operator+=(MVector<dtype> &v1, const MVector<dtype> &v) {
    assert(v1.size() == v.size());
    for (int i = 0; i < v1.size(); i++) {
        v1[i] += v[i];
    }
    return v1;
}

// MVector -= MVector
template <typename dtype>
inline MVector<dtype> operator-=(MVector<dtype> &v1, const MVector<dtype> &v) {
    assert(v1.size() == v.size());
    for (unsigned idx{0}; idx < v1.size(); ++idx) {
        v1[idx] -= v[idx];
    }
    return v1;
}

// Output function for MVector
template <typename dtype>
inline std::ostream &operator<<(std::ostream &os, const MVector<dtype> &rhs) {
    std::size_t n = rhs.size();
    os << "[";
    for (std::size_t i = 0; i < n; i++) {
        os << rhs[i];
        if (i != (n - 1)) os << ", ";
    }
    os << "]";
    return os;
}

// initialize_normal function initialises the vector to a random k-sparse vector
// where each non-zero component has zero mean and unit variance.
template <typename dtype>
MVector<dtype> initialize_normal(const int &n, const int &k) {
    assert(k < n);  // k must be less than n
    // Create norm distribution random numbers
    static std::default_random_engine e(time(0));
    static std::normal_distribution<double> norm(0.0, 1.0);
    std::vector<dtype> v(n);
    for (int i = 0; i < k; ++i) {
        v[i] = static_cast<dtype>(norm(e));
    }
    // Shuffle the vector to make the order of 0 element random
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
    MVector<dtype> v_new(n);  // Create a MVector variable to store the result
    for (int j = 0; j < n; ++j) {
        v_new[j] = v[j];
    }
    return v_new;
}
}  // namespace mmath

#endif