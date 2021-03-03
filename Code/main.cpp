#include <iostream>
#include <fstream>
#include <iomanip>
#include "MVector.h"
#include "MMatrix.h"

// function of Steepest Descent Algorithm
std::pair<MMatrix, int> SDLS(const MMatrix& A, const MVector& b, MVector& x0, int MaxIt = 1000, double Tol = 1e-6){
    int n = x0.size(), iter;
    MMatrix X(n, MaxIt + 1, 0.0); // Use MMatrix X to store all vector solution from the algorithm
    MVector r = A.transpose()*(b - A*x0);
    MVector x = x0;
    X.vector_evalue(0, x0);
    for (iter = 1; iter < MaxIt; iter++){
        double alpha = (Transpose(r)*r)[0]/(Transpose(r)*A.transpose()*A*r)[0];
        x += alpha*r;
        r -= alpha*A.transpose()*A*r;
        for (int i = 0; i < n; ++i) {
            X(i, iter) = x[i]; // Store the vector at the i^th stage of SDLS
        }
        if (r.LTwoNorm() < Tol) break; // Check if solution is accurate enough
    }
    X = X.alter_size(n, iter + 1); // Delete 0 element
    return std::make_pair(X, iter - 1);
}

// function of  Normalised Iterative Thresholding Algorithm
// Note that NIHT doesn't always converge, the viable iter is used to judge whether the algorithm converges or not.
// If R.size() > 1, that is to say, iter > 1, the algorithm converges. Whereas if R.size() = 1, the algorithm does not
// converge.
std::pair<MVector, MVector> NIHT(const MMatrix& A, const MVector& b, int k, int MaxIt = 1e3, double Tol = 1e-6){
    int iter;
    // Initialize x, r
    // Use R to store the residual in each step
    MVector x = A.transpose()*b, r = A.transpose()*(b - A*x), R(MaxIt, 0.0);
    x.threshold(k);
    R[0] = r.LTwoNorm();
    for (iter = 0; iter < MaxIt; iter++){
        double alpha = (Transpose(r)*r)[0]/(Transpose(r)*A.transpose()*A*r)[0];
        x += alpha*A.transpose()*(b - A*x);
        x.threshold(k);
        r = A.transpose()*(b - A*x);
        R[iter + 1] = r.LTwoNorm();
        if(R[iter + 1] > R[iter]){ // If the residuals donot descend, the algorithm does not converge
            iter = -1; // Assign -1 to iter and break
            break;
        }
        if(r.LTwoNorm() < Tol) break; // Check if solution is accurate enough
    }
    // When row < 2k, check whether the solution is correct or not
    if(A.Rows() < 2*k && (A*x - b).LTwoNorm() > Tol) iter = -1;
    R.resize(iter + 2); // Delete 0 element
    return std::make_pair(R, x);
}

namespace Exercise1{
    void exercise_main(){
        std::ofstream data1_1, data1_2, data1_3;
        data1_1.open("data1_1.csv");
        data1_2.open("data1_2.csv");
        data1_3.open("data1_3.csv");
        MMatrix A_1(3, 2, 0.0), A_2(3, 2, 0.0), A_3(3, 2, 0.0);
        MVector b(3, 0.0), x0(2, 0);
        // Assign values to matrix and vector
        A_1 = {1, 2, 2, 1, -1, 0};
        A_2 = {1, 2, 2, 1, 1.8, -2};
        A_3 = {1, 2, 2, 1, 3, -3};
        b = {10, -1, 0};
        std::pair<MMatrix, int> sol_1 = SDLS(A_1, b, x0), sol_2 = SDLS(A_2, b, x0), sol_3 = SDLS(A_3, b, x0);
        for (int i = 0; i < sol_1.first.Rows(); ++i) {
            for (int j = 0; j < sol_1.first.Cols(); ++j) {
                data1_1 << sol_1.first(i, j);
                data1_1 << ((j == sol_1.first.Cols() - 1) ? '\n' : ',');
            }
        }
        for (int i = 0; i < sol_2.first.Rows(); ++i) {
            for (int j = 0; j < sol_2.first.Cols(); ++j) {
                data1_2 << sol_2.first(i, j);
                data1_2 << ((j == sol_2.first.Cols() - 1) ? '\n' : ',');
            }
        }
        for (int i = 0; i < sol_3.first.Rows(); ++i) {
            for (int j = 0; j < sol_3.first.Cols(); ++j) {
                data1_3 << sol_3.first(i, j);
                data1_3 << ((j == sol_3.first.Cols() - 1) ? '\n' : ',');
            }
        }
        data1_1.close();
        data1_2.close();
        data1_3.close();
        // Print the outcome
        std::cout << "The solution of the problem is x = " << sol_1.first.last_column() << std::endl
                  << "Iteration step = " << sol_1.second << std::endl;
        std::cout << "last row: [-1 0], iteration step = " << sol_1.second << std::endl
                  << "last row: [1.8 -2], iteration step = " << sol_2.second << std::endl
                  << "last row: [-2 -2], iteration step = " << sol_3.second << std::endl;
    }
}

namespace Exercise2_1{
    void exercise_main(){
        MMatrix A = initialize_normal_matrix(1000, 1);
        std::ofstream data2_1;
        data2_1.open("data2_1.csv");
        for (int i = 0; i < 1000; ++i) {
            data2_1 << A(i, 0) << '\n';
        }
        data2_1.close();
        for (int j = 0; j < 10; ++j) {
            std::cout << initialize_normal(50, 20) << std::endl;
        }
    }
}

namespace Exercise2_2{
    void exercise_main(){
        int n = 90, m = 100, k = 10;
        std::ofstream data2_2;
        data2_2.open("data2_2.csv");
        MVector x = initialize_normal(n, k);
        for (int j = 0; j < 30; ++j) {
            MMatrix A = initialize_normal_matrix(m, n);
            for(int i = 0; i < NIHT(A, A*x, k).first.size(); ++i){
                data2_2 << NIHT(A, A*x, k).first[i] << ',';
            }
            data2_2 << std::endl;
        }
        data2_2.close();
    }
}

namespace Exercise2_3{
    void exercise_main(){
        std::ofstream data2_3, data2_4;
        data2_3.open("data2_3.csv");
        data2_4.open("data2_4.csv");
        MVector test_1 = initialize_normal(100, 99);
        for (int i = 5; i < 99; ++i) {
            MVector temp = test_1;
            const clock_t begin_time = clock();
            temp.threshold(i);
            double seconds = double(clock() - begin_time)/CLOCKS_PER_SEC;
            data2_3 << std::setprecision(8) << seconds << ',';
        }
        data2_3.close();
        for (int j = 50; j < 300; j = j + 10) {
            MVector test_2 = initialize_normal(j, 40);
            const clock_t begin_time = clock();
            test_2.threshold(40);
            double seconds = double(clock() - begin_time)/CLOCKS_PER_SEC;
            data2_4 << std::setprecision(8) << seconds << ',';
        }
        data2_4.close();
    }
}

namespace Exercise3_1{
    void exercise_main(){
        std::ofstream data3_1;
        data3_1.open("data3_1.csv");
        int n = 200, T = 100, successful_try;
        std::vector<int> k = {20, 50};
        for(int sparsity : k){
            for(int rows = 4; rows < 200; rows += 5){
                successful_try = 0;
                MVector x = initialize_normal(n, sparsity);
                for (int j = 0; j < T; ++j) {
                    MMatrix A = initialize_normal_matrix(rows, n);
                    std::pair<MVector, MVector> sol = NIHT(A, A*x, sparsity);
                    if(sol.first.size() > 1) successful_try++;
                }
                data3_1 << rows << ',' << successful_try << std::endl;
            }
        }
        data3_1.close();
    }
}

namespace Exercise3_2{
    void exercise_main(){
        std::ofstream data3_2;
        data3_2.open("data3_2.csv");
        int n = 100, T = 100, successful_try;
        for(int k = 4; k < n; k += 5){
            MVector x = initialize_normal(n, k);
            for(int rows = 4; rows < n; rows += 5){
                successful_try = 0;
                for (int j = 0; j < T; ++j) {
                    MMatrix A = initialize_normal_matrix(rows, n);
                    std::pair<MVector, MVector> sol = NIHT(A, A*x, k);
                    if(sol.first.size() > 1) successful_try++;
                }
                data3_2 << successful_try << ',';
            }
            data3_2 << std::endl;
        }
        data3_2.close();

    }
}

int main() {
//    using namespace Exercise1;
//    using namespace Exercise2_1;
    using namespace Exercise2_2;
//    using namespace Exercise2_3;
//    using namespace Exercise3_1;
//    using namespace Exercise3_2;
    exercise_main();
    return 0;
}
