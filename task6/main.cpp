#include <iostream>
#include <functional>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


double random(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

double fi(double x, double *xi, int i, int len_xi) {
    double result = 1;
    for (int j = 0; j < len_xi; j++) {
        if (j != i) {
            result *= (x - xi[j]) / (xi[i] - xi[j]);
        }
    }
    return result;
}

double func(double x) {
    return sin(x);
}

void f_args(double **args, int len_1, int len_2, double a, double step) {
    double current = a;
    for (int i = 0; i < len_1; i++) {
        for (int j = 0; j < len_2; j++) {
            args[i][j] = current;
            if (j != len_2 - 1) {
                current += step;
            }
        }
    }
}

void f_r_args(double **args, double **random_args, int k, int len2, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            random_args[i][j] = random(args[i][0], args[i][len2 - 1]);
        }
    }
}


void f_A(double** A, double **random_args, double **args, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i * (n - 1) + j1) <= (i * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
                        int i_u = i * (n - 1) + j1;
                        int j_u = i * (n - 1) + j2;
                        A[j_u-i_u][i_u] +=
                                fi(random_args[i][q], args[i], j1, n) * fi(random_args[i][q], args[i], j2, n);
                    }
//                    A[i * (n - 1) + j2][i * (n - 1) + j1] = A[i * (n - 1) + j1][i * (n - 1) + j2];
                }
            }
        }
    }
}

void f_A_test(double** A, double **random_args, double **args, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i * (n - 1) + j1) <= (i * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
//                        int i_h = i * (n - 1) + j1;
//                        int j_h = i * (n - 1) + j2;
                        A[i * (n - 1) + j1][i * (n - 1) + j2] +=
                                fi(random_args[i][q], args[i], j1, n) * fi(random_args[i][q], args[i], j2, n);
                    }
                    A[i * (n - 1) + j2][i * (n - 1) + j1] = A[i * (n - 1) + j1][i * (n - 1) + j2];
                }
            }
        }
    }
}

void f_A_G(double** A,double** A_G,int n,int m){
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++){
                A_G[i][j] = A[i][j];
            }
        }
    for(int i=n;i<2*n-1;i++){
        for(int j=0;j<m;j++){
            A_G[i][j] = A_G[i-n+1][j];
        }
    }
}

void f_B(double* B, function<double(double)> y, double **random_args, double **args, int l, int m, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            for (int q = 0; q < l; q++) {
                B[i * (n - 1) + j] += (y(random_args[i][q]) * fi(random_args[i][q], args[i], j, n));
            }
        }
    }
}

void gauss(double **,int m,int n){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){

        }
    }
}

int main() {

    function<double(double)> y = func;

    double a = -3.14;
    double b = 3.14;

    int k = 4;
    int n = 3;
    int m = (n - 1) * k + 1;
    int l = 15;
    double h = (b - a) / (k * (n - 1));
    int le = 2*n-1;

    double **args = new double *[k];
    for (int i = 0; i < k; i++) {
        args[i] = new double[n];
    }

    double **random_args = new double *[k];
    for (int i = 0; i < k; i++) {
        random_args[i] = new double[l];
    }

    double **A = new double *[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[m];

    }
    double **A_G = new double *[2*n-1];
    for (int i = 0; i < 2*n-1; i++) {
        A_G[i] = new double[m];

    }

    double **A_test = new double *[m];
    for (int i = 0; i < m; i++) {
        A_test[i] = new double[m];

    }

    double B[m];
    for(int i=0;i<m;i++){
//        B[i] = 0.0;
    }

    double X[m];
    for(int i=0;i<m;i++){
        X[i] = 0.0;
    }


    MatrixXd A_E(m, m);
    A_E.setZero();

    VectorXd B_E(m);
    VectorXd X_E(m);
    B_E.setZero();

    f_args(args, k, n, a, h);

    f_r_args(args, random_args, k, n, l);

    f_A(A, random_args, args, l, k, n);
    f_A_G(A,A_G,n,m);
    f_A_test(A_test, random_args, args, l, k, n);
    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            cout << setw(10) << A_test[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            cout << setw(10) << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for(int i=0;i<2*n-1;i++){
        for(int j=0;j<m;j++){
            cout << setw(10) << A_G[i][j] << " ";
        }
        cout << endl;
    }
    f_B(B, y, random_args, args, l, m, k, n);
//    X = A.colPivHouseholderQr().solve(B);

    return 0;
}
