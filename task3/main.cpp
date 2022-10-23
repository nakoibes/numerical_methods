#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>

using namespace std;

double func(double x) {
    return sin(x);
}

double random(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

double fi(double x, double *xi, int i, int len_xi) {
    i = i % len_xi;
    double result = 1;
    for (int j = 0; j < len_xi; j++) {
        if (j != i) {
            result *= (x - xi[j]) / (xi[i] - xi[j]);
        }
    }
    return result;
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

void f_random_args(double **args, double *random_args, int k, int len2, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            random_args[i * l + j] = random(args[i][0], args[i][len2 - 1]);
        }
    }
}

void f_A(double **A, double *random_args, double **args, int l, int m, int k, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i <= j) {
                double value = 0;
                for (int q = 0; q < l * k; q++) {
                    value += (fi(random_args[q], args[i / n], i, n) * fi(random_args[q], args[j / n], j, n));
                }
                A[i][j] = value;
                A[j][i] = value;
            }
        }
    }
}

void f_B(double *B, function<double(double)> y, double *random_args, double **args, int l, int m, int k, int n) {
    for (int i = 0; i < m; i++) {
        double value = 0;
        for (int j = 0; j < l * k; j++) {
            value += (y(random_args[j]) * fi(random_args[j], args[i / n], i, n));
        }
        B[i] = value;
    }
}

int main() {
    srand((unsigned int) time(nullptr));

    function<double(double)> y = func;

    double a = -3;
    double b = 3;

    double delta = 0.0001;


    int k = 3;
    int n = 3;
    int m = (n - 1) * k + 1;
    int l = 5;
    double h = (b - a) / (k * (n - 1));


    double **args = new double *[k];
    for (int i = 0; i < k; i++) {
        args[i] = new double[n];
    }

//    double **random_args = new double *[k];
//    for (int i = 0; i < k; i++) {
//        random_args[i] = new double[l];
//    }

    double random_args[l * k];


    double **A = new double *[m];
    for (int i = 0; i < m; i++) {
        A[i] = new double[m];
    }

    double B[m];

    f_args(args, k, n, a, h);
    f_random_args(args, random_args, k, n, l);
    f_A(A, random_args, args, l, m, k, n);
    f_B(B, y, random_args, args, l, m, l, n);

    for (int i = 0; i < k; i++) {
        delete[] args[i];
    }
    delete[] args;


    for (int i = 0; i < m; i++) {
        delete[] A[i];
    }
    delete[] A;

//    for (int i = 0; i < k; i++) {
//        delete[] random_args[i];
//    }
//    delete[] random_args;

    return 0;
}
