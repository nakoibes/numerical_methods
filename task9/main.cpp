#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


double p(double x) {
    return exp(x);
}

double q(double x) {
    return x * x;
}

double f(double x) {
    return sin(x);
}

void f_even_args(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

void formA(double **A, double *X, double h, int n) {
    for (int i = 0; i < n; i++) {
        A[i][2] = (2.0 + p(X[i + 1]) * h) / (2.0 * h * h);
        A[i][0] = (-4.0 - 2.0*h*h*q(X[i + 1])) / (2.0 * h * h);
        A[i][1] = (2.0 - p(X[i + 1]) * h) / (2.0 * h * h);;
    }
}

void prog(double **M, double *Y, double *F, int n) {
    double A[n];
    double B[n];
    A[0] = -M[0][1] / M[0][0];
    B[0] = F[0] / M[0][0];
    for (int i = 1; i < n - 1; i++) {
        double den = M[i][2] * A[i - 1] + M[i][0];
        A[i] = -M[i][1] / den;
        B[i] = (F[i] - M[i][2] * B[i - 1]) / den;
    }
    A[n - 1] = 0.0;
    B[n - 1] = (F[n - 1] - M[n - 1][2] * B[n - 2]) / (M[n - 1][2] * A[n - 2] + M[n - 1][0]);
    Y[n - 1] = B[n-1];
    for (int i = n - 2; i > -1; i--) {
        Y[i] = A[i] * Y[i + 1] + B[i];
    }
}


void f_diff(double *X, double *Y, double h, int n, double ya, double yb) {

    double F_p[n - 2];
    double Y_p[n - 2];

    double **A = new double *[n - 2];
    for (int i = 0; i < n - 2; i++) {
        A[i] = new double[3];
    }

    formA(A, X, h, n - 2);
    
    for (int i = 1; i < n - 2; i++) {
        F_p[i - 1] = f(X[i]);
    }

    F_p[0] -= A[0][2] * ya;
    F_p[n - 3] -= A[n-3][1] * yb;


    prog(A, Y_p, F_p, n - 2);

    for (int i = 1; i < n - 1; i++) {
        Y[i] = Y_p[i - 1];
    }

}

void write_f(double *Y1_ra, double *X, int n) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n; i++) {
        fout << X[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_ra[i] << " ";
    }


    fout.close();
}


int main() {
    double a = 0.0;
    double b = 1.0;
    double ya = 1.0;
    double yb = 1.0;
    int n1 = 107;
    double h1 = (b - a) / (n1 - 1);

    double X_h[n1];

    double Y1_ra_h[n1];
    Y1_ra_h[0] = ya;
    Y1_ra_h[n1 - 1] = yb;


    f_even_args(X_h, n1, a, b);

    f_diff(X_h, Y1_ra_h, h1, n1, ya, yb);

    write_f(Y1_ra_h, X_h, n1);


    return 0;
}
