#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double f1(double x, double y1, double y2) {
    return y2;
}

double f2(double x, double y1, double y2) {
    return y2 * exp(x) + y1 * x * x + sin(x);
}

void f_even_args(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

void eil(double *Y1, double *Y2, double *X, int n, double h) {
    for (int i = 1; i < n; i++) {
        Y1[i] = Y1[i - 1] + h * f1(X[i-1], Y1[i-1], Y2[i-1]);
        Y2[i] = Y2[i - 1] + h * f2(X[i-1], Y1[i-1], Y2[i-1]);
    }
}

void eil_mod(double *Y1, double *Y2, double *X, int n, double h) {
    double Y1_[n];
    double Y2_[n];
    for (int i = 1; i < n; i++) {
        Y1_[i] = Y1[i - 1] + h * f1(X[i-1], Y1[i-1], Y2[i-1]);
        Y2_[i] = Y2[i - 1] + h * f2(X[i-1], Y1[i-1], Y2[i-1]);
        Y1[i] = Y1[i - 1] + h * (f1(X[i-1], Y1[i-1], Y2[i-1])+f1(X[i], Y1_[i], Y2_[i]))/2.0;
        Y2[i] = Y2[i - 1] + h * (f2(X[i-1], Y1[i-1], Y2[i-1])+f2(X[i], Y1_[i], Y2_[i]))/2.0;
    }
}

void write_f(double *x,double* y, int n) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n; i++) {
        fout << x[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << y[i] << " ";
    }
    fout.close();
}
int main() {
    double a = 0.0;
    double b = 1.0;
    int n = 101;
    double h = (b - a) / (n - 1);

    double X[n];

    double Y1[n];
    double Y2[n];

//    double **Y = new double *[n];
//    for (int i = 0; i < n; i++) {
//        Y[i] = new double[2];
//    }

    Y1[0] = 0;
    Y2[0] = 0;
    f_even_args(X,n,a,b);
    eil_mod(Y1,Y2,X,n,h);
    write_f(X,Y1,n);
    system("python3 vis.py");
    return 0;
}
