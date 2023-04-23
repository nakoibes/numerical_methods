#include <iostream>
#include <cmath>
#include <fstream>

#define PI  3.14159265358979323846

using namespace std;

double f(double t, double x) {
    return t * sin(PI * x);
}

double fi(double x) {
    return -sin(PI * x);
}

double psi0(double t) {
    return t;
}

double psi1(double t) {
    return t;
}

void init_u(double **u, double h, double tau, int n, int m) {
    for (int i = 0; i < m; i++) {
        u[0][i] = fi(i * h);
    }
    for (int i = 0; i < n; i++) {
        u[i][0] = psi0(i * tau);
        u[i][m - 1] = psi1(i * tau);
    }
}

void exp_sch(double **u, double h, double tau, int n, int m) {
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m - 1; j++) {
            u[i][j] = tau * (u[i - 1][j + 1] - 2.0 * u[i - 1][j] + u[i - 1][j - 1]) / (h * h) +
                      f((i - 1) * tau, j * h) * tau + u[i - 1][j];
        }
    }
}

void write_d(double **u, int n, int m) {
    ofstream fout("dots.txt");
    fout << m << endl;
    fout << n << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fout << u[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();
}

int main() {
    int n = 1001;
    int m = 11;

    double a = 0.0;
    double b = 1.0;
    double t0 = 0.0;
    double T = 1.0;

    double h = (b - a) / (m - 1);
    double tau = (T - t0) / (n - 1);

    double **u = new double *[n];
    for (int i = 0; i < n; i++) {
        u[i] = new double[m];
    }

    init_u(u, h, tau, n, m);

//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < m; j++) {
//
//        }
//    }

    exp_sch(u, h, tau, n, m);
    write_d(u, n, m);


    return 0;
}
