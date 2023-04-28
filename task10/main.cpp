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
    Y[n - 1] = B[n - 1];
    for (int i = n - 2; i > -1; i--) {
        Y[i] = A[i] * Y[i + 1] + B[i];
    }
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

void fill_M(double **M, double h, double tau, int n, double sig) {
    M[0][0] = 2.0 * sig / (h * h) + 1.0 / tau;
    M[0][1] = -sig / (h * h);
    M[0][2] = 0.0;
    M[n - 1][0] = 2.0 * sig / (h * h) + 1.0 / tau;
    M[n - 1][1] = 0.0;
    M[n - 1][2] = -sig / (h * h);
    for (int i = 1; i < n - 1; i++) {
        M[i][0] = 2.0 * sig / (h * h) + 1.0 / tau;
        M[i][1] = -sig / (h * h);
        M[i][2] = -sig / (h * h);
    }
}

void fill_F(double **u, double *F, double h, double tau, int n, double sig, int l) {
    for (int i = 1; i < l + 1; i++) {
        F[i - 1] = u[n][i] / tau + (1 - sig) * (u[n][i + 1] - 2.0 * u[n][i] + u[n][i - 1]) / (h * h) +
                   f(n * tau, i * h) * (1 - sig) + f((n + 1) * tau, i * h) * sig;
    }
    F[0] += u[n + 1][0] * sig / (h * h);
    F[l - 1] += u[n + 1][l + 1] * sig / (h * h);
}

void exp_wei(double **u, double h, double tau, int n, int m, double sig) {
    double F[m - 2];
    double X[m - 2];

    double **M = new double *[m - 2];
    for (int i = 0; i < m - 2; i++) {
        M[i] = new double[3];
    }
    fill_M(M, h, tau, m - 2, sig);
    for (int i = 0; i < n - 1; i++) {
        fill_F(u, F, h, tau, i, sig, m - 2);
        prog(M, X, F, m - 2);
        for (int j = 1; j < m - 1; j++) {
            u[i+1][j] = X[j-1];
//            cout << X[j-1] << " ";
        }
        cout << endl;
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

void write_d(double **u, int n, int m,string filename) {
    ofstream fout(filename);
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

    double sig = 0.3;

    double h = (b - a) / (m - 1);
    double tau = (T - t0) / (n - 1);

    double **u_exp = new double *[n];
    for (int i = 0; i < n; i++) {
        u_exp[i] = new double[m];
    }

    double **u_wei = new double *[n];
    for (int i = 0; i < n; i++) {
        u_wei[i] = new double[m];
    }


    init_u(u_exp, h, tau, n, m);
    init_u(u_wei, h, tau, n, m);

    exp_sch(u_exp, h, tau, n, m);
    exp_wei(u_wei, h, tau, n, m, sig);


    write_d(u_exp, n, m,"dots1.txt");
    write_d(u_wei, n, m,"dots2.txt");


    return 0;
}
