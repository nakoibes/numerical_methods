#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

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
            u[i + 1][j] = X[j - 1];
//            cout << X[j-1] << " ";
        }
//        cout << endl;
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

void write_d(double **u, int n, int m, string filename) {
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

double calc_abs_err_1_M(double **M1, double **M2, int n, int m) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result += abs(M1[i][j] - M2[4*i][2 * j]);
        }
    }
    return result;
}

double calc_abs_err_2_M(double **M1, double **M2, int n, int m) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result += pow(abs(M1[i][j] - M2[4*i][2 * j]), 2);
        }
    }
    return sqrt(result);
}

double calc_abs_err_cheb_M(double **M1, double **M2, int n, int m) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (abs(M1[i][j] - M2[4*i][2 * j]) > result) {
                result = abs(M1[i][j] - M2[4*i][2 * j]);
            }
        }
    }
    return result;
}

double calc_rel_err_1_M(double **M1, double **M2, int n, int m) {
    double eps = 0.000000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            num += abs(M1[i][j] - M2[4*i][2 * j]);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            den += abs(M1[i][j]);
        }
    }
    if (den < eps) {
        den = 1.0;
    }
    return num / den;
}

double calc_rel_err_2_M(double **M1, double **M2, int n, int m) {
    double eps = 0.000000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            num += pow(abs(M1[i][j] - M2[4*i][2 * j]), 2);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            den += pow(abs(M1[i][j]), 2);
        }
    }
    if (den < eps) {
        den = 1.0;
    }
    return sqrt(num) / sqrt(den);
}

double calc_rel_err_cheb_M(double **M1, double **M2, int n, int m) {
    double eps = 0.000000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (abs(M1[i][j] - M2[4*i][2 * j]) > num) {
                num = abs(M1[i][j] - M2[4*i][2 * j]);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (abs(M1[i][j]) > den) {
                den = abs(M1[i][j]);
            }
        }
    }
    if (den < eps) {
        den = 1.0;
    }
    return num / den;
}

void write_e(double **u_e_h, double **u_e_h2, double **u_w_h, double **u_w_h2, int m1, int m2, int n, string filename) {
    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);
    fout << setw(15) << left << "Explicit scheme"
         << endl;

    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1_M(u_e_h, u_e_h2, n, m1) << endl
         << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2_M(u_e_h, u_e_h2, n, m1) << endl
         << setw(15) << left << "abs_err_cheb"
         << setw(15) << left << calc_abs_err_cheb_M(u_e_h, u_e_h2, n, m1) << endl
         << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1_M(u_e_h, u_e_h2, n, m1) << endl
         << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2_M(u_e_h, u_e_h2, n, m1) << endl
         << setw(15) << left << "rel_err_cheb"
         << setw(15) << left << calc_rel_err_cheb_M(u_e_h, u_e_h2, n, m1) << endl
         << setw(15) << left << "Wight scheme" << endl
         << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1_M(u_w_h, u_w_h2, n, m1) << endl
         << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2_M(u_w_h, u_w_h2, n, m1) << endl
         << setw(15) << left << "abs_err_cheb"
         << setw(15) << left << calc_abs_err_cheb_M(u_w_h, u_w_h2, n, m1) << endl
         << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1_M(u_w_h, u_w_h2, n, m1) << endl
         << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2_M(u_w_h, u_w_h2, n, m1) << endl
         << setw(15) << left << "rel_err_cheb"
         << setw(15) << left << calc_rel_err_cheb_M(u_w_h, u_w_h2, n, m1) << endl
         << endl;
}

int main() {
    int n1 = 1001;
    int m1 = 11;
    int m2 = 2 * m1 - 1;
    int n2 = 4*n1-3;

    double a = 0.0;
    double b = 1.0;
    double t0 = 0.0;
    double T = 1.0;

    double sig = 0.3;

    double h = (b - a) / (m1 - 1);
    double h2 = (b - a) / (m2 - 1);
    double tau = (T - t0) / (n1 - 1);
    double tau2 = (T - t0) / (n2 - 1);

    double **u_exp_h = new double *[n1];
    for (int i = 0; i < n1; i++) {
        u_exp_h[i] = new double[m1];
    }
    double **u_exp_h2 = new double *[n2];
    for (int i = 0; i < n2; i++) {
        u_exp_h2[i] = new double[m2];
    }

    double **u_wei_h = new double *[n1];
    for (int i = 0; i < n1; i++) {
        u_wei_h[i] = new double[m1];
    }
    double **u_wei_h2 = new double *[n2];
    for (int i = 0; i < n2; i++) {
        u_wei_h2[i] = new double[m2];
    }


    init_u(u_exp_h, h, tau, n1, m1);
    init_u(u_wei_h, h, tau, n1, m1);

    init_u(u_exp_h2, h2, tau2, n2, m2);
    init_u(u_wei_h2, h2, tau2, n2, m2);

    exp_sch(u_exp_h, h, tau, n1, m1);
    exp_wei(u_wei_h, h, tau, n1, m1, sig);

    exp_sch(u_exp_h2, h2, tau2, n2, m2);
    exp_wei(u_wei_h2, h2, tau2, n2, m2, sig);

    write_d(u_exp_h, n1, m1, "dots1.txt");
    write_d(u_wei_h, n1, m1, "dots2.txt");

    write_e(u_exp_h, u_exp_h2, u_wei_h, u_wei_h2, m1, m2, n1, "errs");

    system("python3 vis.py");

    return 0;
}
