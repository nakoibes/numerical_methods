#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

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
        Y1[i] = Y1[i - 1] + h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        Y2[i] = Y2[i - 1] + h * f2(X[i - 1], Y1[i - 1], Y2[i - 1]);
    }
}

void eil_mod(double *Y1, double *Y2, double *X, int n, double h) {
    double Y1_[n];
    double Y2_[n];
    for (int i = 1; i < n; i++) {
        Y1_[i] = Y1[i - 1] + h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        Y2_[i] = Y2[i - 1] + h * f2(X[i - 1], Y1[i - 1], Y2[i - 1]);
        Y1[i] = Y1[i - 1] + h * (f1(X[i - 1], Y1[i - 1], Y2[i - 1]) + f1(X[i], Y1_[i], Y2_[i])) / 2.0;
        Y2[i] = Y2[i - 1] + h * (f2(X[i - 1], Y1[i - 1], Y2[i - 1]) + f2(X[i], Y1_[i], Y2_[i])) / 2.0;
    }
}

void write_f(double *Y1_e, double *Y1_e_m, double *Y1_rk1, double *Y1_rk4, double *Y1_a,
             double *Y2_e, double *Y2_e_m, double *Y2_rk1, double *Y2_rk4, double *Y2_a, double *X, int n) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n; i++) {
        fout << X[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_e[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_e_m[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_rk1[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_rk4[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_a[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y2_e[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y2_e_m[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y2_rk1[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y2_rk4[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y2_a[i] << " ";
    }

    fout.close();
}

void ru_ku2(double *Y1, double *Y2, double *X, int n, double h, double al) {
    for (int i = 1; i < n; i++) {
        double f1_ = f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double f2_ = f2(X[i - 1], Y1[i - 1], Y2[i - 1]);
        Y1[i] = Y1[i - 1] + h * ((1 - al) * f1_ + al * f1(X[i - 1] + h / (2 * al), Y1[i - 1] + h / (2 * al) * f1_,
                                                          Y2[i - 1] + h / (2 * al) * f1_));
        Y2[i] = Y2[i - 1] + h * ((1 - al) * f2_ + al * f2(X[i - 1] + h / (2 * al), Y1[i - 1] + h / (2 * al) * f2_,
                                                          Y2[i - 1] + h / (2 * al) * f2_));
    }
}

void ru_ku4(double *Y1, double *Y2, double *X, int n, double h) {
    for (int i = 1; i < n; i++) {
        double k1 = h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double k2 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k1 / 2, Y2[i - 1] + k1 / 2);
        double k3 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k2 / 2, Y2[i - 1] + k2 / 2);
        double k4 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k3, Y2[i - 1] + k3);
        Y1[i] = Y1[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        k1 = h * f2(X[i - 1], Y1[i - 1], Y2[i - 1]);
        k2 = h * f2(X[i - 1] + h / 2, Y1[i - 1] + k1 / 2, Y2[i - 1] + k1 / 2);
        k3 = h * f2(X[i - 1] + h / 2, Y1[i - 1] + k2 / 2, Y2[i - 1] + k2 / 2);
        k4 = h * f2(X[i - 1] + h / 2, Y1[i - 1] + k3, Y2[i - 1] + k3);
        Y2[i] = Y2[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }
}

void adams(double *Y1, double *Y2, double *X, int n, double h) {
    for (int i = 1; i < 4; i++) {
        double k1 = h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double k2 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k1 / 2, Y2[i - 1] + k1 / 2);
        double k3 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k2 / 2, Y2[i - 1] + k2 / 2);
        double k4 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k3, Y2[i - 1] + k3);
        Y1[i] = Y1[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        k1 = h * f2(X[i - 1], Y1[i - 1], Y2[i - 1]);
        k2 = h * f2(X[i - 1] + h / 2, Y1[i - 1] + k1 / 2, Y2[i - 1] + k1 / 2);
        k3 = h * f2(X[i - 1] + h / 2, Y1[i - 1] + k2 / 2, Y2[i - 1] + k2 / 2);
        k4 = h * f2(X[i - 1] + h / 2, Y1[i - 1] + k3, Y2[i - 1] + k3);
        Y2[i] = Y2[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }
    for (int i = 4; i < n; i++) {
        double k1 = h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double k2 = h * f1(X[i - 2], Y1[i - 2], Y2[i - 2]);
        double k3 = h * f1(X[i - 3], Y1[i - 3], Y2[i - 3]);
        Y1[i] = Y1[i - 1] + (23 * k1 - 16 * k2 + 5 * k3) / 12.0;
        k1 = h * f2(X[i - 1], Y1[i - 1], Y2[i - 1]);
        k2 = h * f2(X[i - 2], Y1[i - 2], Y2[i - 2]);
        k3 = h * f2(X[i - 3], Y1[i - 3], Y2[i - 3]);
        Y2[i] = Y2[i - 1] + (23 * k1 - 16 * k2 + 5 * k3) / 12.0;
    }
}

double calc_abs_err_1(double *func, double *inter, int n_err) {
    double result = 0.0;
    for (int i = 0; i < n_err; i++) {
        result += abs(func[i] - inter[i]);
    }
    return result;
}

double calc_abs_err_2(double *func, double *inter, int n_err) {
    double result = 0.0;
    for (int i = 0; i < n_err; i++) {
        result += pow(abs(func[i] - inter[i]), 2);
    }
    return sqrt(result);
}

double calc_abs_err_cheb(double *func, double *inter, int n_err) {
    double result = 0.0;
    for (int i = 0; i < n_err; i++) {
        if (abs(func[i] - inter[i]) > result) {
            result = abs(func[i] - inter[i]);
        }
    }
    return result;
}

double calc_rel_err_1(double *func, double *inter, int n_err) {
    double eps = 0.000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n_err; i++) {
        num += (abs(func[i] - inter[i]));
    }
    for (int i = 0; i < n_err; i++) {
        den += abs(func[i]);
    }
    if (abs(den - 0.0) < eps) {
        den = 1.0;
    }
    return num / den;
}

double calc_rel_err_2(double *func, double *inter, int n_err) {
    double eps = 0.000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n_err; i++) {
        num += (pow(abs(func[i] - inter[i]), 2));
    }
    for (int i = 0; i < n_err; i++) {
        den += pow(func[i], 2);
    }
    if (abs(den - 0.0) < eps) {
        den = 1.0;
    }
    return sqrt(num) / sqrt(den);
}

double calc_rel_err_cheb(double *func, double *inter, int n_err) {
    double eps = 0.000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n_err; i++) {
        if (abs(func[i] - inter[i]) > num) {
            num = (abs(func[i] - inter[i]));
        }
    }
    for (int i = 0; i < n_err; i++) {
        if (abs(func[i]) > den) {
            den = (abs(func[i]));
        }
    }
    if (abs(den - 0.0) < eps) {
        den = 1.0;
    }
    return num / den;
}

void sh(double *Y, double *Y_c, int n) {
    for (int i = 0; i < n; i++) {
        Y_c[i] = Y[2 * i];
    }
}

void write_errs1(double *Y1_e_h, double *Y2_e_h, double *Y1_e_m_h, double *Y2_e_m_h,
                 double *Y1_rk2_h, double *Y2_rk2_h, double *Y1_rk4_h, double *Y2_rk4_h,
                 double *Y1_a_h, double *Y2_a_h, double *Y1_e_h2, double *Y2_e_h2,
                 double *Y1_e_m_h2, double *Y2_e_m_h2, double *Y1_rk2_h2, double *Y2_rk2_h2,
                 double *Y1_rk4_h2, double *Y2_rk4_h2, double *Y1_a_h2, double *Y2_a_h2,
                 double *Y1_e_2h, double *Y2_e_2h, double *Y1_e_m_2h, double *Y2_e_m_2h,
                 double *Y1_rk2_2h, double *Y2_rk2_2h, double *Y1_rk4_2h, double *Y2_rk4_2h,
                 double *Y1_a_2h, double *Y2_a_2h, int n1, int n2, int n3, string filename) {

    double Y1_e_h2_s[n1];
    double Y1_e_m_h2_s[n1];
    double Y1_rk2_h2_s[n1];
    double Y1_rk4_h2_s[n1];
    double Y1_a_h2_s[n1];
    sh(Y1_e_h2, Y1_e_h2_s, n1);
    sh(Y1_e_m_h2, Y1_e_m_h2_s, n1);
    sh(Y1_rk2_h2, Y1_rk2_h2_s, n1);
    sh(Y1_rk4_h2, Y1_rk4_h2_s, n1);
    sh(Y1_a_h2, Y1_a_h2_s, n1);

    double Y2_e_h2_s[n1];
    double Y2_e_m_h2_s[n1];
    double Y2_rk2_h2_s[n1];
    double Y2_rk4_h2_s[n1];
    double Y2_a_h2_s[n1];
    sh(Y2_e_h2, Y2_e_h2_s, n1);
    sh(Y2_e_m_h2, Y2_e_m_h2_s, n1);
    sh(Y2_rk2_h2, Y2_rk2_h2_s, n1);
    sh(Y2_rk4_h2, Y2_rk4_h2_s, n1);
    sh(Y2_a_h2, Y2_a_h2_s, n1);

    double Y1_e_h_s[n1];
    double Y1_e_m_h_s[n1];
    double Y1_rk2_h_s[n1];
    double Y1_rk4_h_s[n1];
    double Y1_a_h_s[n1];
    sh(Y1_e_h, Y1_e_h_s, n2);
    sh(Y1_e_m_h, Y1_e_m_h_s, n2);
    sh(Y1_rk2_h, Y1_rk2_h_s, n2);
    sh(Y1_rk4_h, Y1_rk4_h_s, n2);
    sh(Y1_a_h, Y1_a_h_s, n2);

    double Y2_e_h_s[n1];
    double Y2_e_m_h_s[n1];
    double Y2_rk2_h_s[n1];
    double Y2_rk4_h_s[n1];
    double Y2_a_h_s[n1];
    sh(Y2_e_h, Y2_e_h_s, n2);
    sh(Y2_e_m_h, Y2_e_m_h_s, n2);
    sh(Y2_rk2_h, Y2_rk2_h_s, n2);
    sh(Y2_rk4_h, Y2_rk4_h_s, n2);
    sh(Y2_a_h, Y2_a_h_s, n2);

    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);

    fout << setw(15) << left << "error"
         << setw(15) << left << "y step h;h/2"
         << setw(15) << left << "y step h;2h"
         << setw(15) << left << "y' step h;h/2"
         << setw(15) << left << "y' step h;2h"
         << endl;

    fout << "eiler" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1(Y1_e_h, Y1_e_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y1_e_2h, Y1_e_h_s, n2)
         << setw(15) << left << calc_abs_err_1(Y2_e_h, Y2_e_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y2_e_2h, Y2_e_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y1_e_h, Y1_e_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y1_e_2h, Y1_e_h_s, n2)
         << setw(15) << left << calc_abs_err_2(Y2_e_h, Y2_e_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y2_e_2h, Y2_e_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y1_e_h, Y1_e_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y1_e_2h, Y1_e_h_s, n2)
         << setw(15) << left << calc_abs_err_cheb(Y2_e_h, Y2_e_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y2_e_2h, Y2_e_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y1_e_h, Y1_e_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y1_e_2h, Y1_e_h_s, n2)
         << setw(15) << left << calc_rel_err_1(Y2_e_h, Y2_e_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y2_e_2h, Y2_e_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y1_e_h, Y1_e_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y1_e_2h, Y1_e_h_s, n2)
         << setw(15) << left << calc_rel_err_2(Y2_e_h, Y2_e_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y2_e_2h, Y2_e_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y1_e_h, Y1_e_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y1_e_2h, Y1_e_h_s, n2)
         << setw(15) << left << calc_rel_err_cheb(Y2_e_h, Y2_e_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y2_e_2h, Y2_e_h_s, n2)
         << endl;

    fout << "eiler_mod" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << calc_abs_err_1(Y1_e_m_h, Y1_e_m_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y1_e_m_2h, Y1_e_m_h_s, n2)
         << setw(15) << left << calc_abs_err_1(Y2_e_m_h, Y2_e_m_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y2_e_m_2h, Y2_e_m_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y1_e_m_h, Y1_e_m_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y1_e_m_2h, Y1_e_m_h_s, n2)
         << setw(15) << left << calc_abs_err_2(Y2_e_m_h, Y2_e_m_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y2_e_m_2h, Y2_e_m_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y1_e_m_h, Y1_e_m_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y1_e_m_2h, Y1_e_m_h_s, n2)
         << setw(15) << left << calc_abs_err_cheb(Y2_e_m_h, Y2_e_m_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y2_e_m_2h, Y2_e_m_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y1_e_m_h, Y1_e_m_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y1_e_m_2h, Y1_e_m_h_s, n2)
         << setw(15) << left << calc_rel_err_1(Y2_e_m_h, Y2_e_m_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y2_e_m_2h, Y2_e_m_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y1_e_m_h, Y1_e_m_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y1_e_m_2h, Y1_e_m_h_s, n2)
         << setw(15) << left << calc_rel_err_2(Y2_e_m_h, Y2_e_m_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y2_e_m_2h, Y2_e_m_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y1_e_m_h, Y1_e_m_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y1_e_m_2h, Y1_e_m_h_s, n2)
         << setw(15) << left << calc_rel_err_cheb(Y2_e_m_h, Y2_e_m_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y2_e_m_2h, Y2_e_m_h_s, n2)
         << endl;

    fout << "runge 2" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << calc_abs_err_1(Y1_rk2_h, Y1_rk2_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y1_rk2_2h, Y1_rk2_h_s, n2)
         << setw(15) << left << calc_abs_err_1(Y2_rk2_h, Y2_rk2_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y2_rk2_2h, Y2_rk2_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y1_rk2_h, Y1_rk2_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y1_rk2_2h, Y1_rk2_h_s, n2)
         << setw(15) << left << calc_abs_err_2(Y2_rk2_h, Y2_rk2_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y2_rk2_2h, Y2_rk2_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y1_rk2_h, Y1_rk2_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y1_rk2_2h, Y1_rk2_h_s, n2)
         << setw(15) << left << calc_abs_err_cheb(Y2_rk2_h, Y2_rk2_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y2_rk2_2h, Y2_rk2_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y1_rk2_h, Y1_rk2_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y1_rk2_2h, Y1_rk2_h_s, n2)
         << setw(15) << left << calc_rel_err_1(Y2_rk2_h, Y2_rk2_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y2_rk2_2h, Y2_rk2_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y1_rk2_h, Y1_rk2_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y1_rk2_2h, Y1_rk2_h_s, n2)
         << setw(15) << left << calc_rel_err_2(Y2_rk2_h, Y2_rk2_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y2_rk2_2h, Y2_rk2_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y1_rk2_h, Y1_rk2_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y1_rk2_2h, Y1_rk2_h_s, n2)
         << setw(15) << left << calc_rel_err_cheb(Y2_rk2_h, Y2_rk2_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y2_rk2_2h, Y2_rk2_h_s, n2)
         << endl;

    fout << "runge 4" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << calc_abs_err_1(Y1_rk4_h, Y1_rk4_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y1_rk4_2h, Y1_rk4_h_s, n2)
         << setw(15) << left << calc_abs_err_1(Y2_rk4_h, Y2_rk4_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y2_rk4_2h, Y2_rk4_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y1_rk4_h, Y1_rk4_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y1_rk4_2h, Y1_rk4_h_s, n2)
         << setw(15) << left << calc_abs_err_2(Y2_rk4_h, Y2_rk4_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y2_rk4_2h, Y2_rk4_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y1_rk4_h, Y1_rk4_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y1_rk4_2h, Y1_rk4_h_s, n2)
         << setw(15) << left << calc_abs_err_cheb(Y2_rk4_h, Y2_rk4_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y2_rk4_2h, Y2_rk4_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y1_rk4_h, Y1_rk4_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y1_rk4_2h, Y1_rk4_h_s, n2)
         << setw(15) << left << calc_rel_err_1(Y2_rk4_h, Y2_rk4_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y2_rk4_2h, Y2_rk4_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y1_rk4_h, Y1_rk4_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y1_rk4_2h, Y1_rk4_h_s, n2)
         << setw(15) << left << calc_rel_err_2(Y2_rk4_h, Y2_rk4_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y2_rk4_2h, Y2_rk4_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y1_rk4_h, Y1_rk4_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y1_rk4_2h, Y1_rk4_h_s, n2)
         << setw(15) << left << calc_rel_err_cheb(Y2_rk4_h, Y2_rk4_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y2_rk4_2h, Y2_rk4_h_s, n2)
         << endl;

    fout << "adams" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << calc_abs_err_1(Y1_a_h, Y1_a_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y1_a_2h, Y1_a_h_s, n2)
         << setw(15) << left << calc_abs_err_1(Y2_a_h, Y2_a_h2_s, n1)
         << setw(15) << left << calc_abs_err_1(Y2_a_2h, Y2_a_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y1_a_h, Y1_a_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y1_a_2h, Y1_a_h_s, n2)
         << setw(15) << left << calc_abs_err_2(Y2_a_h, Y2_a_h2_s, n1)
         << setw(15) << left << calc_abs_err_2(Y2_a_2h, Y2_a_h_s, n2)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y1_a_h, Y1_a_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y1_a_2h, Y1_a_h_s, n2)
         << setw(15) << left << calc_abs_err_cheb(Y2_a_h, Y2_a_h2_s, n1)
         << setw(15) << left << calc_abs_err_cheb(Y2_a_2h, Y2_a_h_s, n2)

         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y1_a_h, Y1_a_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y1_a_2h, Y1_a_h_s, n2)
         << setw(15) << left << calc_rel_err_1(Y2_a_h, Y2_a_h2_s, n1)
         << setw(15) << left << calc_rel_err_1(Y2_a_2h, Y2_a_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y1_a_h, Y1_a_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y1_a_2h, Y1_a_h_s, n2)
         << setw(15) << left << calc_rel_err_2(Y2_a_h, Y2_a_h2_s, n1)
         << setw(15) << left << calc_rel_err_2(Y2_a_2h, Y2_a_h_s, n2)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y1_a_h, Y1_a_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y1_a_2h, Y1_a_h_s, n2)
         << setw(15) << left << calc_rel_err_cheb(Y2_a_h, Y2_a_h2_s, n1)
         << setw(15) << left << calc_rel_err_cheb(Y2_a_2h, Y2_a_h_s, n2)
         << endl;

    fout.close();
}

int main() {
    double a = 0.0;
    double b = 1.0;
    int n1 = 101;
    int n2 = (n1 - 1) / 2 + 1;
    int n3 = (n1 - 1) * 2 + 1;
    double h1 = (b - a) / (n1 - 1);
    double h2 = (b - a) / (n2 - 1);
    double h3 = (b - a) / (n3 - 1);

    double X_h[n1];
    double X_h2[n3];
    double X_2h[n2];

    double Y1_e_h[n1];
    double Y2_e_h[n1];
    Y1_e_h[0] = 0;
    Y2_e_h[0] = 0;

    double Y1_e_m_h[n1];
    double Y2_e_m_h[n1];
    Y1_e_m_h[0] = 0;
    Y2_e_m_h[0] = 0;

    double Y1_rk2_h[n1];
    double Y2_rk2_h[n1];
    Y1_rk2_h[0] = 0;
    Y2_rk2_h[0] = 0;

    double Y1_rk4_h[n1];
    double Y2_rk4_h[n1];
    Y1_rk4_h[0] = 0;
    Y2_rk4_h[0] = 0;

    double Y1_a_h[n1];
    double Y2_a_h[n1];
    Y1_a_h[0] = 0;
    Y2_a_h[0] = 0;


    double Y1_e_h2[n3];
    double Y2_e_h2[n3];
    Y1_e_h2[0] = 0;
    Y2_e_h2[0] = 0;

    double Y1_e_m_h2[n3];
    double Y2_e_m_h2[n3];
    Y1_e_m_h2[0] = 0;
    Y2_e_m_h2[0] = 0;

    double Y1_rk2_h2[n3];
    double Y2_rk2_h2[n3];
    Y1_rk2_h2[0] = 0;
    Y2_rk2_h2[0] = 0;

    double Y1_rk4_h2[n3];
    double Y2_rk4_h2[n3];
    Y1_rk4_h2[0] = 0;
    Y2_rk4_h2[0] = 0;

    double Y1_a_h2[n3];
    double Y2_a_h2[n3];
    Y1_a_h2[0] = 0;
    Y2_a_h2[0] = 0;

    double Y1_e_2h[n2];
    double Y2_e_2h[n2];
    Y1_e_2h[0] = 0;
    Y2_e_2h[0] = 0;

    double Y1_e_m_2h[n2];
    double Y2_e_m_2h[n2];
    Y1_e_m_2h[0] = 0;
    Y2_e_m_2h[0] = 0;

    double Y1_rk2_2h[n2];
    double Y2_rk2_2h[n2];
    Y1_rk2_2h[0] = 0;
    Y2_rk2_2h[0] = 0;

    double Y1_rk4_2h[n2];
    double Y2_rk4_2h[n2];
    Y1_rk4_2h[0] = 0;
    Y2_rk4_2h[0] = 0;

    double Y1_a_2h[n2];
    double Y2_a_2h[n2];
    Y1_a_2h[0] = 0;
    Y2_a_2h[0] = 0;

    f_even_args(X_h, n1, a, b);
    f_even_args(X_h2, n3, a, b);
    f_even_args(X_2h, n2, a, b);

    eil(Y1_e_h, Y2_e_h, X_h, n1, h1);
    eil_mod(Y1_e_m_h, Y2_e_m_h, X_h, n1, h1);
    ru_ku2(Y1_rk2_h, Y2_rk2_h, X_h, n1, h1, 0.7);
    ru_ku4(Y1_rk4_h, Y2_rk4_h, X_h, n1, h1);
    adams(Y1_a_h, Y2_a_h, X_h, n1, h1);

    eil(Y1_e_h2, Y2_e_h2, X_h2, n3, h3);
    eil_mod(Y1_e_m_h2, Y2_e_m_h2, X_h2, n3, h3);
    ru_ku2(Y1_rk2_h2, Y2_rk2_h2, X_h2, n3, h3, 0.7);
    ru_ku4(Y1_rk4_h2, Y2_rk4_h2, X_h2, n3, h3);
    adams(Y1_a_h2, Y2_a_h2, X_h2, n3, h3);

    eil(Y1_e_2h, Y2_e_2h, X_2h, n2, h2);
    eil_mod(Y1_e_m_2h, Y2_e_m_2h, X_2h, n2, h2);
    ru_ku2(Y1_rk2_2h, Y2_rk2_2h, X_2h, n2, h2, 0.7);
    ru_ku4(Y1_rk4_2h, Y2_rk4_2h, X_2h, n2, h2);
    adams(Y1_a_2h, Y2_a_2h, X_2h, n2, h2);


    write_errs1(Y1_e_h, Y2_e_h, Y1_e_m_h, Y2_e_m_h, Y1_rk2_h, Y2_rk2_h, Y1_rk4_h, Y2_rk4_h, Y1_a_h, Y2_a_h, Y1_e_h2,
                Y2_e_h2, Y1_e_m_h2, Y2_e_m_h2, Y1_rk2_h2, Y2_rk2_h2, Y1_rk4_h2, Y2_rk4_h2, Y1_a_h2, Y2_a_h2, Y1_e_2h,
                Y2_e_2h, Y1_e_m_2h, Y2_e_m_2h, Y1_rk2_2h, Y2_rk2_2h, Y1_rk4_2h, Y2_rk4_2h, Y1_a_2h, Y2_a_2h, n1, n2, n3,
                "errors");

    write_f(Y1_e_h, Y1_e_m_h, Y1_rk2_h, Y1_rk4_h, Y1_a_h,
            Y2_e_h, Y2_e_m_h, Y2_rk2_h, Y2_rk4_h, Y2_a_h, X_h, n1);
//    system("python3 vis.py");
    return 0;
}
