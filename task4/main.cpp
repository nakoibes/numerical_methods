#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double f(double x) {
    return sin(x);
}

double df(double x) {
    return cos(x);
}

void f_even(double *x, int n, double a, double b) {
    x[0] = a;
    x[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        x[i] = x[i - 1] + step;
    }
}

void fill_func(double *func, double *x, int n) {
    for (int i = 0; i < n; i++) {
        func[i] = f(x[i]);
    }
}

void fill_dfunc(double *func, double *x, int n) {
    for (int i = 0; i < n; i++) {
        func[i] = df(x[i]);
    }
}

void differentiate(double *y, double *dy, int n, double h) {
    dy[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);
    dy[n - 1] = (y[n - 3] - 4 * y[n - 2] + 3 * y[n - 1]) / (2 * h);
    for (int i = 1; i < n - 1; i++) {
        dy[i] = (y[i + 1] - y[i - 1]) / (2 * h);
    }
}

void calc_main_err(double *main_err, double *dy1, double *dy2, int m1) {
    for (int i = 0; i < m1; i += 1) {
        main_err[i] = (dy1[i] - dy2[2 * i]) / (0.5 * 0.5 - 1.0);
    }
}

void calc_runge(double *runge, double *main_err, double *dy, int n) {
    for (int i = 0; i < n; i++) {
        runge[i] = dy[i] + main_err[i];
    }
}

double calc_abs_err_1(double *v1, double *v2, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += abs(v1[i] - v2[i]);
    }
    return result;
}

double calc_abs_err_2(double *v1, double *v2, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += pow(abs(v1[i] - v2[i]), 2);
    }
    return sqrt(result);
}

double calc_abs_err_cheb(double *v1, double *v2, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        if (abs(v1[i] - v2[i]) > result) {
            result = abs(v1[i] - v2[i]);
        }
    }
    return result;
}

double calc_rel_err_1(double *v1, double *v2, int n) {
    double eps = 0.000000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        num += (abs(v1[i] - v2[i]));
    }
    for (int i = 0; i < n; i++) {
        den += abs(v1[i]);
    }
    if (den < eps) {
        den = 1.0;
    }
    return num / den;
}

double calc_rel_err_2(double *v1, double *v2, int n) {
    double eps = 0.000000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        num += (pow(abs(v1[i] - v2[i]), 2));
    }
    for (int i = 0; i < n; i++) {
        den += pow(v1[i], 2);
    }
    if (den < eps) {
        den = 1.0;
    }
    return sqrt(num) / sqrt(den);
}

double calc_rel_err_cheb(double *v1, double *v2, int n) {
    double eps = 0.000000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        if (abs(v1[i] - v2[i]) > num) {
            num = (abs(v1[i] - v2[i]));
        }
    }
    for (int i = 0; i < n; i++) {
        if (abs(v1[i]) > den) {
            den = (abs(v1[i]));
        }
    }
    if (den < eps) {
        den = 1.0;
    }
    return num / den;
}

void write_errs1(double *dy_h, double *der1, int n1, double *dy_h2, double *der2, int n2,
                 double *runge, double *main_err, string filename) {
    double zeros[n1];
    for (int i = 0; i < n1; i++) { zeros[i] = 0.0; }
    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);
    fout << setw(15) << left << "error"
         << setw(15) << left << "step h"
         << setw(15) << left << "step h/2"
         << setw(15) << left << "runge"
         << setw(15) << left << "main_err"
         << endl;

    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1(dy_h, der1, n1)
         << setw(15) << left << calc_abs_err_1(dy_h2, der2, n2)
         << setw(15) << left << calc_abs_err_1(dy_h, runge, n1)
         << setw(15) << left << calc_abs_err_1(zeros, main_err, n1)
         << endl;

    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(dy_h, der1, n1)
         << setw(15) << left << calc_abs_err_2(dy_h2, der2, n2)
         << setw(15) << left << calc_abs_err_2(dy_h, runge, n1)
         << setw(15) << left << calc_abs_err_2(zeros, main_err, n1)
         << endl;

    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(dy_h, der1, n1)
         << setw(15) << left << calc_abs_err_cheb(dy_h2, der2, n2)
         << setw(15) << left << calc_abs_err_cheb(dy_h, runge, n1)
         << setw(15) << left << calc_abs_err_cheb(zeros, main_err, n1)
         << endl;

    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(dy_h, der1, n1)
         << setw(15) << left << calc_rel_err_1(dy_h2, der2, n2)
         << setw(15) << left << calc_rel_err_1(dy_h, runge, n1)
         << endl;

    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(dy_h, der1, n1)
         << setw(15) << left << calc_rel_err_2(dy_h2, der2, n2)
         << setw(15) << left << calc_rel_err_2(dy_h, runge, n1)
         << endl;

    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(dy_h, der1, n1)
         << setw(15) << left << calc_rel_err_cheb(dy_h2, der2, n2)
         << setw(15) << left << calc_rel_err_cheb(dy_h, runge, n1)
         << endl;

    fout.close();
}

void write_f(double *x_v, double *y_v, double *x_h, double *der1, double *x_h2, double *der2,
             double *runge, int n1, int n2, int n_v) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n_v; i++) {
        fout << x_v[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_v; i++) {
        fout << y_v[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << x_h[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << der1[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n2; i++) {
        fout << x_h2[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n2; i++) {
        fout << der2[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << x_h[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << runge[i] << " ";
    }
    fout << endl;
    fout.close();
}

int main() {

    double a = -3.14;
    double b = 3.14;

    int m_v = 2000;
    int m1 = 11;
    int m2 = 2 * m1 - 1;
    double h1 = (b - a) / (m1 - 1);
    double h2 = (b - a) / (m2 - 1);

    double x_h[m1];
    double y_h[m1];
    double dy_h[m1];
    double main_err[m1];
    double runge[m1];
    double der1[m1];
    double x_h2[m2];
    double y_h2[m2];
    double dy_h2[m2];
    double der2[m2];
    double x_v[m_v];
    double y_v[m_v];


    f_even(x_v, m_v, a, b);
    f_even(x_h, m1, a, b);
    f_even(x_h2, m2, a, b);
    fill_func(y_h, x_h, m1);
    fill_func(y_h2, x_h2, m2);
    fill_dfunc(y_v, x_v, m_v);

    fill_dfunc(dy_h, x_h, m1);
    fill_dfunc(dy_h2, x_h2, m2);

    differentiate(y_h, der1, m1, h1);
    differentiate(y_h2, der2, m2, h2);
    calc_main_err(main_err, der1, der2, m1);
    calc_runge(runge, main_err, der1, m1);

    write_errs1(dy_h, der1, m1, dy_h2, der2, m2, runge, main_err, "errs1");

    write_f(x_v, y_v, x_h, der1, x_h2, der2, runge, m1, m2, m_v);

    system("python3 vis.py");

    return 0;
}
