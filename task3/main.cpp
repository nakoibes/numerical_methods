#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double f(double x) {
    return sin(x);
}

double random(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

double fi(double x, double *xi, int i, int len_xi) {
    double result = 1.0;
    for (int j = 0; j < len_xi; j++) {
        if (j != i) {
            result *= (x - xi[j]) / (xi[i] - xi[j]);
        }
    }
    return result;
}

void f_even2(double **x, int k, int n, double a, double h) {
    double current = a;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            x[i][j] = current;
            if (j != n - 1) {
                current += h;
            }
        }
    }
}

void f_r_x(double **x, double **x_r, int k, int n, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            x_r[i][j] = random(x[i][0], x[i][n - 1]);
        }
    }
}

MatrixXd f_A(MatrixXd A, double **x_r, double **x, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i * (n - 1) + j1) <= (i * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
                        A(i * (n - 1) + j1, i * (n - 1) + j2) +=
                                fi(x_r[i][q], x[i], j1, n) * fi(x_r[i][q], x[i], j2, n);
                    }
                    A(i * (n - 1) + j2, i * (n - 1) + j1) = A(i * (n - 1) + j1, i * (n - 1) + j2);
                }
            }
        }
    }
    return A;
}

VectorXd f_B(VectorXd B, double **x_r, double **x, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            for (int q = 0; q < l; q++) {
                B[i * (n - 1) + j] += (f(x_r[i][q]) * fi(x_r[i][q], x[i], j, n));
            }
        }
    }
    return B;
}


void f_even(double *x, int n, double a, double b) {
    x[0] = a;
    x[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        x[i] = x[i - 1] + step;
    }
}

void f_y(MatrixXd X, double *x_v, double *y_v, double **x, int n_v, int n, int k) {
    double eps = 0.0000001;
    for (int i = 0; i < n_v; i++) {
        for (int j = 0; j < k; j++) {
            if (x_v[i] - x[j][n - 1] <= eps && x_v[i] - x[j][0] >= -eps) {
                for (int q = 0; q < n; q++) {
                    y_v[i] += X(j * (n - 1) + q) * fi(x_v[i], x[j], q, n);
                }
            }
        }
    }
}

void fill_func(double *func, double *x, int n) {
    for (int i = 0; i < n; i++) {
        func[i] = f(x[i]);
    }
}

void write_f(double *x_v, double *y_v, double *x_k, double *y_k, double *func, double *d_knots,
             int m, int n_v, int k) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n_v; i++) {
        fout << x_v[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_v; i++) {
        fout << y_v[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < m; i++) {
        fout << x_k[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < m; i++) {
        fout << y_k[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < k + 1; i++) {
        fout << d_knots[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_v; i++) {
        fout << func[i] << " ";
    }
    fout.close();
}

void f_div_k(double *div_k, double **all_k, int k, int n) {
    div_k[0] = all_k[0][0];
    for (int i = 0; i < k; i++) {
        div_k[i + 1] = all_k[i][n - 1];
    }
}

void zeronize(double *data, int n) {
    for (int i = 0; i < n; i++) {
        data[i] = 0.0;
    }
}

void f_x(double **X_r, double *x_r, int k, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            x_r[i * l + j] = X_r[i][j];
        }
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

void write_errs(double *f, double *y, int n_err, string filename) {
    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);
    fout << "abs_err_1 " << scientific << calc_abs_err_1(f, y, n_err) << endl;
    fout << "abs_err_2 " << calc_abs_err_2(f, y, n_err) << endl;
    fout << "abs_err_c " << calc_abs_err_cheb(f, y, n_err) << endl;
    fout << "rel_err_1 " << calc_rel_err_1(f, y, n_err) << endl;
    fout << "rel_err_2 " << calc_rel_err_2(f, y, n_err) << endl;
    fout << "rel_err_c " << calc_rel_err_cheb(f, y, n_err) << endl;
    fout.close();
}

int main() {
    srand((unsigned int) time(nullptr));

    double a = -3.14;
    double b = 3.14;

    int n_v = 2000;
    int k = 4;
    int n = 3;
    int m = (n - 1) * k + 1;
    int l = 15;
    double h = (b - a) / (m - 1);
    int n_err = 100 * (m - 1) + 1;


    double **X_k = new double *[k];
    for (int i = 0; i < k; i++) {
        X_k[i] = new double[n];
    }

    double d_knots[k + 1];

    double x_k[m];
    double y_k[m];
    double func[n_v];
    double x_v[n_v];
    double y_v[n_v];

    double **X_r = new double *[k];
    for (int i = 0; i < k; i++) {
        X_r[i] = new double[l];
    }
    double x_r[l * k];
    double f_r[l * k];
    double y_r[l * k];

    double x100[n_err];
    double f100[n_err];
    double y100[n_err];

    MatrixXd A(m, m);
    A.setZero();

    VectorXd B(m);
    VectorXd X(m);
    B.setZero();

    f_even2(X_k, k, n, a, h);
    f_div_k(d_knots, X_k, k, n);

    f_r_x(X_k, X_r, k, n, l);

    A = f_A(A, X_r, X_k, l, k, n);
//    cout << A << endl;
    B = f_B(B, X_r, X_k, l, k, n);
    X = A.colPivHouseholderQr().solve(B);

    f_even(x_v, n_v, a, b);
    f_y(X, x_v, y_v, X_k, n_v, n, k);

    fill_func(func, x_v, n_v);

    f_even(x_k, m, a, b);
    zeronize(y_k, m);

    f_x(X_r, x_r, k, l);
    fill_func(f_r, x_r, l * k);
    f_y(X, x_r, y_r, X_k, l * k, n, k);

    f_even(x100, n_err, a, b);
    fill_func(f100, x100, n_err);
    f_y(X, x100, y100, X_k, n_err, n, k);

    write_errs(f_r, y_r, l * k, "rand_errs.txt");
    write_errs(f100, y100, n_err, "errs.txt");

    write_f(x_v, y_v, x_k, y_k, func, d_knots, m, n_v, k);

    system("python3 vis.py");

    for (int i = 0; i < k; i++) {
        delete[] X_k[i];
        delete[] X_r[i];
    }
    delete[] X_k;
    delete[] X_r;

    return 0;
}