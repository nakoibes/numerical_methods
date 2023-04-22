#include <iostream>
#include <functional>
#include <cmath>
#include <fstream>

using namespace std;

double f(double x) {
    return sin(x);
}

void f_even2(double **x, int k, int n1, double a, double h) {
    double current = a;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n1; j++) {
            x[i][j] = current;
            if (j != n1 - 1) {
                current += h;
            }
        }
    }
}

void f_func2(double **X, double **Y, int k, int n1) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n1; j++) {
            Y[i][j] = f(X[i][j]);
        }
    }
}

void f_even(double *x, int n, double a, double b) {
    x[0] = a;
    x[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        x[i] = x[i - 1] + step;
    }
}

void interp(double *y, double *x, double *y_k, double *x_k, int n, int n_k) {
    for (int i = 0; i < n; i++) {
        y[i] = 0.0;
    }
    double mult[n_k];
    for (int i = 0; i < n_k; i++) {
        double prod = 1.0;
        for (int j = 0; j < n_k; j++) {
            if (j != i) {
                prod *= (x_k[i] - x_k[j]);
            }
        }
        mult[i] = y_k[i] / prod;
    }
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n_k; k++) {
            double prod = 1;
            for (int j = 0; j < n_k; j++) {
                if (k != j) {
                    prod *= (x[i] - x_k[j]);
                }
            }
            y[i] += mult[k] * prod;
        }
    }
}

void f_y(double *y, double *x, double **Y, double **X, int k, int n1) {
    double eps = 0.000001;
    int start = 0;
    int end = 0;
    int l;
    for (int i = 0; i < k; i++) {
        double current = X[i][n1 - 1];
        while (current - x[end] > eps) {
            end++;
        }
        l = end - start + 1;
        interp(&y[start], &x[start], Y[i], X[i], l, n1);
        start = end;
    }
}


void f_knots(double *x_k, double *y_k, double **X_k, double **Y_k, int k, int n1) {
    int l = 0;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n1; j++) {
            if (i != k - 1) {
                if (j != n1 - 1) {
                    x_k[l] = X_k[i][j];
                    y_k[l] = Y_k[i][j];
                    l++;
                }
            } else {
                x_k[l] = X_k[i][j];
                y_k[l] = Y_k[i][j];
                l++;
            }
        }
    }
}

void write_f(double *x_v, double *y_v, double *x_k, double *y_k, double *func, double *el_k_x, double *el_k_y,
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
        fout << el_k_x[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < k + 1; i++) {
        fout << el_k_y[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_v; i++) {
        fout << func[i] << " ";
    }
    fout.close();
}

void f_func(double *y, double *x, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = f(x[i]);
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

void write_errs(double *err_f, double *err_y, int n) {
    ofstream fout("errs.txt");
    fout << fixed;
    fout.precision(6);
    fout << "abs_err_1   " << scientific << calc_abs_err_1(err_f, err_y, n) << endl;
    fout << "abs_err_2   " << calc_abs_err_2(err_f, err_y, n) << endl;
    fout << "abs_err_c   " << calc_abs_err_cheb(err_f, err_y, n) << endl;
    fout << "rel_err_1   " << calc_rel_err_1(err_f, err_y, n) << endl;
    fout << "rel_err_2   " << calc_rel_err_2(err_f, err_y, n) << endl;
    fout << "rel_err_c   " << calc_rel_err_cheb(err_f, err_y, n) << endl;
    fout.close();
}

void f_div_k(double *div_k, double **all_k, int k, int n) {
    div_k[0] = all_k[0][0];
    for (int i = 0; i < k; i++) {
        div_k[i + 1] = all_k[i][n];
    }
}

int main() {
    double a = -3.14;
    double b = 3.14;

    int k = 4;
    int n = 2;  //степень
    int n_v = 2000;

    int m = n * k + 1;

    double int_l = (b - a) / k;
    double h = int_l / n;

    int n_err = 100 * m - 99;

    double **X_k = new double *[k];
    for (int i = 0; i < k; i++) {
        X_k[i] = new double[n + 1];
    }

    double **Y_k = new double *[k];
    for (int i = 0; i < k; i++) {
        Y_k[i] = new double[n + 1];
    }

    double x_k[m];
    double y_k[m];
    double x_v[n_v];
    double y_v[n_v];
    double func[n_v];
    double func100[n_err];
    double x100[n_err];
    double y100[n_err];

    double el_k_x[k + 1];
    double el_k_y[k + 1];

    f_even2(X_k, k, n + 1, a, h);
    f_div_k(el_k_x, X_k, k, n);
    f_func2(X_k, Y_k, k, n + 1);
    f_div_k(el_k_y, Y_k, k, n);
    f_even(x_v, n_v, a, b);
    f_y(y_v, x_v, Y_k, X_k, k, n + 1);
    f_knots(x_k, y_k, X_k, Y_k, k, n + 1);

    f_func(func, x_v, n_v);

    f_even(x100, n_err, a, b);
    f_func(func100, x100, n_err);
    f_y(y100, x100, Y_k, X_k, k, n + 1);

    write_errs(func100, y100, n_err);
    write_f(x_v, y_v, x_k, y_k, func, el_k_x, el_k_y, m, n_v, k);

    system("python3 vis.py");

    for (int i = 0; i < k; i++) {
        delete[] X_k[i];
    }
    delete[] X_k;

    for (int i = 0; i < k; i++) {
        delete[] Y_k[i];
    }
    delete[] Y_k;

    return 0;
}
