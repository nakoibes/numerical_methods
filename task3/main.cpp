#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double func(double x) {
    return sin(x);
}

double random(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

double fi(double x, double *xi, int i, int len_xi) {
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

void f_random_args(double **args, double **random_args, int k, int len2, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            random_args[i][j] = random(args[i][0], args[i][len2 - 1]);
        }
    }
}

MatrixXd f_A(MatrixXd A, double **random_args, double **args, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i * (n - 1) + j1) <= (i * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
                        A(i * (n - 1) + j1, i * (n - 1) + j2) +=
                                fi(random_args[i][q], args[i], j1, n) * fi(random_args[i][q], args[i], j2, n);
                    }
                    A(i * (n - 1) + j2, i * (n - 1) + j1) = A(i * (n - 1) + j1, i * (n - 1) + j2);
                }
            }
        }
    }
    return A;
}

VectorXd f_B(VectorXd B, function<double(double)> y, double **random_args, double **args, int l, int m, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            for (int q = 0; q < l; q++) {
                B[i * (n - 1) + j] += (y(random_args[i][q]) * fi(random_args[i][q], args[i], j, n));
            }
        }
    }
    return B;
}


void f_repr_args(double *repr_args, int n_repr, double a, double b) {
    repr_args[0] = a;
    repr_args[n_repr - 1] = b;
    double step = (b - a) / (n_repr - 1);
    for (int i = 1; i < n_repr - 1; i++) {
        repr_args[i] = repr_args[i - 1] + step;
    }
}

void f_repr_vals(MatrixXd X, double *repr_args, double *repr_vals, double **args, int n_repr, int n, int k) {
    double eps = 0.00001;
    for (int i = 0; i < n_repr; i++) {
        for (int j = 0; j < k; j++) {
            if (repr_args[i] - args[j][n - 1] <= eps && repr_args[i] - args[j][0] >= -eps) {
                for (int q = 0; q < n; q++) {
                    repr_vals[i] += X(j * (n - 1) + q) * fi(repr_args[i], args[j], q, n);
                }
            }
        }
    }
}

void fill_func(double *func, function<double(double)> y, double *res_args, int n_repr) {
    for (int i = 0; i < n_repr; i++) {
        func[i] = y(res_args[i]);
    }
}

void
write_f(double *res_args, double *res_vals, double *inter_args, double *inter_vals, double *func, double *div_knots,
        int n_inter,
        int n_res, int k) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n_res; i++) {
        fout << res_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_res; i++) {
        fout << res_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_inter; i++) {
        fout << inter_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_inter; i++) {
        fout << inter_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < k + 1; i++) {
        fout << div_knots[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_res; i++) {
        fout << func[i] << " ";
    }
    fout.close();
}

void f_div_knots(double *div_knots, double **args, int k, int n) {
    div_knots[0] = args[0][0];
    for (int i = 0; i < k; i++) {
        div_knots[i + 1] = args[i][n - 1];
    }
}

void zeronize(double *array, int n) {
    for (int i = 0; i < n; i++) {
        array[i] = 0;
    }
}

void f_random_args_err(double **args, double *random_args, int k, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            random_args[i * l + j] = args[i][j];
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

void write_errs(double *err_func, double *err_vals, int n_err, string filename) {
    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);
    fout << "abs_err_1 " << scientific << calc_abs_err_1(err_func, err_vals, n_err) << endl;
    fout << "abs_err_2 " << calc_abs_err_2(err_func, err_vals, n_err) << endl;
    fout << "abs_err_c " << calc_abs_err_cheb(err_func, err_vals, n_err) << endl;
    fout << "rel_err_1 " << calc_rel_err_1(err_func, err_vals, n_err) << endl;
    fout << "rel_err_2 " << calc_rel_err_2(err_func, err_vals, n_err) << endl;
    fout << "rel_err_c " << calc_rel_err_cheb(err_func, err_vals, n_err) << endl;
    fout.close();
}

void f_even_args(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

int main() {
    srand((unsigned int) time(nullptr));

    function<double(double)> y = func;

    double a = -7;
    double b = 7;

    int n_repr = 2000;
    int k = 3;
    int n = 4;
    int m = (n - 1) * k + 1;
    int l = 15;
    double h = (b - a) / (k * (n - 1));
    int n_err = 100 * (m - 1) + 1;


    double **args = new double *[k];
    for (int i = 0; i < k; i++) {
        args[i] = new double[n];
    }

    double div_knots[k + 1];

    double knots[m];
    double knots_vals[m];
    double func[n_repr];
    double repr_args[n_repr];
    double repr_vals[n_repr];

    double **random_args = new double *[k];
    for (int i = 0; i < k; i++) {
        random_args[i] = new double[l];
    }
    double err_args_rand[l * k];
    double err_func_rand[l * k];
    double err_vals_rand[l * k];

    double err_args[n_err];
    double err_func[n_err];
    double err_vals[n_err];

    MatrixXd A(m, m);
    A.setZero();

    VectorXd B(m);
    VectorXd X(m);
    B.setZero();

    f_args(args, k, n, a, h);
    f_div_knots(div_knots, args, k, n);

    f_random_args(args, random_args, k, n, l);

    A = f_A(A, random_args, args, l, k, n);
    B = f_B(B, y, random_args, args, l, m, k, n);
    X = A.colPivHouseholderQr().solve(B);

    f_repr_args(repr_args, n_repr, a, b);
    f_repr_vals(X, repr_args, repr_vals, args, n_repr, n, k);

    fill_func(func, y, repr_args, n_repr);

    f_repr_args(knots, m, a, b);
    zeronize(knots_vals, m);

    f_random_args_err(random_args, err_args_rand, k, l);
    fill_func(err_func_rand, y, err_args_rand, l * k);
    f_repr_vals(X, err_args_rand, err_vals_rand, args, l * k, n, k);

    f_even_args(err_args, n_err, a, b);
    fill_func(err_func, y, err_args, n_err);
    f_repr_vals(X, err_args, err_vals, args, n_err, n, k);

    write_errs(err_func_rand, err_vals_rand, l * k, "rand_errs.txt");
    write_errs(err_func, err_vals, n_err, "errs.txt");

    write_f(repr_args, repr_vals, knots, knots_vals, func, div_knots, m, n_repr, k);

    system("python3 repr.py");

    for (int i = 0; i < k; i++) {
        delete[] args[i];
    }
    delete[] args;

    return 0;
}
