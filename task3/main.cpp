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

//    i = i % len_xi;
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

MatrixXd f_A(MatrixXd A, double **random_args, double **args, int l, int m, int k, int n) {
    for (int i1 = 0; i1 < k; i1++) {
        for (int j1 = 0; j1 < n; j1++) {
//            for (int i2 = 0; i2 < k; i2++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i1 * (n - 1) + j1) <= (i1 * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
                        A(i1 * (n - 1) + j1, i1 * (n - 1) + j2) +=
                                fi(random_args[i1][q], args[i1], j1, n) * fi(random_args[i1][q], args[i1], j2, n);
                    }
                    A(i1 * (n - 1) + j2, i1 * (n - 1) + j1) = A(i1 * (n - 1) + j1, i1 * (n - 1) + j2);
//                        double value = 0;
//                        for (int q = 0; q < l * k; q++) {
//                            value += (fi(random_args[q], args[i / n], i, n) * fi(random_args[q], args[j / n], j, n));
//                        }
//                        A(i, j) = value;
//                        A(j, i) = value;
                }

            }
        }
    }
    return A;
}

VectorXd f_B(VectorXd B, function<double(double)> y, double **random_args, double **args, int l, int m, int k, int n) {
    for (int i = 0; i < k; i++) {
//        double value = 0;
        for (int j = 0; j < n; j++) {
            for (int q = 0; q < l; q++) {
                B[i * (n - 1) + j] += (y(random_args[i][q]) * fi(random_args[i][q], args[i], j, n));
//                value += (y(random_args[j]) * fi(random_args[j], args[i / n], i, n));
            }
        }
//        B[i] = value;
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

void f_repr_vals(MatrixXd X, double *repr_args, double *repr_vals, double **args, int n_repr, int m, int n, int k) {
    for (int i = 0; i < n_repr; i++) {
//        double value = 0;
        for (int j = 0; j < k; j++) {
            if (repr_args[i] < args[j][n - 1] && repr_args[i] > args[j][0]) {
                for (int q = 0; q < n; q++) {
                    repr_vals[i] += X(j * (n - 1) + q) * fi(repr_args[i], args[j], q, n);
                }
            }
        }
//        repr_vals[i] = value;
    }
}

void fill_func(double *func, function<double(double)> y, double *res_args, int n_repr) {
    for (int i = 0; i < n_repr; i++) {
        func[i] = y(res_args[i]);
    }
}

void write_f(double *res_args, double *res_vals, double *inter_args, double *inter_vals, double *func, int n_inter,
             int n_res) {
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
    for (int i = 0; i < n_res; i++) {
        fout << func[i] << " ";
    }
    fout.close();
}

void f_knots(double *knots, double **args, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n - 1; j++) {
            knots[i * n + j] = args[i][j];
        }
    }
}

void zeronize(double *array, int n) {
    for (int i = 0; i < n; i++) {
        array[i] = 0;
    }
}

int main() {
    srand((unsigned int) time(nullptr));

    function<double(double)> y = func;

    double a = -4;
    double b = 4;

    double delta = 0.0001;

    int n_repr = 2000;
    int k = 3;
    int n = 4;
    int m = (n - 1) * k + 1;
    int l = 5;
    double h = (b - a) / (k * (n - 1));
//    double h_repr = (b-a)/(n_repr-1)


    double **args = new double *[k];
    for (int i = 0; i < k; i++) {
        args[i] = new double[n];
    }

    double knots[m];
    double knots_vals[m];
    double func[n_repr];
    double repr_args[n_repr];
    double repr_vals[n_repr];
//    double **random_args = new double *[k];
//    for (int i = 0; i < k; i++) {
//        random_args[i] = new double[l];
//    }

//    double random_args[l * k];
    double **random_args = new double *[k];
    for (int i = 0; i < k; i++) {
        random_args[i] = new double[l];
    }


    MatrixXd A(m, m);
    A.setZero();

//    double **A = new double *[m];
//    for (int i = 0; i < m; i++) {
//        A[i] = new double[m];
//    }

    VectorXd B(m);
    VectorXd X(m);
    B.setZero();

    f_args(args, k, n, a, h);
    f_random_args(args, random_args, k, n, l);
    A = f_A(A, random_args, args, l, m, k, n);
//    cout << A << endl;
    B = f_B(B, y, random_args, args, l, m, k, n);
    X = A.colPivHouseholderQr().solve(B);
    f_repr_args(repr_args, n_repr, a, b);
    f_repr_vals(X, repr_args, repr_vals, args, n_repr, m, n, k);
    fill_func(func, y, repr_args, n_repr);
    f_repr_args(knots, m, a, b);
    zeronize(knots_vals, m);
    write_f(repr_args, repr_vals, knots, knots_vals, func, m, n_repr);
    system("python3 repr.py");

//    cout << X;
//    MatrixXd T(2, 2);
//    T(0,0) = 2;
//    T(0,1) = 2;
//    T(1,0) = 1;
//    T(1,1) = 2;
//    VectorXd U(2);
//    VectorXd R(2);
//    U(0) = 2;
//    U(1) = 3;
//    R = T.colPivHouseholderQr().solve(U);
//    cout << R;

    for (int i = 0; i < k; i++) {
        delete[] args[i];
    }
    delete[] args;


//    for (int i = 0; i < m; i++) {
//        delete[] A[i];
//    }
//    delete[] A;

//    for (int i = 0; i < k; i++) {
//        delete[] random_args[i];
//    }
//    delete[] random_args;

    return 0;
}
