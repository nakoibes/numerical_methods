#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>

using namespace std;


double p(double x) {
    return exp(x);
}

double q(double x) {
    return x * x;
}

double f(double x) {
    return sin(x);
}

void f_even_args(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

void formA(double **A, double *X, double h, int n) {
    for (int i = 0; i < n; i++) {
        A[i][2] = (2.0 + p(X[i + 1]) * h) / (2.0 * h * h);
        A[i][0] = (-4.0 - 2.0 * h * h * q(X[i + 1])) / (2.0 * h * h);
        A[i][1] = (2.0 - p(X[i + 1]) * h) / (2.0 * h * h);;
    }
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


void f_diff(double *Y, double *X, double h, int n) {

    double F_p[n - 2];
    double Y_p[n - 2];

    double **A = new double *[n - 2];
    for (int i = 0; i < n - 2; i++) {
        A[i] = new double[3];
    }

    formA(A, X, h, n - 2);

    for (int i = 1; i < n - 2; i++) {
        F_p[i - 1] = f(X[i]);
    }

    F_p[0] -= A[0][2] * Y[0];
    F_p[n - 3] -= A[n - 3][1] * Y[n - 1];


    prog(A, Y_p, F_p, n - 2);

    for (int i = 1; i < n - 1; i++) {
        Y[i] = Y_p[i - 1];
    }

}

void write_f(double *Y1, double *Y2, double *Y3, double *Y4, double *X, int n) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n; i++) {
        fout << X[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y2[i] << " ";
    }

    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y3[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y4[i] << " ";
    }
    fout << endl;


    fout.close();
}

double f1(double x, double y1, double y2) {
    return y2;
}

double f2(double x, double y1, double y2) {
    return y2 * exp(x) + y1 * x * x + sin(x);
}

double f2_2(double x, double y1, double y2) {
    return y2 * y2 * exp(x) + y1 * y1 * x * x + sin(x);
}

double g(double x, double y1, double y2, double Y, double U) {
    return 2 * y1 * x * x * Y + 2 * y2 * exp(x) * U;
}


void ru_ku4(double *Y1, double *Y2, double *X, int n, double h, function<double(double, double, double)> f2_) {

    for (int i = 1; i < n; i++) {
        double k1 = h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double k2 = h * f1(X[i - 1] + h / 2.0, Y1[i - 1] + h * k1 / 2.0, Y2[i - 1] + h * k1 / 2.0);
        double k3 = h * f1(X[i - 1] + h / 2.0, Y1[i - 1] + h * k2 / 2.0, Y2[i - 1] + h * k2 / 2.0);
        double k4 = h * f1(X[i - 1] + h, Y1[i - 1] + h * k3, Y2[i - 1] + h * k3);
        Y1[i] = Y1[i - 1] + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        k1 = h * f2_(X[i - 1], Y1[i - 1], Y2[i - 1]);
        k2 = h * f2_(X[i - 1] + h / 2.0, Y1[i - 1] + h * k1 / 2.0, Y2[i - 1] + h * k1 / 2.0);
        k3 = h * f2_(X[i - 1] + h / 2.0, Y1[i - 1] + h * k2 / 2.0, Y2[i - 1] + h * k2 / 2.0);
        k4 = h * f2_(X[i - 1] + h, Y1[i - 1] + h * k3, Y2[i - 1] + h * k3);
        Y2[i] = Y2[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }
}

void ru_ku4_g(double *Y1, double *Y2, double *X, double *y, double *u, int n, double h) {

    for (int i = 1; i < n; i++) {
        double k1 = h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double k2 = h * f1(X[i - 1] + h / 2.0, Y1[i - 1] + h * k1 / 2.0, Y2[i - 1] + h * k1 / 2.0);
        double k3 = h * f1(X[i - 1] + h / 2.0, Y1[i - 1] + h * k2 / 2.0, Y2[i - 1] + h * k2 / 2.0);
        double k4 = h * f1(X[i - 1] + h, Y1[i - 1] + h * k3, Y2[i - 1] + h * k3);
        Y1[i] = Y1[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        k1 = h * g(X[i - 1], y[i - 1], u[i - 1], Y1[i - 1], Y2[i - 1]);
        k2 = h * g(X[i - 1] + h / 2.0, y[i - 1] + h * k1 / 2.0, u[i - 1] + h * k1 / 2.0, Y1[i - 1] + h * k1 / 2.0,
                   Y2[i - 1] + h * k1 / 2.0);
        k3 = h * g(X[i - 1] + h / 2.0, y[i - 1] + h * k2 / 2.0, u[i - 1] + h * k2 / 2.0, Y1[i - 1] + h * k2 / 2.0,
                   Y2[i - 1] + h * k2 / 2.0);
        k4 = h * g(X[i - 1] + h, y[i - 1] + h * k3, u[i - 1] + h * k3, Y1[i - 1] + h * k3, Y2[i - 1] + h * k3);
        Y2[i] = Y2[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }
}

void shoot(double *Y, double *X, int n, double h) {
    double eps = 0.000001;
    double t1 = -5.0;
    double t2 = 5.0;
    double yb = Y[n - 1];
    double U[n];
    double fi;

    U[0] = (t1 + t2) / 2.0;
    ru_ku4(Y, U, X, n, h, f2);
    fi = Y[n - 1] - yb;
    if (fi < 0.0) {
        t1 = (t1 + t2) / 2.0;
    } else {
        t2 = (t1 + t2) / 2.0;
    }

    int i = 2;
    while (true) {
        U[0] = (t1 + t2) / 2.0;
        ru_ku4(Y, U, X, n, h, f2);

        if (abs(fi - (Y[n - 1] - yb)) < eps) {
            cout << i << " shooting iter" << endl;
//            cout << t1 << endl;
            break;
        }
        fi = Y[n - 1] - yb;
        if (fi < 0.0) {
            t1 = (t1 + t2) / 2.0;
        } else {
            t2 = (t1 + t2) / 2.0;
        }

        if (i == 5000) {
            cout << "INFINITE LOOP" << endl;
            break;
        }
        i++;
    }
//    int m =10;//выбор области в которой будем решать уравнение
//    double T[m];
//    for (int i = 1; i < m; i++) {
//        T[i] = -5 + i;
//    }
//    for (int i = 1; i < m; i++) {
//        U[0] = T[i];
//        ru_ku4(Y,U,X,n,h);
//        cout << Y[n-1] << " " << endl;
//    }

}

void newton(double *Y, double *X, int n, double h, bool sec) {
    double eps = 0.000001;
    double s = 0.0;
    double yb = Y[n - 1];
    double fi_s;
    double fi_s_h;
    double h_d = 0.01;
    double U[n];

    double Y1[n];
    double U1[n];
    Y1[0] = 0.0;
    U1[0] = 1.0;

    int i = 2;
    double s_n;
    double dfi;

    if (sec) {
        U[0] = s + h_d;
        ru_ku4(Y, U, X, n, h, f2_2);
        fi_s_h = Y[n - 1] - yb;
        U[0] = s;
        ru_ku4(Y, U, X, n, h, f2_2);
        fi_s = Y[n - 1] - yb;
        dfi = ((fi_s_h - fi_s) / h_d);
    } else {
        U[0] = s;
        ru_ku4(Y, U, X, n, h, f2_2);
        fi_s = Y[n - 1] - yb;

        ru_ku4_g(Y1, U1, X, Y, U, n, h);
        dfi = Y1[n - 1];

    }
    s = s - fi_s / dfi;
    while (true) {
        if (sec) {
            U[0] = s + h_d;
            ru_ku4(Y, U, X, n, h, f2_2);
            fi_s_h = Y[n - 1] - yb;
            U[0] = s;
            ru_ku4(Y, U, X, n, h, f2_2);
            fi_s = Y[n - 1] - yb;
            dfi = ((fi_s_h - fi_s) / h_d);
        } else {
            U[0] = s;
            ru_ku4(Y, U, X, n, h, f2_2);
            fi_s = Y[n - 1] - yb;

            ru_ku4_g(Y1, U1, X, Y, U, n, h);
            dfi = Y1[n - 1];

        }
        s_n = s - fi_s / dfi;
        if (abs(s_n - s) < eps) {
            cout << i << " newton iter" << endl;
//            cout << s_n << endl;
            break;
        }
        s = s_n;
        if (i == 5000) {
            cout << "INFINITE LOOP" << endl;
            break;
        }
        i++;
    }

//    for (int j = 0; j < 20; j++) {
//        U[0] = -1 + 0.1 * j;
//        ru_ku4(Y, U, X, n, h, f2_2);
//        cout << Y[n - 1] - yb << " ";
//    }
//    cout << endl;

}

void sh(double *Y, double *Y_c, int n) {
    for (int i = 0; i < n; i++) {
        Y_c[i] = Y[2 * i];
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

void write_errs(double *Y_diff_h, double *Y_shoo_h, double *Y_sec_h, double *Y_new_h,
                double *Y_diff_h2, double *Y_shoo_h2, double *Y_sec_h2, double *Y_new_h2,
                int n1, string filename) {

    double Y_diff_h2_s[n1];
    double Y_shoo_h2_s[n1];
    double Y_sec_h2_s[n1];
    double Y_new_h2_s[n1];
    sh(Y_diff_h2, Y_diff_h2_s, n1);
    sh(Y_shoo_h2, Y_shoo_h2_s, n1);
    sh(Y_sec_h2, Y_sec_h2_s, n1);
    sh(Y_new_h2, Y_new_h2_s, n1);

//    double U_diff_h2_s[n1];
//    double U_shoo_h2_s[n1];
//    double U_sec_h2_s[n1];
//    double U_new_h2_s[n1];
//    sh(U_diff_h2, U_diff_h2_s, n1);
//    sh(U_shoo_h2, U_shoo_h2_s, n1);
//    sh(U_sec_h2, U_sec_h2_s, n1);
//    sh(U_new_h2, U_new_h2_s, n1);

    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);

    fout << setw(15) << left << "error"
         << setw(15) << left << "y step h;h/2"
         //         << setw(15) << left << "y' step h;h/2"
         << endl;

    fout << "Finite difference" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1(Y_diff_h, Y_diff_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_1(U_diff_h, U_diff_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y_diff_h, Y_diff_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_2(U_diff_h, U_diff_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y_diff_h, Y_diff_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_cheb(U_diff_h, U_diff_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y_diff_h, Y_diff_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_1(U_diff_h, U_diff_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y_diff_h, Y_diff_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_2(U_diff_h, U_diff_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y_diff_h, Y_diff_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_cheb(U_diff_h, U_diff_h2_s, n1)
         << endl;

    fout << "Shooting" << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y_shoo_h, Y_shoo_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_2(U_shoo_h, U_shoo_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y_shoo_h, Y_shoo_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_cheb(U_shoo_h, U_shoo_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y_shoo_h, Y_shoo_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_1(U_shoo_h, U_shoo_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y_shoo_h, Y_shoo_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_2(U_shoo_h, U_shoo_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y_shoo_h, Y_shoo_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_cheb(U_shoo_h, U_shoo_h2_s, n1)
         << endl;

    fout << "Secant" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << calc_abs_err_1(Y_sec_h, Y_sec_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_1(U_sec_h, U_sec_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y_sec_h, Y_sec_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_2(U_sec_h, U_sec_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y_sec_h, Y_sec_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_cheb(U_sec_h, U_sec_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y_sec_h, Y_sec_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_1(U_sec_h, U_sec_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y_sec_h, Y_sec_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_2(U_sec_h, U_sec_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y_sec_h, Y_sec_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_cheb(U_sec_h, U_sec_h2_s, n1)
         << endl;

    fout << "Newton" << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << calc_abs_err_1(Y_new_h, Y_new_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_1(U_new_h2, U_new_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(Y_new_h, Y_new_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_2(U_new_h2, U_new_h2_s, n1)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(Y_new_h, Y_new_h2_s, n1)
         //         << setw(15) << left << calc_abs_err_cheb(U_new_h2, U_new_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(Y_new_h, Y_new_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_1(U_new_h2, U_new_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(Y_new_h, Y_new_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_2(U_new_h2, U_new_h2_s, n1)
         << endl;
    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(Y_new_h, Y_new_h2_s, n1)
         //         << setw(15) << left << calc_rel_err_cheb(U_new_h2, U_new_h2_s, n1)
         << endl;

    fout.close();
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double ya = 0.0;
    double yb = 0.0;
    int n1 = 51;
    int n2 = 2 * n1 - 1;
    double h1 = (b - a) / (n1 - 1);
    double h2 = (b - a) / (n2 - 1);

    double X_h[n1];
    double X_h2[n2];

    double Y_ra_h[n1];
    Y_ra_h[0] = ya;
    Y_ra_h[n1 - 1] = yb;
    double Y_ra_h2[n2];
    Y_ra_h2[0] = ya;
    Y_ra_h2[n1 - 1] = yb;

    double Y_shoo_h[n1];
    Y_shoo_h[0] = ya;
    Y_shoo_h[n1 - 1] = yb;
    double Y_shoo_h2[n2];
    Y_shoo_h2[0] = ya;
    Y_shoo_h2[n1 - 1] = yb;

    double Y_sec_h[n1];
    Y_sec_h[0] = ya;
    Y_sec_h[n1 - 1] = yb;
    double Y_sec_h2[n2];
    Y_sec_h2[0] = ya;
    Y_sec_h2[n1 - 1] = yb;

    double Y_new_h[n1];
    Y_new_h[0] = ya;
    Y_new_h[n1 - 1] = yb;
    double Y_new_h2[n2];
    Y_new_h2[0] = ya;
    Y_new_h2[n1 - 1] = yb;

    f_even_args(X_h, n1, a, b);
    f_even_args(X_h2, n1, a, b);

    f_diff(Y_ra_h, X_h, h1, n1);
    shoot(Y_shoo_h, X_h, n1, h1);

    newton(Y_sec_h, X_h, n1, h1, true);
    newton(Y_new_h, X_h, n1, h1, false);


    f_diff(Y_ra_h2, X_h2, h2, n2);
    shoot(Y_shoo_h2, X_h2, n2, h2);

    newton(Y_sec_h2, X_h2, n2, h2, true);
    newton(Y_new_h2, X_h2, n2, h2, false);

    write_errs(Y_ra_h, Y_shoo_h, Y_sec_h, Y_new_h, Y_ra_h2, Y_shoo_h2, Y_sec_h2, Y_new_h2, n1, "errors");

    write_f(Y_ra_h, Y_shoo_h, Y_sec_h, Y_new_h, X_h, n1);

    system("python3 vis.py");

    return 0;
}
