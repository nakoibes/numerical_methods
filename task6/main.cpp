#include <iostream>
#include <functional>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


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

double func(double x) {
    return sin(x);
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

void f_r_args(double **args, double **random_args, int k, int len2, int l) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < l; j++) {
            random_args[i][j] = random(args[i][0], args[i][len2 - 1]);
        }
    }
}


void f_A(double **A, double **random_args, double **args, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i * (n - 1) + j1) <= (i * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
                        int i_u = i * (n - 1) + j1;
                        int j_u = i * (n - 1) + j2;
                        A[j_u - i_u][i_u] +=
                                fi(random_args[i][q], args[i], j1, n) * fi(random_args[i][q], args[i], j2, n);
                    }
                }
            }
        }
    }
}

void f_A_G(double **A, double **A_G, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A_G[i][j] = A[i][j];
        }
    }
    for (int i = n; i < 2 * n - 1; i++) {
        for (int j = 0; j < m; j++) {
            A_G[i][j] = A_G[i - n + 1][j];
        }
    }
}

void f_B(double *B, function<double(double)> y, double **random_args, double **args, int l, int m, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            for (int q = 0; q < l; q++) {
                B[i * (n - 1) + j] += (y(random_args[i][q]) * fi(random_args[i][q], args[i], j, n));
            }
        }
    }
}

VectorXd f_B_E(double *B, int m) {
    VectorXd B_E(m);
    for (int i = 0; i < m; i++) {
        B_E(i) = B[i];
    }
    return B_E;
}

void gauss_s(double **A, double *B, int m, int n, double eps) {
    for (int al = 0; al < m - 1; al++) {
        for (int be = 1; be < n; be++) {
            int i = be + n - 1;
            int j = al;
            if (abs(A[i][j]) > eps) {
                double koef = A[i][j] / A[0][j];
                A[i][j] = 0;
                B[al + be] = B[al + be] - B[al] * koef;
                for (int k = 1; k < n; k++) {
                    if (al + k < m) {
                        int i2;
                        int j2;
                        int i1 = k;
                        int j1 = al;
                        if (al + be <= al + k) {
                            i2 = k - be;
                            j2 = al + be;
                        } else {
                            i2 = be - k + n - 1;
                            j2 = al + k;
                        }
                        A[i2][j2] = A[i2][j2] - A[i1][j1] * koef;

                    } else {
                        break;
                    }
                }
            } else {
                continue;
            }

        }
    }

}

//void gauss_s_(double **Aup, double **Adown, double *B, int m, int n, double eps) {
//    for (int i_h = 0; i_h < m - 1; i_h++) {
//        for (int j_h = 1; i_h < n; j_h++) {
//            int fi = i_h + n;
//            if (fi > m) {
//                fi = m;
//            }
//            int i1 = j_h;
//            int j1 = i_h;
//            double koef = -Adown[i1][j1] / Adown[0][j1];
//            Adown[i1][j1] = 0;
//            for (int k = i_h + 1; k < fi; k++) {
//                if (i_h + j_h >= k) {
//                    i1 = i_h+j_h-k;
//                    j1 = k;
//                    Adown[i1][j1] += A[i][k] * koef;
//                }
//                else{
//                    i1 = k-i_h-j_h;
//                    j1= i_h+j_h;
//                }
//
//                B[i + j] += B[i] * koef;
//            }
//        }
//    }
//}

void gauss_r(double **A, double *B, int m, int n, double eps) {
    for (int al = m - 1; al > 0; al--) {
        for (int be = 1; be < n; be++) {
            int i = be;
            int j = al - be;
            if (abs(A[i][j]) > eps) {
                double koef = A[i][j] / A[0][al];
                A[i][j] = 0;
                B[al - be] = B[al - be] - B[al] * koef;
            } else {
                continue;
            }
        }
    }
}

void gauss(double **A, double *B, double *X, int m, int n) {
    double eps = 0.0000001;
    gauss_s(A, B, m, n, eps);
    gauss_r(A, B, m, n, eps);
    for (int i = 0; i < m; i++) {
        X[i] = B[i] / A[0][i];
    }
}

void g_s_LU(double **A, double *B, int m, int n) {
    for (int i = 0; i < m - 1; i++) {
        for (int j = 1; j < n; j++) {
            if (i + j < m) {
                B[i + j] -= B[i] * A[j][i];
                A[j][i] = 0.0;
            }
        }
    }
}

void g_s_H(double **A, double *B, int m, int n) {
    for (int i = 0; i < m - 1; i++) {
        for (int j = 1; j < n; j++) {
            if (i + j < m) {
                B[i + j] -= B[i] * A[j][i] / A[0][i];
                A[j][i] = 0.0;
            }
        }
    }
}

void gauss_LU(double **L, double **U, double *B, double *X, int m, int n) {
    double Y[m];
    double eps = 0.0000001;
    g_s_LU(L, B, m, n);
    for (int i = 0; i < m; i++) {
        Y[i] = B[i];
    }

    gauss_r(U, Y, m, n, eps);
    for (int i = 0; i < m; i++) {
        X[i] = Y[i] / U[0][i];
    }
}

void gauss_H(double **H, double **H_t, double *B, double *X, int m, int n) {
    double Y[m];
    double eps = 0.0000001;
    g_s_H(H, B, m, n);
    for (int i = 0; i < m; i++) {
        Y[i] = B[i] / H[0][i];
    }

    gauss_r(H_t, Y, m, n, eps);
    for (int i = 0; i < m; i++) {
        X[i] = Y[i] / H_t[0][i];
    }
}

MatrixXd f_A_E(double **A, int m, int n) {
    MatrixXd A_E(m, m);
    A_E.setZero();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (i + j < m) {
                A_E(j, i + j) = A[i][j];
                A_E(i + j, j) = A[i][j];
            }
        }
    }
    return A_E;
}

void prep_L(double **L, int m) {
    for (int i = 0; i < m; i++) {
        L[0][i] = 1;
    }
}

void LU(double **A, double **L, double **U, int m, int n) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            U[i][j] = A[i][j];
            int i_h = j;
            int j_h = j + i;
            int M;
            if (i_h > j_h) {
                M = i_h;
            } else {
                M = j_h;
            }
            int st = M - n + 1;
            if (st < 0) {
                st = 0;
            }

            for (int k = st; k < i_h; k++) {
                int i_l = i_h - k;
                int j_l = k;
                int i_u = j_h - k;
                int j_u = k;
                U[i][j] -= L[i_l][j_l] * U[i_u][j_u];
            }
        }
        for (int i = 1; i < n; i++) {
            L[i][j] = A[i][j];
            int i_h = j + i;
            int j_h = j;
            int M;
            if (i_h > j_h) {
                M = i_h;
            } else {
                M = j_h;
            }
            int st = M - n + 1;
            if (st < 0) {
                st = 0;
            }
            for (int k = st; k < j_h; k++) {
                int i_l = i_h - k;
                int j_l = k;
                int i_u = j_h - k;
                int j_u = k;
                L[i][j] -= L[i_l][j_l] * U[i_u][j_u];
            }
            int j_u = j_h;
            L[i][j] = L[i][j] / U[0][j_u];
        }

    }
}

void copy_vec(double *A, double *B, int m) {
    for (int i = 0; i < m; i++) {
        B[i] = A[i];
    }
}

void holec(double **A, double **H, int m, int n) {
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            int i_h = j;
            int j_h = j + i;
            int M;
            if (i_h > j_h) {
                M = i_h;
            } else {
                M = j_h;
            }
            int st = M - n + 1;
            if (st < 0) {
                st = 0;
            }
            if (i_h != j_h) {
                H[j_h - i_h][i_h] = A[i][j];
                for (int k = st; k < i_h; k++) {
                    H[j_h - i_h][i_h] -= H[i_h - k][k] * H[j_h - k][k];
                }
                H[j_h - i_h][i_h] = H[j_h - i_h][i_h] / H[0][i_h];
            } else {
                H[0][i_h] = A[i][j];
                for (int k = st; k < i_h; k++) {
                    H[0][i_h] -= H[i_h - k][k] * H[i_h - k][k];
                }
                H[0][i_h] = sqrt(H[0][i_h]);
            }
        }
    }
}

void copy_mat(double **H, double **H_t, int m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            H_t[i][j] = H[i][j];
        }
    }
}

double c_norm_1(double *X, double *Y, int m) {
    double norm = 0;
    for (int i = 0; i < m; i++) {
        norm += abs(X[i] - Y[i]);
    }
    return norm;
}

double scal_pr(double *A, double *B, int m) {
    double product = 0;
    for (int i = 0; i < m; i++) {
        product += A[i] * B[i];
    }
    return product;
}

void mat_vec(double **A, double *B, double *res, int m, int n) {
    for (int i = 0; i < m; i++) {
        int st = i - n + 1;
        if (st < 0)
            st = 0;
        int fi = i + n;
        if (fi > m) {
            fi = m;
        }
        for (int j = st; j < fi; j++) {
            if (i >= j) {
                res[i] += A[i - j][j] * B[j];
            } else {
                res[i] += A[j - i][i] * B[j];
            }
        }
    }
}

void vec_dif(double *A, double *B, double *res, int m) {
    for (int i = 0; i < m; i++) {
        res[i] = A[i] - B[i];
    }
}

void rel_comp(double **A, double *Y, double *X, double *B, int q, int m, int n) {
    X[q] = B[q];
    int st = q - n + 1;
    if (st < 0)
        st = 0;
    for (int i = st; i < q; i++) {
        X[q] -= A[q - i][i] * X[i];
    }
    int fi = q + n;
    if (fi > m) {
        fi = m;
    }
    for (int i = q + 1; i < fi; i++) {
        X[q] -= A[i - q][q] * Y[i];
    }
    X[q] = X[q] / A[0][q];
}

void rel_iter(double **A, double *Y, double *X, double *B, int m, int n, double w) {
    for (int i = 0; i < m; i++) {
        rel_comp(A, Y, X, B, i, m, n);
    }
    for (int i = 0; i < m; i++) {
        X[i] = (X[i] * w + (1.0 - w) * Y[i]);
    }
}

int over_rel(double **A, double *Y, double *X, double *B, double w, int m, int n, double eps, bool silent) {
    int i = 1;
    while (true) {
        rel_iter(A, Y, X, B, m, n, w);
        if (c_norm_1(X, Y, m) < eps) {
            if (!silent) {
                cout << "relaxation converged for " << i << " steps" << endl;
            }
            break;
        }

        for (int j = 0; j < m; j++) {
            Y[j] = X[j];
        }
        if (i == 10000) {
            if (!silent) {
                cout << "infinite loop" << endl;
            }

            break;
        }
        i++;
    }
    return i;
}

void conjug(double **A, double *Y, double *X, double *B, int m, int n, double eps) {
    int i = 1;
    int count = 0;
    double r_o[m];
    double r_n[m];
    double Ax[m];
    for (int j = 0; j < m; j++) {
        Ax[j] = 0.0;
    }
    mat_vec(A, Y, Ax, m, n);
    for (int j = 0; j < m; j++) {
        r_o[j] = B[j] - Ax[j];
    }
    double p[m];
    for (int j = 0; j < m; j++) {
        p[j] = r_o[j];
    }
    double al;
    double be;
    while (true) {
        double Apk[m];
        for (int j = 0; j < m; j++) {
            Apk[j] = 0.0;
        }
        mat_vec(A, p, Apk, m, n);
        double scalrr = scal_pr(r_o, r_o, m);
        al = scalrr / scal_pr(Apk, p, m);
        for (int j = 0; j < m; j++) {
            X[j] = Y[j] + al * p[j];
        }
        for (int j = 0; j < m; j++) {
            r_n[j] = r_o[j] - al * Apk[j];
        }
        be = scal_pr(r_n, r_n, m) / scalrr;
        for (int j = 0; j < m; j++) {
            p[j] = r_n[j] + be * p[j];
        }

        if (c_norm_1(X, Y, m) < eps) {
            count += 1;
            if (count == 3) {
                cout << "conjugate converged for " << i << " steps" << endl;
                break;
            }
        } else {
            count = 0;
        }
        for (int j = 0; j < m; j++) {
            Y[j] = X[j];
            X[j] = 0.0;
            r_o[j] = r_n[j];
        }

        if (i == 10000) {
            cout << "infinite loop" << endl;
            break;
        }
        i++;
    }
}

double calc_abs_err_1(double *A, double *B, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += abs(A[i] - B[i]);
    }
    return result;
}

double calc_abs_err_2(double *A, double *B, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += pow(abs(A[i] - B[i]), 2);
    }
    return sqrt(result);
}

double calc_abs_err_cheb(double *A, double *B, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        if (abs(A[i] - B[i]) > result) {
            result = abs(A[i] - B[i]);
        }
    }
    return result;
}

double calc_rel_err_1(double *A, double *B, int n) {
    double eps = 0.000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        num += (abs(A[i] - B[i]));
    }
    for (int i = 0; i < n; i++) {
        den += abs(A[i]);
    }
    if (abs(den - 0.0) < eps) {
        den = 1.0;
    }
    return num / den;
}

double calc_rel_err_2(double *A, double *B, int n) {
    double eps = 0.000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        num += (pow(abs(A[i] - B[i]), 2));
    }
    for (int i = 0; i < n; i++) {
        den += pow(A[i], 2);
    }
    if (abs(den - 0.0) < eps) {
        den = 1.0;
    }
    return sqrt(num) / sqrt(den);
}

double calc_rel_err_cheb(double *A, double *B, int n) {
    double eps = 0.000001;
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; i++) {
        if (abs(A[i] - B[i]) > num) {
            num = (abs(A[i] - B[i]));
        }
    }
    for (int i = 0; i < n; i++) {
        if (abs(A[i]) > den) {
            den = (abs(A[i]));
        }
    }
    if (abs(den - 0.0) < eps) {
        den = 1.0;
    }
    return num / den;
}

double opt_w(double **A, double *Y, double *X, double *B, int m, int n, double eps) {
    double X0[m];
    copy_vec(Y, X0, m);
    int best_i = 10000;
    int cur;
    double best_w;
    for (int i = 1; i < 99; i++) {
        double w = 1 + 0.01 * i;
        cur = over_rel(A, X0, X, B, w, m, n, eps, true);
        if (cur < best_i) {
            best_i = cur;
            best_w = w;
        }
        copy_vec(Y, X0, m);
    }
    return best_w;
}

void
write_e(double *X_e, double *X_g, double *X_LU, double *X_hol, double *X_rel, double *X_conj, int m, string filename) {
    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);
    fout << setw(15) << left << "error"
         << setw(15) << left << "gauss"
         << setw(15) << left << "LU"
         << setw(15) << left << "cholesky"
         << setw(15) << left << "relax"
         << setw(15) << left << "conjugate"
         << endl;
    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1(X_e, X_g, m)
         << setw(15) << left << calc_abs_err_1(X_e, X_LU, m)
         << setw(15) << left << calc_abs_err_1(X_e, X_hol, m)
         << setw(15) << left << calc_abs_err_1(X_e, X_rel, m)
         << setw(15) << left << calc_abs_err_1(X_e, X_conj, m)
         << endl;
    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(X_e, X_g, m)
         << setw(15) << left << calc_abs_err_2(X_e, X_LU, m)
         << setw(15) << left << calc_abs_err_2(X_e, X_hol, m)
         << setw(15) << left << calc_abs_err_2(X_e, X_rel, m)
         << setw(15) << left << calc_abs_err_2(X_e, X_conj, m)
         << endl;
    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(X_e, X_g, m)
         << setw(15) << left << calc_abs_err_cheb(X_e, X_LU, m)
         << setw(15) << left << calc_abs_err_cheb(X_e, X_hol, m)
         << setw(15) << left << calc_abs_err_cheb(X_e, X_rel, m)
         << setw(15) << left << calc_abs_err_cheb(X_e, X_conj, m)
         << endl;

    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(X_e, X_g, m)
         << setw(15) << left << calc_rel_err_1(X_e, X_LU, m)
         << setw(15) << left << calc_rel_err_1(X_e, X_hol, m)
         << setw(15) << left << calc_rel_err_1(X_e, X_rel, m)
         << setw(15) << left << calc_rel_err_1(X_e, X_conj, m)
         << endl;

    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(X_e, X_g, m)
         << setw(15) << left << calc_rel_err_2(X_e, X_LU, m)
         << setw(15) << left << calc_rel_err_2(X_e, X_hol, m)
         << setw(15) << left << calc_rel_err_2(X_e, X_rel, m)
         << setw(15) << left << calc_rel_err_2(X_e, X_conj, m)
         << endl;

    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(X_e, X_g, m)
         << setw(15) << left << calc_rel_err_cheb(X_e, X_LU, m)
         << setw(15) << left << calc_rel_err_cheb(X_e, X_hol, m)
         << setw(15) << left << calc_rel_err_cheb(X_e, X_rel, m)
         << setw(15) << left << calc_rel_err_cheb(X_e, X_conj, m)
         << endl;
    fout.close();
}

int main() {
    srand((unsigned int) time(nullptr));
    function<double(double)> y = func;

    double a = -3.14;
    double b = 3.14;

    double iter_par = 0.000001;

    int k = 3;
    double w = 1.1;
    int n = 5;
    int m = (n - 1) * k + 1;
    int l = 13;
    double h = (b - a) / (k * (n - 1));

    double **args = new double *[k];
    for (int i = 0; i < k; i++) {
        args[i] = new double[n];
    }

    double **random_args = new double *[k];
    for (int i = 0; i < k; i++) {
        random_args[i] = new double[l];
    }

    double **A = new double *[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[m];
    }

    double **L = new double *[n];
    for (int i = 0; i < n; i++) {
        L[i] = new double[m];
    }

    double **H = new double *[n];
    for (int i = 0; i < n; i++) {
        H[i] = new double[m];
    }

    double **H_t = new double *[n];
    for (int i = 0; i < n; i++) {
        H_t[i] = new double[m];
    }

    double **U = new double *[n];
    for (int i = 0; i < n; i++) {
        U[i] = new double[m];
    }

    double **A_g = new double *[2 * n - 1];
    for (int i = 0; i < 2 * n - 1; i++) {
        A_g[i] = new double[m];

    }

    double B[m];

    double X_e[m];
    double X_g[m];
    for (int i = 0; i < m; i++) {
        X_g[i] = 0.0;
    }
    double X_LU[m];
    for (int i = 0; i < m; i++) {
        X_LU[i] = 0.0;
    }
    double X_hol[m];
    for (int i = 0; i < m; i++) {
        X_hol[i] = 0.0;
    }
    double X_rel[m];
    for (int i = 0; i < m; i++) {
        X_rel[i] = 0.0;
    }
    double X_conj[m];
    for (int i = 0; i < m; i++) {
        X_conj[i] = 0.0;
    }

    double Y_rel[m];
    for (int i = 0; i < m; i++) {
        Y_rel[i] = 0.0;
    }

    double Y_conj[m];
    for (int i = 0; i < m; i++) {
        Y_conj[i] = 0.0;
    }

    double B_g[m];
    double B_LU[m];
    double B_hol[m];

    double AX_g[m];
    double AX_LU[m];
    double AX_hol[m];
    double AX_rel[m];
    double AX_conj[m];

    double resi_g[m];
    double resi_LU[m];
    double resi_hol[m];
    double resi_rel[m];
    double resi_conj[m];

    MatrixXd A_E(m, m);
    A_E.setZero();

    VectorXd B_E(m);
    VectorXd X_E(m);
    B_E.setZero();

    f_args(args, k, n, a, h);

    f_r_args(args, random_args, k, n, l);

    f_A(A, random_args, args, l, k, n);

    f_A_G(A, A_g, n, m);

    A_E = f_A_E(A_g, m, n);

    f_B(B, y, random_args, args, l, m, k, n);
    copy_vec(B, B_g, m);
    copy_vec(B, B_LU, m);
    copy_vec(B, B_hol, m);
    B_E = f_B_E(B, m);

    f_A_G(A, A_g, n, m);
    gauss(A_g, B_g, X_g, m, n);

    prep_L(L, m);
    LU(A, L, U, m, n);
    gauss_LU(L, U, B_LU, X_LU, m, n);

    holec(A, H, m, n);
    copy_mat(H, H_t, m, n);
    gauss_H(H, H_t, B_hol, X_hol, m, n);

    w = opt_w(A, Y_rel, X_rel, B, m, n, iter_par);
    over_rel(A, Y_rel, X_rel, B, w, m, n, iter_par, false);

    conjug(A, Y_conj, X_conj, B, m, n, iter_par);

    X_E = A_E.colPivHouseholderQr().solve(B_E);
    for (int i = 0; i < m; i++) {
        X_e[i] = X_E(i);
    }

    write_e(X_e, X_g, X_LU, X_hol, X_rel, X_conj, m, "Xerrors");

    mat_vec(A, X_g, AX_g, m, n);
    mat_vec(A, X_LU, AX_LU, m, n);
    mat_vec(A, X_hol, AX_hol, m, n);
    mat_vec(A, X_rel, AX_rel, m, n);
    mat_vec(A, X_conj, AX_conj, m, n);

    write_e(B, AX_g, AX_LU, AX_hol, AX_rel, AX_conj, m, "Residuals");

    return 0;
}
