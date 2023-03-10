#include <iostream>
#include <functional>
#include <cmath>
#include <iomanip>
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
//                    A[i * (n - 1) + j2][i * (n - 1) + j1] = A[i * (n - 1) + j1][i * (n - 1) + j2];
                }
            }
        }
    }
}

void f_A_test(double **A, double **random_args, double **args, int l, int k, int n) {
    for (int i = 0; i < k; i++) {
        for (int j1 = 0; j1 < n; j1++) {
            for (int j2 = 0; j2 < n; j2++) {
                if ((i * (n - 1) + j1) <= (i * (n - 1) + j2)) {
                    for (int q = 0; q < l; q++) {
//                        int i_h = i * (n - 1) + j1;
//                        int j_h = i * (n - 1) + j2;
                        A[i * (n - 1) + j1][i * (n - 1) + j2] +=
                                fi(random_args[i][q], args[i], j1, n) * fi(random_args[i][q], args[i], j2, n);
                    }
                    A[i * (n - 1) + j2][i * (n - 1) + j1] = A[i * (n - 1) + j1][i * (n - 1) + j2];
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





//for square matrix
//    double eps = 0.0000001;
//    for (int i = 0; i < m - 1; i++) {
//        for (int j = 0; j < n; j++) {
//            if (i + j + 1 < m) {
//                if (abs(A[i + j + 1][i]) > eps) {
//                    double koef = 1.0 / A[i + j + 1][i];
//                    A[i + j + 1][i] = 0;
//                    for (int k = 1; k < n; k++) {
//                        if (i + k < m) {
//                            A[i + j + 1][i + k] = A[i + j + 1][i + k] - A[i + j][i + k] * koef;
//                        }
//                    }
//                } else {
//                    break;
//                }
//            } else {
//                break;
//            }
//        }
//    }
}

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

void m_v_p(double **A, double *B, double *C, int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            C[i] += A[i][j] * B[j];
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
                B[i + j] -= B[i] * A[j][i]/A[0][i];
                A[j][i] = 0.0;
            }
        }
    }
}

//void g_r_LU(double **A, double *B, int m, int n){
//    for (int i = m-1; i > 0; i--) {
//        for (int j = 1; j < n; j++) {
//
//                B[j] -= B[i] * A[j][i];
//                A[j][i] = 0.0;
//
//        }
//    }
//}

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
        Y[i] = B[i]/H[0][i];
    }

    gauss_r(H_t, Y, m, n, eps);
    for (int i = 0; i < m; i++) {
        X[i] = Y[i] / H_t[0][i];
    }
}

MatrixXd f_A_E(double **A, int m, int n) {
    MatrixXd A_E(m, m);
    A_E.setZero();
    for (int i = 0; i < 2 * n - 1; i++) {
        for (int j = 0; j < m; j++) {
//            cout << A_E << endl;
            if (i < n) {
                if (i + j < m) {
                    A_E(j, i + j) = A[i][j];
                } else {
                    continue;
                }
            } else {
                if (i + j - n + 1 < m) {
                    A_E(i + j - n + 1, j) = A[i][j];
                } else {
                    continue;
                }
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

void prep_L_test(double **L, int m) {
    for (int i = 0; i < m; i++) {
        L[i][i] = 1;
    }
}

void LU_test(double **A, double **L, double **U, int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (i <= j) {
                U[i][j] = A[i][j];
                for (int k = 0; k < i; k++) {
                    U[i][j] -= (L[i][k] * U[k][j]);
                }
            } else {
                L[i][j] = A[i][j];
                for (int k = 0; k < j; k++) {
                    L[i][j] -= (L[i][k] * U[k][j]);
                }
                L[i][j] = L[i][j] / U[j][j];
            }
        }
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


void mat_pro(double **A, double **B, double **C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void mat_dif(double **A, double **B, double **C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] - B[i][j];
        }
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

int main() {

    function<double(double)> y = func;

    double a = -3.14;
    double b = 3.14;

    int k = 3;
    int n = 3;
    int m = (n - 1) * k + 1;
    int l = 13;
    double h = (b - a) / (k * (n - 1));
    int le = 2 * n - 1;

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

    double **L_H = new double *[n];
    for (int i = 0; i < n; i++) {
        L_H[i] = new double[m];
    }

    double **L_test = new double *[m];
    for (int i = 0; i < m; i++) {
        L_test[i] = new double[m];
    }

    double **H_test = new double *[m];
    for (int i = 0; i < m; i++) {
        H_test[i] = new double[m];
    }

    double **H_test_t = new double *[m];
    for (int i = 0; i < m; i++) {
        H_test_t[i] = new double[m];
    }


    double **U = new double *[n];
    for (int i = 0; i < n; i++) {
        U[i] = new double[m];
    }

    double **U_test = new double *[m];
    for (int i = 0; i < m; i++) {
        U_test[i] = new double[m];
    }


    double **A_G = new double *[2 * n - 1];
    for (int i = 0; i < 2 * n - 1; i++) {
        A_G[i] = new double[m];

    }

    double **A_test = new double *[m];
    for (int i = 0; i < m; i++) {
        A_test[i] = new double[m];

    }

    double **C = new double *[m];
    for (int i = 0; i < m; i++) {
        C[i] = new double[m];

    }

    double **D = new double *[m];
    for (int i = 0; i < m; i++) {
        D[i] = new double[m];

    }

    double B[m];
//    for (int i = 0; i < m; i++) {
//    }

    double X[m];
    for (int i = 0; i < m; i++) {
        X[i] = 0.0;
    }


    MatrixXd A_E(m, m);
    A_E.setZero();

    VectorXd B_E(m);
    VectorXd X_E(m);
    B_E.setZero();

    f_args(args, k, n, a, h);

    f_r_args(args, random_args, k, n, l);

    f_A(A, random_args, args, l, k, n);

    f_A_G(A, A_G, n, m);

    A_E = f_A_E(A_G, m, n);
//    cout << A_E << endl;


    f_A_test(A_test, random_args, args, l, k, n);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            cout << setw(10) << A_test[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    //half
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;


//    for (int i = 0; i < 2 * n - 1; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A_G[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

    f_B(B, y, random_args, args, l, m, k, n);
    B_E = f_B_E(B, m);

//    double **C = new double *[5];
//    for (int i = 0; i < 5; i++) {
//        C[i] = new double[5];
//    }
//
//    C[0][0] = 1.0;
//    C[0][1] = 2.0;
//    C[0][2] = 3.0;
//    C[0][3] = 0.0;
//    C[0][4] = 0.0;
//    C[1][0] = 7.0;
//    C[1][1] = 8.0;
//    C[1][2] = 9.0;
//    C[1][3] = 0.0;
//    C[1][4] = 0.0;
//    C[2][0] = 2.0;
//    C[2][1] = 3.0;
//    C[2][2] = 4.0;
//    C[2][3] = 5.0;
//    C[2][4] = 1.0;
//    C[3][0] = 0.0;
//    C[3][1] = 0.0;
//    C[3][2] = 4.0;
//    C[3][3] = 5.0;
//    C[3][4] = 1.0;
//    C[4][0] = 0.0;
//    C[4][1] = 0.0;
//    C[4][2] = 7.0;
//    C[4][3] = 7.0;
//    C[4][4] = 2.0;


//    n = 2;
//    m = 3;

//    double **P_test = new double *[3];
//    for (int i = 0; i < 3; i++) {
//        P_test[i] = new double[3];
//    }
//
//    P_test[0][0] = 1.0;
//    P_test[0][1] = 0.0;
//    P_test[0][2] = 2.0;
//    P_test[1][0] = 0.0;
//    P_test[1][1] = 1.0;
//    P_test[1][2] = 0.0;
//    P_test[2][0] = 0.0;
//    P_test[2][1] = 0.0;
//    P_test[2][2] = 4.0;
//
//    double B_test[3] = {1.0, 1.0, 1.0};
//    double C_test[3] = {0.0, 0.0, 0.0};


//    double **A_G1 = new double *[2 * n - 1];
//    for (int i = 0; i < 2 * n - 1; i++) {
//        A_G1[i] = new double[m];
//
//    }

//    double B1[3]{1, 1, 1};

//    f_A_G(A, A_G, n, m);

//    gauss(A_G, B, X, m, n);


    double **T = new double *[2];
    for (int i = 0; i < 2; i++) {
        T[i] = new double[3];
    }
    T[0][0] = 1.0;
    T[0][1] = 5;
    T[0][2] = 10;
    T[1][0] = 2;
    T[1][1] = 3;
    T[1][2] = 0;
//    T[2][0] = 5;
//    T[2][1] = -1;
//    T[2][2] = 5;

//    for (int i = 0; i < 2 * n - 1; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A_G[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;


//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            A_test[i][j] = 0.0;
//        }
//    }

//    for (int i = 0; i < 2 * n - 1; i++) {
//        for (int j = 0; j < m; j++) {
//            if (i <= n) {
//                if (i + j < m) {
//                    A_test[j][i + j] = A_G[i][j];
//                } else {
//                    continue;
//                }
//            } else {
//                if (i + j - 2 < m) {
//                    A_test[i + j - 2][j] = A_G[i][j];
//                } else {
//                    continue;
//                }
//            }
//        }
//    }



//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A_test[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A_test[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;





//    prep_L_test(L_test, m);
//    LU_test(A_test, L_test, U_test, m);

//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << L_test[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;
//
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << U_test[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

//    mat_pro(L_test, U_test, C, m);
//
//    mat_dif(A_test, C, D, m);

//    m_v_p(P_test,B_test,C_test,3);

//        for (int i = 0; i < 3; i++) {
//            cout << setw(10) << C_test[i] << " ";
//        }
//        cout << endl;

    prep_L(L, m);
    LU(A, L, U, m, n);
    gauss_LU(L, U, B, X, m, n);


//    double **T_t = new double *[2];
//    for (int i = 0; i < 2; i++) {
//        T_t[i] = new double[3];
//    }

// double B_test[3]={1,1,1};

//    holec(A, H, m, n);
////        for (int i = 0; i < 2; i++) {
////        for (int j = 0; j < 3; j++) {
////            H_t[i][j] = H[i][j];
////        }
////    }
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < m; j++) {
//            H_t[i][j] = H[i][j];
//        }
//    }
//    gauss_H(H, H_t, B, X, m, n);


    X_E = A_E.colPivHouseholderQr().solve(B_E);
    cout << endl;

    for (int i = 0; i < m; i++) {
        cout << X_E(i) << " ";
    }
    cout << endl;

    for (int i = 0; i < m; i++) {
        cout << X[i] << " ";
    }
    cout << endl;


//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << U[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            L_test[i][j] = 0;
//            U_test[i][j] = 0;
//        }
//    }

//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << U_test[i][j] << " ";
//        }
//        cout << endl;
//    }

//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < m; j++) {
//            if (i + j < m) {
//                H_test[i + j][j] = H[i][j];
//                H_test_t[j][i + j] = H_t[i][j];
//            }
//        }
//    }
//
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << H_test[i][j] << " ";
//        }
//        cout << endl;
//    }
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << H_test_t[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//
//    mat_pro(H_test, H_test_t, C, m);
//
//    mat_dif(A_test, C, D, m);
//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << D[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;


    return 0;
}
