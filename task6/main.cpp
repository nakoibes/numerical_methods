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

VectorXd f_B_E(double* B,int m){
    VectorXd B_E(m);
    for(int i =0;i<m;i++){
        B_E(i) = B[i];
    }
    return B_E;
}

void gauss(double **A, double *B,double *X, int m, int n) {
    double eps = 0.0000001;
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
    for (int al = m - 1; al > 0; al--) {
        for (int be = 1; be < n; be++) {
            int i = be;
            int j = al - be;
            if (abs(A[i][j]) > eps) {
                double koef = A[i][j] / A[0][al];
                A[i][j] = 0;
//                A[i][j] = A[i][j] - A[0][al] * koef;
                B[al - be] = B[al - be] - B[al] * koef;
            } else {
                continue;
            }
        }
    }

    for(int i =0;i<m;i++){
        X[i] = B[i]/A[0][i];
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

MatrixXd f_A_E(double** A,int m,int n){
    MatrixXd A_E(m,m);
    A_E.setZero();
    for (int i = 0; i < 2 * n - 1; i++) {
        for (int j = 0; j < m; j++) {
//            cout << A_E << endl;
            if (i < n) {
                if (i + j < m) {
                    A_E(j,i + j) = A[i][j];
                } else {
                    continue;
                }
            } else {
                if (i + j - n+1 < m) {
                    A_E(i + j - n+1,j) = A[i][j];
                } else {
                    continue;
                }
            }
        }
    }
    return A_E;
}

void prep_L(double** L,int n,int m){
    for(int i =0;i<m;i++){
        L[0][i] = 1;
    }
}

int main() {

    function<double(double)> y = func;

    double a = -3.14;
    double b = 3.14;

    int k = 3;
    int n = 4;
    int m = (n - 1) * k + 1;
    int l = 15;
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


    double **U = new double *[n];
    for (int i = 0; i < n; i++) {
        U[i] = new double[m];
    }


    double **A_G = new double *[2 * n - 1];
    for (int i = 0; i < 2 * n - 1; i++) {
        A_G[i] = new double[m];

    }

    double **A_test = new double *[m];
    for (int i = 0; i < m; i++) {
        A_test[i] = new double[m];

    }

    double B[m];
    for (int i = 0; i < m; i++) {
    }

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
    A_E = f_A_E(A_G,m,n);
    cout << A_E << endl;
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


    for (int i = 0; i < 2 * n - 1; i++) {
        for (int j = 0; j < m; j++) {
            cout << setw(10) << A_G[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    f_B(B, y, random_args, args, l, m, k, n);
    B_E = f_B_E(B,m);

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
//
//    double **D = new double *[m];
//    for (int i = 0; i < m; i++) {
//        D[i] = new double[m];
//    }
//
//    D[0][0] = 1.0;
//    D[0][1] = 5.0;
//    D[0][2] = 4.0;
//    D[1][0] = 2.0;
//    D[1][1] = 6.0;
//    D[1][2] = 0.0;
//    D[2][0] = 2.0;
//    D[2][1] = 6.0;
//    D[2][2] = 0.0;


//    double **A_G1 = new double *[2 * n - 1];
//    for (int i = 0; i < 2 * n - 1; i++) {
//        A_G1[i] = new double[m];
//
//    }

//    double B1[3]{1, 1, 1};

//    f_A_G(A, A_G, n, m);

    gauss(A_G, B,X, m, n);


//    for (int i = 0; i < 2 * n - 1; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A_G[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A_test[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 2 * n - 1; i++) {
        for (int j = 0; j < m; j++) {
            if (i <= n) {
                if (i + j < m) {
                    A_test[j][i + j] = A_G[i][j];
                } else {
                    continue;
                }
            } else {
                if (i + j - 2 < m) {
                    A_test[i + j - 2][j] = A_G[i][j];
                } else {
                    continue;
                }
            }
        }
    }



//    for (int i = 0; i < m; i++) {
//        for (int j = 0; j < m; j++) {
//            cout << setw(10) << A_test[i][j] << " ";
//        }
//        cout << endl;
//    }
//    cout << endl;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            cout << setw(10) << A_test[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    for (int i = 0; i < m; i++) {
        cout << X[i] << " ";
    }

    X_E = A_E.colPivHouseholderQr().solve(B_E);
    cout << endl;
    for (int i = 0; i < m; i++) {
        cout << X_E(i) << " ";
    }

    return 0;
}
