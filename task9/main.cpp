#include <iostream>
#include <cmath>
#include <fstream>

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


void f_diff(double *X, double *Y, double h, int n) {

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

void write_f(double *Y1_ra, double *X, int n) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n; i++) {
        fout << X[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n; i++) {
        fout << Y1_ra[i] << " ";
    }


    fout.close();
}

double f1(double x, double y1, double y2) {
    return y2;
}

double f2(double x, double y1, double y2) {
    return y2 * exp(x) + y1 * x * x + sin(x);
}

double g(double x, double y1, double y2, double Y, double U) {
    return 2 * y1 * x * x * Y + 2 * y2 * exp(x) * U;
//    return y2 * exp(x) + y1 * x * x + sin(x);
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

void ru_ku4_g(double *Y1, double *Y2, double *X, double *y, double *u,double* y2,double*u2, int n, double h) {

    for (int i = 1; i < n; i++) {
        double k1 = h * f1(X[i - 1], Y1[i - 1], Y2[i - 1]);
        double k2 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k1 / 2, Y2[i - 1] + k1 / 2);
        double k3 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k2 / 2, Y2[i - 1] + k2 / 2);
        double k4 = h * f1(X[i - 1] + h / 2, Y1[i - 1] + k3, Y2[i - 1] + k3);
        Y1[i] = Y1[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        k1 = h * g(X[i - 1], y[i - 1], u[i - 1], Y1[i - 1], Y2[i - 1]);
        k2 = h * g(X[i - 1] + h / 2, y2[i - 1], u2[i - 1], Y1[i - 1] + k1 / 2, Y2[i - 1] + k1 / 2);
        k3 = h * g(X[i - 1] + h / 2, y2[i - 1], u2[i - 1], Y1[i - 1] + k2 / 2, Y2[i - 1] + k2 / 2);
        k4 = h * g(X[i - 1] + h / 2, y2[i - 1], u2[i - 1], Y1[i - 1] + k3, Y2[i - 1] + k3);
        Y2[i] = Y2[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }
}

void shoot(double *Y, double *X, int n, double h) {
    double eps = 0.001;
    double t1 = -5.0;
    double t2 = 5.0;
    double yb = Y[n - 1];
    double Y1[n];
    double Y2[n];
    double U1[n];
    double U[n];
    double U2[n];
    double res = 0.0;

    U1[0] = t1;
    ru_ku4(Y1, U1, X, n, h);

    U2[0] = t2;
    ru_ku4(Y2, U2, X, n, h);

    double l_v = Y1[n - 1] - yb;
    double r_v = Y2[n - 1] - yb;

    int i = 1;
    while (true) {
        U[0] = (t1 + t2) / 2.0;
        ru_ku4(Y, U, X, n, h);

        if (abs(res - (Y[n - 1] - yb)) < eps) {
            cout << i << " iter" << endl;
            cout << t1 << endl;
            break;
        }
        res = Y[n - 1] - yb;
        if (res < 0.0) {
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


//    int m =10;
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
    double eps = 0.0001;
    double s = 0;
    double yb = Y[n - 1];
    double ya = Y[0];
    double fi_s;
    double fi_s_h;
    double h_d = 0.01;
    double U1[n];
    double Y2[n];

    double U4[n];
    double Y4[n];

    double U2[n];

    Y2[0] = 0.0;
    U2[0] = 1.0;

    double Y3[2*n-1];
    double Xh2[2*n-1];
    f_even_args(Xh2,2*n-1,0.0,1.0);//todo
    double U3[2*n-1];
    double h2 = h/2.0;
    Y3[0] = ya;

    int i = 1;
    double s_n;
    double dfi;
    while (true) {
        if (sec) {
            U1[0] = s + h_d;
            ru_ku4(Y, U1, X, n, h);
            fi_s_h = Y[n - 1] - yb;
            U1[0] = s;
            ru_ku4(Y, U1, X, n, h);
            fi_s = Y[n - 1] - yb;
            dfi = ((fi_s_h - fi_s) / h_d);
        } else {
            U3[0] = s;
            ru_ku4(Y3, U3, Xh2, 2*n-1, h2);

            U1[0] = s;
            ru_ku4(Y, U1, X, n, h);
            for (int j = 0; j < n; j++) {
                Y[i] = Y3[2*i];
            }
            for (int j = 0; j < n; j++) {
                U1[i] = U3[2*i];
            }

            for (int j = 0; j < n-1; j++) {
                Y4[i] = Y3[2*i+1];
            }
            for (int j = 0; j < n-1; j++) {
                U4[i] = U3[2*i+1];
            }

            fi_s = Y[n - 1] - yb;
            ru_ku4_g(Y2,U2,X,Y,U1,Y4,U4,n,h);
            dfi = Y2[n-1];
        }
        s_n = s - fi_s / dfi;
        if (abs(s_n - s) < eps) {
            cout << i << " iter" << endl;
            cout << s_n << endl;
            break;
        }
        s = s_n;
        if (i == 5000) {
            cout << "INFINITE LOOP" << endl;
            break;
        }
        i++;
    }
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double ya = 0.0;
    double yb = 0.0;
    int n1 = 107;
    double h1 = (b - a) / (n1 - 1);

    double X_h[n1];

    double Y_ra_h[n1];
    Y_ra_h[0] = ya;
    Y_ra_h[n1 - 1] = yb;

    double Y_shoo_h[n1];
    Y_shoo_h[0] = ya;
    Y_shoo_h[n1 - 1] = yb;

    double Y_sec_h[n1];
    Y_sec_h[0] = ya;
    Y_sec_h[n1 - 1] = yb;


    f_even_args(X_h, n1, a, b);

//    f_diff(X_h, Y_ra_h, h1, n1);
//    shoot(Y_shoo_h, X_h, n1, h1);
    newton(Y_sec_h, X_h, n1, h1, false);

    write_f(Y_sec_h, X_h, n1);


    return 0;
}
