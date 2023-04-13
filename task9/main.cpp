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
    A[0][0] = 1.0;
    A[0][1] = 0.0;
    A[n - 1][0] = 1.0;
    A[n - 1][2] = 0.0;
    for (int i = 1; i < n-1; i++) {
        A[i][2] = (2.0 - p(X[i]) * h) / (2 * h * h);
        A[i][0] = (-4.0 - q(X[i])) / (2 * h * h);
        A[i][1] = (2.0 - p(X[i]) * h) / (2 * h * h);;
    }
}

void prog(double** M,double* Y,double* F,int n){
    double A_p = -M[0][1]/M[0][0];
    double B_p = F[0]/M[0][0];
    double A_c;
    double B_c;
    for (int i = 1; i < n-1; i++) {
        double den = M[i][2]*A_p+M[i][0];
        A_c = -M[i][1]/den;
        B_c = (F[i]-M[i][2]*B_p)/den;
        Y[i] = (-B_c + F[i]/M[i][2])/(A_c+M[i][0]/M[i][2]);
        A_p = A_c;
        B_p = B_c;
    }
}

void formF(double *F,double* X, int n) {
    for (int i = 1; i < n - 1; i++) {
        F[i] = f(X[i]);
    }
}
void write_f(double *Y1_ra,  double *X, int n) {
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



    int main() {
        double a = 0.0;
        double b = 1.0;
        double ya = 0.0;
        double yb = 0.0;
        int n1 = 11;
        double h1 = (b - a) / (n1 - 1);

        double X_h[n1];

        double F[n1];
        F[0] = ya;
        F[n1 - 1] = yb;

        double Y1_ra_h[n1];
        Y1_ra_h[0] = ya;
        Y1_ra_h[n1 - 1] = yb;

        double **A = new double *[n1];
        for (int i = 0; i < n1; i++) {
            A[i] = new double[3];
        }

        f_even_args(X_h, n1, a, b);

        formA(A,X_h,h1,n1);

        formF(F,X_h,n1);

        prog(A,Y1_ra_h,F,n1);
        write_f(Y1_ra_h,X_h,n1);

        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < 3; j++) {
                cout << A[i][j] << " ";
            }
            cout << endl;
        }


        return 0;
    }
