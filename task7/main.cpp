#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double f1(double x, double y) {
    return y - x * x;
}

double f2(double x, double y) {
    return x - exp(y) + 2;
}

double fi1(double y) {
    return -sqrt(y);
}

double fi2(double x) {
    return log(x + 2);
}

void simple(double *X0, double *X, double eps, bool left) {
    int i = 1;
    double x_p = X0[0];
    double y_p = X0[1];
    while (true) {
        double x_c;
        if (left) {
            x_c = fi1(y_p);
        } else {
            x_c = -fi1(y_p);
        };
        double y_c = fi2(x_p);
        if (abs(x_p - x_c) + abs(y_p - y_c) < eps) {
            X[0] = x_c;
            X[1] = y_c;
            cout << "simple converged for " << i << " steps, start: " << X0[0] << " ; " << X0[1] << endl;
            break;
        }
        if (i == 10000) {
            cout << "INFINITE LOOP" << endl;
            break;
        }
        x_p = x_c;
        y_p = y_c;
        i++;
    }
}

double df1dx(double x, double y) {
    return -2 * x;
}

double df1dy(double x, double y) {
    return 1;
}

double df2dx(double x, double y) {
    return 1;
}

double df2dy(double x, double y) {
    return -exp(y);
}

double df1dx_d(double x, double y, double h) {
    return (f1(x + h, y) - f1(x, y)) / h;
}

double df1dy_d(double x, double y, double h) {
    return (f1(x, y + h) - f1(x, y)) / h;
}

double df2dx_d(double x, double y, double h) {
    return (f2(x + h, y) - f1(x, y)) / h;
}

double df2dy_d(double x, double y, double h) {
    return (f2(x, y + h) - f1(x, y)) / h;
}

double newton(double *X0, double *X, double eps, bool con, bool dis, string o_str,bool sile) {
    int i = 1;
    double h = 0.01;
    double x_p = X0[0];
    double y_p = X0[1];
    double d1;
    double d2;
    double d3;
    double d4;
    if (!dis) {
        d1 = df1dx(x_p, y_p);
        d2 = df1dy(x_p, y_p);
        d3 = df2dx(x_p, y_p);
        d4 = df2dy(x_p, y_p);
    } else {
        d1 = df1dx_d(x_p, y_p, h);
        d2 = df1dy_d(x_p, y_p, h);
        d3 = df2dx_d(x_p, y_p, h);
        d4 = df2dy_d(x_p, y_p, h);
    }

    while (true) {
        if (!con) {
            if (!dis) {
                d1 = df1dx(x_p, y_p);
                d2 = df1dy(x_p, y_p);
                d3 = df2dx(x_p, y_p);
                d4 = df2dy(x_p, y_p);
            } else {
                d1 = df1dx_d(x_p, y_p, h);
                d2 = df1dy_d(x_p, y_p, h);
                d3 = df2dx_d(x_p, y_p, h);
                d4 = df2dy_d(x_p, y_p, h);
            }
        }
        double x_c = x_p + (-f1(x_p, y_p) * d4 + f2(x_p, y_p) * d2) /
                           (d1 * d4 - d2 * d3);
        double y_c = y_p + (-d1 * f2(x_p, y_p) + f1(x_p, y_p) * d3) /
                           (d1 * d4 - d2 * d3);

        if (abs(x_p - x_c) + abs(y_p - y_c) < eps) {
            X[0] = x_c;
            X[1] = y_c;
            if(!sile) {
                cout << o_str << " converged for " << i << " steps, start: " << X0[0] << " ; " << X0[1] << endl;
            }
            return i;
        }
        if (i == 5000) {
            if(!sile) {
                cout << "INFINITE LOOP" << endl << endl;
            }
            return 5000;
        }
        x_p = x_c;
        y_p = y_c;
        i++;

    }
}

void root_a(double** XY,double* X01,double* X,double eps,int k,double s_x,double s_y){
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            X01[0] = s_x + 1.0/k*i;
            X01[1] = s_y + 1.0/k*j;
            XY[i][j] = newton(X01, X, eps, false, true, "newton_d",true);
        }
    }
}
void write_d(double **XY1,double **XY2,int k) {
    ofstream fout("dots.txt");
    fout << k << endl;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            fout << XY1[i][j] << " ";
        }
        fout << endl;
    }
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            fout << XY2[i][j] << " ";
        }
        fout << endl;
    }
    fout.close();
}
int main() {
    double eps = 0.000001;
    double X01[2] = {-0.5, 0.3};
    double X[2];

    int k=100;

    double **XY1 = new double *[k];
    for (int i = 0; i < k; i++) {
        XY1[i] = new double[k];
    }
    double **XY2 = new double *[k];
    for (int i = 0; i < k; i++) {
        XY2[i] = new double[k];
    }


//    system("python3 vis.py");

    simple(X01, X, eps, true);
    cout << X[0] << endl << X[1] << endl;
    X01[0] = 1.1;
    X01[1] = 1.1;
    simple(X01, X, eps, false);
    cout << X[0] << endl << X[1] << endl;

    X01[0] = -0.5;
    X01[1] = 0.3;
    newton(X01, X, eps, false, false, "newton",false);
    cout << X[0] << endl << X[1] << endl;
    X01[0] = 1.1;
    X01[1] = 1.1;
    newton(X01, X, eps, false, false, "newton",false);
    cout << X[0] << endl << X[1] << endl;

    X01[0] = -0.5;
    X01[1] = 0.3;
    newton(X01, X, eps, true, false, "newton_c",false);
    cout << X[0] << endl << X[1] << endl;
    X01[0] = 1.1;
    X01[1] = 1.1;
    newton(X01, X, eps, true, false, "newton_c",false);
    cout << X[0] << endl << X[1] << endl;

    X01[0] = -0.5;
    X01[1] = 0.3;
    newton(X01, X, eps, false, true, "newton_d",false);
    cout << X[0] << endl << X[1] << endl;
    X01[0] = 1.1;
    X01[1] = 1.1;
    newton(X01, X, eps, false, true, "newton_d",false);
    cout << X[0] << endl << X[1] << endl;

    root_a(XY1,X01,X,eps,k,-1.1,-0.15);
    root_a(XY2,X01,X,eps,k,0.6,0.6);
    write_d(XY1,XY2,k);

    for (int i = 0; i < k; i++) {
        delete[] XY1[i];
        delete[] XY2[i];
    }
    delete[] XY1;
    delete[] XY2;

    return 0;
}
