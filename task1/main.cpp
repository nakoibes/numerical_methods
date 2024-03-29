#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double f(double x) {
    return sin(x);
}

double random(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

void f_i_x(double *x, int n, double a, double b, double eps) {
    x[0] = a;
    x[n - 1] = b;
    int i = 1;
    while (i < n - 1) {
        double value = random(a, b);
        bool good_v = true;
        for (int j = 0; j < i; j++) {
            if (abs(x[j] - value) < eps) {
                good_v = false;
                break;
            }
        }
        if (good_v) {
            x[i] = value;
            i++;
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
        double prod = 1;
        for (int j = 0; j < n_k; j++) {
            if (j != i) {
                prod *= (x_k[i] - x_k[j]);
            }
        }
        mult[i] = y_k[i] / prod;
    }
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n_k; k++) {
            double prod = 1.0;
            for (int j = 0; j < n_k; j++) {
                if (k != j) {
                    prod *= (x[i] - x_k[j]);
                }
            }
            y[i] += mult[k] * prod;
        }
    }
}

void write_d(double *x_v, double *y_v, double *x_k, double *y_k, double *func, int n_i, int n_v) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n_v; i++) {
        fout << x_v[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_v; i++) {
        fout << y_v[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_i; i++) {
        fout << x_k[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_i; i++) {
        fout << y_k[i] << " ";
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

int main() {
    srand((unsigned int) time(nullptr));

    double a = -3.14;
    double b = 3.14;
    double eps = 0.0001;
    int n_k = 20;
    int n_vis = 2000;

    double x_k[n_k];
    double y_k[n_k];
    double x_v[n_vis];
    double y_v[n_vis];

    double func[n_vis];

    f_i_x(x_k, n_k, a, b, eps);
    f_func(y_k, x_k, n_k);
    f_even(x_v, n_vis, a, b);
    interp(y_v, x_v, y_k, x_k, n_vis, n_k);

    f_func(func, x_v, n_vis);

    write_d(x_v, y_v, x_k, y_k, func, n_k, n_vis);

    system("python3 vis.py");
    return 0;
}
