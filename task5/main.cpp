#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>

using namespace std;

double f(double x) {
    return pow(x, 10) + 1;
}

void f_even(double *x, int n, double a, double b) {
    x[0] = a;
    x[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        x[i] = x[i - 1] + step;
    }
}

double i_func(double x) {
    return 1.0 / 11.0 * pow(x, 11) + x;
}

double rectangle(double *x, double h, int n) {
    double sum = 0.0;
    for (int i = 0; i < n - 1; i++) {
        sum += f(x[i] + h / 2.0);
    }
    return sum * h;
}

double trapeze(double *x, int n, double h) {
    double sum = f(x[0]) / 2.0 + f(x[n - 1]) / 2.0;
    for (int i = 1; i < n - 1; i++) {
        sum += (f(x[i]));
    }
    return sum * h;
}

double simps(double *x, int n, double h) {
    double sum = f(x[0]) + f(x[n - 1]);
    for (int i = 1; i < n - 1; i++) {
        sum += 2.0 * f(x[i]);
    }
    for (int i = 0; i < n - 1; i++) {
        sum += 4.0 * f(x[i] + h / 2.0);
    }
    return sum * h / 6.0;
}

double newton_cotes(double *x, int n, double len) {
    double h = len / 4.0;
    double sum = 0.0;
    for (int i = 0; i < n - 1; i++) {
        sum += (7.0 * f(x[i]) + 32.0 * f(x[i] + h) + 12.0 * f(x[i] + 2.0 * h) + 32.0 * f(x[i] + 3.0 * h) +
                7.0 * f(x[i + 1]));
    }
    return sum * 2.0 * h / 45.0;
}

double gauss(double *x, int n, double h) {
    double sum = 0.0;
    double c1 = 5.0 / 9.0;
    double c2 = 8.0 / 9.0;
    double c3 = 5.0 / 9.0;
    double ksi1 = -sqrt(3.0 / 5.0);
    double ksi3 = sqrt(3.0 / 5.0);
    for (int i = 0; i < n - 1; i++) {
        sum += (c1 * f(h / 2.0 * ksi1 + (x[i] + x[i + 1]) / 2.0) + c2 * f((x[i] + x[i + 1]) / 2.0) +
                c3 * f(h / 2.0 * ksi3 + (x[i] + x[i + 1]) / 2.0));
    }
    return sum * h / 2.0;
}

void write_errs(double rectangle_h, double trapeze_h, double simps_h, double newton_cotes_h, double gauss_h,
                double rectangle_h2, double trapeze_h2, double simps_h2, double newton_cotes_h2, double gauss_h2,
                double integral) {
    ofstream fout("errs.txt");
    fout << fixed;
    fout.precision(6);
    fout << setw(15) << left << scientific << "scheme"
         << setw(15) << left << "h_integral"
         << setw(15) << left << "h/2_integral"
         << setw(15) << left << "error h"
         << setw(15) << left << "error h/2"
         << setw(15) << left << "integral"
         << endl;
    fout << setw(15) << left << "rectangle"
         << setw(15) << left << rectangle_h
         << setw(15) << left << rectangle_h2
         << setw(15) << left << abs(rectangle_h - integral) / integral
         << setw(15) << left << abs(rectangle_h2 - integral) / integral
         << setw(15) << left << integral
         << endl;
    fout << setw(15) << left << "trapeze"
         << setw(15) << left << trapeze_h
         << setw(15) << left << trapeze_h2
         << setw(15) << left << abs(trapeze_h - integral) / integral
         << setw(15) << left << abs(trapeze_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "simpson"
         << setw(15) << left << simps_h
         << setw(15) << left << simps_h2
         << setw(15) << left << abs(simps_h - integral) / integral
         << setw(15) << left << abs(simps_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "newton_cotes"
         << setw(15) << left << newton_cotes_h
         << setw(15) << left << newton_cotes_h2
         << setw(15) << left << abs(newton_cotes_h - integral) / integral
         << setw(15) << left << abs(newton_cotes_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "gauss"
         << setw(15) << left << gauss_h
         << setw(15) << left << gauss_h2
         << setw(15) << left << abs(gauss_h - integral) / integral
         << setw(15) << left << abs(gauss_h2 - integral) / integral
         << endl;
}

int main() {

    double a = 0.0;
    double b = 3.0;

    int k = 5;
    int m1 = k + 1;
    int m2 = 2 * k + 1;

    double h1 = (b - a) / k;
    double h2 = (b - a) / (2.0 * k);


    double x_h[m1];
    double x_h2[m2];

    f_even(x_h, m1, a, b);
    f_even(x_h2, m2, a, b);

    double rectangle_h = rectangle(x_h, h1, m1);
    double trapeze_h = trapeze(x_h, m1, h1);
    double simps_h = simps(x_h, m1, h1);
    double newton_cotes_h = newton_cotes(x_h, m1, h1);
    double gauss_h = gauss(x_h, m1, h1);

    double rectangle_h2 = rectangle(x_h2, h2, m2);
    double trapeze_h2 = trapeze(x_h2, m2, h2);
    double simps_h2 = simps(x_h2, m2, h2);
    double newton_cotes_h2 = newton_cotes(x_h2, m2, h2);
    double gauss_h2 = gauss(x_h2, m2, h2);

    double integral = i_func(b) - i_func(a);

    write_errs(rectangle_h, trapeze_h, simps_h, newton_cotes_h, gauss_h,
               rectangle_h2, trapeze_h2, simps_h2, newton_cotes_h2, gauss_h2, integral);


    return 0;
}
