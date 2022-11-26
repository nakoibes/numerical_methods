#include <iostream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>

using namespace std;

double func(double x) {
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

double i_func(double a, double b) {
    return -cos(b) + cos(a);
}

//double rectangle(double* args, function<double(double)> y, double h, int n){
//    double sum = y(args[0])/2+y(args[n-1])/2;
//    for(int i=1;i<n-1;i++){
//        sum+=y(args[i]);
//    }
//    return sum*h;
//}
double rectangle(double *args, function<double(double)> y, double h, int n) {
    double sum = 0;
    for (int i = 0; i < n - 1; i++) {
        sum += y(args[i] + h / 2);
    }
    return sum * h;
}

double trapeze(double *args, function<double(double)> y, int n, double h) {
    double sum = y(args[0]) / 2 + y(args[n - 1]) / 2;
    for (int i = 1; i < n - 2; i++) {
        sum += (y(args[i]));
    }
    return sum * h;
}

double simps(double *args, function<double(double)> y, int n, double len) {
    double sum = y(args[0]) + y(args[n - 1]);
    for (int i = 1; i < n - 1; i++) {
        sum += 2 * y(args[i]);
    }
    for (int i = 0; i < n - 1; i++) {
        sum += 4 * y(args[i] + len / 2);
    }
    return sum * len / 6;
}

double newton_cotes(double *args, function<double(double)> y, int n, double len) {
    double h = len / 4;
    double sum = 0;
    for (int i = 0; i < n - 1; i++) {
        sum += (7 * y(args[i]) + 32 * y(args[i] + h) + 12 * y(args[i] + 2 * h) + 32 * y(args[i] + 3 * h) +
                7 * y(args[i + 1]));
    }
    return sum * 2 * h / 45;
}

double gauss(double *args, function<double(double)> y, int n, double len) {
    double sum = 0;
    double c1 = 5.0 / 9.0;
    double c2 = 8.0 / 9.0;
    double c3 = 5.0 / 9.0;
    double ksi1 = -sqrt(3 / 5);
    double ksi3 = sqrt(3 / 5);
    for (int i = 0; i < n - 1; i++) {
//        cout << c1*y(len/2*ksi1+(args[i]+args[i+1])/2) << endl;
//        cout << c2*y((args[i]+args[i+1])/2) << endl;
//        cout << c3*y(len/2*ksi3+(args[i]+args[i+1])/2) << endl;
        sum += (c1 * y(len / 2 * ksi1 + (args[i] + args[i + 1]) / 2) + c2 * y((args[i] + args[i + 1]) / 2) +
                c3 * y(len / 2 * ksi3 + (args[i] + args[i + 1]) / 2));
    }
    return sum * len / 2;
}

void write_errs(double rectangle_h, double trapeze_h, double simps_h, double newton_cotes_h, double gauss_h,
                double rectangle_h2, double trapeze_h2, double simps_h2, double newton_cotes_h2, double gauss_h2,
                double integral) {
    ofstream fout("errs");
    fout << fixed;
    fout.precision(6);
    fout << setw(15) << left << scientific << "scheme"
         << setw(15) << left << "error h"
         << setw(15) << left << "error h/2"
         << endl;
    fout << setw(15) << left << "rectangle"
         << setw(15) << left << (rectangle_h - integral) / integral
         << setw(15) << left << (rectangle_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "trapeze"
         << setw(15) << left << (trapeze_h - integral) / integral
         << setw(15) << left << (trapeze_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "simpson"
         << setw(15) << left << (simps_h - integral) / integral
         << setw(15) << left << (simps_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "newton_cotes"
         << setw(15) << left << (newton_cotes_h - integral) / integral
         << setw(15) << left << (newton_cotes_h2 - integral) / integral
         << endl;
    fout << setw(15) << left << "gauss"
         << setw(15) << left << (gauss_h - integral) / integral
         << setw(15) << left << (gauss_h2 - integral) / integral
         << endl;
}

int main() {
    function<double(double)> y = func;

    double a = -3.14;
    double b = 0;

    int k = 3;
    int m1 = k + 1;
    int m2 = 2 * k + 1;

    double h1 = (b - a) / k;
    double h2 = (b - a) / (2 * k);


    double args1[m1];
    double args2[m2];

    f_even_args(args1, m1, a, b);
    f_even_args(args2, m2, a, b);

    double rectangle_h = rectangle(args1, y, h1, m1);
    double trapeze_h = trapeze(args1, y, m1, h1);
    double simps_h = simps(args1, y, m1, h1);
    double newton_cotes_h = newton_cotes(args1, y, m1, h1);
    double gauss_h = gauss(args1, y, m1, h1);

    double rectangle_h2 = rectangle(args2, y, h2, m2);
    double trapeze_h2 = trapeze(args2, y, m2, h2);
    double simps_h2 = simps(args2, y, m2, h2);
    double newton_cotes_h2 = newton_cotes(args2, y, m2, h2);
    double gauss_h2 = gauss(args2, y, m2, h2);

    double internal = i_func(a, b);

    write_errs(rectangle_h, trapeze_h, simps_h, newton_cotes_h, gauss_h,
               rectangle_h2, trapeze_h2, simps_h2, newton_cotes_h2, gauss_h2, internal);


    return 0;
}
