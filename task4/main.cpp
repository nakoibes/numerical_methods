#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>

using namespace std;

double func(double x) {
    return x*x*x;
}

double d_func(double x) {
    return 3*x*x;
}

void f_even_args(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

void fill_func(double *func, function<double(double)> y, double *args, int n) {
    for (int i = 0; i < n; i++) {
        func[i] = y(args[i]);
    }
}

void differentiate(double *vals, double *der, int n, double h) {
    der[0] = (-3 * vals[0] + 4 * vals[1] - vals[2]) / (2 * h);
    der[n - 1] = (vals[n - 3] - 4 * vals[n - 2] + 3 * vals[n - 1]) / (2 * h);
    for (int i = 1; i < n - 1; i++) {
        der[i] = (vals[i + 1] - vals[i - 1]) / (2 * h);
    }
}

void calc_main_err(double *main_err, double *der1, double *der2, int m1) {
    for (int i = 0; i < m1; i+=1) {
        main_err[i] = (der1[i] - der2[2*i]) / (0.5 * 0.5 - 1);
    }
}

void calc_runge(double *runge, double *main_err, double *der, int n) {
    for (int i = 0; i < n; i++) {
        runge[i] = der[i] + main_err[i];
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

void write_errs1(double *err_func1, double *err_vals1, int n_err1,
                 double *err_func2, double *err_vals2, int n_err2,
                 double *err_vals3,
                 double* main_err,
                 string filename) {
    double zeros[n_err1];
    for(int i=0;i<n_err1;i++){zeros[i] = 0.0;}
    ofstream fout(filename);
    fout << fixed;
    fout.precision(6);
    fout << setw(15) << left << "error"
         << setw(15) << left << "step h"
         << setw(15) << left << "step h/2"
         << setw(15) << left << "runge"
         << setw(15) << left << "main_err"
         << endl;

    fout << setw(15) << left << "abs_err_1"
         << setw(15) << left << scientific << calc_abs_err_1(err_func1, err_vals1, n_err1)
         << setw(15) << left << calc_abs_err_1(err_func2, err_vals2, n_err2)
         << setw(15) << left << calc_abs_err_1(err_func1, err_vals3, n_err1)
         << setw(15) << left << calc_abs_err_1(zeros, main_err, n_err1)
         << endl;

    fout << setw(15) << left << "abs_err_2"
         << setw(15) << left << calc_abs_err_2(err_func1, err_vals1, n_err1)
         << setw(15) << left << calc_abs_err_2(err_func2, err_vals2, n_err2)
         << setw(15) << left << calc_abs_err_2(err_func1, err_vals3, n_err1)
         << setw(15) << left << calc_abs_err_2(zeros, main_err, n_err1)
         << endl;

    fout << setw(15) << left << "abs_err_c"
         << setw(15) << left << calc_abs_err_cheb(err_func1, err_vals1, n_err1)
         << setw(15) << left << calc_abs_err_cheb(err_func2, err_vals2, n_err2)
         << setw(15) << left << calc_abs_err_cheb(err_func1, err_vals3, n_err1)
         << setw(15) << left << calc_abs_err_cheb(zeros, main_err, n_err1)
         << endl;

    fout << setw(15) << left << "rel_err_1"
         << setw(15) << left << calc_rel_err_1(err_func1, err_vals1, n_err1)
         << setw(15) << left << calc_rel_err_1(err_func2, err_vals2, n_err2)
         << setw(15) << left << calc_rel_err_1(err_func1, err_vals3, n_err1)
         << endl;

    fout << setw(15) << left << "rel_err_2"
         << setw(15) << left << calc_rel_err_2(err_func1, err_vals1, n_err1)
         << setw(15) << left << calc_rel_err_2(err_func2, err_vals2, n_err2)
         << setw(15) << left << calc_rel_err_2(err_func1, err_vals3, n_err1)
         << endl;

    fout << setw(15) << left << "rel_err_c"
         << setw(15) << left << calc_rel_err_cheb(err_func1, err_vals1, n_err1)
         << setw(15) << left << calc_rel_err_cheb(err_func2, err_vals2, n_err2)
         << setw(15) << left << calc_rel_err_cheb(err_func1, err_vals3, n_err1)
         << endl;

    fout.close();
}

void write_f(double *f_args, double *f_vals, double *a1_args, double *a1_vals, double *a2_args, double *a2_vals,
        double *ar_args, double *ar_vals, int n1, int n2, int n_viz) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n_viz; i++) {
        fout << f_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_viz; i++) {
        fout << f_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << a1_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << a1_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n2; i++) {
        fout << a2_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n2; i++) {
        fout << a2_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << ar_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n1; i++) {
        fout << ar_vals[i] << " ";
    }
    fout << endl;
    fout.close();
}

int main() {
    function<double(double)> y = func;
    function<double(double)> d_y = d_func;

    double a = 0;
    double b = 5;

    int m_viz = 2000;
    int m1 = 11;
    int m2 = 2 * m1 - 1;
    double h1 = (b - a) / (m1 - 1);
    double h2 = (b - a) / (m2 - 1);

    double args1[m1];
    double vals1[m1];
    double d_vals1[m1];
    double main_err[m1];
    double runge[m1];
    double der1[m1];
    double args2[m2];
    double vals2[m2];
    double d_vals2[m2];
    double der2[m2];
    double viz_args[m_viz];
    double viz_vals[m_viz];


    f_even_args(viz_args, m_viz, a, b);
    f_even_args(args1, m1, a, b);
    f_even_args(args2, m2, a, b);
    fill_func(vals1, y, args1, m1);
    fill_func(vals2, y, args2, m2);
    fill_func(viz_vals, d_y, viz_args, m_viz);

    fill_func(d_vals1, d_y, args1, m1);
    fill_func(d_vals2, d_y, args2, m2);

    differentiate(vals1, der1, m1, h1);
    differentiate(vals2, der2, m2, h2);
    calc_main_err(main_err, der1, der2, m1);
    calc_runge(runge, main_err, der1, m1);

    write_errs1(d_vals1, der1, m1, d_vals2, der2, m2, runge,main_err, "errs1");

    write_f(viz_args, viz_vals, args1, der1, args2, der2, args1, runge, m1, m2, m_viz);

    system("python3 repr.py");


    return 0;
}
