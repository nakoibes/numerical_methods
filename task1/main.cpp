#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>

using namespace std;

double func(double x) {
    return sin(x);
}

double random(double min, double max) {
    return (double) (rand()) / RAND_MAX * (max - min) + min;
}

void f_inter_args(double *args, int length, double a, double b, double delta) {
    args[0] = a;
    args[length - 1] = b;
    int i = 1;
    while (i < length - 1) {
        double value = random(a, b);
        bool good_value = true;
        for (int j = 0; j < i; j++) {
            if (abs(args[j] - value) < delta) {
                good_value = false;
                break;
            }
        }
        if (good_value) {
            args[i] = value;
            i++;
        }
    }
}

void f_inter_vals(double *vals, double *args, function<double(double)> y, int length) {
    for (int i = 0; i < length; i++) {
        vals[i] = y(args[i]);
    }

}

void f_repr_args(double *repr_args, int n_repr, double a, double b) {
    repr_args[0] = a;
    repr_args[n_repr - 1] = b;
    double step = (b - a) / (n_repr - 1);
    for (int i = 1; i < n_repr - 1; i++) {
        repr_args[i] = repr_args[i - 1] + step;
    }
}

void f_repr_vals(double *repr_vals, double *repr_args, double *inter_vals, double *inter_args, int n_repr, int n_inter) {
    for(int i =0;i<n_repr;i++){
        repr_vals[i] = 0;
    }
    double multipliers[n_inter];
    for (int i = 0; i < n_inter; i++) {
        double product = 1;
        for (int j = 0; j < n_inter; j++) {
            if (j != i) {
                product = product * (inter_args[i] - inter_args[j]);
            }
        }
        multipliers[i] = inter_vals[i] / product;
    }
    for (int i = 0; i < n_repr; i++) {
        for (int k = 0; k < n_inter; k++) {
            double product = 1;
            for (int j = 0; j < n_inter; j++) {
                if (k != j) {
                    product = product * (repr_args[i] - inter_args[j]);
                }
            }
            repr_vals[i] += multipliers[k] * product;
        }
    }
}

void write_f(double *res_args, double *res_vals, double *inter_args, double *inter_vals, double *func, int n_inter,
             int n_res) {
    ofstream fout("dots.txt");
    for (int i = 0; i < n_res; i++) {
        fout << res_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_res; i++) {
        fout << res_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_inter; i++) {
        fout << inter_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_inter; i++) {
        fout << inter_vals[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < n_res; i++) {
        fout << func[i] << " ";
    }
    fout.close();
}

void fill_func(double *func, function<double(double)> y, double *res_args, int n_repr) {
    for (int i = 0; i < n_repr; i++) {
        func[i] = y(res_args[i]);
    }
}

int main() {
    srand((unsigned int) time(nullptr));

    function<double(double)> y = func;

    double a = -5;
    double b = 5;
    double delta = 0.0001;
    int n_inter = 15;
    int n_repr = 2000;

    double inter_args[n_inter];
    double inter_vals[n_inter];
    double repr_args[n_repr];
    double repr_vals[n_repr];

    double func[n_repr];

    f_inter_args(inter_args, n_inter, a, b, delta);
    f_inter_vals(inter_vals, inter_args, y, n_inter);
    f_repr_args(repr_args, n_repr, a, b);
    f_repr_vals(repr_vals, repr_args, inter_vals, inter_args, n_repr, n_inter);

    fill_func(func, y, repr_args, n_repr);

    write_f(repr_args, repr_vals, inter_args, inter_vals, func, n_inter, n_repr);
    system("python3 repr.py");
    return 0;
}
