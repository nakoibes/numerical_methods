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

void f_i_args(double *args, int n, double a, double b, double eps) {
    args[0] = a;
    args[n - 1] = b;
    int i = 1;
    while (i < n - 1) {
        double value = random(a, b);
        bool good_v = true;
        for (int j = 0; j < i; j++) {
            if (abs(args[j] - value) < eps) {
                good_v = false;
                break;
            }
        }
        if (good_v) {
            args[i] = value;
            i++;
        }
    }
}

//void f_i_vals(double *vals, double *args, function<double(double)> y, int n) {
//    for (int i = 0; i < n; i++) {
//        vals[i] = y(args[i]);
//    }
//
//}

void f_even(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

void f_vis_vals(double *vis_vals, double *vis_args, double *i_vals, double *i_args, int n_v, int n_i) {
    for(int i =0; i < n_v; i++){
        vis_vals[i] = 0;
    }
    double mult[n_i];
    for (int i = 0; i < n_i; i++) {
        double product = 1;
        for (int j = 0; j < n_i; j++) {
            if (j != i) {
                product = product * (i_args[i] - i_args[j]);
            }
        }
        mult[i] = i_vals[i] / product;
    }
    for (int i = 0; i < n_v; i++) {
        for (int k = 0; k < n_i; k++) {
            double prod = 1;
            for (int j = 0; j < n_i; j++) {
                if (k != j) {
                    prod *= (vis_args[i] - i_args[j]);
                }
            }
            vis_vals[i] += mult[k] * prod;
        }
    }
}

void write_d(double *res_args, double *res_vals, double *inter_args, double *inter_vals, double *func, int n_inter,
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

void fill_func(double *vals, function<double(double)> y, double *iter, int n) {
    for (int i = 0; i < n; i++) {
        vals[i] = y(iter[i]);
    }
}

int main() {
    srand((unsigned int) time(nullptr));

    function<double(double)> y = func;

    double a = -3.14;
    double b = 3.14;
    double eps = 0.0001;
    int n_i = 20;
    int n_v = 2000;

    double i_args[n_i];
    double i_vals[n_i];
    double v_args[n_v];
    double v_vals[n_v];

    double func[n_v];

    f_i_args(i_args, n_i, a, b, eps);
    fill_func(i_vals,y, i_args, n_i);
    f_even(v_args, n_v, a, b);
    f_vis_vals(v_vals, v_args, i_vals, i_args, n_v, n_i);

    fill_func(func, y, v_args, n_v);

    write_d(v_args, v_vals, i_args, i_vals, func, n_i, n_v);

    system("python3 vis.py");
    return 0;
}
