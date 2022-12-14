#include <iostream>
#include <functional>
#include <cmath>
#include <fstream>

using namespace std;

double func(double x) {
    return sin(x);
}

void f_inter_args(double **args, int len_1, int len_2, double a, double step) {
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

void f_inter_vals(double **args, double **vals, function<double(double)> y, int len_1, int len_2) {
    for (int i = 0; i < len_1; i++) {
        for (int j = 0; j < len_2; j++) {
            vals[i][j] = y(args[i][j]);
        }
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

void interpolate(double *repr_vals, double *repr_args, double *inter_vals, double *inter_args, int n_repr, int n_inter) {
    for (int i = 0; i < n_repr; i++) {
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

void f_repr_vals(double *repr_vals, double *repr_args, double **inter_vals, double **inter_args, int n_repr, int len_1, int len_2) {
    double eps = 0.000001;
    int start = 0;
    int end = 0;
    int len;
    for (int i = 0; i < len_1; i++) {
        double current = inter_args[i][len_2 - 1];
        while (current-repr_args[end] > eps) {
            end++;
        }
        len = end - start + 1;
        interpolate(&repr_vals[start], &repr_args[start], inter_vals[i], inter_args[i], len, len_2);
        start = end;
    }
}

void f_res_int(double *res_int_args, double *res_int_vals, double **int_args, double **int_vals, int len_1, int len_2) {
    int k = 0;
    for (int i = 0; i < len_1; i++) {
        for (int j = 0; j < len_2; j++) {
            if (i != len_1 - 1) {
                if (j != len_2 - 1) {
                    res_int_args[k] = int_args[i][j];
                    res_int_vals[k] = int_vals[i][j];
                    k++;
                }
            }
            else {
                res_int_args[k] = int_args[i][j];
                res_int_vals[k] = int_vals[i][j];
                k++;
            }
        }
    }
}

void write_f(double *res_args, double *res_vals, double *inter_args, double *inter_vals, double *func,double* div_knots_args,double* div_knots_vals,
             int n_inter, int n_res, int k) {
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
    for (int i = 0; i < k+1; i++) {
        fout << div_knots_args[i] << " ";
    }
    fout << endl;
    for (int i = 0; i < k+1; i++) {
        fout << div_knots_vals[i] << " ";
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

double calc_abs_err_1(double *func, double *inter, int n_err) {
    double result = 0;
    for (int i = 0; i < n_err; i++) {
        result += abs(func[i] - inter[i]);
    }
    return result;
}

double calc_abs_err_2(double *func, double *inter, int n_err) {
    double result = 0;
    for (int i = 0; i < n_err; i++) {
        result += pow(abs(func[i] - inter[i]), 2);
    }
    return sqrt(result);
}

double calc_abs_err_cheb(double *func, double *inter, int n_err) {
    double result = 0;
    for (int i = 0; i < n_err; i++) {
        if (abs(func[i] - inter[i]) > result) {
            result = abs(func[i] - inter[i]);
        }
    }
    return result;
}

double calc_rel_err_1(double *func, double *inter, int n_err) {
    double num = 0;
    double den = 0;
    for (int i = 0; i < n_err; i++) {
        num += (abs(func[i] - inter[i]));
    }
    for (int i = 0; i < n_err; i++) {
        den += abs(func[i]);
    }
    if(den == 0.0){
        den = 1;
    }
    return num / den;
}

double calc_rel_err_2(double *func, double *inter, int n_err) {
    double num = 0;
    double den = 0;
    for (int i = 0; i < n_err; i++) {
        num += (pow(abs(func[i] - inter[i]), 2));
    }
    for (int i = 0; i < n_err; i++) {
        den += pow(func[i], 2);
    }
    if(den == 0.0){
        den = 1;
    }
    return sqrt(num)/sqrt(den);
}

double calc_rel_err_cheb(double *func, double *inter, int n_err) {
    double num = 0;
    double den = 0;
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
    if(den == 0.0){
        den = 1;
    }
    return num / den;
}

void write_errs(double* err_func,double* err_vals,int n_err) {
    ofstream fout("errs.txt");
    fout << fixed;
    fout.precision(6);
    fout << "abs_err_1 " << scientific << calc_abs_err_1(err_func, err_vals, n_err) << endl;
    fout << "abs_err_2 " << calc_abs_err_2(err_func, err_vals, n_err) << endl;
    fout << "abs_err_c " << calc_abs_err_cheb(err_func, err_vals, n_err) << endl;
    fout << "rel_err_1 " << calc_rel_err_1(err_func, err_vals, n_err) << endl;
    fout << "rel_err_2 " << calc_rel_err_2(err_func, err_vals, n_err) << endl;
    fout << "rel_err_c " << calc_rel_err_cheb(err_func, err_vals, n_err) << endl;
    fout.close();
}

void f_div_knots_args(double* div_knots,double** args,int k,int n){
    div_knots[0] = args[0][0];
    for(int i = 0;i<k;i++){
        div_knots[i+1] = args[i][n];
    }
}
void f_div_knots_vals(double* div_knots,double** vals,int k,int n){
    div_knots[0] = vals[0][0];
    for(int i = 0;i<k;i++){
        div_knots[i+1] = vals[i][n];
    }
}

int main() {
    function<double(double)> y = func;

    double a = 0.0;
    double b = 8.0;

    int k = 2;
    int n = 5;//??????????????
    int n_repr = 2000;

    int m = n * k + 1;

    double int_len = (b - a) / k;
    double h = int_len / n;

    int n_err = 100*m-99;//right

    double **int_args = new double *[k];
    for (int i = 0; i < k; i++) {
        int_args[i] = new double[n + 1];
    }

    double **int_vals = new double *[k];
    for (int i = 0; i < k; i++) {
        int_vals[i] = new double[n + 1];
    }

    double res_int_args[m];
    double res_int_vals[m];
    double repr_args[n_repr];
    double repr_vals[n_repr];
    double func[n_repr];
    double err_func[n_err];
    double err_args[n_err];
    double err_vals[n_err];

    double div_knots_args[k+1];
    double div_knots_vals[k+1];

    f_inter_args(int_args, k, n + 1, a, h);
    f_div_knots_args(div_knots_args, int_args,k,n);
    f_inter_vals(int_args, int_vals, y, k, n + 1);
    f_div_knots_vals(div_knots_vals, int_vals,k,n);
    f_repr_args(repr_args, n_repr, a, b);
    f_repr_vals(repr_vals, repr_args, int_vals, int_args, n_repr, k, n + 1);
    f_res_int(res_int_args, res_int_vals, int_args, int_vals, k, n + 1);

    fill_func(func, y, repr_args, n_repr);

    f_repr_args(err_args, n_err, a, b);
    fill_func(err_func, y, err_args, n_err);
    f_repr_vals(err_vals, err_args, int_vals, int_args, n_err, k, n + 1);

    write_errs(err_func,err_vals,n_err);
    write_f(repr_args, repr_vals, res_int_args, res_int_vals, func,div_knots_args,div_knots_vals, m, n_repr,k);

    system("python3 repr.py");

    for (int i = 0; i < k; i++) {
        delete[] int_args[i];
    }
    delete[] int_args;

    for (int i = 0; i < k; i++) {
        delete[] int_vals[i];
    }
    delete[] int_vals;

    return 0;
}
