#include <iostream>
#include <functional>
#include <cmath>

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
        int good_value = 1;
        for (int j = 0; j < i; j++) {
            if (abs(args[j] - value) < delta) {
                good_value = 0;
                break;
            }
        }
        if (good_value) {
            args[i] = value;
            i++;
        }
    }
}

void f_args(double *args, int n, double a, double b) {
    args[0] = a;
    args[n - 1] = b;
    double step = (b - a) / (n - 1);
    for (int i = 1; i < n - 1; i++) {
        args[i] = args[i - 1] + step;
    }
}

int main() {
    srand((unsigned int) time(nullptr));

    function<double(double)> y = func;

    double a = -3;
    double b = 3;

    int k = 3;
    int n = 3;

    double h = 1;
    double int_len = (b - a) / k;

    int m_loc = (int) (int_len / h) + 1;
    int m = (int) ((b - a) / h) + 1;

    double int_args[k][m_loc];
    double int_vals[k][m_loc];

    int i = 0;
//    cout << (int)30/0.3 << endl;

//    system("python3 1.py");
    return 0;
}
