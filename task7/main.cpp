#include <iostream>
#include <cmath>

using namespace std;

double fi1(double x, double y, double l) {
//    l = 1/(exp(y));
    return x - l * (x - exp(y) + 2);
}

double fi2(double x, double y, double l) {
//    l = 1/(2*x);
    return y - l * (y - x * x);
}

double fi1_(double y) {
    return -sqrt(y);
}

double fi2_(double x) {
    return  log(x+2);
}

void s_iter(double *X0, double *X, double eps,bool left) {
    int i = 1;
    double x_p = X0[0];
    double y_p = X0[1];
    while (true) {
        double x_c;
        if(left) {
            x_c = fi1_(y_p);
        }
        else{
            x_c = -fi1_(y_p);
        };
        double y_c = fi2_(x_p);
        if (abs(x_p - x_c) + abs(y_p - y_c) < eps) {
            X[0] = x_c;
            X[1] = y_c;
            cout << "simple converged for " << i << " steps" << endl;
            break;
        }
        if (i == 10000) {
            cout << "infinite loop" << endl;
            break;
        }
        x_p = x_c;
        y_p = y_c;
        i++;
    }
}

int main() {
    double eps = 0.0001;
    double X01[2] = {-0.5, 0.3};
    double X[2];
    s_iter(X01, X, eps,true);
    cout << X[0] << endl << X[1] << endl;
    X01[0] = 1.1;
    X01[1] = 1.1;
    s_iter(X01, X, eps,false);
    cout << X[0] << endl << X[1] << endl;
    return 0;
}
