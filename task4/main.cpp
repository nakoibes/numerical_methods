#include <iostream>
#include <cmath>
#include <fstream>
#include <functional>
using namespace std;

double func(double x) {
    return sin(x);
}

double d_func(double x) {
    return cos(x);
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
void differentiate(double* args,double* vals,double* der ,int n,double h){
    der[0] = (-3*vals[0]+4*vals[1]-vals[2])/(2*h);
    der[n-1] = (vals[n-3]-4*vals[n-2]+3*vals[n-1])/(2*h);
    for(int i=1;i<n-2;i++){
        der[i] = (vals[i+1]-vals[i-1])/(2*h);
    }
}

int main() {
    function<double(double)> y = func;
    function<double(double)> d_y = d_func;

    double a = -5.0;
    double b = 5.0;

    int m1 = 30;
    int m2 = 2*m1-1;
    double h1 = (b - a) / (m1-1);

    double args1[m1];
    double vals1[m1];
    double d_vals1[m1];
    double der1[m1];
    double args2[m2];
    double vals2[m2];
    double der2[m1];

    f_even_args(args1,m1,a,b);
    f_even_args(args2,m2,a,b);
    fill_func(vals1,y,args1,m1);
    fill_func(vals2,y,args2,m2);
    fill_func(d_vals1,y,args1,m1);

    fill_func(d_vals1,d_y,args1,m1);

    differentiate(args1,vals1,der1,m1,h1);


    return 0;
}
