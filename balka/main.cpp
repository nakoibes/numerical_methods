#include <iostream>
#include <fstream>
#include <string>
//#include <Eigen/Dense>
#include "eigen/Eigen/Dense"

using namespace std;


double *loc_force(double q1, double q2,double l) {
    double *f = new double[4];
    f[0] = (q2-q1)*3*l/20.0 + q1*l/2.0;
    f[1] = (q2-q1)*l*l/30.0 + q1*l*l/12.0;
    f[2] = (q2-q1)*7*l/20.0 + q1*l/2.0;
    f[3] = -(q2-q1)*l*l/20.0 - q1*l*l/12.0;
    return f;
}

double **calc_A_loc(double E, double J,double l) {
    double **a = new double *[4];
    for (int i = 0; i < 4; i++) {
        a[i] = new double[4];
    }
    a[0][0] = a[2][2] = 12 * E * J/(l*l*l);
    a[0][1] = a[1][0] = 6 * E * J/(l*l);
    a[0][2] = a[2][0] = -12 * E * J/(l*l*l);
    a[0][3] = a[3][0] = 6 * E * J/(l*l);
    a[1][1] = a[3][3] = 4 * E * J/l;
    a[1][2] = a[2][1] = -6 * E * J/(l*l);
    a[1][3] = a[3][1] = 2 * E * J/l;
    a[2][3] = a[3][2] = -6 * E * J/(l*l);
    return a;
}


void add_hinge(int el, double **A, double *F) {
    int p = el * 2;
    for (int i = 0; i < 52; i++) {
        if (p == i) {
            A[p][i] = 1;
            A[i][p] = 1;
        } else {
            A[p][i] = 0;
            A[i][p] = 0;
        }
    }
    F[p] = 0;
}


void add_seal(int el, double **A, double *F) {
    int p = el * 2;
    for (int i = 0; i < 52; i++) {
        A[p][i] = 0;
        A[i][p] = 0;
        A[p + 1][i] = 0;
        A[i][p + 1] = 0;
    }
    A[p][p] = 1;
    A[p + 1][p + 1] = 1;
    F[p] = 0;
    F[p + 1] = 0;
}


void add_force(int el, double *F, double P) {
    int p = el * 2;
    F[p] += P;
}

void add_elst_sup(int el, double **K, double c) {
    int p = el * 2;
    K[p][p] += c;
}

void add_moment(int el, double *F, double M) {
    int p = el * 2 + 1;
    F[p] -= M;
}

void add_distr_load(int first, int last, double q_1, double q_2, double *F,double l) {
    double k = (q_2 - q_1) / (last - first);
    for (int el = first; el < last; el++) {
        double q_1_loc = q_1 + k * (el - first);
        double q_2_loc = q_1 + k * (el + 1 - first);
        double *f = loc_force(q_1_loc, q_2_loc,l);
        int p = el * 2;
        for (int i = 0; i < 4; i++) {
            F[p + i] += f[i];
        }
        delete[] f;
    }

}


double calc_w(int i, int N, int el, double *koefs,double l) {
    double w = 0.0;
    double z = (static_cast<double>(i)*l) / N;
    int p = el * 2;
    w += koefs[p] * (1.0 - 3.0 * z * z/(l*l) + 2.0 * z * z * z/(l*l*l));
    w += koefs[p + 1] * (z - 2.0 * z * z/l + z * z * z/(l*l));
    w += koefs[p + 2] * (3.0 * z * z/(l*l) - 2.0 * z * z * z/(l*l*l));
    w += koefs[p + 3] * (-z * z/l + z * z * z/(l*l));
    return w;
}


double calc_th(int i, int N, int el, double *koefs,double l) {
    double th = 0.0;
    double z = (static_cast<double>(i)*l) / N;
    int p = el * 2;
    th += koefs[p] * (-6.0 * z/(l*l) + 6.0 * z * z/(l*l*l));
    th += koefs[p + 1] * (1 - 4.0 * z/l + 3.0 * z * z/(l*l));
    th += koefs[p + 2] * (6.0 * z/(l*l) - 6.0 * z * z/(l*l*l));
    th += koefs[p + 3] * (-2.0 * z/l + 3.0 * z * z/(l*l));
    return th;
}


double calc_M(int i, int N, int el, double *koefs, double E, double J,double l) {
    double M = 0.0;
    double z = (static_cast<double>(i)*l) / N;
    int p = el * 2;
    M += koefs[p] * (-6.0/(l*l) + 12.0 * z/(l*l*l));
    M += koefs[p + 1] * (-4.0/l + 6.0 * z/(l*l));
    M += koefs[p + 2] * (6.0/(l*l) - 12.0 * z/(l*l*l));
    M += koefs[p + 3] * (-2.0/l + 6.0 * z/(l*l));
    return M * E * J;
}


double calc_Q(int i, int N, int el, double *koefs, double E, double J,double l) {
    double Q = 0.0;
    int p = el * 2;
    Q += koefs[p] * 12.0/(l*l*l);
    Q += koefs[p + 1] * 6.0/(l*l);
    Q += koefs[p + 2] * (-12.0)/(l*l*l);
    Q += koefs[p + 3] * 6.0/(l*l);
    return Q * E * J;
}

Eigen::MatrixXd create_A_E(double **A, int n) {
    Eigen::MatrixXd A_E(n, n);
    A_E.setZero();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_E(i, j) = A[i][j];
        }
    }
    return A_E;
}

Eigen::VectorXd create_B_E(double *B, int n) {
    Eigen::VectorXd B_E(n);
    B_E.setZero();
    for (int i = 0; i < n; i++) {
        B_E(i) = B[i];
    }
    return B_E;
}

double *create_koefs(Eigen::VectorXd X, int n) {
    double *koefs = new double[n];
    for (int i = 0; i < n; i++) {
        koefs[i] = X(i);
    }
    return koefs;
}


void write(double *koefs, int N, double E, double J,double l) {
    ofstream stream("/home/nakoibes/prak_balka/cpp_output.txt");
    for (int el = 0; el < 25; el++) {
        for (int i = 0; i < N; i++) {
            stream << el*l + (static_cast<double>(i)*l) / N << ','
                   << calc_w(i, N, el, koefs,l) << ','
                   << calc_th(i, N, el, koefs,l) << ','
                   << calc_M(i, N, el, koefs, E, J,l) << ','
                   << calc_Q(i, N, el, koefs, E, J,l) << endl;
        }
    }
    stream.close();
}


int main() {
    int N = 50;

    string s2;
    ifstream file2("/home/nakoibes/prak_balka/params.txt");
    getline(file2, s2);
    double J = stod(s2);
    getline(file2, s2);
    double E = stod(s2);
    getline(file2, s2);
    double q_1 = stod(s2);
    getline(file2, s2);
    double q_2 = stod(s2);
    getline(file2, s2);
    double c = stod(s2);
    getline(file2, s2);
    double M = stod(s2);
    getline(file2, s2);
    double P = stod(s2);


    string s;
    ifstream file("/home/nakoibes/prak_balka/dots.txt");
    getline(file, s);
    int Pk = stoi(s);
    getline(file, s);
    int R1k = stoi(s);
    getline(file, s);
    int R2k = stoi(s);
    getline(file, s);
    int R3k = stoi(s);
    getline(file, s);
    int q1k = stoi(s);
    getline(file, s);
    int q2k = stoi(s);
    getline(file, s);
    int qlenk = stoi(s);
    getline(file, s);
    int Mk = stoi(s);
    getline(file, s);
    int kk = stoi(s);
    getline(file, s);
    int sealk = stoi(s);
    getline(file, s);
    double l = stod(s);


    double **A = new double *[52];
    for (int i = 0; i < 52; i++) {
        A[i] = new double[52];
        for (int j = 0; j < 52; j++) {
            A[i][j] = 0;
        }
    }
    double **a = calc_A_loc(E, J,l);

    for (int el = 0; el < 25; el++) {
        int p = el * 2;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                A[p + i][p + j] += a[i][j];
            }
        }
    }
//    for (int i = 0; i < 52; i++) {
//        for (int j = 0; j < 52; j++) {
//            cout <<  A[i][j] << " ";
//        }
//        cout << endl;
//    }


    double *F = new double[52];
    for (int i = 0; i < 52; i++) {
        F[i] = 0;
    }


    add_seal(sealk, A, F);
    add_hinge(R1k, A, F);
    add_distr_load(q1k, q2k, q_1, q_2, F,l);
    add_moment(Mk, F, M);
    add_hinge(R2k, A, F);
    add_hinge(R3k, A, F);
    add_force(Pk, F, P);
    add_elst_sup(kk, A, c);

    Eigen::MatrixXd A_E = create_A_E(A, 52);

    Eigen::VectorXd koefs_E(52);
    koefs_E.setZero();

    Eigen::MatrixXd B_E = create_B_E(F, 52);

    koefs_E = A_E.colPivHouseholderQr().solve(B_E);

    double *koefs = create_koefs(koefs_E, 52);

    write(koefs, N, E, J,l);


    for (int i = 0; i < 52; i++) { delete[] A[i]; }
    delete[] A;
    for (int i = 0; i < 4; i++) { delete[] a[i]; }
    delete[] a;
    delete[] F;
    delete[] koefs;
    return 0;
}