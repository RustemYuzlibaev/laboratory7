#include <iostream>
#include <cmath>

using namespace std;

double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};

double K, tau0;
double F, F1, F2, F3, Fa, Fb;
int Nexp = 10, Nit = 1000;
double eps = 1e-12, alpha, beta, n, a = 0, b = 2;
double a11, a12, a21, a22, b1, b2;
double d, d1, d2;

void Ktau0(int Nexp, double *x, double *y, double n) {
    a11 = 0, a12 = 0, a21 = 0, a22 = 0, b1 = 0, b2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        a11 = a11 + pow(x[i], 2 * n);
        a12 = a12 + pow(x[i], n);
        b1 = b1 + pow(x[i], n) * y[i];
        b2 = b2 + y[i];
    }
    a21 = a12;
    a22 = Nexp;
    d = a11*a22-a12*a21;
    d1 = b1*a22-b2*a12;
    d2 = b2*a11-b1*a21;

    K = d1/d;
    tau0 = d2/d;
}

void FF(int Nexp, double *x, double *y, double n) {
    F = 0, F1 = 0, F2 = 0, F3 = 0;
    for (int i = 0; i < Nexp; ++i) {
        F = F + pow(tau0 + K * pow(x[i], n) - y[i], 2) / Nexp;
        F1 = F1 + 2 * (tau0 + K * pow(x[i], n) - y[i]) / Nexp;
        F2 = F2 + 2 * (tau0 + K * pow(x[i], n) - y[i]) * pow(x[i], n) / Nexp;
        F3 = F3 + 2 * (tau0 + K * pow(x[i], n) - y[i]) * K * pow(x[i], n) * log(x[i]) / Nexp;
    }
}

int main() {


    for (int i = 0; i < Nit; ++i) {
        alpha = a + 2 / (3 + pow(5, 0.5))*(b-a);
        n = alpha;

        Ktau0(Nexp, x, y, n);
        FF(Nexp, x, y, n);
        Fa = F;

        beta = a + 2 / (1 + pow(5, 0.5))*(b-a);
        n = beta;

        Ktau0(Nexp, x, y, n);
        FF(Nexp, x, y, n);
        Fb = F;

        if (Fa <= Fb) {
            b = beta;
            n = alpha;
        } else {
            a = alpha;
            n = beta;
        }

        if ((b - a) < eps) {
            Ktau0(Nexp, x, y, n);
            FF(Nexp, x, y, n);
            break;
        }
    }

    cout << "tau0: " << tau0 << endl;
    cout << "K: " << K << endl;
    cout << "n: " << n << endl;
    cout << "F: " << F << endl;
    cout << "F1: " << F1 << endl;
    cout << "F2: " << F2 << endl;
    cout << "F3: " << F3 << endl;

    return 0;
}
