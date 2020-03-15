#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

double xb[10], yb[10];
double s1, s2, tau0b, L, L1, L2;
int Nexp = 10, M = 1000;
double a, b, eps, delta, alpha, beta, n, F1, F2, tau0, K;

void Func(double n) {
    s1 = s2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        s1 = s1 + (yb[i] - pow(xb[i], n)) * (1 - pow(xb[i], n));
        s2 = s2 + pow(1 - pow(xb[i], n), 2);
    }
    tau0b = s1 / s2;
}

void LFunc(double n) {
    L = 0;
    for (int i = 0; i < Nexp; ++i) {
        L = L + pow(yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n), 2);
    }
}

void L1L2Func(double n) {
    L1 = L2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        L1 = L1 + 2 * yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n) * (-1 * pow(xb[i], n));
        L2 = L2 + 2 * yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n) * (-1 + tau0b) * pow(xb[i], n) * log(xb[i]);
    }
}

int main() {
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};


    cout << "Outer cycle method. NAME LASTNAME group BPOi-16" << endl;
    cout << "r" << setw(8) << "tau0" << setw(8) << "K" << setw(8) << "n" << setw(8) << "L" << setw(16) << "L1" << setw(8) << "L2" << endl;

    for (int r = 0; r < Nexp; ++r) {
        for (int i = 0; i < Nexp; ++i) {
            xb[i] = x[i] / x[r];
            yb[i] = y[i] / y[r];
        }
        a = 0, b = 3;
        eps = 1e-6;
        delta = eps / 3;
        for (int j = 0; j < M; ++j) {
            alpha = (a + b) / 2 - delta;
            n = alpha;
            Func(n);
            LFunc(n);
            F1 = L;

            beta = (a + b) / 2 + delta;
            n = beta;
            Func(n);
            LFunc(n);
            F2 = L;

            if (F1 < F2) {
                b = beta;
                n = alpha;
            } else {
                a = alpha;
                n = beta;
            }

            if (b - a < eps) {
                tau0 = tau0b * y[r];
                K = (y[r] - tau0) / pow(x[r], n);
                LFunc(n);
                L1L2Func(n);
            }

        }
        cout << r << "  " << tau0 << "  " << K << "  " << n << "  " << L << "  " << L1 << "  "  << L2 << endl;
    }
    return 0;
}
