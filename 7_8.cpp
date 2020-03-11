#include <iostream>
#include <cmath>
using namespace std;

double z1, z2, z3, Nexp = 10, za, zb, K, K1, K2, ta, tb;
double z11, z12, t1, t2;

void Func(double *x, double *y, double n) {
    z1 = z2 = z3 = z11 = z12 = t1 = t2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        z1 = z1+y[i]*pow(x[i], n)*log(x[i]);
        z2 = z2 + y[i];
        z3 = z3 + pow(x[i], n)*log(x[i]);
        z11 = z11 + pow(x[i], 2*n)*log(x[i]);
        z12 = z12 + pow(x[i], n);
        t1 = t1 + y[i]*pow(x[i], n);
        t2 = t2 + pow(x[i], 2*n);
    }

    za = z1 - z2 * z3 / Nexp;
    zb = z11 - z12 * z3 / Nexp;
    K1 = za/zb;
    ta = t1 - z2 * z12 / Nexp;
    tb = t2 - z12 * z12 / Nexp;
    K2 = ta / tb;
}

int main() {
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};

    double tau0, a = 0.2, b = 2, eps = 1e-12, alpha, beta, n;
    double z1, z2, z3, z11, z12, t1, t2;
    int M = 1000, Nexp = 10;
    double F, F1, F2, F3;
    for (int i = 0; i < M; ++i) {
        alpha = a + 2 / (3 + pow(5, 0.5)) * (b-a);
        beta = a + 2 / (1 + pow(5, 0.5)) * (b-a);
        n = alpha;
        Func(x, y, n);
        F1 = abs(K1 - K2);
        n = beta;
        Func(x, y, n);
        F2 = abs(K1 - K2);

        if (F1 <= F2) {
            b = beta;
            n = alpha;
        } else {
            a = alpha;
            n = beta;
        }

        if (b - a < eps) {
            Func(x, y, n);
            K = (K1 + K2) / 2;
            tau0 = 0;
            for (int i = 0; i < Nexp; ++i) {
                tau0 += y[i]-K*pow(x[i], n);
            }
            tau0 /= Nexp;
            F = F1 = F2 = F3 = 0;
            for (int i = 0; i < Nexp; ++i) {
                F += pow(y[i] - tau0 - K * pow(x[i], n), 2) / Nexp;
                F1 = F1 + 2 * (y[i] - tau0 - K * pow(x[i], n)) / Nexp;
                F2 = F2 + 2 * (y[i] - tau0 - K * pow(x[i], n) * pow(x[i], n)) / Nexp;
                F3 = F3 + 2 * (y[i] - tau0 - K * pow(x[i], n)* K * pow(x[i], n) * log(x[i])) / Nexp;
            }
        }
    }
    cout << "tau0: " << tau0 << endl;
    cout << "K: " << K << endl;
    cout << "n: " << n << endl;
    cout << "Function: " << F << endl;
    cout << "Derivatives: " << abs(F1) << "  " << abs(F2) << "  " << abs(F3) << endl;

    return 0;
}
