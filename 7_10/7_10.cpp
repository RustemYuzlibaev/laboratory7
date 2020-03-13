#include <iostream>
#include <cmath>

int main() {
    int Nexp = 10;
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};
    double xb[10], yb[10];

    double a = 0, b = 2, eps = 1e-6, delta = eps / 3;
    double alpha, beta, n, xx;

    for (int r = 0; r < Nexp; ++r) {
        for (int i = 0; i < Nexp; ++i) {
            xb[i] = x[i] / x[r];
            yb[i] = y[i] / y[r];
        }

        for (int j = 0; j < 1000; ++j) {
            alpha = (a + b) / 2 - delta;
            n = alpha;
            // ????? I don't know where the code
            beta = (a + b) / 2 + delta;
            n = alpha;
            // ????? I don't know where the code
            if (T1 < T2) {
                b = beta;
                xx = alpha;
            } else {
                a = alpha;
                xx = beta;
            }

            if (b - a < eps) {
                tau0 = tau0b * y[r];
                KK = (y[r] - tau0) / pow(x[r], n);
            }
        }

    }
    return 0;
}
