#include <iostream>
#include <cmath>
using namespace std;


int main() {
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};

    double tau0 = 1, K = 10, n = 1, eps = 1e-6, lambda = -0.0001, dev;
    int Iter = 0, Nit = 1000000, Nexp = 10;
    double F1, F2, F3;

    for (int i = 0; i < Nit; ++i) {
        Iter += 1;

        F1 = 0;
        for (int i = 0; i < Nexp; ++i) {
            F1 += F1 += (tau0 + K * pow(x[i], n) - y[i]);
        }

        tau0 += lambda * F1;

        F2 = 0;
        for (int i = 0; i < Nexp; ++i) {
            F2 += (tau0 + K * pow(x[i], n) - y[i]) * pow(x[i], n);
        }

        K += lambda * F2;

        F3 = 0;
        for (int i = 0; i < Nexp; ++i) {
            F3 += (tau0 + K * pow(x[i], n) - y[i]) * K * pow(x[i], n) * log(x[i]);
        }

        n += lambda * F3;

        dev = sqrt(pow(F1, 2) + pow(F2, 2) + pow(F3, 2));
        if (dev < eps) {
            break;
        }
    }

    cout << "Simple iteration method for solving systems of nonlinear equations" << endl;
    cout << "NAME LASTNAME group BPOi-16" << endl;
    cout << "tau0: " << tau0 << endl;
    cout << "K: " << K << endl;
    cout << "n: " << n << endl;
    cout << "Iteration: " << Iter << endl;
    cout << "Nevyzki: " << abs(F1) << "  " << abs(F2) << "  " << abs(F3);
    return 0;
}
