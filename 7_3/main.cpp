#include <iostream>
#include <cmath>
using namespace std;


int main() {
    int Nexp = 10;
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};
    double tau0 = 1.0, eps = 1e-9, K = 10, n = 1;
    double F, F1, F2, F3;
    double F11, F12, F13, F21, F22, F23, F31, F32, F33;
    double delt, delta, delta1, delta2, delta3, delt1, delt2, delt3, dtau0, dK, dn, dev;
    double alphak = 0.1;
    int kiter = 0;
    do {
        F = 0;
        for (int i = 0; i < Nexp; ++i) {
            F = F + pow(tau0 + K * pow(x[i], n) - y[i], 2);
        }
        F1 = 0;
        F2 = 0;
        F3 = 0;
        for (int i = 0; i < Nexp; ++i) {
            F1 = F1 + 2 * (tau0 + K * pow(x[i], n) - y[i]);
            F2 = F2 + 2 * (tau0 + K * pow(x[i], n) - y[i]) * pow(x[i], n);
            F3 = F3 + 2 * (tau0 + K * pow(x[i], n) - y[i]) * K * pow(x[i], n) * log(x[i]);
        }
        F12 = F13 = F22 = F23 = F33 = 0;
        for (int i = 0; i < Nexp; ++i) {
            F12 = F12 + 2 * pow(x[i], n);
            F13 = F13 + 2 * K * pow(x[i], n) * log(x[i]);
            F22 = F22 + 2 * pow(x[i], 2 * n);
            F23 = F23 + 2 * (tau0 + 2 * K * pow(x[i], n) - y[i]) * pow(x[i], n) * log(x[i]);
            F33 = F33 + 2 * (tau0 + 2 * K * pow(x[i], n) - y[i]) * K * pow(x[i], n) * log(x[i]) * log(x[i]);
        }
        F11 = 2 * Nexp;
        F21 = F12;
        F32 = F23;
        F31 = F13;

        delt = F11 * F22 * F33 + F12 * F23 * F31 + F13 * F21 * F32;
        delta = delt - F13 * F22 * F31 - F12 * F21 * F33 - F11 * F32 * F23;
        delt1 = F1 * F22 * F33 + F12 * F23 * F3 + F13 * F2 * F32;
        delta1 = delt1 - F13 * F22 * F3 - F23 * F32 * F1 - F2 * F12 * F33;
        delt2 = F11 * F2 * F33 + F1 * F23 * F31 + F13 * F3 * F21;
        delta2 = delt2 - F13 * F2 * F31 - F23 * F3 * F11 - F33 * F21 * F1;
        delt3 = F11 * F22 * F3 + F12 * F2 * F31 + F1 * F21 * F32;
        delta3 = delt3 - F1 * F22 * F31 - F2 * F32 * F11 - F3 * F21 * F12;
        dtau0 = -delta1 / delta;
        dK = -delta2 / delta;
        dn = -delta3 / delta;


        tau0 = tau0 + alphak * dtau0;
        K = K + alphak * dK;
        n = n + alphak * dn;
        dev = sqrt(pow(dtau0, 2) + pow(dK, 2) + pow(dn, 2));
        kiter += 1;

    } while (dev < eps);

    cout << "Newton method. NAME LASTNAME group BPOi-16" << endl;
    cout << "tau0: " << tau0 << endl;
    cout << "K: " << K << endl;
    cout << "n: " << n << endl;
    cout << "Function: " << F << endl;
    cout << "Derivatives: " << abs(F1) << "  " << abs(F2) << "  " << abs(F3) << endl;
    return 0;
}