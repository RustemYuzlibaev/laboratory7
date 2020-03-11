#include <iostream>
#include <cmath>
using namespace std;


int main()
{
    int Nexp = 10;
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203,  26.270, 27.233, 28.109};

    double F1, F2, F3, F, F0;
    double tau0 = 1, n = 1, n0, tau00;
    double K = 10, K0;
    double eps1 = 1e-9, eps2 = 1e-13;
    double H = 0.001;
    int nit = 50000;
    bool wasReached = false;

    for (int i = 0; i < nit; i++) {
        tau00 = tau0;
        K0 = K;
        n0 = n;

        if (wasReached) {
            break;
        }

        F1 = 0, F2 = 0, F3 = 0;
        for (int i = 0; i < Nexp; i++) {
            F1 = F1 + 2 * (tau0 + K * pow(x[i], n) - y[i]);
            F2 = F2 + 2 * (tau0 + K * pow(x[i], n) - y[i]) * pow(x[i], n);
            F3 = F3 + 2 * (tau0 + K * pow(x[i], n) - y[i]) * K * pow(x[i], n) * log(x[i]);
        }

        do {
            F = 0;
            for (int i = 0; i < Nexp; ++i) {
                F+=pow(tau0+K*pow(x[i], n)-y[i], 2);
            }

            F0 = F;
            tau0 = tau0-H*F1;
            K=K-H*F2;
            n=n-H*F3;

            F = 0;
            for (int i = 1; i < Nexp; ++i) {
                F+=pow(tau0+K*pow(x[i], n)-y[i], 2);
            }
        } while(F > F0);

        if ( pow(pow(tau0-tau00, 2)+pow(K-K0, 2)+pow(n-n0, 2), 0.5) < eps1) {
            cout << "tauo: " << tau0 << endl;
            cout << "K: " << K << endl;
            cout << "n: " << n << endl;
            cout << "F: " << F << endl;
            cout << "F1: " << F1 << endl;
            cout << "F2: " << F2 << endl;
            cout << "F3: " << F3 << endl;
            wasReached = true;
            break;
        }

        if (abs(F-F0)<eps2) {
            cout << "tauo: " << tau0 << endl;
            cout << "K: " << K << endl;
            cout << "n: " << n << endl;
            cout << "F: " << F << endl;
            cout << "F1: " << F1 << endl;
            cout << "F2: " << F2 << endl;
            cout << "F3: " << F3 << endl;
            wasReached = true;
            break;
        }
    }
    return 0;
}