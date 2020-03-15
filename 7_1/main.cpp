#include <iostream>
#include <cmath>

using namespace std;

double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};

double tau00 = 1, K0 = 10, n0 = 1, H = 0.0001, eps = 1e-6, F, F0;
int Nit = 100, Nexp = 10;
double tau0 = tau00, K = K0, n = n0;

void Func() {
    F = 0;
    for (int i = 0; i < Nexp; ++i) {
        F = F + pow((tau0 + K * pow(x[i], n) - y[i]), 2);
    }
}

int main() {

    Func();

    for (int i = 0; i < Nit ; ++i) {
        do {
            F0 = F;
            tau0 += H;
            Func();
        } while (F - F0 > 0);

        do {
            F0 = F;
            K += H;
            Func();
        } while (F - F0 > 0);

        do {
            F0 = F;
            n += H;
            Func();
        } while (F - F0 > 0);

        if (abs(H) > eps / 2) {
            continue;
        } else {
            H = (-H) / 2;
        }
    }
    cout << "Coordinate descent method. NAME LASTNAME group BPOi-16" << endl;
    cout << "tau0:  " << tau0 << endl;
    cout << "K:  " << K << endl;
    cout << "n:  " << n << endl;
    cout << "F:  " << F << endl;

    return 0;
}
