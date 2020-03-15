#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double xb[10], yb[10];
double tau0, tau0b = 0.5, tau0b0, eps, H = 0.001, n = 0.5, n0, K, a;
int Nmax = 1000, Nexp = 10;
double L, L0, L1, L2;

void F100() {
    L = 0;
    for (int i = 0; i < Nexp; ++i) {
        L += pow(yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n), 2);
    }
}

void F200() {
    L1 = 0, L2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        L1 = L1 + 2 * (yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n)) * (-1 + pow(xb[i], n));
        L2 = L2 + 2 * (yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n)) * (-1 + tau0b) * pow(xb[i], n) * log(xb[i]);
    }
}

int main() {
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};
    cout << "Reference point input method" << endl;
    cout << "NAME LASTNAME" << endl;
    cout << "r\ttau0\tK\tn\tL\tL1\tL2" << endl;

    for (int r = 0; r < Nexp; ++r) {
        for (int i = 0; i < Nexp; ++i) {
            xb[i] = x[i]/x[r];
            yb[i] = y[i]/y[r];
        }

        tau0b = 0.5, n = 0.5;
        eps = 1e-9, H = 0.001, Nmax = 1000;

        for (int i = 0; i < Nmax; ++i) {
            a = 0;
            tau0b0 = tau0b;
            n0 = n;
            F100();
            F200();

            do {
                L0 = L;
                a += H;
                tau0b = tau0b - a * L1;
                n -= a * L2;
                F100();
                F200();

            } while( L >  L0);


            if ((pow(tau0b - tau0b0, 2) + pow(n - n0, 2)) < eps) {
                tau0 = tau0b*y[r];
                K = (y[r]-tau0)/pow(x[r], n);
                break;
            }
        }

        cout.precision(4);
        cout << r+1 << "\t" << fixed << tau0 << "\t" << K << "\t" << n << "\t" << L << "\t" << abs(L1) << "\t" << abs(L2) << endl;
    }
    return 0;
}
