#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double a1, b1, eps1, delta1, alpha1, beta1;
double tau0, tau0b, F, F1, F2, F11, F21, n, K;
int Nexp = 10;
double xb[10], yb[10];
double a, b, eps, delta, alpha, beta;

void Func() {
    F = 0, F1 = 0, F2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        F = F + pow(yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n), 2);
        F1 = F1 + 2 * (yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n)) * (pow(xb[i], n) - 1);
        F2 = F2 + 2 * (yb[i] - tau0b - (1 - tau0b) * pow(xb[i], n)) * (tau0b - 1) * pow(xb[i], n) * log(xb[i]);
    }
}

void MPD() {
    a = 0, b = 10, eps = 1e-9, delta = eps / 3;
    for (int j = 0; j < 1000; ++j) {
        alpha = (a + b) / 2 - delta;
        n = alpha;
        Func();
        F1 = F;

        beta = (a + b) / 2 + delta;
        n = beta;
        Func();
        F2 = F;

        if (F1 < F2) {
            b = beta;
            n = alpha;
        } else {
            a = alpha;
            n = beta;
        }

        if ((b - a) < eps) {
            return;
        }
    }
}

int main() {
    double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
    double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};

    cout << "Nested outer cycle method. NAME LASTNAME group BPOi-16" << endl;
    cout << "Ref\ttau0\tK\tn\tF\t\tabs(F1)\t\tabs(F2)" << endl;

    for (int r = 0; r < Nexp; ++r) {
        for (int i = 0; i < Nexp; ++i) {
            xb[i] = x[i] / x[r];
            yb[i] = y[i] / y[r];
        }

        a1 = 0;
        b1 = 10;
        eps1 = 1e-6;
        delta1 = eps1 / 3;

        for (int j = 0; j < 1000; ++j) {
            alpha1 = (a1 + b1) / 2 - delta1;
            tau0b = alpha1;
            MPD();
            F11 = F;

            beta1 = (a1 + b1) / 2 + delta1;
            tau0b = beta1;
            MPD();
            F21 = F;

            if (F11 < F21) {
                b1 = beta1;
                tau0b = alpha1;
            } else {
                a1 = alpha1;
                tau0b = beta1;
            }

            if ((b1 - a1) < eps1) {
                tau0 = tau0b * y[r];
                K = (y[r] - tau0)/pow(x[r], n);
                break;
            }
        }
        cout.precision(4);
        cout << r+1 << "\t" << fixed << tau0 << "\t" << K << "\t" << n << "\t" << scientific << setprecision(2) << F << "\t" << abs(F1) << "\t" << abs(F2) << endl;
    }
    return 0;
}

