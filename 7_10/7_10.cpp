#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double s1, s2, s3,tau0, KK, tau0b, T1, T2;
double f, f1, f2;
int Nexp = 10;
double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};
double xb[10], yb[10];

double a = 0, b = 2, eps = 1e-6, delta = eps / 3;
double alpha, beta, n, xx;

void tau() {
    s1 = 0, s2 = 0, s3 = 0;
    for (int i = 0; i < Nexp; ++i) {
        s1 += (yb[i]-pow(xb[i], n))*(1-pow(xb[i], n));
        s2 += pow(1-pow(xb[i], n), 2);
    }
    tau0b = s1/s2;
    for (int i = 0; i < Nexp; ++i) {
        s3 += 2*(yb[i]-tau0b-(1-tau0b)*pow(xb[i], n))*(tau0b-1)*pow(xb[i], n)*log(xb[i]);
    }
}

void derivative() {
    f = 0, f1 = 0, f2 = 0;
    for (int i = 0; i < Nexp; ++i) {
        f += pow(yb[i]-tau0b-(1-tau0b)*pow(xb[i], n), 2);
        f1 += 2*(yb[i]-tau0b-(1-tau0b)*pow(xb[i], n))*(pow(xb[i], n-1));
        f2 += 2*(yb[i]-tau0b-(1-tau0b)*pow(xb[i], n))*(tau0b-1)*pow(xb[i], n)*log(xb[i]);
    }
}

int main() {
    cout << "Ref\ttau0\tK\tn\tL\t\tabs(L1)\t\tabs(L2)" << endl;

    for (int r = 0; r < Nexp; ++r) {
        for (int i = 0; i < Nexp; ++i) {
            xb[i] = x[i] / x[r];
            yb[i] = y[i] / y[r];
        }

        a = 0, b = 2, eps = 1e-6, delta = eps/3;

        for (int j = 0; j < 1000; ++j) {
            alpha = (a + b) / 2 - delta;
            n = alpha;
            tau();
            T1 = abs(s3);
            beta = (a + b) / 2 + delta;
            n = beta;
            tau();
            T2 = abs(s3);
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
                derivative();
                break;
            }
        }
        cout.precision(4);
        cout << r+1 << "\t" << fixed << tau0 << "\t" << KK << "\t" << n << "\t" << scientific << setprecision(2) << f << "\t" << abs(f1) << "\t" << abs(f2) << endl;
    }

    return 0;
}
