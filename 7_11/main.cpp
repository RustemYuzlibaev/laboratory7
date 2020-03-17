#include <iostream>
#include <cmath>
#include <iomanip>



using namespace std;

// Albina Yakhina
double x[10] = {2.9524, 3.5545, 3.9067, 4.1565, 4.3504, 4.5087, 4.6426, 4.7586, 4.8609, 4.9524};
double y[10] = {18.3344, 20.3390, 21.4480, 22.2105, 22.7894, 23.2546, 23.6429, 23.9755, 24.2663, 24.5243};

// Rustem Yuzlibaev
//double x[10] = {1.8808, 2.4829, 2.8351, 3.0849, 3.2788, 3.4371, 3.5710, 3.6870, 3.7893, 3.8808};
//double y[10] = {14.9390, 16.6878, 17.6048, 18.2185, 18.6761, 19.0392, 19.3392, 19.5941, 19.8155, 20.0107};

// Murat Tlyaumbetov
//double x[10] = {1.1468, 1.7488 ,2.1010, 2.3509, 2.5447, 2.7031, 2.8370, 2.9530, 3.0553, 3.1468};
//double y[10] = {9.4820, 12.6591, 14.3814, 15.5576, 16.4475, 17.1616, 17.7569, 18.2667, 18.7120, 19.1071};

// Danil Kudryashov
//double x[10] = {2.1630, 2.7651, 3.1172, 3.3671, 3.5609, 3.7193, 3.8532, 3.9692, 4.0715, 4.1630};
//double y[10] = {15.3495, 17.5595, 18.7677, 19.5940, 20.2192, 20.7205, 21.1381, 21.4955, 21.8075, 22.0840};

// Kamil Kamalov
//double x[10] = {1.0350, 1.6371, 1.9893, 2.2392, 2.4330, 2.5913, 2.7252, 2.8412, 2.9435, 3.0350};
//double y[10] = {9.0897, 11.9319, 13.4246, 14.4290, 15.1818, 15.7817, 16.2791, 16.7033, 17.0726, 17.3992};

// Azamat Ishmuratov
//double x[10] = {2.6401, 3.2422, 3.5944, 3.8443, 4.0381, 4.1964, 4.3303, 4.4463, 4.5486, 4.6401};
//double y[10] = {15.0165, 16.6413, 17.5256, 18.1283, 18.5832, 18.9472, 19.2500, 19.5086, 19.7342, 19.9339};

// Alina Zharmukhametova
//double x[10] = {2.4322, 3.0342, 3.3864, 3.6363, 3.8301, 3.9885, 4.1224, 4.2384, 4.3407, 4.4322};
//double y[10] = {14.8590, 15.8524, 16.3739, 16.7224, 16.9820, 17.1876, 17.3574, 17.5014, 17.6264, 17.7365};

// Tatyana Vorobyova
//double x[10] = {2.4000, 3.0020, 3.3542, 3.6041, 3.7979, 3.9563, 4.0902, 4.2062, 4.3085, 4.4000};
//double y[10] = {12.9691, 13.9516, 14.4672, 14.8118, 15.0684, 15.2717, 15.4395, 15.5819, 15.7054, 15.8143};

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
    cout << "Outer cycle method. Albina Yakhina group BPOi-16" << endl;
    cout << "Ref\ttau0\tK\tn\tL\t\tL1\t\tL2" << endl;
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
        cout.precision(4);
        cout << r+1 << "\t" << fixed << tau0 << "\t" << K << "\t" << n << "\t" << scientific << setprecision(2) << L << "\t" << L1 << "\t" << L2 << endl;
    }
    return 0;
}

