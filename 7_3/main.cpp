#include <iostream>
#include <cmath>
using namespace std;

//double x[10] = {1.0986, 1.3863, 1.6094, 1.7918, 1.9459, 2.0794, 2.1972, 2.3026, 2.3979, 2.4849};
//double y[10] = {13.811, 16.850, 19.176, 21.060, 22.642, 24.005, 25.203, 26.270, 27.233, 28.109};

// Albina Yakhina
//double x[10] = {2.9524, 3.5545, 3.9067, 4.1565, 4.3504, 4.5087, 4.6426, 4.7586, 4.8609, 4.9524};
//double y[10] = {18.3344, 20.3390, 21.4480, 22.2105, 22.7894, 23.2546, 23.6429, 23.9755, 24.2663, 24.5243};

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
double x[10] = {2.4000, 3.0020, 3.3542, 3.6041, 3.7979, 3.9563, 4.0902, 4.2062, 4.3085, 4.4000};
double y[10] = {12.9691, 13.9516, 14.4672, 14.8118, 15.0684, 15.2717, 15.4395, 15.5819, 15.7054, 15.8143};


int main() {
    int Nexp = 10;
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

    } while (dev > eps);

    cout << "Newton method. Tatyana Vorobyova group BPOi-16" << endl;
    cout << "tau0: " << tau0 << endl;
    cout << "K: " << K << endl;
    cout << "n: " << n << endl;
    cout << "Function: " << F << endl;
    cout << "Derivatives: " << abs(F1) << "  " << abs(F2) << "  " << abs(F3) << endl;
    return 0;
}
