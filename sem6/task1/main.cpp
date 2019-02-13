#include <omp.h>
#include <iostream>
#include <vector>
#include <complex>

using namespace std;

vector<complex<double>> rand_vec_norm(int n) {
    vector<complex<double>> vec(n);
    double sum = 0;
    long seed = omp_get_wtime() * 100;
    printf("%ld\n", seed);
    srand48(seed);
#pragma omp parallel for reduction(+: sum)
    for (int i = 0; i < n; ++i) {
        //vec[i] = complex<double>(drand48(), drand48());
        vec[i] = i;
        sum += norm(vec[i]);
    }
    sum = sqrt(sum);
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        vec[i] /= sum;
    }
    return vec;
}

vector<complex<double>> transform(vector<complex<double>> a, int k, vector<complex<double>> u) {
    vector<complex<double>> b(a.size());
    unsigned long long bit = 1 << k;
#pragma omp parallel for
    for (unsigned long long i = 0; i < a.size(); ++i) {
        b[i] = a[i & ~bit] * u[(i & bit) >> (k - 1)] + a[i | bit] * u[((i & bit) >> (k - 1)) + 1];
    }
    return b;
}

int main(int argc, char **argv) {
    int n, k, num_threads;
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &k);
    sscanf(argv[3], "%d", &num_threads);
    omp_set_num_threads(num_threads);
    k = n - k;

    vector<complex<double>> u = {
            M_SQRT1_2, M_SQRT1_2,
            M_SQRT1_2, -M_SQRT1_2 };

    vector<complex<double>> a = rand_vec_norm(1 << n);
    double sum = 0;
    for (int i = 0; i < 1 << n; ++i) {
        sum += norm(a[i]);
        printf("%f %f\n", a[i].real(), a[i].imag());
    }
    printf("%f\n", sum);
    vector<complex<double>> b = transform(a, k, u);
    for (int i = 0; i < 1 << n; ++i) {
        printf("%f %f\n", b[i].real(), b[i].imag());
    }
}
