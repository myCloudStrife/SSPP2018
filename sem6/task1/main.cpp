#include <omp.h>
#include <iostream>
#include <vector>
#include <complex>

using namespace std;

vector<complex<double>> rand_vec_norm(unsigned long long n) {
    vector<complex<double>> vec(n);
    double sum = 0;
    unsigned time = omp_get_wtime() * 100;
#pragma omp parallel
    {
        unsigned seed = time + omp_get_thread_num();
#pragma omp for reduction(+: sum)
        for (unsigned long long i = 0; i < n; ++i) {
            vec[i] = complex<double>(rand_r(&seed), rand_r(&seed));
            sum += norm(vec[i]);
        }
    }
    sum = sqrt(sum);
#pragma omp parallel for
    for (unsigned long long i = 0; i < n; ++i) {
        vec[i] /= sum;
    }
    return vec;
}

vector<complex<double>> transform(const vector<complex<double>> & a, int k,
        const vector<complex<double>> & u) {
    vector<complex<double>> b(a.size());
    unsigned long long bit = 1ull << k;
#pragma omp parallel for
    for (unsigned long long i = 0; i < a.size(); ++i) {
        b[i] = a[i & ~bit] * u[((i & bit) >> k) << 1] + a[i | bit] * u[(((i & bit) >> k) << 1) + 1];
    }
    return b;
}

void print_vec(const vector<complex<double>> & a) {
    for (unsigned long long i = 0; i < a.size(); ++i) {
        printf("%f %f\n", a[i].real(), a[i].imag());
    }
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

    double start_time = omp_get_wtime();
    vector<complex<double>> a = rand_vec_norm(1ull << n);
    print_vec(a);
    print_vec(transform(a, k, u));
    printf("%.10f\n", omp_get_wtime() - start_time);
}
