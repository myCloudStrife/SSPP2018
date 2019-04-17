#include <iostream>
#include <complex>

void usage(int argc, char **argv) {
    printf("Usage: %s <q> <output_file>\n"
            "   q - number of qubits\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        usage(argc, argv);
    }
    int q;
    sscanf(argv[1], "%d", &q);
    unsigned seed = clock();
    uint64_t n = 1ull << q;
    std::complex<double> *vec = new std::complex<double>[n];
    double sum = 0;
    for (uint64_t i = 0; i < n; ++i) {
        vec[i] = std::complex<double>(1.0 / (rand_r(&seed) - rand_r(&seed)),
                1.0 / (rand_r(&seed) - rand_r(&seed)));
        sum += norm(vec[i]);
    }
    sum = sqrt(sum);
    for (uint64_t i = 0; i < n; ++i) {
        vec[i] /= sum;
    }
    FILE *f = fopen(argv[2], "wb");
    fwrite(&n, sizeof(n), 1, f);
    fwrite(vec, sizeof(vec[0]), n, f);
    fclose(f);
    delete [] vec;
}
