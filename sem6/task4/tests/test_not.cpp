#include <iostream>
#include <complex>

void usage(int argc, char **argv) {
    printf("Usage: %s <input_file> <k> <output_file>\n"
            "   k - qubit index\n"
            "File format:\n"
            "   uint64_t n, //vector size, not number of qubits\n"
            "   complex<double> data[n]\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 4) {
        usage(argc, argv);
    }
    int k;
    sscanf(argv[2], "%d", &k);
    uint64_t n;
    FILE *f = fopen(argv[1], "rb");
    if (fread(&n, sizeof(n), 1, f) != 1) {
        exit(1);
    }
    std::complex<double> *source = new std::complex<double>[n];
    std::complex<double> *transformed = new std::complex<double>[n];
    if (fread(source, sizeof(source[0]), n, f) != n) {
        exit(1);
    }
    fclose(f);

    int q = -1;
    uint64_t tmpn = n;
    while (tmpn) {
        q++;
        tmpn >>= 1;
    }
    k = q - k;
    uint64_t bit = 1ull << k;
    for (uint64_t i = 0; i < n; ++i) {
        transformed[i] = source[i ^ bit];
    }

    f = fopen(argv[3], "wb");
    fwrite(&n, sizeof(n), 1, f);
    fwrite(transformed, sizeof(transformed[0]), n, f);
    fclose(f);

    delete [] source;
    delete [] transformed;
}
