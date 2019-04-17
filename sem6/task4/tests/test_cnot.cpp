#include <iostream>
#include <complex>

void usage(int argc, char **argv) {
    printf("Usage: %s <input_file> <k> <l> <output_file>\n"
            "   k, l - qubit indexes\n"
            "File format:\n"
            "   uint64_t n, //vector size, not number of qubits\n"
            "   complex<double> data[n]\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 5) {
        usage(argc, argv);
    }
    int k, l;
    sscanf(argv[2], "%d", &k);
    sscanf(argv[3], "%d", &l);
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
    l = q - l;
    uint64_t bitk = 1ull << k;
    uint64_t bitl = 1ull << l;
    for (uint64_t i = 0; i < n; ++i) {
        if (i & bitk) {
            transformed[i] = source[i ^ bitl];
        } else {
            transformed[i] = source[i];
        }
    }

    f = fopen(argv[4], "wb");
    fwrite(&n, sizeof(n), 1, f);
    fwrite(transformed, sizeof(transformed[0]), n, f);
    fclose(f);

    delete [] source;
    delete [] transformed;
}
