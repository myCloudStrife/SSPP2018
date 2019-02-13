#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdint.h>

using namespace std;

double rand_double() {
    return (double) rand() / rand() - (double) rand() / rand();
}

void usage(int argc, char **argv) {
    printf("Usage: %s <type> <M> <N> <filename>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 5 || (argv[1][0] != 'f' && argv[1][0] != 'd')) {
        usage(argc, argv);
    }
    ofstream out(argv[4], ios::binary);
    uint64_t m, n;
    sscanf(argv[2], "%llu", (long long unsigned *) &m);
    sscanf(argv[3], "%llu", (long long unsigned *) &n);
    srand(time(NULL));
    out.write(argv[1], 1);
    out.write((char *)&m, sizeof(m));
    out.write((char *)&n, sizeof(n));
    for (uint64_t i = 0; i < m; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            if (argv[1][0] == 'd') {
                double x = rand_double();
                out.write((char *) &x, sizeof(x));
            } else {
                float x = rand_double();
                out.write((char *) &x, sizeof(x));
            }
        }
    }
    out.close();
    return 0;
}
