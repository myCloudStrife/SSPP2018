#include <stdlib.h>
#include <stdio.h>

int main(int argc, char **argv) {
    unsigned long long n;
    FILE *f = fopen(argv[1], "rb");
    if (fread(&n, sizeof(n), 1, f) != 1) {
        exit(1);
    }
    double d[2];
    for (unsigned long long i = 0; i < n; ++i) {
        if (fread(d, sizeof(d[0]), 2, f) != 2) {
            exit(1);
        }
        printf("%f + %fi\n", d[0], d[1]);
    }
    fclose(f);
    return 0;
}
