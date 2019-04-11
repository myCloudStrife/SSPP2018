#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>

using namespace std;

void usage(int argc, char **argv) {
    printf("Usage: %s <matrix1> <matrix2>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        usage(argc, argv);
    }
    char type1, type2;
    uint64_t m1, m2, n1, n2;
    ifstream in1, in2;
    in1.open(argv[1], ios::binary);
    in2.open(argv[2], ios::binary);
    in1.read(&type1, 1);
    in2.read(&type2, 1);
    in1.read((char *) &m1, sizeof(m1));
    in2.read((char *) &m2, sizeof(m2));
    in1.read((char *) &n1, sizeof(n1));
    in2.read((char *) &n2, sizeof(n2));
    if (type1 != type2 || m1 != m2 || n1 != n2) {
        printf("Matrix %s and %s not equal\n", argv[1], argv[2]);
        return 1;
    }
    for (uint64_t i = 0; i < m1 * n1; ++i) {
        if (type1 == 'd') {
            double x1, x2;
            in1.read((char *) &x1, sizeof(x1));
            in2.read((char *) &x2, sizeof(x2));
            if (fabs(x1 - x2) > 10e-10) {
                printf("Matrix %s and %s not equal, %f != %f\n", argv[1], argv[2], x1, x2);
                return 1;
            }
        } else {
            float x1, x2;
            in1.read((char *) &x1, sizeof(x1));
            in2.read((char *) &x2, sizeof(x2));
            if (fabs(x1 - x2) > 10e-5) {
                printf("Matrix %s and %s not equal, %f != %f\n", argv[1], argv[2], x1, x2);
                return 1;
            }
        }
    }
    printf("Matrix %s and %s are equal\n", argv[1], argv[2]);
    in1.close();
    in2.close();
    return 0;
}
