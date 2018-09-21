#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void usage(int argc, char **argv) {
    printf("Usage: %s <input_binary_file> <output_text_file>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 3) {
        usage(argc, argv);
        exit(1);
    }
    ifstream in(argv[1], ios::binary);
    ofstream out(argv[2]);
    uint64_t m, n;
    char type;
    in.read(&type, sizeof(type));
    in.read((char *) &m, sizeof(m));
    in.read((char *) &n, sizeof(n));
    for (uint64_t i = 0; i < m; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            if (type == 'd') {
                double x;
                in.read((char *) &x, sizeof(x));
                out << setw(10) << setprecision(4) << x << (j == n - 1 ? '\n' : ' ');
            } else {
                float x;
                in.read((char *) &x, sizeof(x));
                out << setw(10) << setprecision(4) << x << (j == n - 1 ? '\n' : ' ');
            }
        }
    }
    in.close();
    out.close();
}
