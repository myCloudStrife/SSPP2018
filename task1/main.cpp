#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

using namespace std;

int repeat = 1;
double work_time = 0;

struct Info {
    uint64_t n, m;
    char type;
};

void usage(int argc, char **argv) {
    printf("Usage: %s <matrix1> <matrix2> <output_filename> <mode> [repeat_count]>\n", argv[0]);
    exit(1);
}

template <typename T>
T ** scanMatrix(ifstream & in, Info & matrix_info) {
    in.read((char *) &matrix_info.m, sizeof(matrix_info.m));
    in.read((char *) &matrix_info.n, sizeof(matrix_info.n));
    T *data = new T [matrix_info.m * matrix_info.n];
    T **matrix = new T * [matrix_info.m];
    for (uint64_t i = 0; i < matrix_info.m; ++i) {
        matrix[i] = data + i * matrix_info.n;
    }
    for (uint64_t i = 0; i < matrix_info.m; ++i) {
        for (uint64_t j = 0; j < matrix_info.n; ++j) {
            in.read((char *) &matrix[i][j], sizeof(T));
        }
    }
    return matrix;
}

template <typename T>
void writeMatrix(T **matrix, char type, uint64_t m, uint64_t n, const char *file) {
    ofstream out(file, ios::binary);
    out.write(&type, 1);
    out.write((char *)&m, sizeof(m));
    out.write((char *)&n, sizeof(n));
    for (uint64_t i = 0; i < m; ++i) {
        for (uint64_t j = 0; j < n; ++j) {
            out.write((char *) &matrix[i][j], sizeof(T));
        }
    }
    out.close();
}

template <typename T>
T ** mulMatrix(T **A, Info & A_info, T **B, Info & B_info, int mode) {
    T *data = new T[A_info.m * B_info.n];
    T **C = new T *[A_info.m];
    for (uint64_t i = 0; i < A_info.m; ++i) {
        C[i] = data + i * B_info.n;
    }
    chrono::steady_clock::time_point time1, time2;
    switch(mode) {
    case 0:
        time1 = chrono::steady_clock::now();
        for (uint64_t i = 0; i < A_info.m; ++i) {
            for (uint64_t j = 0; j < B_info.n; ++j) {
                for (uint64_t k = 0; k < A_info.n; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        time2 = chrono::steady_clock::now();
        break;
    case 1:
        time1 = chrono::steady_clock::now();
        for (uint64_t i = 0; i < A_info.m; ++i) {
            for (uint64_t k = 0; k < A_info.n; ++k) {
                for (uint64_t j = 0; j < B_info.n; ++j) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        time2 = chrono::steady_clock::now();
        break;
    case 2:
        time1 = chrono::steady_clock::now();
        for (uint64_t k = 0; k < A_info.n; ++k) {
            for (uint64_t i = 0; i < A_info.m; ++i) {
                for (uint64_t j = 0; j < B_info.n; ++j) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        time2 = chrono::steady_clock::now();
        break;
    case 3:
        time1 = chrono::steady_clock::now();
        for (uint64_t j = 0; j < B_info.n; ++j) {
            for (uint64_t i = 0; i < A_info.m; ++i) {
                for (uint64_t k = 0; k < A_info.n; ++k) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        time2 = chrono::steady_clock::now();
        break;
    case 4:
        time1 = chrono::steady_clock::now();
        for (uint64_t j = 0; j < B_info.n; ++j) {
            for (uint64_t k = 0; k < A_info.n; ++k) {
                for (uint64_t i = 0; i < A_info.m; ++i) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        time2 = chrono::steady_clock::now();
        break;
    case 5:
        time1 = chrono::steady_clock::now();
        for (uint64_t k = 0; k < A_info.n; ++k) {
            for (uint64_t j = 0; j < B_info.n; ++j) {
                for (uint64_t i = 0; i < A_info.m; ++i) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        time2 = chrono::steady_clock::now();
        break;
    }
    work_time += chrono::duration_cast<chrono::duration<double>>(time2 - time1).count();
    return C;
}

template <typename T>
void deleteMatrix(T **matrix) {
    delete [] matrix[0];
    delete [] matrix;
}

int main(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        usage(argc, argv);
    }
    int mode;
    if (argc == 6) {
        sscanf(argv[5], "%d", &repeat);
    }
    sscanf(argv[4], "%d", &mode);
    ifstream in(argv[1], ios::binary);
    Info A_info, B_info;
    void *A, *B, *C;
    in.read(&A_info.type, sizeof(A_info.type));
    if (A_info.type == 'd') {
        A = scanMatrix<double>(in, A_info);
    } else {
        A = scanMatrix<float>(in, A_info);
    }
    in.close();
    in.open(argv[2], ios::binary);
    in.read(&B_info.type, sizeof(B_info.type));
    if (B_info.type == 'd') {
        B = scanMatrix<double>(in, B_info);
    } else {
        B = scanMatrix<float>(in, B_info);
    }
    in.close();
    if (A_info.type != B_info.type || A_info.n != B_info.m) {
        printf("Incorrect input matrixes\n");
        exit(1);
    }
    if (A_info.type == 'd') {
        for (int i = 0; i < repeat; ++i) {
            C = mulMatrix((double **) A, A_info, (double **) B, B_info, mode);
        }
        cout << mode << ' ' << work_time / repeat << endl;
        writeMatrix((double **) C, A_info.type, A_info.m, B_info.n, argv[3]);
        deleteMatrix((double **) A);
        deleteMatrix((double **) B);
        deleteMatrix((double **) C);
    } else {
        for (int i = 0; i < repeat; ++i) {
            C = mulMatrix((float **) A, A_info, (float **) B, B_info, mode);
        }
        cout << mode << ' ' << work_time / repeat << endl;
        writeMatrix((float **) C, A_info.type, A_info.m, B_info.n, argv[3]);
        deleteMatrix((float **) A);
        deleteMatrix((float **) B);
        deleteMatrix((float **) C);
    }
    return 0;
}
