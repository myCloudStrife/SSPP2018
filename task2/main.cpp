#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <chrono>
#include <papi.h>

constexpr int EFFECTIVE_CACHE_SIZE = 45;
constexpr int NUM_EVENTS = 2;
const char *modes[] = { "time", "L1_cache_misses", "L2_cache_misses", "Total_cycles" };
//int events[NUM_EVENTS] = { PAPI_L1_TCM, PAPI_L2_TCM };//, PAPI_TOT_CYC };

int repeat = 1;
double work_time[3] = { 0, 0, 0 };

using namespace std;

struct Matrix {
    float *data = nullptr;
    uint64_t rows = 0, columns = 0;
    Matrix(uint64_t rows, uint64_t columns) {
        this->rows = rows;
        this->columns = columns;
        data = new float [rows * columns];
        memset(data, 0, rows * columns * sizeof(data[0]));
    }
    explicit Matrix(const char *filename) {
        ifstream in(filename, ios::binary);
        in.read((char*) &rows, sizeof(char));
        in.read((char*) &rows, sizeof(rows));
        in.read((char*) &columns, sizeof(columns));
        data = new float [rows * columns];
        for (uint64_t i = 0; i < rows * columns; ++i) {
            in.read((char*) &data[i], sizeof(data[0]));
        }
        in.close();
    }
    Matrix(const Matrix & matrix) {
        rows = matrix.rows;
        columns = matrix.columns;
        data = new float [rows * columns];
        memcpy(data, matrix.data, rows * columns * sizeof(data[0]));
    }
    Matrix & operator=(const Matrix & matrix) {
        delete [] data;
        rows = matrix.rows;
        columns = matrix.columns;
        data = new float [rows * columns];
        memcpy(data, matrix.data, rows * columns * sizeof(data[0]));
        return *this;
    }
    float * operator[](uint64_t index) {
        return data + index * columns;
    }
    ~Matrix() {
        delete [] data;
    }
};

void writeMatrix(const char *filename, Matrix & matrix) {
    char type = 'f';
    ofstream out(filename, ios::binary);
    out.write(&type, 1);
    out.write((char *) &matrix.rows, sizeof(matrix.rows));
    out.write((char *) &matrix.columns, sizeof(matrix.columns));
    for (uint64_t i = 0; i < matrix.rows; ++i) {
        for (uint64_t j = 0; j < matrix.columns; ++j) {
            out.write((char *) &matrix[i][j], sizeof(matrix[0][0]));
        }
    }
    out.close();
}

Matrix mulMatrix(Matrix & A, Matrix & B, int mode, long_long *values) {
    Matrix C = Matrix(A.rows, B.columns);
    int event[1];
    long_long tmpvalue;
    switch (mode) {
    case 1:
        event[0] = PAPI_L1_TCM;
        break;
    case 2:
        event[0] = PAPI_L2_TCM;
        break;
    case 3:
        event[0] = PAPI_TOT_CYC;
        break;
    }
    chrono::steady_clock::time_point time1, time2;
    if (mode) {
        if (PAPI_start_counters(event, 1) != PAPI_OK) {
            fprintf(stderr, "PAPI_start_counters failed\n");
            exit(1);
        }
    } else {
        time1 = chrono::steady_clock::now();
    }
    for (uint64_t i = 0; i < A.rows; i += 32) {
        for (uint64_t j = 0; j < B.columns; j += 32) {
            for (uint64_t k = 0; k < A.columns; k += 32) {
                for (uint64_t i2 = i; i2 < min(i + 32, A.rows); ++i2) {
                    for (uint64_t j2 = j; j2 < min(j + 32, B.columns); ++j2) {
                        for (uint64_t k2 = k; k2 < min(k + 32, A.columns); ++k2) {
                            C[i2][j2] += A[i2][k2] * B[k2][j2];
                        }
                    }
                }
            }
        }
    }
    if (mode) {
        if (PAPI_read_counters(&tmpvalue, 1) != PAPI_OK) {
            fprintf(stderr, "PAPI_read_counters failed\n");
            exit(1);
        }
        values[0] += tmpvalue;
        /*if (PAPI_start_counters(event, 1) != PAPI_OK) {
            fprintf(stderr, "PAPI_start_counters failed\n");
            exit(1);
        }*/
    } else {
        time2 = chrono::steady_clock::now();
        work_time[0] += chrono::duration_cast<chrono::duration<double>>(time2 - time1).count();
        time1 = chrono::steady_clock::now();
    }
    for (uint64_t i = 0; i < A.rows; i += 32) {
        for (uint64_t k = 0; k < A.columns; k += 32) {
            for (uint64_t j = 0; j < B.columns; j += 32) {
                for (uint64_t i2 = i; i2 < min(i + 32, A.rows); ++i2) {
                    for (uint64_t k2 = k; k2 < min(k + 32, A.columns); ++k2) {
                        for (uint64_t j2 = j; j2 < min(j + 32, B.columns); ++j2) {
                            C[i2][j2] += A[i2][k2] * B[k2][j2];
                        }
                    }
                }
            }
        }
    }
    if (mode) {
        if (PAPI_read_counters(&tmpvalue, 1) != PAPI_OK) {
            fprintf(stderr, "PAPI_read_counters failed\n");
            exit(1);
        }
        values[1] += tmpvalue;
        /*if (PAPI_start_counters(event, 1) != PAPI_OK) {
            fprintf(stderr, "PAPI_start_counters failed\n");
            exit(1);
        }*/
    } else {
        time2 = chrono::steady_clock::now();
        work_time[1] += chrono::duration_cast<chrono::duration<double>>(time2 - time1).count();
        time1 = chrono::steady_clock::now();
    }
    for (uint64_t i = 0; i < A.rows; i += EFFECTIVE_CACHE_SIZE) {
        for (uint64_t k = 0; k < A.columns; k += EFFECTIVE_CACHE_SIZE) {
            for (uint64_t j = 0; j < B.columns; j += EFFECTIVE_CACHE_SIZE) {
                for (uint64_t i2 = i; i2 < min(i + EFFECTIVE_CACHE_SIZE, A.rows); ++i2) {
                    for (uint64_t k2 = k; k2 < min(k + EFFECTIVE_CACHE_SIZE, A.columns); ++k2) {
                        for (uint64_t j2 = j; j2 < min(j + EFFECTIVE_CACHE_SIZE, B.columns); ++j2) {
                            C[i2][j2] += A[i2][k2] * B[k2][j2];
                        }
                    }
                }
            }
        }
    }
    if (mode) {
        if (PAPI_stop_counters(&tmpvalue, 1) != PAPI_OK) {
            fprintf(stderr, "PAPI_stop_counters failed\n");
            exit(1);
        }
        values[2] += tmpvalue;
        //if (PAPI_)
    } else {
        time2 = chrono::steady_clock::now();
        work_time[2] += chrono::duration_cast<chrono::duration<double>>(time2 - time1).count();
    }
    return C;
}

void usage(int argc, char **argv) {
    printf("Usage: %s <matrix1> <matrix2> <output_filename> <mode> [repeat_count]>\n", argv[0]);
    exit(1);
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
    PAPI_library_init(PAPI_VER_CURRENT);
    long_long values[3] = { 0, 0, 0 };
    Matrix A = Matrix(argv[1]);
    Matrix B = Matrix(argv[2]);
    Matrix C = mulMatrix(A, B, mode, values);
    printf("# mode %15s\n", modes[mode]);
    for (int i = 1; i < repeat; ++i) {
        C = mulMatrix(A, B, mode, values);
    }
    writeMatrix(argv[3], C);
    if (mode) {
        cout << setw(6) << 0 << ' ' << setw(15) << values[0] / repeat << endl;
        cout << setw(6) << 1 << ' ' << setw(15) << values[1] / repeat << endl;
        cout << setw(6) << 2 << ' ' << setw(15) << values[2] / repeat << endl;
    } else {
        cout << setw(6) << 0 << ' ' << setw(15) << work_time[0] / repeat << endl;
        cout << setw(6) << 1 << ' ' << setw(15) << work_time[1] / repeat << endl;
        cout << setw(6) << 2 << ' ' << setw(15) << work_time[2] / repeat << endl;
    }
}
