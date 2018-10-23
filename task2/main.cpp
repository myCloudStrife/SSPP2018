#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <chrono>
#include <papi.h>

constexpr int EFFECTIVE_CACHE_SIZE = 72;
const char *events_name[] = { "Time", "L1_cache_misses", "L2_cache_misses", "Total_cycles",
        "Flop", "TLB_misses" };

int repeat = 1;

using namespace std;

struct Papi_data {
    long long time = 0;
    long_long l1_misses = 0;
    long_long l2_misses = 0;
    long_long total_cycles = 0;
    long_long flop = 0;
    long_long tlb_misses = 0;
    void operator+=(const Papi_data & x) {
        time += x.time;
        l1_misses += x.l1_misses;
        l2_misses += x.l2_misses;
        total_cycles += x.total_cycles;
        flop += x.flop;
        tlb_misses += x.tlb_misses;
    }
};

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

Matrix mulMatrix(Matrix & A, Matrix & B, int mode, Papi_data & value) {
    Matrix C = Matrix(A.rows, B.columns);
    Papi_data tmpvalue;
    long long time1, time2;
    for (int i = 0; i < 2; ++i) {
        if (!i) {
            int events[] = { PAPI_L1_TCM, PAPI_L2_TCM, PAPI_TOT_CYC };
            if (PAPI_start_counters(events, 3) != PAPI_OK) {
                fprintf(stderr, "PAPI_start_counters failed\n");
                exit(1);
            }
        } else {
            int events[] = { PAPI_FP_OPS, PAPI_TLB_DM };
            if (PAPI_start_counters(events, 2) != PAPI_OK) {
                fprintf(stderr, "PAPI_start_counters failed\n");
                exit(1);
            }
        }
        switch (mode) {
        case 0:
            time1 = PAPI_get_real_usec();
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
            time2 = PAPI_get_real_usec();
            break;
        case 1:
            time1 = PAPI_get_real_usec();
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
            time2 = PAPI_get_real_usec();
            break;
        case 2:
            time1 = PAPI_get_real_usec();
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
            time2 = PAPI_get_real_usec();
            break;
        }
        if (!i) {
            if (PAPI_stop_counters(&tmpvalue.l1_misses, 3) != PAPI_OK) {
                fprintf(stderr, "PAPI_stop_counters failed\n");
                exit(1);
            }
            if (time2 < time1) {
                fprintf(stderr, "time error\n");
                exit(1);
            }
            tmpvalue.time = time2 - time1;
        } else {
            if (PAPI_stop_counters(&tmpvalue.flop, 2) != PAPI_OK) {
                fprintf(stderr, "PAPI_stop_counters failed\n");
                exit(1);
            }
            if (time2 < time1) {
                fprintf(stderr, "time error\n");
                exit(1);
            }
            tmpvalue.time = (time2 - time1 + tmpvalue.time) / 2;
        }
    }
    value += tmpvalue;
    return C;
}

void usage(int argc, char **argv) {
    printf("Usage: %s <matrix1> <matrix2> <output_file> <index_mode> [repeat_count]>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 5 && argc != 6) {
        usage(argc, argv);
    }
    int index_mode;
    if (argc == 6) {
        sscanf(argv[5], "%d", &repeat);
    }
    sscanf(argv[4], "%d", &index_mode);
    PAPI_library_init(PAPI_VER_CURRENT);
    Papi_data value;
    Matrix A = Matrix(argv[1]);
    Matrix B = Matrix(argv[2]);
    if (A.columns != B.rows) {
        printf("Incorrect input matrixes\n");
        exit(1);
    }
    Matrix C = mulMatrix(A, B, index_mode, value);
    printf("# size");
    for (int i = 0; i < 6; ++i) {
        printf(" %15s", events_name[i]);
    }
    printf("\n");
    for (int i = 1; i < repeat; ++i) {
        C = mulMatrix(A, B, index_mode, value);
    }
    writeMatrix(argv[3], C);
    printf("%6lu %15lld %15lld %15lld %15lld %15lld %15lld\n", A.rows, value.time / repeat,
            value.l1_misses / repeat, value.l2_misses / repeat, value.total_cycles / repeat,
            value.flop / repeat, value.tlb_misses / repeat);
}
