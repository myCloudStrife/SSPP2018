#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <stdint.h>

using namespace std;

struct Matrix {
    double *data;
    uint64_t rows, columns;
    Matrix(uint64_t rows, uint64_t columns) {
        this->rows = rows;
        this->columns = columns;
        data = new double [rows * columns];
        memset(data, 0, rows * columns * sizeof(data[0]));
    }
    Matrix(const Matrix & matrix) {
        rows = matrix.rows;
        columns = matrix.columns;
        data = new double [rows * columns];
        memcpy(data, matrix.data, rows * columns * sizeof(data[0]));
    }
    Matrix & operator=(const Matrix & matrix) {
        delete [] data;
        rows = matrix.rows;
        columns = matrix.columns;
        data = new double [rows * columns];
        memcpy(data, matrix.data, rows * columns * sizeof(data[0]));
        return *this;
    }
    double * operator[](uint64_t index) {
        return data + index * columns;
    }
    ~Matrix() {
        delete [] data;
    }
};

void readmn(uint64_t & m, uint64_t & n, FILE *matrix_file, FILE *vector_file) {
    char type;
    int rank;
    uint64_t tmp;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (fread(&type, sizeof(type), 1, matrix_file) != 1) {
        exit(1);
    }
    if (!rank && type != 'd') {
        printf("Incorrect input\n");
        exit(1);
    }
    if (fread(&type, sizeof(type), 1, vector_file) != 1) {
        exit(1);
    };
    if (!rank && type != 'd') {
        printf("Incorrect input\n");
        exit(1);
    }
    if (fread(&m, sizeof(m), 1, matrix_file) != 1) {
        exit(1);
    }
    if (fread(&n, sizeof(n), 1, matrix_file) != 1) {
        exit(1);
    }
    if (fread(&tmp, sizeof(tmp), 1, vector_file) != 1) {
        exit(1);
    }
    if (!rank && tmp != n) {
        printf("Incorrect input\n");
        exit(1);
    }
    if (fread(&tmp, sizeof(tmp), 1, vector_file) != 1) {
        exit(1);
    }
    if (!rank && tmp != 1) {
        printf("Incorrect input\n");
        exit(1);
    }
}

void write_vector(double *vector, uint64_t size, const char *filename) {
    FILE *f = fopen(filename, "wb");
    char type = 'd';
    fwrite(&type, sizeof(type), 1, f);
    fwrite(&size, sizeof(size), 1, f);
    uint64_t one = 1;
    fwrite(&one, sizeof(one), 1, f);
    fwrite(vector, sizeof(vector[0]), size, f);
    fclose(f);
}

void usage(int argc, char **argv) {
    printf("Usage: %s <matrix> <vector> <output>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    if (argc != 4) {
        usage(argc, argv);
    }
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    FILE *matrix_file = fopen(argv[1], "rb");
    FILE *vector_file = fopen(argv[2], "rb");
    uint64_t n, m;
    readmn(m, n, matrix_file, vector_file);
    double time = 0;
    double *b, *c;
    c = new double [m];
    if (m >= n) {
        b = new double [n];
        if (fread(b, sizeof(b[0]), n, vector_file) != n) {
            exit(1);
        }
        uint64_t left = rank * m / size, right = (rank + 1) * m / size;
        fseek(matrix_file, left * n * sizeof(double), SEEK_CUR);
        Matrix A(right - left, n);
        if (fread(A[0], sizeof(A[0][0]), n * (right - left), matrix_file) != n * (right - left)) {
            exit(1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        time = MPI_Wtime();
        for (uint64_t i = 0; i < right - left; ++i) {
            c[i] = 0;
            for (uint64_t j = 0; j < n; ++j) {
                c[i] += A[i][j] * b[j];
            }
        }
        if (m % size) {
            int *ranges = NULL, *displs = NULL;
            if (!rank) {
                ranges = new int [size];
                displs = new int [size];
                for (int i = 0; i < size; ++i) {
                    ranges[i] = (i + 1) * m / size - i * m / size;
                    displs[i] = i ? ranges[i - 1] + displs[i - 1] : 0;
                }
            }
            MPI_Gatherv(rank ? c : MPI_IN_PLACE, right - left, MPI_DOUBLE, c, ranges,
                    displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (!rank) {
                delete [] ranges;
                delete [] displs;
            }
        } else {
            MPI_Gather(rank ? c : MPI_IN_PLACE, right - left, MPI_DOUBLE, c, right - left,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    } else {
        uint64_t left = rank * n / size, right = (rank + 1) * n / size;
        uint64_t range = right - left;
        b = new double [range];
        Matrix A(m, range);
        fseek(vector_file, left * sizeof(double), SEEK_CUR);
        if (fread(b, sizeof(b[0]), range, vector_file) != range) {
            exit(1);
        }
        fseek(matrix_file, left * sizeof(double), SEEK_CUR);
        for (uint64_t i = 0; i < m; ++i) {
            if (fread(A[i], sizeof(A[0][0]), range, matrix_file) != range) {
                exit(1);
            }
            fseek(matrix_file, (n - range) * sizeof(double), SEEK_CUR);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        time = MPI_Wtime();
        for (uint64_t i = 0; i < m; ++i) {
            c[i] = 0;
            for (uint64_t j = 0; j < range; ++j) {
                c[i] += A[i][j] * b[j];
            }
        }
        MPI_Reduce(rank ? c : MPI_IN_PLACE, rank ? NULL : c, m, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD);
    }
    time = MPI_Wtime() - time;
    if (!rank) {
        write_vector(c, m, argv[3]);
    }
    double max_time, all_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time, &all_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("Max time: %fs\n", max_time);
        printf("Sum time: %fs\n", all_time);
    }
    fclose(vector_file);
    fclose(matrix_file);
    delete [] b;
    delete [] c;
    MPI_Finalize();
    return 0;
}
