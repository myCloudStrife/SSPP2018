#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <complex>

using namespace std;

complex<double> * rand_vec_norm(unsigned long long n, unsigned seed) {
    complex<double> * vec = new complex<double>[n];
    double sum = 0;
    for (unsigned long long i = 0; i < n; ++i) {
        vec[i] = complex<double>(1.0 / rand_r(&seed), 1.0 / rand_r(&seed));
        sum += norm(vec[i]);
    }
    sum = sqrt(sum);
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (unsigned long long i = 0; i < n; ++i) {
        vec[i] /= sum;
    }
    return vec;
}

complex<double> * read_vec(const char *filename, unsigned long long & n, int size, int rank) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    unsigned long long vec_fullsize;
    MPI_File_read_all(file, &vec_fullsize, 1, MPI_UNSIGNED_LONG_LONG, NULL);
    n = vec_fullsize / size;
    if (vec_fullsize % size) {
        fprintf(stderr, "Can't divide a vector by %d processors\n", size);
        exit(1);
    }
    complex<double> *vec = new complex<double>[n];
    MPI_File_seek(file, rank * n * sizeof(vec[0]), MPI_SEEK_CUR);
    MPI_File_read_all(file, vec, n * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&file);
    return vec;
}

void write_vec(const char *filename, unsigned long long n, int size, int rank,
        complex<double> *vec) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE,
            MPI_INFO_NULL, &file);
    unsigned long long vec_fullsize = n * size;
    if (!rank) {
        MPI_File_write(file, &vec_fullsize, 1, MPI_UNSIGNED_LONG_LONG, NULL);
    } else {
        MPI_File_seek(file, sizeof(vec_fullsize) + rank * n * sizeof(vec[0]), MPI_SEEK_CUR);
    }
    MPI_File_write_all(file, vec, n * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&file);
}

complex<double> * transform(complex<double> *a, unsigned long long n, int k,
        const complex<double> *u, int size, int rank) {
    complex<double> *b = new complex<double>[n];
    unsigned long long vec_fullsize = n * size;
    int target = (rank + size / (1ull << k)) % size;
    if (target != rank) {
        MPI_Sendrecv(a, n * 2, MPI_DOUBLE, target, 7, b, n * 2, MPI_DOUBLE, target, 7,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (target < rank) {
            for (unsigned long long i = 0; i < n; ++i) {
                b[i] = b[i] * u[0] + a[i] * u[1];
            }
        }
    } else {
        int q = -1;
        while (vec_fullsize) {
            q++;
            vec_fullsize >>= 1;
        }
        k = q - k;
        unsigned long long bit = 1ull << k;
        for (unsigned long long i = 0; i < n; ++i) {
            int u_row = ((i & bit) >> k) << 1;
            b[i] = a[i & ~bit] * u[u_row] + a[i | bit] * u[u_row + 1];
        }
    }
    return b;
}

void usage(int argc, char **argv) {
    printf("Usage: %s <input_file or n> <k> <output_file(optional)>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    unsigned long long n, k;
    if (!rank && argc != 3 && argc != 4) {
        usage(argc, argv);
    }
    complex<double> *a, *b;
    FILE *f = fopen(argv[1], "rb");
    if (f) {
        fclose(f);
        a = read_vec(argv[1], n, size, rank);
    } else {
        sscanf(argv[1], "%llu", &n);
        if ((1ull << n) % size) {
            fprintf(stderr, "Can't divide a vector by %d processors\n", size);
            exit(1);
        }
        n = (1ull << n) / size;
        a = rand_vec_norm(n, MPI_Wtime() + rank);
    }
    sscanf(argv[2], "%llu", &k);

    complex<double> u[4] = {
            M_SQRT1_2, M_SQRT1_2,
            M_SQRT1_2, -M_SQRT1_2
    };

    double time, max_time, start_time = MPI_Wtime();
    b = transform(a, n, k, u, size, rank);
    time = MPI_Wtime() - start_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("%f\n", max_time);
    }
    if (argc == 4) {
        write_vec(argv[3], n, size, rank, a);
    }
    delete [] a;
    delete [] b;
    MPI_Finalize();
    return 0;
}
