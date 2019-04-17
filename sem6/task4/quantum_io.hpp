#ifndef SSPP_SEM6_TASK4_QUANTUM_IO_HPP_
#define SSPP_SEM6_TASK4_QUANTUM_IO_HPP_

#include <mpi.h>
#include <omp.h>
#include <complex>

constexpr int MESSAGETAG = 7345;

std::complex<double> * rand_vec_norm(uint64_t n) {
    static int rank = -1;
    if (rank == -1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    unsigned seed;
    if (!rank) {
        seed = MPI_Wtime();
    }
    MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    std::complex<double> * vec = new std::complex<double>[n];
    double sum = 0;
#pragma omp parallel
    {
        seed += rank * omp_get_num_threads() + omp_get_thread_num();
#pragma omp for reduction(+: sum)
        for (uint64_t i = 0; i < n; ++i) {
            vec[i] = std::complex<double>(1.0 / rand_r(&seed),
                    1.0 / rand_r(&seed));
            sum += norm(vec[i]);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum = sqrt(sum);
#pragma omp parallel for
    for (uint64_t i = 0; i < n; ++i) {
        vec[i] /= sum;
    }
    return vec;
}

std::complex<double> * read_vec(const char *filename, uint64_t *n, int size,
        int rank) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL,
            &file);
    uint64_t vec_fullsize;
    MPI_File_read_all(file, &vec_fullsize, 1, MPI_UNSIGNED_LONG_LONG, NULL);
    *n = vec_fullsize / size;
    if (vec_fullsize % size) {
        fprintf(stderr, "Can't divide a vector by %d processors\n", size);
        exit(1);
    }
    std::complex<double> *vec = new std::complex<double>[*n];
    MPI_File_seek(file, rank * *n * sizeof(vec[0]) + sizeof(vec_fullsize),
            MPI_SEEK_SET);
    MPI_File_read_all(file, vec, *n * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&file);
    return vec;
}

void write_vec(const char *filename, uint64_t n, int size, int rank,
        std::complex<double> *vec) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE,
            MPI_INFO_NULL, &file);
    uint64_t vec_fullsize = n * size;
    if (!rank) {
        MPI_File_write(file, &vec_fullsize, 1, MPI_UNSIGNED_LONG_LONG, NULL);
    } else {
        MPI_File_seek(file, sizeof(vec_fullsize) + rank * n * sizeof(vec[0]),
                MPI_SEEK_CUR);
    }
    MPI_File_write_all(file, vec, n * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&file);
}

void print_vec(std::complex<double> *vec, uint64_t n, int size, int rank) {
    char ready = 0;
    if (rank) {
        MPI_Recv(&ready, 1, MPI_CHAR, rank - 1, MESSAGETAG, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    }
    for (uint64_t i = 0; i < n; ++i) {
        printf("%f + %fi\n", vec[i].real(), vec[i].imag());
    }
    if (rank != size - 1) {
        MPI_Send(&ready, 1, MPI_CHAR, rank + 1, MESSAGETAG, MPI_COMM_WORLD);
    }
}

#endif /* SSPP_SEM6_TASK4_QUANTUM_IO_HPP_ */
