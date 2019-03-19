#include <mpi.h>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <climits>

constexpr int EXPERIMENT_COUNT = 10;

using namespace std;

double normal_dis_gen(unsigned & seed) {
    double S = 0;
    for (int i = 0; i < 12; ++i) {
        S += (double) rand_r(&seed) / RAND_MAX;
    }
    return S - 6;
}

complex<double> * rand_vec_norm(unsigned long long n, unsigned seed) {
    complex<double> * vec = new complex<double>[n];
    double sum = 0;
#pragma omp parallel
    {
        seed += omp_get_thread_num();
#pragma omp for reduction(+: sum)
        for (unsigned long long i = 0; i < n; ++i) {
            vec[i] = complex<double>(1.0 / rand_r(&seed), 1.0 / rand_r(&seed));
            sum += norm(vec[i]);
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sum = sqrt(sum);
#pragma omp parallel for
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
        unsigned long long need_count = n * 2;
        int sendrecvcount = INT_MAX - 1;
        unsigned long long i = 0;
        while (need_count > INT_MAX) {
            MPI_Sendrecv(a + i * sendrecvcount / 2, sendrecvcount, MPI_DOUBLE, target, 7,
                    b + i * sendrecvcount / 2, sendrecvcount, MPI_DOUBLE, target, 7,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            i++;
            need_count -= sendrecvcount;
        }
        MPI_Sendrecv(a + i * sendrecvcount / 2, need_count, MPI_DOUBLE, target, 7,
                b + i * sendrecvcount / 2, need_count, MPI_DOUBLE, target, 7,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (target < rank) {
#pragma omp parallel for
            for (unsigned long long i = 0; i < n; ++i) {
                b[i] = b[i] * u[2] + a[i] * u[3];
            }
        } else {
#pragma omp parallel for
            for (unsigned long long i = 0; i < n; ++i) {
                b[i] = a[i] * u[0] + b[i] * u[1];
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
#pragma omp parallel for
        for (unsigned long long i = 0; i < n; ++i) {
            int u_row = ((i & bit) >> k) * 2;
            b[i] = a[i & ~bit] * u[u_row] + a[i | bit] * u[u_row + 1];
        }
    }
    return b;
}

void usage(int argc, char **argv) {
    printf("Usage: %s <input_file or n> <eps> <thread_num> <output_file(optional)>\n", argv[0]);
    exit(1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    unsigned long long n;
    int num_threads, q;
    double epsilon;
    if (!rank && argc != 4 && argc != 5) {
        usage(argc, argv);
    }
    sscanf(argv[2], "%lf", &epsilon);
    sscanf(argv[3], "%d", &num_threads);
    omp_set_num_threads(num_threads);
    complex<double> *ideal, *noise, *tmp;
    FILE *f = fopen(argv[1], "rb");
    if (f) {
        fclose(f);
        ideal = read_vec(argv[1], n, size, rank);
        unsigned long long vec_fullsize = n * size;
        q = -1;
        while (vec_fullsize) {
            q++;
            vec_fullsize >>= 1;
        }
    } else {
        sscanf(argv[1], "%d", &q);
        if ((1ull << q) % size) {
            fprintf(stderr, "Can't divide a vector by %d processors\n", size);
            exit(1);
        }
        n = (1ull << q) / size;
        double time;
        if (!rank) {
            time = MPI_Wtime();
        }
        MPI_Bcast(&time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        ideal = rand_vec_norm(n, time + rank * num_threads);
    }

    complex<double> u[4] = {
            M_SQRT1_2, M_SQRT1_2,
            M_SQRT1_2, -M_SQRT1_2
    };
    complex<double> u_noise[4];

    unsigned seed;
    if (!rank) {
        seed = MPI_Wtime() * 100;
    }
    MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    double time = 0, max_time, start_time;

    for (int i = 0; i < EXPERIMENT_COUNT; ++i) {
        if (i) {
            delete [] ideal;
            delete [] noise;
            ideal = rand_vec_norm(n, seed + rank * num_threads);
            if (!rank) {
                seed = MPI_Wtime() * 100000;
            }
            MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        }
        noise = ideal;
        for (int k = 1; k <= q; ++k) {
            tmp = transform(ideal, n, k, u, size, rank);
            if (k > 1) {
                delete [] ideal;
            }
            ideal = tmp;

            double phi = epsilon * normal_dis_gen(seed);
            u_noise[0] = u[0] * cos(phi) - u[1] * sin(phi);
            u_noise[1] = u[0] * sin(phi) + u[1] * cos(phi);
            u_noise[2] = u[2] * cos(phi) - u[3] * sin(phi);
            u_noise[3] = u[2] * sin(phi) + u[3] * cos(phi);

            start_time = MPI_Wtime();
            tmp = transform(noise, n, k, u_noise, size, rank);
            time += MPI_Wtime() - start_time;

            delete [] noise;
            noise = tmp;
        }
        double fidelity_loc = 0, fidelity;
        for (unsigned long long i = 0; i < n; ++i) {
            fidelity_loc += abs(ideal[i] * conj(noise[i]));
        }
        MPI_Reduce(&fidelity_loc, &fidelity, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (!rank) {
            printf("Loss: %f\n", 1.0 - fidelity);
        }
    }
    time /= EXPERIMENT_COUNT;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("Spend time: %f\n", max_time);
    }

    if (argc == 5) {
        write_vec(argv[4], n, size, rank, noise);
    }
    delete [] ideal;
    delete [] noise;
    MPI_Finalize();
    return 0;
}
