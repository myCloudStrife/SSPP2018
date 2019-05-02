#include "../task4/quantum_transform.hpp"
#include "../task4/quantum_io.hpp"

typedef std::complex<double> complexd;

unsigned index_reverse(unsigned index, unsigned bit_count) {
    unsigned ret = 0;
    if (!bit_count) {
        return index;
    }
    unsigned bit = 1 << (bit_count - 1);
    while (index) {
        if (index & 1) {
            ret ^= bit;
        }
        index >>= 1;
        bit >>= 1;
    }
    return ret;
}

void writeFourier(const char *filename, uint64_t n, int size, int rank,
        complexd *vec) {
    MPI_File file;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE,
            MPI_INFO_NULL, &file);
    uint64_t vec_fullsize = n * size;
    if (!rank) {
        MPI_File_write(file, &vec_fullsize, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_STATUS_IGNORE);
    }
    int logsize = -1, tmpsize = size;
    while (tmpsize) {
        logsize++;
        tmpsize >>= 1;
    }
    int q = -1;
    while (vec_fullsize) {
        q++;
        vec_fullsize >>= 1;
    }
    for (unsigned i = 0; i < n; ++i) {
        unsigned ri = index_reverse(i, q - logsize);
        if (ri > i) {
            std::swap(vec[i], vec[ri]);
        }
    }
    MPI_Datatype filetype;
    int array_size = size * 2;
    int start_array = index_reverse(rank, logsize) * 2;
    int subsize_array = 2;
    MPI_Type_create_subarray(1, &array_size, &subsize_array, &start_array,
            MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(file, sizeof(n), MPI_DOUBLE, filetype, "native",
            MPI_INFO_NULL);
    MPI_File_write_all(file, vec, n * 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Type_free(&filetype);
    MPI_File_close(&file);
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t n;
    int num_threads;
    sscanf(argv[2], "%d", &num_threads);
    omp_set_num_threads(num_threads);
    complexd *a;
    FILE *f = fopen(argv[1], "rb");
    int q;
    if (f) {
        fclose(f);
        a = read_vec(argv[1], &n, size, rank);
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
        a = rand_vec_norm(n);
    }
    double start_time, time, max_time;
    complexd *tmp;
    start_time = MPI_Wtime();
    for (int i = 1; i <= q; ++i) {
        tmp = transform_adamar(a, n, i);
        delete [] a;
        a = tmp;
        for (int j = i + 1; j <= q; ++j) {
            double phi = M_PI / (1 << (j - i));
            tmp = transform_cRw(a, n, i, j, phi);
            delete [] a;
            a = tmp;
        }
    }
    time = MPI_Wtime() - start_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("Time: %f\n", max_time);
    }
    if (argc == 4) {
        writeFourier(argv[3], n, size, rank, a);
    }
    delete [] a;
    MPI_Finalize();
}
