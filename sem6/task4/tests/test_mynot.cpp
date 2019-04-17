#include "../quantum_transform.hpp"
#include "../quantum_io.hpp"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t n;
    int k;
    int num_threads;
    sscanf(argv[2], "%d", &k);
    sscanf(argv[4], "%d", &num_threads);
    omp_set_num_threads(num_threads);
    std::complex<double> *a = read_vec(argv[1], &n, size, rank);
    std::complex<double> *b = transform_not(a, n, k);
    write_vec(argv[3], n, size, rank, b);
    delete [] a;
    delete [] b;
    MPI_Finalize();
}
