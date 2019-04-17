#include "../quantum_transform.hpp"
#include "../quantum_io.hpp"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    uint64_t n;
    int k, l;
    int num_threads;
    sscanf(argv[2], "%d", &k);
    sscanf(argv[3], "%d", &l);
    sscanf(argv[5], "%d", &num_threads);
    omp_set_num_threads(num_threads);
    std::complex<double> *a = read_vec(argv[1], &n, size, rank);
    std::complex<double> *b = transform_cnot(a, n, k, l);
    write_vec(argv[4], n, size, rank, b);
    delete [] a;
    delete [] b;
    MPI_Finalize();
}
