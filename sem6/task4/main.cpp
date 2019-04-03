#include "quantum_transform.hpp"
#include "quantum_io.hpp"

using namespace std;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    unsigned long long n;
    int num_threads;
    sscanf(argv[3], "%d", &num_threads);
    omp_set_num_threads(num_threads);
    complex<double> *a = read_vec(argv[1], n, size, rank);
    print_vec(a, n, size, rank);
    complex<double> *b = transform_cnot(a, n, 2, 1);
    if (!rank) {
        printf("---   ---   ---\n");
    }
    print_vec(b, n, size, rank);
    write_vec(argv[2], n, size, rank, b);
    delete [] a;
    delete [] b;
    MPI_Finalize();
}
