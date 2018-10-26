#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>

using namespace std;

void usage(int argc, char **argv) {
    printf("Usage: %s <first_number> <last_number> <output_file>\n", argv[0]);
    exit(1);
}

int find_primes(int left, int right, vector<int> & known_primes, MPI_File & output, double & time) {
    MPI_Status status;
    int len = (right - left + 7) >> 3;
    unsigned char *numbers = new unsigned char [len];
    memset(numbers, 0, len);
    for (int i = 0; i < known_primes.size(); ++i) {
        int first = left % known_primes[i] ? (left / known_primes[i] + 1) * known_primes[i] : left;
        for (int j = first; j < right; j += known_primes[i]) {
            numbers[(j - left) >> 3] |= (1 << ((j - left) & 7));
        }
    }
    time = MPI_Wtime();
    int count = 0;
    string prime_numbers;
    for (int i = left; i < right; ++i) {
        if (!(numbers[(i - left) >> 3] & (1 << ((i - left) & 7)))) {
            count++;
            prime_numbers += to_string(i) + ' ';
        }
    }
    delete [] numbers;
    MPI_File_write_shared(output, prime_numbers.data(), prime_numbers.length(), MPI_CHAR, &status);
    return count;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        usage(argc, argv);
    }
    MPI_Init(&argc, &argv);
    int left, right, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status status;
    MPI_File output;
    string prime_numbers;
    MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
            &output);
    sscanf(argv[1], "%d", &left);
    sscanf(argv[2], "%d", &right);
    if (left < 2) {
        left = 2;
    }
    if (right < 2) {
        MPI_Finalize();
        return 0;
    }
    double start_time, finish_time;
    int count = 0;
    vector<int> common_prime_numbers;
    start_time = MPI_Wtime();
    int n = sqrt(right);
    int len = (n + 8) >> 3;
    unsigned char *numbers = new unsigned char [len];
    memset(numbers, 0, len);
    for (int i = 2; i <= n; ++i) {
        if (!(numbers[i >> 3] & (1 << (i & 7)))) {
            if (!rank && i >= left) {
                count++;
                prime_numbers += to_string(i) + ' ';
            }
            common_prime_numbers.push_back(i);
            for (int j = i * 2; j <= n; j += i) {
                numbers[j >> 3] |= (1 << (j & 7));
            }
        }
    }
    if (!rank) {
        MPI_File_write_shared(output, prime_numbers.data(), prime_numbers.length(), MPI_CHAR,
                &status);
        prime_numbers.clear();
    }
    delete [] numbers;
    left = max(n + 1, left);
    int range = right - left + 1;
    count += find_primes(left + range * rank / size, left + range * (rank + 1) / size,
            common_prime_numbers, output, finish_time);
    int total;
    MPI_Reduce(&count, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double max_time, all_time, time = finish_time - start_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time, &all_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("Found %d prime numbers\n", total);
        printf("Total time: %fs\n", all_time);
        printf("Max time: %fs\n", max_time);
    }

    MPI_File_close(&output);
    MPI_Finalize();
    return 0;
}
