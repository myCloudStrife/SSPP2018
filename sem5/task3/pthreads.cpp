#include <pthread.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <ctime>

using namespace std;

void usage(int argc, char **argv) {
    printf("Usage: %s <first_number> <last_number> <output_file> <threads_num>\n", argv[0]);
    exit(1);
}

struct Data {
    int left, right;
    vector<int> & known_primes, & new_primes;
    double time;
};


void * find_primes(void *args) {
    timespec time1, time2;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time1);
    Data *d = (Data *) args;
    vector<bool> numbers(d->right - d->left, true);
    for (int i = 0; i < d->known_primes.size(); ++i) {
        int prime = d->known_primes[i];
        int first = d->left % prime ? (d->left / prime + 1) * d->known_primes[i] : d->left;
        for (int j = first; j < d->right; j += prime) {
            numbers[j - d->left] = false;
        }
    }
    for (int i = d->left; i < d->right; ++i) {
        if (numbers[i - d->left]) {
            d->new_primes.push_back(i);
        }
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time2);
    d->time = time2.tv_sec - time1.tv_sec + 1e-9 * (time2.tv_nsec - time1.tv_nsec);
    return nullptr;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        usage(argc, argv);
    }
    int left, right, size;
    sscanf(argv[1], "%d", &left);
    sscanf(argv[2], "%d", &right);
    sscanf(argv[4], "%d", &size);
    if (left < 2) {
        left = 2;
    }
    if (right < 2) {
        return 0;
    }
    int total = 0;
    vector<int> common_prime_numbers;
    int n = (int) sqrt(right);
    vector<bool> numbers(n + 1, true);
    for (int i = 2; i <= n; ++i) {
        if (numbers[i]) {
            if (i >= left) {
                total++;
            }
            common_prime_numbers.push_back(i);
            for (int j = i * 2; j <= n; j += i) {
                numbers[j] = false;
            }
        }
    }
    int range = right - max(n + 1, left) + 1;
    vector<int> extern_primes[size];
    double times[size];
    vector<Data> data;
    pthread_t tids[size];
    for (int i = 0; i < size; ++i) {
        data.push_back({
                (int) (max(n + 1, left) + (long long) range * i / size),
                (int) (max(n + 1, left) + (long long) range * (i + 1) / size),
                common_prime_numbers,
                extern_primes[i],
                times[i]
        });
    }
    for (int i = 0; i < size; ++i) {
        pthread_create(&tids[i], nullptr, find_primes, &data[i]);
    }
    for (int i = 0; i < size; ++i) {
        pthread_join(tids[i], nullptr);
    }
    FILE *f = fopen(argv[3], "w");
    for (int i = common_prime_numbers.size() - 1; i >= 0; --i) {
        if (common_prime_numbers[i] < left) {
            break;
        }
        fprintf(f, "%d\n", common_prime_numbers[i]);
    }
    for (int i = 0; i < size; ++i) {
        total += extern_primes[i].size();
        for (auto num : extern_primes[i]) {
            fprintf(f, "%d\n", num);
        }
    }
    double all_time = 0, max_time = 0;
    for (int i = 0; i < size; ++i) {
        all_time += data[i].time;
        if (max_time < data[i].time) {
            max_time = data[i].time;
        }
    }
    printf("Found %d prime numbers\n", total);
    printf("Total time: %fs\n", all_time);
    printf("Max time: %fs\n", max_time);
    fclose(f);
    return 0;
}
