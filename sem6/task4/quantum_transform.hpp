#ifndef SSPP_SEM6_TASK4_QUANTUM_TRANSFORM_HPP_
#define SSPP_SEM6_TASK4_QUANTUM_TRANSFORM_HPP_

#include <mpi.h>
#include <omp.h>
#include <complex>
#include <climits>
#include <utility>

constexpr int SENDRECVTAG = 7321;

std::complex<double> * quant_transform(std::complex<double> *a, uint64_t n,
        int k, std::complex<double> *u) {
    static int size = -1, rank = -1;
    if (size == -1 || rank == -1) {
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }

    std::complex<double> *b = new std::complex<double>[n];
    uint64_t vec_fullsize = n * size;
    int target = rank ^ (size >> k);
    if (target != rank) {
        uint64_t need_count = n * 2;
        int sendrecvcount = INT_MAX - 1;
        uint64_t i = 0;
        while (need_count > INT_MAX) {
            MPI_Sendrecv(a + i * sendrecvcount / 2, sendrecvcount, MPI_DOUBLE,
                    target, SENDRECVTAG, b + i * sendrecvcount / 2,
                    sendrecvcount, MPI_DOUBLE, target, SENDRECVTAG,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            i++;
            need_count -= sendrecvcount;
        }
        MPI_Sendrecv(a + i * sendrecvcount / 2, need_count, MPI_DOUBLE, target,
                SENDRECVTAG, b + i * sendrecvcount / 2, need_count, MPI_DOUBLE,
                target, SENDRECVTAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (target < rank) {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; ++i) {
                b[i] = b[i] * u[1] + a[i] * u[3];
            }
        } else {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; ++i) {
                b[i] = a[i] * u[2] + b[i] * u[0];
            }
        }
    } else {
        int q = -1;
        while (vec_fullsize) {
            q++;
            vec_fullsize >>= 1;
        }
        k = q - k;
        uint64_t bit = 1ull << k;
#pragma omp parallel for
        for (uint64_t i = 0; i < n; ++i) {
            int u_row = ((i & bit) >> k) * 2;
            b[i] = a[i & ~bit] * u[u_row] + a[i | bit] * u[u_row + 1];
        }
    }
    return b;
}

std::complex<double> * quant_transform(std::complex<double> *a, uint64_t n,
        int k1, int k2, std::complex<double> *u) {
    static int size = -1, rank = -1;
    if (size == -1 || rank == -1) {
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }

    bool needswap = false;
    if (k2 < k1) {
        needswap = true;
        std::swap(k1, k2);
    }
    std::complex<double> *a00, *a01, *a10, *a11;
    int owner00, owner01, owner10, owner11;
    owner00 = rank & ~(size >> k1) & ~(size >> k2);
    owner01 = (rank & ~(size >> k1)) | (size >> k2);
    owner10 = (rank | (size >> k1)) & ~(size >> k2);
    owner11 = rank | (size >> k1) | (size >> k2);

    if (rank != owner00) {
        a00 = new std::complex<double>[n];
        MPI_Sendrecv(a, n * 2, MPI_DOUBLE, owner00, SENDRECVTAG, a00, n * 2,
                MPI_DOUBLE, owner00, SENDRECVTAG, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    } else {
        a00 = a;
    }

    if (rank != owner10) {
        a10 = new std::complex<double>[n];
        MPI_Sendrecv(a, n * 2, MPI_DOUBLE, owner10, SENDRECVTAG, a10, n * 2,
                MPI_DOUBLE, owner10, SENDRECVTAG, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    } else {
        a10 = a;
    }

    if (rank != owner01) {
        if (owner01 != owner00) {
            a01 = new std::complex<double>[n];
            MPI_Sendrecv(a, n * 2, MPI_DOUBLE, owner01, SENDRECVTAG, a01, n * 2,
                    MPI_DOUBLE, owner01, SENDRECVTAG, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
        } else {
            a01 = a00;
        }
    } else {
        a01 = a;
    }

    if (rank != owner11) {
        if (owner11 != owner10) {
            a11 = new std::complex<double>[n];
            MPI_Sendrecv(a, n * 2, MPI_DOUBLE, owner11, SENDRECVTAG, a11, n * 2,
                    MPI_DOUBLE, owner11, SENDRECVTAG, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
        } else {
            a11 = a10;
        }
    } else {
        a11 = a;
    }
    std::complex<double> *b = new std::complex<double>[n];

    if (a00 != a01) {
        int me = 0;
        if (rank == owner00) {
            me = 0;
        } else if (rank == owner01) {
            me = needswap ? 0b1000 : 0b100;
        } else if (rank == owner10) {
            me = needswap ? 0b100 : 0b1000;
        } else {
            me = 0b1100;
        }
        if (needswap) {
            swap(a01, a10);
        }
#pragma omp parallel for
        for (uint64_t i = 0; i < n; ++i) {
            b[i] = u[me + 0] * a00[i] + u[me + 1] * a01[i] + u[me + 2] * a10[i]
                    + u[me + 3] * a11[i];
        }
    } else if (a00 != a10) {
        int me = 0;
        if (rank == owner10) {
            me = needswap ? 0b100 : 0b1000;
        }

        int q = -1;
        uint64_t vec_fullsize = n * size;
        while (vec_fullsize) {
            q++;
            vec_fullsize >>= 1;
        }
        int k = q - k2;
        uint64_t bit = 1ull << k;

        if (!needswap) {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; ++i) {
                int u_row = me + ((i & bit) >> k) * 0b100;
                b[i] = u[u_row + 0] * a00[i & ~bit]
                        + u[u_row + 1] * a00[i | bit]
                        + u[u_row + 2] * a10[i & ~bit]
                        + u[u_row + 3] * a10[i | bit];
            }
        } else {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; ++i) {
                int u_row = me + ((i & bit) >> k) * 0b1000;
                b[i] = u[u_row + 0] * a00[i & ~bit]
                        + u[u_row + 1] * a10[i & ~bit]
                        + u[u_row + 2] * a00[i | bit]
                        + u[u_row + 3] * a10[i | bit];
            }
        }
    } else {
        int q = -1;
        uint64_t vec_fullsize = n * size;
        while (vec_fullsize) {
            q++;
            vec_fullsize >>= 1;
        }
        k1 = q - k1;
        k2 = q - k2;
        uint64_t bit1 = 1ull << k1;
        uint64_t bit2 = 1ull << k2;

        if (!needswap) {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; ++i) {
                int u_row = ((i & bit1) >> k1) * 0b1000
                        + ((i & bit2) >> k2) * 0b100;
                b[i] = u[u_row + 0] * a[i & ~bit1 & ~bit2]
                        + u[u_row + 1] * a[(i & ~bit1) | bit2]
                        + u[u_row + 2] * a[(i | bit1) & ~bit2]
                        + u[u_row + 3] * a[i | bit1 | bit2];
            }
        } else {
#pragma omp parallel for
            for (uint64_t i = 0; i < n; ++i) {
                int u_row = ((i & bit1) >> k1) * 0b100
                        + ((i & bit2) >> k2) * 0b1000;
                b[i] = u[u_row + 0] * a[i & ~bit1 & ~bit2]
                        + u[u_row + 1] * a[(i | bit1) & ~bit2]
                        + u[u_row + 2] * a[(i & ~bit1) | bit2]
                        + u[u_row + 3] * a[i | bit1 | bit2];
            }
        }
    }

    if (a != a00) {
        delete [] a00;
    }
    if (a != a01 && a00 != a01) {
        delete [] a01;
    }
    if (a != a10 && a00 != a10 && a01 != a10) {
        delete [] a10;
    }
    if (a != a11 && a00 != a11 && a01 != a11 && a10 != a11) {
        delete [] a11;
    }

    return b;
}

std::complex<double> * transform_adamar(std::complex<double> *a, uint64_t n,
        int k) {
    static std::complex<double> u[4] = {
            M_SQRT1_2, M_SQRT1_2,
            M_SQRT1_2, -M_SQRT1_2
    };
    return quant_transform(a, n, k, u);
}

std::complex<double> * transform_n_adamar(std::complex<double> *a, uint64_t n) {
    static int size = -1;
    if (size == -1) {
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }
    uint64_t vec_fullsize = n * size;
    int q = -1;
    while (vec_fullsize) {
        q++;
        vec_fullsize >>= 1;
    }
    std::complex<double> *b = a, *tmp;
    for (int k = 1; k <= q; ++k) {
        tmp = transform_adamar(b, n, k);
        if (k > 1) {
            delete [] b;
        }
        b = tmp;
    }
    return b;
}

std::complex<double> * transform_not(std::complex<double> *a, uint64_t n,
        int k) {
    static std::complex<double> u[4] = {
            0, 1,
            1, 0
    };
    return quant_transform(a, n, k, u);
}

std::complex<double> * transform_Rw(std::complex<double> *a, uint64_t n, int k,
        double phi) {
    std::complex<double> u[4] = {
            1, 0,
            0, std::complex<double>(cos(phi), sin(phi))
    };
    return quant_transform(a, n, k, u);
}

std::complex<double> * transform_cnot(std::complex<double> *a, uint64_t n,
        int k1, int k2) {
    static std::complex<double> u[16] = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0
    };
    return quant_transform(a, n, k1, k2, u);
}

std::complex<double> * transform_cRw(std::complex<double> *a, uint64_t n,
        int k1, int k2, double phi) {
    static std::complex<double> u[16] = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, std::complex<double>(cos(phi), sin(phi))
    };
    return quant_transform(a, n, k1, k2, u);
}

#endif /* SSPP_SEM6_TASK4_QUANTUM_TRANSFORM_HPP_ */
