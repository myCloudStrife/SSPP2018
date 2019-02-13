#include <mpi.h>
#include <stdint.h>
#include <cstdlib>
#include <iostream>


struct Matrix {
    double *data;
    uint64_t rows, columns;
    Matrix() {
        rows = 0;
        columns = 0;
        data = NULL;
    }
    Matrix(uint64_t rows, uint64_t columns) {
        this->rows = rows;
        this->columns = columns;
        data = new double [rows * columns];
        memset(data, 0, rows * columns * sizeof(data[0]));
    }
    Matrix(const Matrix & matrix) {
        rows = matrix.rows;
        columns = matrix.columns;
        data = new double [rows * columns];
        memcpy(data, matrix.data, rows * columns * sizeof(data[0]));
    }
    Matrix & operator=(const Matrix & matrix) {
        delete [] data;
        rows = matrix.rows;
        columns = matrix.columns;
        data = new double [rows * columns];
        memcpy(data, matrix.data, rows * columns * sizeof(data[0]));
        return *this;
    }
    double * operator[](uint64_t index) {
        return data + index * columns;
    }
    ~Matrix() {
        delete [] data;
    }
};

void usage(int argc, char **argv) {
    printf("Usage: %s <matrix1> <matrix2> <output>\n", argv[0]);
    exit(1);
}

void readMatrix(Matrix & matrix, MPI_Comm & square, int *cube_dims, int *cube_coords,
        char *filename, uint64_t & m, uint64_t & n) {
    MPI_File file;
    MPI_File_open(square, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    char type;
    MPI_File_read_all(file, &type, 1, MPI_CHAR, NULL);
    MPI_File_read_all(file, &m, 1, MPI_UNSIGNED_LONG_LONG, NULL);
    MPI_File_read_all(file, &n, 1, MPI_UNSIGNED_LONG_LONG, NULL);
    int block_m, block_n;
    block_m = (cube_coords[0] + 1) * m / cube_dims[0] - cube_coords[0] * m / cube_dims[0];
    MPI_Offset offset = sizeof(type) + sizeof(m) + sizeof(n) +
            (cube_coords[0] * m / cube_dims[0]) * n * sizeof(double);
    MPI_Datatype filetype;
    int array_size = n;
    int start_array = cube_coords[1] * n / cube_dims[1];
    block_n = (cube_coords[1] + 1) * n / cube_dims[1] - start_array;
    MPI_Type_create_subarray(1, &array_size, &block_n, &start_array, MPI_ORDER_C,
            MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    matrix = Matrix(block_m, block_n);
    MPI_File_read_all(file, matrix[0], block_m * block_n, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Type_free(&filetype);
    MPI_File_close(&file);
}

void writeMatrix(Matrix & matrix, MPI_Comm & square, int *cube_dims, int *cube_coords,
        char *filename, uint64_t m, uint64_t n) {
    MPI_File file;
    MPI_File_open(square, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
    char type = 'd';
    if (!cube_coords[0] && !cube_coords[1]) {
        MPI_File_write(file, &type, 1, MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write(file, &m, 1, MPI_UNSIGNED_LONG_LONG, MPI_STATUS_IGNORE);
        MPI_File_write(file, &n, 1, MPI_UNSIGNED_LONG_LONG, MPI_STATUS_IGNORE);
    }
    MPI_Offset offset = sizeof(type) + sizeof(m) + sizeof(n) +
            (cube_coords[0] * m / cube_dims[0]) * n * sizeof(double);
    MPI_Datatype filetype;
    int array_size = n;
    int start_array = cube_coords[1] * n / cube_dims[1];
    int block_m = matrix.rows, block_n = matrix.columns;
    MPI_Type_create_subarray(1, &array_size, &block_n, &start_array, MPI_ORDER_C,
            MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(file, offset, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);
    MPI_File_write_all(file, matrix[0], block_m * block_n, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Type_free(&filetype);
    MPI_File_close(&file);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm cube, square, line_x, line_y, line_z;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int cube_rank;
    int cube_dims[] = { 0, 0, 0 };
    int cube_periods[] = { 0, 0, 0 };
    int cube_coords[] = { 0, 0, 0 };
    MPI_Dims_create(size, 3, cube_dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, cube_dims, cube_periods, 0, &cube);
    MPI_Comm_rank(cube, &cube_rank);
    if (!cube_rank) {
        if (argc != 4) {
            usage(argc, argv);
        }
        if (cube_dims[0] != cube_dims[1] || cube_dims[1] != cube_dims[2]) {
            fprintf(stderr, "Wrong topology, dimensions: %d, %d, %d\n", cube_dims[0],
                    cube_dims[1], cube_dims[2]);
            exit(1);
        }
    }
    double load_data_time, calculation_time, print_matrix_time, max_time;
    MPI_Cart_coords(cube, cube_rank, 3, cube_coords);
    int square_dims[] = { 1, 1, 0 };
    MPI_Cart_sub(cube, square_dims, &square);
    int line_dims[] = { 1, 0, 0 };
    MPI_Cart_sub(cube, line_dims, &line_x);
    line_dims[0] = 0;
    line_dims[1] = 1;
    MPI_Cart_sub(cube, line_dims, &line_y);
    line_dims[1] = 0;
    line_dims[2] = 1;
    MPI_Cart_sub(cube, line_dims, &line_z);
    Matrix A, B, C;
    uint64_t C_m = 0, C_n = 0;

    MPI_Barrier(cube);
    load_data_time = MPI_Wtime();
    //load A matrix
    if (!cube_coords[2]) { //bottom square
        uint64_t m, n;
        readMatrix(A, square, cube_dims, cube_coords, argv[1], m, n);
        C_m = m;
        if (cube_coords[1]) { //firstly send block up, then receive broadcast
            int dest_coord = cube_coords[1];
            int dest_rank;
            MPI_Cart_rank(line_z, &dest_coord, &dest_rank);
            int block_sizes[] = { (int) A.rows, (int) A.columns };
            MPI_Send(block_sizes, 2, MPI_INT, dest_rank, 2, line_z);
            MPI_Send(A[0], A.rows * A.columns, MPI_DOUBLE, dest_rank, 1, line_z);
            A.columns = n / cube_dims[1];
            dest_coord = 0;
            MPI_Cart_rank(line_y, &dest_coord, &dest_rank);
            MPI_Bcast(A[0], A.rows * A.columns, MPI_DOUBLE, dest_rank, line_y);
        } else { // only send broadcast
            int root_rank;
            MPI_Comm_rank(line_y, &root_rank);
            MPI_Bcast(A[0], A.rows * A.columns, MPI_DOUBLE, root_rank, line_y);
        }
    } else if (cube_coords[1] == cube_coords[2]) {
        int block_sizes[2];
        int source_coord = 0;
        int source_rank;
        MPI_Cart_rank(line_z, &source_coord, &source_rank);
        MPI_Recv(block_sizes, 2, MPI_INT, source_rank, 2, line_z, MPI_STATUS_IGNORE);
        A = Matrix(block_sizes[0], block_sizes[1]);
        MPI_Recv(A[0], A.rows * A.columns, MPI_DOUBLE, source_rank, 1, line_z, MPI_STATUS_IGNORE);
        MPI_Comm_rank(line_y, &source_rank);
        MPI_Bcast(block_sizes, 2, MPI_INT, source_rank, line_y);
        MPI_Bcast(A[0], A.rows * A.columns, MPI_DOUBLE, source_rank, line_y);
    } else {
        int block_sizes[2];
        int root_coord = cube_coords[2];
        int root_rank;
        MPI_Cart_rank(line_y, &root_coord, &root_rank);
        MPI_Bcast(block_sizes, 2, MPI_INT, root_rank, line_y);
        A = Matrix(block_sizes[0], block_sizes[1]);
        MPI_Bcast(A[0], A.rows * A.columns, MPI_DOUBLE, root_rank, line_y);
    }

    //load B matrix
    if (!cube_coords[2]) { //bottom square
        uint64_t m, n;
        readMatrix(B, square, cube_dims, cube_coords, argv[2], m, n);
        C_n = n;
        if (cube_coords[0]) { //firstly send block up, then receive broadcast
            int dest_coord = cube_coords[0];
            int dest_rank;
            MPI_Cart_rank(line_z, &dest_coord, &dest_rank);
            int block_sizes[] = { (int) B.rows, (int) B.columns };
            MPI_Send(block_sizes, 2, MPI_INT, dest_rank, 2, line_z);
            MPI_Send(B[0], B.rows * B.columns, MPI_DOUBLE, dest_rank, 1, line_z);
            B.rows = m / cube_dims[0];
            dest_coord = 0;
            MPI_Cart_rank(line_x, &dest_coord, &dest_rank);
            MPI_Bcast(B[0], B.rows * B.columns, MPI_DOUBLE, dest_rank, line_x);
        } else { // only send broadcast
            int root_rank;
            MPI_Comm_rank(line_x, &root_rank);
            MPI_Bcast(B[0], B.rows * B.columns, MPI_DOUBLE, root_rank, line_x);
        }
    } else if (cube_coords[0] == cube_coords[2]) {
        int block_sizes[2];
        int source_coord = 0;
        int source_rank;
        MPI_Cart_rank(line_z, &source_coord, &source_rank);
        MPI_Recv(block_sizes, 2, MPI_INT, source_rank, 2, line_z, MPI_STATUS_IGNORE);
        B = Matrix(block_sizes[0], block_sizes[1]);
        MPI_Recv(B[0], B.rows * B.columns, MPI_DOUBLE, source_rank, 1, line_z, MPI_STATUS_IGNORE);
        MPI_Comm_rank(line_x, &source_rank);
        MPI_Bcast(block_sizes, 2, MPI_INT, source_rank, line_x);
        MPI_Bcast(B[0], B.rows * B.columns, MPI_DOUBLE, source_rank, line_x);
    } else {
        int block_sizes[2];
        int root_coord = cube_coords[2];
        int root_rank;
        MPI_Cart_rank(line_x, &root_coord, &root_rank);
        MPI_Bcast(block_sizes, 2, MPI_INT, root_rank, line_x);
        B = Matrix(block_sizes[0], block_sizes[1]);
        MPI_Bcast(B[0], B.rows * B.columns, MPI_DOUBLE, root_rank, line_x);
    }
    load_data_time = MPI_Wtime() - load_data_time;

    if (A.columns != B.rows) {
        fprintf(stderr, "Error on (%d,%d,%d) node, wrong sizes\n", cube_coords[0],
                cube_coords[1], cube_coords[2]);
        exit(1);
    }

    MPI_Barrier(cube);
    calculation_time = MPI_Wtime();
    C = Matrix(A.rows, B.columns);
    for (uint64_t i = 0; i < A.rows; ++i) {
        for (uint64_t k = 0; k < A.columns; ++k) {
            for (uint64_t j = 0; j < B.columns; ++j) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    int root_coord = 0;
    int root_rank;
    MPI_Cart_rank(line_z, &root_coord, &root_rank);
    if (cube_coords[2]) {
        MPI_Reduce(C[0], NULL, C.rows * C.columns, MPI_DOUBLE, MPI_SUM, root_rank, line_z);
    } else {
        MPI_Reduce(MPI_IN_PLACE, C[0], C.rows * C.columns, MPI_DOUBLE, MPI_SUM, root_rank, line_z);
    }
    calculation_time = MPI_Wtime() - calculation_time;

    if (!cube_coords[2]) {
        MPI_Barrier(square);
        print_matrix_time = MPI_Wtime();
        writeMatrix(C, square, cube_dims, cube_coords, argv[3], C_m, C_n);
        print_matrix_time = MPI_Wtime() - print_matrix_time;
        int coords[] = { 0, 0 };
        int root;
        MPI_Cart_rank(square, coords, &root);
        MPI_Reduce(&print_matrix_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, root, square);
        if (!cube_coords[0] && !cube_coords[1]) {
            printf("Print matrix max time: %f\n", max_time);
        }
    }
    MPI_Reduce(&calculation_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, cube);
    if (!cube_rank) {
        printf("Calculation max time: %f\n", max_time);
    }
    MPI_Reduce(&load_data_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, cube);
    if (!cube_rank) {
        printf("Load data max time: %f\n", max_time);
    }

    MPI_Finalize();
}
