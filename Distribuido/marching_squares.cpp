#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>

// Función objetivo - círculo de radio 200
double f(double x, double y) {
    return x * x + y * y - 200.0 * 200.0;
}

bool out(double x, double y) {
    return f(x, y) > 0;
}

// Calcula el punto de interpolación en un borde.
double interpolation(double val0, double val1) {
    // Evita la división por cero si los valores son iguales.
    if (std::abs(val1 - val0) < 1e-9) return 0.5;
    return val0 / (val0 - val1);
}

struct Point2D {
    double x, y;
};

struct LineSegment {
    Point2D p1, p2;
};

class ParallelMarchingSquaresEPS {
private:
    int rank, size;
    MPI_Comm cart_comm;
    int dims[2], coords[2];
    int neighbors[4]; // 0:left, 1:right, 2:up, 3:down

    int global_M, global_N;
    double x_min, x_max, y_min, y_max;
    double dx, dy;

    int local_M_no_ghost, local_N_no_ghost;
    int local_M_with_ghost, local_N_with_ghost;
    int global_start_i, global_start_j;

    std::vector<double> local_values; 
    std::vector<LineSegment> local_segments;

    MPI_Datatype col_type;
    bool mpi_types_created;

public:
    ParallelMarchingSquaresEPS(double xmin, double xmax, double ymin, double ymax,
                               int M, int N)
        : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax),
          global_M(M), global_N(N), mpi_types_created(false) {

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Crear topología cartesiana 2D
        dims[0] = dims[1] = 0;
        MPI_Dims_create(size, 2, dims);

        if (rank == 0) {
            std::cout << "Creating " << dims[0] << "x" << dims[1]
                      << " process grid for " << size << " processes\n";
        }

        int periods[2] = {0, 0}; 
        int reorder = 1;
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
        MPI_Cart_coords(cart_comm, rank, 2, coords);

        MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[1]); // left, right
        MPI_Cart_shift(cart_comm, 1, 1, &neighbors[2], &neighbors[3]); // down, up (MPI Cartesiano es como una matriz)

        local_M_no_ghost = global_M / dims[0] + (coords[0] < global_M % dims[0] ? 1 : 0);
        local_N_no_ghost = global_N / dims[1] + (coords[1] < global_N % dims[1] ? 1 : 0);

        global_start_i = (coords[0] * (global_M / dims[0])) + std::min(coords[0], global_M % dims[0]);
        global_start_j = (coords[1] * (global_N / dims[1])) + std::min(coords[1], global_N % dims[1]);

        local_M_with_ghost = local_M_no_ghost + 2;
        local_N_with_ghost = local_N_no_ghost + 2;

        dx = (x_max - x_min) / global_M;
        dy = (y_max - y_min) / global_N;

        local_values.resize((local_M_with_ghost + 1) * (local_N_with_ghost + 1));

        
        MPI_Type_vector(local_N_no_ghost + 1, 1, local_M_with_ghost + 1, MPI_DOUBLE, &col_type);
        MPI_Type_commit(&col_type);
        mpi_types_created = true;
    }

    ~ParallelMarchingSquaresEPS() {
        if (mpi_types_created) {
            MPI_Type_free(&col_type);
        }
        if (cart_comm != MPI_COMM_NULL) {
            MPI_Comm_free(&cart_comm);
        }
    }

    void initialize_values() {
        for (int j = 0; j <= local_N_no_ghost; ++j) {
            for (int i = 0; i <= local_M_no_ghost; ++i) {
                double x = x_min + (global_start_i + i) * dx;
                double y = y_min + (global_start_j + j) * dy;
                local_values[(j + 1) * (local_M_with_ghost + 1) + (i + 1)] = f(x, y);
            }
        }
    }

    void exchange_ghost_values() {
        MPI_Request requests[8];
        int req_count = 0;

        if (neighbors[0] != MPI_PROC_NULL) {
            MPI_Isend(&local_values[1 * (local_M_with_ghost + 1) + 1], 1, col_type, neighbors[0], 0, cart_comm, &requests[req_count++]);
            MPI_Irecv(&local_values[1 * (local_M_with_ghost + 1) + 0], 1, col_type, neighbors[0], 1, cart_comm, &requests[req_count++]);
        }
        if (neighbors[1] != MPI_PROC_NULL) {
            MPI_Isend(&local_values[1 * (local_M_with_ghost + 1) + local_M_no_ghost], 1, col_type, neighbors[1], 1, cart_comm, &requests[req_count++]);
            MPI_Irecv(&local_values[1 * (local_M_with_ghost + 1) + local_M_with_ghost], 1, col_type, neighbors[1], 0, cart_comm, &requests[req_count++]);
        }

        if (neighbors[2] != MPI_PROC_NULL) {
            MPI_Isend(&local_values[1 * (local_M_with_ghost + 1) + 1], local_M_no_ghost + 1, MPI_DOUBLE, neighbors[2], 2, cart_comm, &requests[req_count++]);
            MPI_Irecv(&local_values[0 * (local_M_with_ghost + 1) + 1], local_M_no_ghost + 1, MPI_DOUBLE, neighbors[2], 3, cart_comm, &requests[req_count++]);
        }
        if (neighbors[3] != MPI_PROC_NULL) {
            MPI_Isend(&local_values[local_N_no_ghost * (local_M_with_ghost + 1) + 1], local_M_no_ghost + 1, MPI_DOUBLE, neighbors[3], 3, cart_comm, &requests[req_count++]);
            MPI_Irecv(&local_values[(local_N_no_ghost + 1) * (local_M_with_ghost + 1) + 1], local_M_no_ghost + 1, MPI_DOUBLE, neighbors[3], 2, cart_comm, &requests[req_count++]);
        }

        MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
    }
    
    void process_cell(int i_local, int j_local) {
        double x0 = x_min + (global_start_i + i_local) * dx;
        double y0 = y_min + (global_start_j + j_local) * dy;
        int i = i_local + 1;
        int j = j_local + 1;
        double v[4];
        v[0] = local_values[j * (local_M_with_ghost + 1) + i];
        v[1] = local_values[j * (local_M_with_ghost + 1) + (i + 1)];
        v[2] = local_values[(j + 1) * (local_M_with_ghost + 1) + (i + 1)];
        v[3] = local_values[(j + 1) * (local_M_with_ghost + 1) + i];
        int code = (v[0] > 0) | ((v[1] > 0) << 1) | ((v[2] > 0) << 2) | ((v[3] > 0) << 3);

        if (code == 0 || code == 15) return;

        Point2D m[4];
        m[0] = {x0 + dx * interpolation(v[0], v[1]), y0};
        m[1] = {x0 + dx, y0 + dy * interpolation(v[1], v[2])};
        m[2] = {x0 + dx * interpolation(v[3], v[2]), y0 + dy};
        m[3] = {x0, y0 + dy * interpolation(v[0], v[3])};

        switch (code) {
            case 1: case 14: local_segments.push_back({m[0], m[3]}); break;
            case 2: case 13: local_segments.push_back({m[0], m[1]}); break;
            case 3: case 12: local_segments.push_back({m[1], m[3]}); break;
            case 4: case 11: local_segments.push_back({m[1], m[2]}); break;
            case 5: case 10:
                if ((f(x0 + dx/2, y0 + dy/2) > 0) != (v[0] > 0)) {
                    local_segments.push_back({m[0], m[3]});
                    local_segments.push_back({m[1], m[2]});
                } else {
                    local_segments.push_back({m[0], m[1]});
                    local_segments.push_back({m[2], m[3]});
                }
                break;
            case 6: case 9:  local_segments.push_back({m[0], m[2]}); break;
            case 7: case 8:  local_segments.push_back({m[2], m[3]}); break;
        }
    }

    void process_all_local_cells() {
        for (int j = 0; j < local_N_no_ghost; ++j) {
            for (int i = 0; i < local_M_no_ghost; ++i) {
                process_cell(i, j);
            }
        }
    }

    void write_eps_gather_and_write(const std::string& filename) {
        
        MPI_Datatype segment_type;
        MPI_Type_contiguous(4, MPI_DOUBLE, &segment_type);
        MPI_Type_commit(&segment_type);

        int local_segment_count = local_segments.size();
        std::vector<int> recv_counts;
        std::vector<int> displacements;
        std::vector<LineSegment> all_segments;
        int total_segments = 0;

        if (rank == 0) {
            recv_counts.resize(size);
        }

        MPI_Gather(&local_segment_count, 1, MPI_INT, 
                   recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            displacements.resize(size);
            displacements[0] = 0;
            total_segments += recv_counts[0];
            for (int i = 1; i < size; ++i) {
                total_segments += recv_counts[i];
                displacements[i] = displacements[i-1] + recv_counts[i-1];
            }
            all_segments.resize(total_segments);
        }

       
        MPI_Gatherv(local_segments.data(), local_segment_count, segment_type,
                    all_segments.data(), recv_counts.data(), displacements.data(),
                    segment_type, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::ofstream eps_file(filename);
            eps_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
            eps_file << "%%BoundingBox: " << x_min << " " << y_min << " " << x_max << " " << y_max << "\n";
            eps_file << "%%Title: Marching Squares Contour (Gather-and-Write)\n";
            eps_file << "0.5 setlinewidth\n";
            eps_file << std::fixed << std::setprecision(6);

            for(const auto& seg : all_segments) {
                eps_file << seg.p1.x << " " << seg.p1.y << " moveto "
                         << seg.p2.x << " " << seg.p2.y << " lineto stroke\n";
            }

            eps_file << "showpage\n%%EOF\n";
            eps_file.close();

            std::cout << "Written " << total_segments << " segments to " << filename << std::endl;
        }
        
        MPI_Type_free(&segment_type);
    }

    void run() {
        double start_time, init_time, compute_time, end_time;

        MPI_Barrier(cart_comm);
        start_time = MPI_Wtime();

        initialize_values();
        
        MPI_Barrier(cart_comm);
        init_time = MPI_Wtime();

        local_segments.clear();
        local_segments.reserve(local_M_no_ghost * local_N_no_ghost / 4);

        exchange_ghost_values();
        process_all_local_cells();

        MPI_Barrier(cart_comm);
        compute_time = MPI_Wtime();

        write_eps_gather_and_write("output_parallel_optimized.eps");
        
        MPI_Barrier(cart_comm);
        end_time = MPI_Wtime();

        double times[3] = {init_time - start_time,
                           compute_time - init_time,
                           end_time - compute_time};
        double max_times[3];
        MPI_Reduce(times, max_times, 3, MPI_DOUBLE, MPI_MAX, 0, cart_comm);

        if (rank == 0) {
            std::cout << "Total time: " << (end_time - start_time)*1000 << "ms\n";
            std::cout << "  Initialization: " << (max_times[0])*1000 << "ms\n";
            std::cout << "  Computation: " << (max_times[1])*1000 << "ms\n";
            std::cout << "  I/O: " << (max_times[2])*1000 << "ms\n";
        }
    }
};

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);


    {
        double x_min = -300.0, x_max = 300.0;
        double y_min = -300.0, y_max = 300.0;
        int M = 1000, N = 1000;

        if (argc > 1) M = std::atoi(argv[1]);
        if (argc > 2) N = std::atoi(argv[2]);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank == 0) {
            std::cout << "Running Parallel Marching Squares with PURE MPI output\n";
        }

        ParallelMarchingSquaresEPS pms(x_min, x_max, y_min, y_max, M, N);
        pms.run();
    } 

    MPI_Finalize();
    return 0;
}
