#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>

double f(double x, double y) {
    return x * x + y * y - 200.0 * 200.0;
}


double interpolation(double v0, double v1) {
    if (std::abs(v1 - v0) < 1e-9) return 0.5;
    return v0 / (v0 - v1);
}

struct Point2D {
    double x, y;
};

struct LineSegment {
    Point2D p1, p2;
};

class MarchingSquaresMPI {
private:
    int rank, size;
    MPI_Comm cart_comm;
    int dims[2], coords[2];

    int global_M, global_N;
    double x_min, x_max, y_min, y_max;
    double dx, dy;

    int local_M, local_N;
    int global_start_i, global_start_j;

    std::vector<LineSegment> all_local_segments;

public:
    MarchingSquaresMPI(double xmin, double xmax, double ymin, double ymax, int M, int N)
        : x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax),
          global_M(M), global_N(N) {

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        dims[0] = dims[1] = 0;
        MPI_Dims_create(size, 2, dims);

        int periods[2] = {0, 0}; 
        int reorder = 1;
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
        MPI_Cart_coords(cart_comm, rank, 2, coords);

        local_M = global_M / dims[0] + (coords[0] < global_M % dims[0] ? 1 : 0);
        local_N = global_N / dims[1] + (coords[1] < global_N % dims[1] ? 1 : 0);

        global_start_i = (coords[0] * (global_M / dims[0])) + std::min(coords[0], global_M % dims[0]);
        global_start_j = (coords[1] * (global_N / dims[1])) + std::min(coords[1], global_N % dims[1]);

        dx = (x_max - x_min) / global_M;
        dy = (y_max - y_min) / global_N;

        all_local_segments.reserve(local_M * local_N / 10);
    }

    ~MarchingSquaresMPI() {
        if (cart_comm != MPI_COMM_NULL) {
            MPI_Comm_free(&cart_comm);
        }
    }

    void process_cell_independent(int i_local, int j_local, std::vector<LineSegment>& segments) {
        double x0 = x_min + (global_start_i + i_local) * dx;
        double y0 = y_min + (global_start_j + j_local) * dy;
        
        double v[4];
        v[0] = f(x0,      y0);      
        v[1] = f(x0 + dx, y0);      
        v[2] = f(x0 + dx, y0 + dy); 
        v[3] = f(x0,      y0 + dy); 
        
        int code = (v[0] > 0) | ((v[1] > 0) << 1) | ((v[2] > 0) << 2) | ((v[3] > 0) << 3);

        if (code == 0 || code == 15) return;

        Point2D m[4];
        m[0] = {x0 + dx * interpolation(v[0], v[1]), y0};      
        m[1] = {x0 + dx, y0 + dy * interpolation(v[1], v[2])}; 
        m[2] = {x0 + dx * interpolation(v[3], v[2]), y0 + dy}; 
        m[3] = {x0, y0 + dy * interpolation(v[0], v[3])};      
        switch (code) {
            case 1: case 14: segments.push_back({m[0], m[3]}); break;
            case 2: case 13: segments.push_back({m[0], m[1]}); break;
            case 3: case 12: segments.push_back({m[1], m[3]}); break;
            case 4: case 11: segments.push_back({m[1], m[2]}); break;
            case 6: case 9:  segments.push_back({m[0], m[2]}); break;
            case 7: case 8:  segments.push_back({m[2], m[3]}); break;
            case 5:
                segments.push_back({m[0], m[3]});
                segments.push_back({m[1], m[2]});
                break;
            case 10: 
                segments.push_back({m[0], m[1]});
                segments.push_back({m[2], m[3]});
                break;
        }
    }

    void process_all_local_cells() {
        all_local_segments.clear();
        for (int j = 0; j < local_N; ++j) {
            for (int i = 0; i < local_M; ++i) {
                process_cell_independent(i, j, all_local_segments);
            }
        }
    }

 
    void write_eps_parallel(const std::string& filename) {
        std::ostringstream local_eps_stream;
        local_eps_stream << std::fixed << std::setprecision(6);
        for (const auto& seg : all_local_segments) {
            local_eps_stream << seg.p1.x << " " << seg.p1.y << " moveto "
                             << seg.p2.x << " " << seg.p2.y << " lineto stroke\n";
        }
        std::string local_data = local_eps_stream.str();
        
        if (rank == 0) {
            std::ofstream eps_file(filename);
            eps_file << "%!PS-Adobe-3.0 EPSF-3.0\n";
            eps_file << "%%BoundingBox: " << x_min << " " << y_min << " " 
                     << x_max << " " << y_max << "\n";
            eps_file << "0.5 setlinewidth\n";
        }
        
        MPI_Barrier(cart_comm);
        
        MPI_File fh;
        MPI_File_open(cart_comm, filename.c_str(), 
                      MPI_MODE_WRONLY | MPI_MODE_APPEND, 
                      MPI_INFO_NULL, &fh);
        
        MPI_File_write_ordered(fh, local_data.c_str(), local_data.size(), 
                               MPI_CHAR, MPI_STATUS_IGNORE);
        
        MPI_File_close(&fh);
        
        if (rank == 0) {
            std::ofstream eps_file(filename, std::ios::app);
            eps_file << "showpage\n";
        }
    }

    void run() {
        double start_time, compute_end_time, total_end_time;

        MPI_Barrier(cart_comm);
        start_time = MPI_Wtime();

        process_all_local_cells();

        MPI_Barrier(cart_comm);
        compute_end_time = MPI_Wtime();

        write_eps_parallel("output_mpi_optimized.eps");
        
        MPI_Barrier(cart_comm);
        total_end_time = MPI_Wtime();

        if (rank == 0) {
            double compute_time = compute_end_time - start_time;
            double io_time = total_end_time - compute_end_time;
            double total_time = total_end_time - start_time;
            double cells_per_second = (static_cast<double>(global_M) * global_N) / compute_time;

            std::cout << "--- Resultados del Enfoque Optimizado ---\n";
            std::cout << "Grid: " << global_M << "x" << global_N 
                      << " en " << size << " procesos (" << dims[0] << "x" << dims[1] << ")\n";
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "Tiempo de Cómputo: " << compute_time * 1000 << " ms\n";
            std::cout << "Tiempo de I/O Paralelo: " << io_time * 1000 << " ms\n";
            std::cout << "Tiempo Total: " << total_time * 1000 << " ms\n";
            std::cout << std::scientific;
            std::cout << "Rendimiento: " << cells_per_second / 1e6 << " MCells/s\n";
            std::cout << std::defaultfloat;
            std::cout << "Salida guardada en: output_mpi_optimized.eps\n";
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

        MarchingSquaresMPI ms(x_min, x_max, y_min, y_max, M, N);
        ms.run();
    }

    MPI_Finalize();
    return 0;
}   