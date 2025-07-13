#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <chrono>
#include <memory>

using namespace std;
using namespace std::chrono;

struct Point2D {
    double x, y;
};

struct Segment {
    Point2D p1, p2;
};

struct alignas(64) ThreadData {
    vector<Segment> segments;
    char padding[64 - sizeof(vector<Segment>)]; 
};

inline double f(double x, double y) {
    return x * x + y * y - 200.0 * 200.0;
}

inline double interpolation(double f0, double f1) {
    double diff = f0 - f1;
    return (diff == 0.0) ? 0.5 : f0 / diff;
}

void draw_segment(stringstream& ss, double x1, double y1, double x2, double y2) {
    ss << x1 << " " << y1 << " moveto\n";
    ss << x2 << " " << y2 << " lineto\n";
    ss << "stroke\n";
}

void marching_squares(double x_min, double x_max, double y_min, double y_max, int M, int N, stringstream& final_buffer) {
    const double dx = (x_max - x_min) / M;
    const double dy = (y_max - y_min) / N;

    int max_threads = omp_get_max_threads();
    
    vector<unique_ptr<ThreadData>> thread_data(max_threads);
    for (int i = 0; i < max_threads; i++) {
        thread_data[i] = make_unique<ThreadData>();
    }

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        
        int estimated_segments = (M * N) / num_threads / 25; 
        thread_data[tid]->segments.reserve(estimated_segments);

        #pragma omp for schedule(static) nowait
        for (int idx = 0; idx < M * N; idx++) {
            int i = idx / N;
            int j = idx % N;

            const double x0 = x_min + i * dx;
            const double x1 = x0 + dx;
            const double y0 = y_min + j * dy;
            const double y1 = y0 + dy;
            
            const double f_vals[4] = {
                f(x0, y0), f(x1, y0), f(x1, y1), f(x0, y1)
            };
            
            const int code = (f_vals[0] > 0) | 
                           ((f_vals[1] > 0) << 1) | 
                           ((f_vals[2] > 0) << 2) | 
                           ((f_vals[3] > 0) << 3);

            if (code == 0 || code == 15) continue;

            Point2D m[4];
            
            auto calc_m0 = [&]() { m[0] = {x0 + dx * interpolation(f_vals[0], f_vals[1]), y0}; };
            auto calc_m1 = [&]() { m[1] = {x1, y0 + dy * interpolation(f_vals[1], f_vals[2])}; };
            auto calc_m2 = [&]() { m[2] = {x0 + dx * interpolation(f_vals[3], f_vals[2]), y1}; };
            auto calc_m3 = [&]() { m[3] = {x0, y0 + dy * interpolation(f_vals[0], f_vals[3])}; };
            
            auto& segments = thread_data[tid]->segments;
            
            switch (code) {
                case 1: case 14: 
                    calc_m0(); calc_m3();
                    segments.push_back({m[0], m[3]}); 
                    break;
                case 2: case 13: 
                    calc_m0(); calc_m1();
                    segments.push_back({m[0], m[1]}); 
                    break;
                case 3: case 12: 
                    calc_m1(); calc_m3();
                    segments.push_back({m[1], m[3]}); 
                    break;
                case 4: case 11: 
                    calc_m1(); calc_m2();
                    segments.push_back({m[1], m[2]}); 
                    break;
                case 5:
                    calc_m0(); calc_m1(); calc_m2(); calc_m3();
                    segments.push_back({m[0], m[3]});
                    segments.push_back({m[1], m[2]});
                    break;
                case 6: case 9:  
                    calc_m0(); calc_m2();
                    segments.push_back({m[0], m[2]}); 
                    break;
                case 7: case 8:  
                    calc_m2(); calc_m3();
                    segments.push_back({m[2], m[3]}); 
                    break;
                case 10:
                    calc_m0(); calc_m1(); calc_m2(); calc_m3();
                    segments.push_back({m[0], m[1]});
                    segments.push_back({m[2], m[3]});
                    break;
            }
        }
    } 

    size_t total_size = 0;
    for (int i = 0; i < max_threads; i++) {
        total_size += thread_data[i]->segments.size();
    }
    
    final_buffer.str("");
    final_buffer.clear();
    
    for (int i = 0; i < max_threads; i++) {
        for (const auto& seg : thread_data[i]->segments) {
            draw_segment(final_buffer, seg.p1.x, seg.p1.y, seg.p2.x, seg.p2.y);
        }
    }
}

int main(int argc, char* argv[]) {
    int M = 1000, N = 1000, num_threads = 0;

    if (argc > 1) M = stoi(argv[1]);
    if (argc > 2) N = stoi(argv[2]);
    if (argc > 3) num_threads = stoi(argv[3]);

    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    
    #pragma omp parallel
    #pragma omp single
    num_threads = omp_get_num_threads();
    
    cout << "Ejecutando OMP (VERSION OPTIMIZADA) con " << num_threads << " hilos y grid " << M << "x" << N << "." << endl;
    
    stringstream final_buffer;
    
    auto start = high_resolution_clock::now();
    marching_squares(-300.0, 300.0, -300.0, 300.0, M, N, final_buffer);
    auto end = high_resolution_clock::now();
    
    auto duration_ms = duration_cast<milliseconds>(end - start);
    
    cout << "FINAL_TIME_MS:" << duration_ms.count() << endl;

    ofstream eps("output_omp_optimized.eps");
    eps << "%!PS-Adobe-3.0 EPSF-3.0\n";
    eps << "%%BoundingBox: -600 -600 600 600\n";
    eps << final_buffer.str();
    eps << "showpage\n";
    eps.close();
    
    return 0;
}