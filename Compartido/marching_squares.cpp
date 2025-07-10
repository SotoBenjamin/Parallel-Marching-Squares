#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <omp.h>

using namespace std;

// ... (las funciones f, out, draw_segment, interpolation no cambian) ...
double f(double x, double y) {
    return x*x + y*y - 200.0*200.0;
}

bool out(double x, double y) {
    return f(x, y) > 0;
}

void draw_segment(ofstream& eps,
                         double x0, double y0,
                         double x1, double y1) {
    eps << x0 << " " << y0 << " moveto\n";
    eps << x1 << " " << y1 << " lineto\n";
    eps << "stroke\n";
}

double interpolation(double x0, double y0,
                     double x1, double y1){
    double f0 = f(x0,y0);
    double f1 = f(x1,y1);
    double df = f1 - f0;
    return f0/(-df);
}

struct Point2D {
    double x, y;
};

// Estructura para guardar un segmento de línea
struct Segment {
    Point2D p1, p2;
};

void marching_squares(double x_min, double x_max,double y_min, double y_max,int M, int N,ofstream& eps) {
    double dx = (x_max - x_min) / M;
    double dy = (y_max - y_min) / N;

    vector<vector<Segment>> thread_segments;
    int num_threads = 1;

    num_threads = omp_get_max_threads();
    thread_segments.resize(num_threads);

    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
           
            int tid = omp_get_thread_num();
            double x0 = x_min + i * dx;
            double x1 = x_min + (i + 1) * dx;
            double y0 = y_min + j * dy;
            double y1 = y_min + (j + 1) * dy;

            bool b0 = out(x0, y0);
            bool b1 = out(x1, y0);
            bool b2 = out(x1, y1);
            bool b3 = out(x0, y1);
            int code = (b0 ? 1 : 0) | (b1 ? 2 : 0) | (b2 ? 4 : 0) | (b3 ? 8 : 0);
            
            if (code == 0 || code == 15) continue;

            Point2D m[4];
            m[0] = { x_min + (i + interpolation(x0, y0, x1, y0)) * dx, y0 };
            m[1] = { x1, y_min + (j + interpolation(x1, y0, x1, y1)) * dy };
            m[2] = { x_min + (i + interpolation(x0, y1, x1, y1)) * dx, y1 };
            m[3] = { x0, y_min + (j + interpolation(x0, y0, x0, y1)) * dy };
            
     
            switch (code) {
                case 1: case 14: thread_segments[tid].push_back({m[0], m[3]}); break;
                case 2: case 13: thread_segments[tid].push_back({m[0], m[1]}); break;
                case 3: case 12: thread_segments[tid].push_back({m[1], m[3]}); break;
                case 4: case 11: thread_segments[tid].push_back({m[1], m[2]}); break;
                case 5: 
                    thread_segments[tid].push_back({m[0], m[3]});
                    thread_segments[tid].push_back({m[1], m[2]});
                    break;
                case 6: case 9:  thread_segments[tid].push_back({m[0], m[2]}); break;
                case 7: case 8:  thread_segments[tid].push_back({m[2], m[3]}); break;
                case 10:
                    thread_segments[tid].push_back({m[0], m[1]});
                    thread_segments[tid].push_back({m[2], m[3]});
                    break;
            }
        }
    }

    for (const auto& segments : thread_segments) {
        for (const auto& seg : segments) {
            draw_segment(eps, seg.p1.x, seg.p1.y, seg.p2.x, seg.p2.y);
        }
    }
}

int main(int argc, char* argv[]) {

    int num_threads = omp_get_max_threads();
    if (argc > 1) {
        try {
            num_threads = stoi(argv[1]);
        } catch (const exception& e) {
            cerr << "Argumento invalido. Se necesita un numero entero para los hilos." << endl;
            return 1;
        }
    }
    omp_set_num_threads(num_threads);
    cout << "Ejecutando en modo PARALELO con " << num_threads << " hilos." << endl;
    
    ofstream eps("output.eps");
    eps << "%!PS-Adobe-3.0 EPSF-3.0\n";
    eps << "%%BoundingBox: -600 -600 600 600\n";

    double start = omp_get_wtime();
    marching_squares(-300.0, 300.0, -300.0, 300.0, 1000, 1000, eps);
    double end = omp_get_wtime();
    
    double elapsed_time = end - start;
    
    cout << "Tiempo de ejecucion: " << (elapsed_time * 1000.0) << " ms" << endl;

    eps << "showpage\n";
    eps.close();
    return 0;
}