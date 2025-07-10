#include <bits/stdc++.h>
#include <chrono>
using namespace std;
using namespace std::chrono;

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

void marching_squares(double x_min, double x_max,double y_min, double y_max,int M, int N,ofstream& eps) {
    double dx = (x_max - x_min) / M;
    double dy = (y_max - y_min) / N;
    
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            double x0 = x_min + i * dx;
            double x1 = x_min + (i + 1) * dx;
            double y0 = y_min + j * dy;
            double y1 = y_min + (j + 1) * dy;
            
            bool b0 = out(x0, y0);
            bool b1 = out(x1, y0);
            bool b2 = out(x1, y1);
            bool b3 = out(x0, y1);
            
            int code = (b0 ? 1 : 0) |
                      (b1 ? 1 << 1: 0) |
                      (b2 ? 1 << 2 : 0) |
                      (b3 ? 1 << 3 : 0);
            
            double t_bottom = interpolation(x0,y0,x1,y0);
            double t_right = interpolation(x1,y0,x1,y1);
            double t_top = interpolation(x1,y1,x0,y1);
            double t_left = interpolation(x0,y1,x0,y0);
            
            double xm_bottom = (x0 + dx*t_bottom);
            double ym_right = (y0 + dy*t_right);
            double xm_top = (x1 - dx*t_top);
            double ym_left = (y1 - dy*t_left);
            
            Point2D m[4] = {
                { xm_bottom , y0 }, // bottom
                { x1, ym_right}, // right
                { xm_top , y1 }, // top
                { x0, ym_left } // left
            };
            
            switch (code) {
                case 0: case 15:
                    break;
                case 1: case 14:
                    draw_segment(eps, m[3].x, m[3].y, m[0].x, m[0].y);
                    break;
                case 2: case 13:
                    draw_segment(eps, m[0].x, m[0].y, m[1].x, m[1].y);
                    break;
                case 3: case 12:
                    draw_segment(eps, m[3].x, m[3].y, m[1].x, m[1].y);
                    break;
                case 4: case 11:
                    draw_segment(eps, m[1].x, m[1].y, m[2].x, m[2].y);
                    break;
                case 5:
                    draw_segment(eps, m[3].x, m[3].y, m[0].x, m[0].y);
                    draw_segment(eps, m[1].x, m[1].y, m[2].x, m[2].y);
                    break;
                case 6: case 9:
                    draw_segment(eps, m[0].x, m[0].y, m[2].x, m[2].y);
                    break;
                case 7: case 8:
                    draw_segment(eps, m[3].x, m[3].y, m[2].x, m[2].y);
                    break;
                case 10:
                    draw_segment(eps, m[0].x, m[0].y, m[1].x, m[1].y);
                    draw_segment(eps, m[2].x, m[2].y, m[3].x, m[3].y);
                    break;
            }
        }
    }
}

int main() {
    ofstream eps("output.eps");
    eps << "%!PS-Adobe-3.0 EPSF-3.0\n";
    eps << "%%BoundingBox: -600 -600 600 600\n";
    
    auto start = high_resolution_clock::now();
    
    marching_squares(-300.0, 300.0, -300.0, 300.0, 1000, 1000, eps);
    
    auto end = high_resolution_clock::now();
    
    auto duration_ms = duration_cast<milliseconds>(end - start);

    
    cout << "Tiempo de la función marching_squares:" << endl;
    cout << duration_ms.count() << " ms" << endl;

    eps << "showpage\n";
    eps.close();
    
    return 0;
}