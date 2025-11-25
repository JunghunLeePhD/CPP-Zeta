#include "Plotter.h"
#include "HardyZ.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <functional>
#include <ranges>
#include <complex> 

PlotCanvas::PlotCanvas(int w, int h) : width(w), height(h) {
    pixels.resize(width * height * 3);
    std::fill(pixels.begin(), pixels.end(), 0); 
}

void PlotCanvas::set_pixel(int x, int y, const Color& c) {
    if (x >= 0 && x < width && y >= 0 && y < height) {
        int idx = (y * width + x) * 3;
        pixels[idx] = c.r; 
        pixels[idx+1] = c.g; 
        pixels[idx+2] = c.b;
    }
}

int PlotCanvas::to_screen_x(double val) const {
    return static_cast<int>(val * (width - 1));
}

int PlotCanvas::to_screen_y(double val) const {
    return height - 1 - static_cast<int>(val * (height - 1));
}

void PlotCanvas::draw_line_raw(int x0, int y0, int x1, int y1, const Color& c) {
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy, e2;
    while (true) {
        set_pixel(x0, y0, c);
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

PlotCanvas& PlotCanvas::fill_background(const Color& c) {
    for (size_t i = 0; i < pixels.size(); i += 3) {
        pixels[i] = c.r; 
        pixels[i+1] = c.g; 
        pixels[i+2] = c.b;
    }
    return *this;
}

PlotCanvas& PlotCanvas::draw_function(std::function<double(double)> func, const Color& c) {
    int prev_px = to_screen_x(0);
    int prev_py = to_screen_y(func(0));
    for (int i = 1; i < width; ++i) {
        double x = static_cast<double>(i) / width;
        double y = func(x);
        if (y > 1.0 || y < 0.0) continue;
        int px = to_screen_x(x);
        int py = to_screen_y(y);
        draw_line_raw(prev_px, prev_py, px, py, c);
        prev_px = px; 
        prev_py = py;
    }
    return *this;
}

PlotCanvas& PlotCanvas::draw_baseline(int y_pos, const Color& c) {
    for(int x=0; x<width; ++x){
        set_pixel(x, y_pos, c);
    }
    return *this;
}


int map_val(double val, double min_val, double max_val, int screen_size) {
    double t = (val - min_val) / (max_val - min_val);
    return static_cast<int>(t * (screen_size - 1));
}

int map_y_val(double val, double min_val, double max_val, int screen_height) {
    double t = (val - min_val) / (max_val - min_val);
    return screen_height - 1 - static_cast<int>(t * (screen_height - 1));
}

void PlotCanvas::animate_function(const std::string& folder, 
                                  std::function<double(double)> func, 
                                  const double start_x, const double end_x, 
                                  const int total_frames, 
                                  const Color& startC, const Color& endC) {
    
    std::cout << "Animating in " << folder << "..." << std::endl;
    
    double view_min_x = start_x;
    double view_max_x = end_x;
    
    // We need to define a Y range that fits the Hardy Z function (usually +/- 5 or 10)
    double view_min_y = -6.0;
    double view_max_y =  6.0;

    double step = (end_x - start_x) / total_frames;
    double curr_x = start_x;
    double curr_y = func(curr_x);

    int prev_px = map_val(curr_x, view_min_x, view_max_x, width);
    int prev_py = map_y_val(curr_y, view_min_y, view_max_y, height);

    for (int i = 0; i < total_frames; ++i) {
        double next_x = curr_x + step;
        double next_y = func(next_x);

        int px = map_val(next_x, view_min_x, view_max_x, width);
        int py = map_y_val(next_y, view_min_y, view_max_y, height);

        double t = (double)i / (total_frames - 1);
        Color c = Color::lerp(startC, endC, t);

        draw_line_raw(prev_px, prev_py, px, py, c);

        curr_x = next_x;
        prev_px = px;
        prev_py = py;

        std::stringstream ss;
        ss << folder << "/frame_" << std::setfill('0') << std::setw(4) << i << ".ppm";
        save(ss.str());
        
        if (i % 50 == 0) std::cout << "Frame " << i << "\r" << std::flush;
    }
    std::cout << "\nDone." << std::endl;
}

void PlotCanvas::animate_complex_zeta(const std::string& folder,
                                      std::function<double(double)> hardy_func,
                                      std::function<double(double)> theta_func,
                                      double t_start, double t_end,
                                      int total_frames,
                                      const Color& startC,
                                      const Color& endC) {
    
    std::cout << "Animating Complex Zeta in " << folder << "..." << std::endl;

    double view_min = -8.0;
    double view_max =  8.0;

    auto to_screen = [&](std::complex<double> z) -> std::pair<int, int> {
        return {
            map_val(z.real(), view_min, view_max, width),
            map_y_val(z.imag(), view_min, view_max, height)
        };
    };

    int center_x = map_val(0, view_min, view_max, width);
    int center_y = map_y_val(0, view_min, view_max, height);
    Color axis_col = {80, 80, 80};
    draw_line_raw(0, center_y, width, center_y, axis_col);
    draw_line_raw(center_x, 0, center_x, height, axis_col);

    // Initial State
    double step = (t_end - t_start) / total_frames;
    
    // RENAMED: Use 'current_t' to avoid confusion
    double current_t = t_start; 

    double z_val = hardy_func(current_t);
    double th_val = theta_func(current_t);
    
    std::complex<double> current_zeta = std::polar(z_val, -th_val);
    auto [prev_px, prev_py] = to_screen(current_zeta);

    for (int i = 0; i < total_frames; ++i) {
        // FIX: Rename loop variable to 'progress' for Color Lerp
        double progress = (double)i / (total_frames - 1);
        Color c = Color::lerp(startC, endC, progress);

        // FIX: Use current_t
        double next_t = current_t + step;
        
        double next_z = hardy_func(next_t);
        double next_th = theta_func(next_t);
        std::complex<double> next_zeta = std::polar(next_z, -next_th);

        auto [px, py] = to_screen(next_zeta);

        draw_line_raw(prev_px, prev_py, px, py, c);

        // Update State
        current_t = next_t;
        prev_px = px;
        prev_py = py;

        std::stringstream ss;
        ss << folder << "/frame_" << std::setfill('0') << std::setw(4) << i << ".ppm";
        save(ss.str());

        if (i % 50 == 0) std::cout << "Frame " << i << " (t=" << current_t << ")\r" << std::flush;
    }
    std::cout << "\nComplex Animation Done." << std::endl;
}

void PlotCanvas::save(const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) { std::cerr << "Error opening " << filename << std::endl; return; }
    file << "P6\n" << width << " " << height << "\n255\n";
    file.write(reinterpret_cast<const char*>(pixels.data()), pixels.size());
}