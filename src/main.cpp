#include <iostream>
#include <vector>
#include <string>
#include "HardyZ.h"
#include "Plotter.h"

int main() {
    Color black(0,0,0);
    Color gray(100,100,100);
    Color white(255,255,255);
    Color blue(170,220,255);
    Color gold(255,215,0);

    double start_x = 10000.0;
    double end_x = 10100.0;  
    int frame = 6000;      
    int height = 300;
    int axis_y = height / 2; // Middle of screen


    std::function<double(double)> hardyEM = [](double x){
        return Zeta::Hardy::compute(x, Zeta::Method::EulerMaclaurin);
    };

    std::function<double(double)> hardyRS = [](double x){
        return Zeta::Hardy::compute(x, Zeta::Method::RiemannSiegel);
    };

    std::function<double(double)> theta = [](double t) {
        return Zeta::theta<double>(t);
    }; 
    
    system("mkdir -p output/frames_hardyEM");
    PlotCanvas(600, height)
        .fill_background(black)
        .draw_baseline(axis_y, gray)
        .animate_function("output/frames_hardyEM", hardyEM, start_x, end_x, frame, blue, gold);


    system("mkdir -p output/frames_hardyRS");
    PlotCanvas(600, height)
        .fill_background(black)
        .draw_baseline(axis_y, gray) 
        .animate_function("output/frames_hardyRS", hardyRS, start_x, end_x, frame, blue, gold);


    system("mkdir -p output/frames_zeta");
    PlotCanvas(600, 600).fill_background(black)
          .animate_complex_zeta(
              "output/frames_zeta",
              hardyEM, 
              theta,
              start_x, 
              end_x, 
              frame, 
              blue,
              gold
          );

    std::cout << "All tasks completed. \nTo create the video, run:" << std::endl;
    std::cout << "ffmpeg -framerate 300 -i output/frames_hardyEM/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output/hardyEM.mp4" << std::endl;
    std::cout << "ffmpeg -framerate 300 -i output/frames_hardyRS/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output/hardyRS.mp4" << std::endl;
    std::cout << "ffmpeg -framerate 300 -i output/frames_zeta/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output/zeta.mp4" << std::endl;
};
