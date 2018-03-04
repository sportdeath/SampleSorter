#ifndef PLOTTING
#define PLOTTING

#include <vector>
#include <string>
#include "Plotting/gnuplot-iostream.h"

class Plotting {
  public:
    static void plotVector(
        std::vector<double> y, 
        std::string title = "Title",
        std::string xaxis = "X Axis",
        std::string yaxis = "Y Axis",
        double xAxisScale = 1, 
        double xAxisOffset = 0,
        bool histogram = false,
        double logscale = 0,
        std::vector<double> markers = {});
    static void plotPair(
        std::vector<std::pair<double, double> > xy, 
        std::string title = "Title",
        std::string xaxis = "X Axis",
        std::string yaxis = "Y Axis",
        bool histogram = false,
        double logscale = false,
        std::vector<double> markers = {});
};

#endif
