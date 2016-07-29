#include <vector>

#include "Plotting/Plotting.hpp"
#include "Plotting/gnuplot-iostream.h"

void Plotting::plotVector(std::vector<double> y, 
                     double xAxisScale, 
                     double xAxisOffset) {

  std::vector<std::pair<double, double> > xyPoints(y.size());
  for (long i = 0; i < y.size(); i++) {
    xyPoints[i] = std::make_pair(i * xAxisScale + xAxisOffset, y[i]);
  }

  plotPair(xyPoints);
}

void Plotting::plotPair(std::vector<std::pair<double, double> > xy) {
  Gnuplot gp;
  gp << "plot" << gp.file1d(xy) << "w l" << std::endl;
  std::cin.get();
}
