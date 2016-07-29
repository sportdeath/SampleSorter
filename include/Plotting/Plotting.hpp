#ifndef PLOTTING
#define PLOTTING

#include <vector>

class Plotting {
  public:
    static void plotVector(std::vector<double> y, double xAxisScale = 1, double xAxisOffset = 0);
    static void plotPair(std::vector<std::pair<double, double> > xy);
};

#endif
