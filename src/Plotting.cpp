#include <vector>

#include "Plotting/Plotting.hpp"
#include "Plotting/gnuplot-iostream.h"

void Plotting::plotVector(
    std::vector<double> y, 
    std::string title,
    std::string xaxis,
    std::string yaxis,
    double xAxisScale, 
    double xAxisOffset,
    bool histogram,
    double logscale,
    std::vector<double> markers) {

  std::vector<std::pair<double, double> > xyPoints(y.size());
  for (long i = 0; i < y.size(); i++) {
    xyPoints[i] = std::make_pair(i * xAxisScale + xAxisOffset, y[i]);
  }

  plotPair(xyPoints, title, xaxis, yaxis, histogram, logscale, markers);
}

void Plotting::plotPair(
    std::vector<std::pair<double, double> > xy, 
    std::string title,
    std::string xaxis,
    std::string yaxis,
    bool histogram,
    double logscale,
    std::vector<double> markers) {

  Gnuplot gp;

  if (markers.size() > 0) {
    for (auto const & marker: markers) {
      gp << "set arrow from " << std::fixed << std::setprecision(3) << marker << ",0 to " << marker << ",8 nohead lc rgb \'#4183c4\' lw 1.5 lt 2 dt 2" << std::endl;
    }
  }

  if (histogram) {
    gp << "set style data histogram" << std::endl;
    gp << "set style histogram cluster gap 1" << std::endl;
    gp << "set style fill solid border -1" << std::endl;
    gp << "set boxwidth 2 absolute" << std::endl;
  } else {
    auto comparison = 
        [](std::pair<double, double> a, std::pair<double, double> b) {
            return a.first < b.first;
        };
    auto min = std::min_element(xy.begin(), xy.end(), comparison);
    auto max = std::max_element(xy.begin(), xy.end(), comparison);
    gp << "set xrange [" << min -> first << ":" << max -> first <<  "]" << std::endl;
  }

  if (logscale > 0) {
    gp  << "set logscale x" << std::endl;
    gp << "set xrange [" << logscale << ":]" << std::endl;
  }

  gp << "set title font \"Helvetica:Bold,18\" tc \"#24292e\"" << std::endl;
  gp << "set xlabel font \"Helvetica:Bold,15\" tc \"#24292e\"" << std::endl;
  gp << "set ylabel font \"Helvetica:Bold,15\" tc \"#24292e\"" << std::endl;
  gp << "set xtics font \"Helvetica,12\" tc \"#24292e\"" << std::endl;
  gp << "set ytics font \"Helvetica,12\" tc \"#24292e\"" << std::endl;
  gp << "set key off"  << std::endl;
  gp << "set title \"" << title << "\"" << std::endl;
  gp << "set xlabel \"" << xaxis << "\"" << std::endl;
  gp << "set ylabel \"" << yaxis << "\"" << std::endl;
  gp << "plot" << gp.file1d(xy);
  if (histogram) {
    gp << "using 2:xtic(1)";
    gp << " lc rgb \"#24292e\"";
  } else {
    gp << "w l";
    gp << " lt rgb \"#24292e\"";
  }
  gp << std::endl;
  std::cin.get();
}
