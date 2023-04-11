#include <iostream>
#include <chrono>
#include "thermal.h"

using Clock = std::chrono::high_resolution_clock;
using seconds_fp = std::chrono::duration<double, std::chrono::seconds::period>;

int main(int argc, char* argv[]) {
  const double totalTime = 5e-5;

  if (argc < 2) {
    std::cerr << "Usage: ./" << argv[0] << " <output filename>" << std::endl;
    return 0;
  }

  std::vector<Line2d*> edges;
  std::fstream file(argv[1], std::fstream::trunc | std::fstream::out);
  if (!file.is_open()) {
    std::cerr << "Error: could not open file " << argv[1] << '\n';
    return 0;
  }

  MFDSolver2d solver;
  size_t nodeIndex;
  auto tick = Clock::now();

  edges.push_back(new StraightLine2d(0, 0, 8, 0));
  edges.push_back(new StraightLine2d(8, 0, 8, 3));
  edges.push_back(new StraightLine2d(8, 3, 5, 6));
  edges.push_back(new StraightLine2d(5, 6, 0, 6));
  edges.push_back(new StraightLine2d(0, 6, 0, 0));

  ThermoShape2d customPlate(edges);
  customPlate.generateMesh(0.001, 0.001);
  if (!customPlate.lookupNode(6.5, 4.5, &nodeIndex))
    return 0;

  customPlate.addBoundary(GM_EDGE, 0, TM_TEMPERATURE, 100);
  customPlate.addBoundary(GM_EDGE, 1, TM_TEMPERATURE, 50);
  customPlate.addBoundary(GM_EDGE, 2, TM_TEMPERATURE, 50);
  customPlate.addBoundary(GM_EDGE, 3, TM_TEMPERATURE, 50);
  customPlate.addBoundary(GM_EDGE, 4, TM_TEMPERATURE, 300);
  customPlate.addBoundary(GM_NODE, nodeIndex, TM_CONVECTION, 0);

  std::cout << "Calculating grid of " << 8 / 0.001 << 'x' << 6 / 0.001 << '\n';

  solver.solve(customPlate, totalTime, file);
  // solver.printTimeLayer(file);
  auto time = Clock::now() - tick;
  std::cout << "The calculation of " << (int)(totalTime / solver.step())
            << " iterations took "
            << std::chrono::duration_cast<seconds_fp>(time).count()
            << " seconds\n";

  for (auto& edge : edges) {
    delete edge;
    edge = nullptr;
  }
  file.close();
  return 0;
}
