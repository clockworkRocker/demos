#include "solver.h"

#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
    // -------- Check arguments --------
    if (argc < 2){
        std::cout << "Usage: " << argv[0] << "<equation coeffs file> [output path]\n";
        return 0;
    }

    // -------- Initialise everything --------
    std::fstream inputFile, outputFileLinear, outputFileCubic;
    std::string outputPath((argc > 2) ? argv[2] : "");
    const std::string outputPathLinear = outputPath + "linear.csv";
    const std::string outputPathCubic = outputPath + "cubic.csv";

    inputFile.open(argv[1], std::ios::in);
    if(!inputFile.is_open()){
        std::cout << "main: could not open file " << argv[1] << '\n';
        return 1;
    }

    LinearODE<2> equation{};
    for(unsigned i = 0; i < equation.order() + 2; ++i) {
        inputFile >> equation.coeffs[i];
        if(inputFile.bad()){
            std::cout << "main: error reading input file\n";
            inputFile.close();
            return 1;
        }
    }
    std::vector<BoundaryCondition> conditions;
    conditions.emplace_back(boundary_t::FIRST, 1, -6);
    //conditions.emplace_back(FIRST, 8, -165.267);
    conditions.emplace_back(boundary_t::SECOND, 8, -6);

    // -------- Calculate linear --------
    outputFileLinear.open(outputPathLinear, std::ios::out);
    Solver1D<LinearODE<2>> solver1(ElemType::Linear1D, outputFileLinear);
    solver1.solve(equation, conditions, {1, 8}, 40);
    outputFileLinear.close();

    // -------- Calculate cubic --------
    outputFileCubic.open(outputPathCubic, std::ios::out);
    Solver1D<LinearODE<2>> solver2(ElemType::Cubic1D, outputFileCubic);
    solver2.solve(equation, conditions, {1, 8}, 40);
    outputFileCubic.close();

    return 0;
}
