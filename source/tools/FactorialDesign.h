#ifndef FACTORIAL_DESIGN_H
#define FACTORIAL_DESIGN_H

#include <vector>
#include <string>
#include <tuple>

class FactorialDesign {
public:
    FactorialDesign();
    void showSignTable();
    void showImpacts();
    void showLinearRegressionModel();
    void showResiduals();

private:
    int k;
    int r;
    int NUM_ROWS;
    std::vector<std::vector<int>> resultsMatrix;
    std::vector<std::vector<int>> indexCombinations;
    std::vector<std::string> factorLabels;
    std::vector<std::string> columnLabels;
    std::vector<std::vector<int>> inputs;
    std::vector<double> resultsMean;
    std::vector<std::tuple<int, std::vector<int>, std::vector<int>, double>> signTable;
    std::vector<double> effectSums;
    std::vector<double> effectAverages;
    std::vector<int> squareSums;
    double squareSumTotal;
    std::vector<double> factorVariations;
    double squareSumError;
    double errorVariation;

    void generateIndexCombinations();
    std::vector<std::vector<int>> generateCombinations(const std::vector<int>& elements, int length);
    void createColumnLabels();
    void generateInputs();
    void calculateResultsMean();
    void createTable();
    void calculateStatistics();
    double linearRegression(const std::vector<int>& xs);
};

#endif
