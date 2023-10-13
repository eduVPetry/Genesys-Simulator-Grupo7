#include <iostream>
#include <iomanip>
#include <sstream>
#include <numeric>
#include <cmath>
#include <tuple>
#include "FactorialDesign.h"

void FactorialDesign::generateIndexCombinations() {
    std::vector<int> indexes(k);
    std::iota(indexes.begin(), indexes.end(), 0);

    for (int length = 1; length <= k; ++length) {
        std::vector<std::vector<int>> combinations = generateCombinations(indexes, length);
        indexCombinations.insert(indexCombinations.end(), combinations.begin(), combinations.end());
    }
}

std::vector<std::vector<int>> FactorialDesign::generateCombinations(const std::vector<int>& elements, int length) {
    std::vector<std::vector<int>> result;
    int n = elements.size();

    if (length <= 0 || length > n) return result;  // Handle invalid input

    // Initialize the first combination (e.g., {0, 1, 2, ..., length-1})
    std::vector<int> combination(length);
    std::iota(combination.begin(), combination.end(), 0);

    while (combination[0] < n - length + 1) {
        std::vector<int> current_combination;
        for (int index : combination) current_combination.push_back(elements[index]);
        result.push_back(current_combination);

        // Find the rightmost element that can be incremented
        int i = length - 1;
        while (i >= 0 && combination[i] == i + n - length) i--;

        if (i < 0) break; // All combinations have been generated

        // Increment the found element and set the following elements
        // to be consecutive values
        combination[i]++;
        for (int j = i + 1; j < length; ++j) combination[j] = combination[j - 1] + 1;
    }

    return result;
}

void FactorialDesign::createColumnLabels() {
    factorLabels.resize(k);
    for (int i = 0; i < k; ++i)
        factorLabels[i] = "F" + std::to_string(i + 1);

    columnLabels = {"Exp.", "I"};
    std::string label;
    for (const std::vector<int>& indexCombination : indexCombinations) {
        label = "";
        for (int index : indexCombination)
            label += factorLabels[index];
        columnLabels.push_back(label);
    }
    for (int i = 1; i <= r; ++i)
        columnLabels.push_back("R" + std::to_string(i));
    columnLabels.push_back("R~");
}

void FactorialDesign::generateInputs() {
    for (int i = 0; i < NUM_ROWS; ++i) {
        std::vector<int> input(k);
        for (int j = 0; j < k; ++j)
            input[j] = (i & (1 << j)) ? 1 : -1;
        inputs.push_back(input);
    }
}

void FactorialDesign::calculateResultsMean() {
    double sum;
    for (int j = 0; j < NUM_ROWS; ++j) {
        sum = 0;
        for (int i = 0; i < r; ++i) {
            sum += resultsMatrix[i][j];
        }
        resultsMean.push_back(sum / r);
    }
}

void FactorialDesign::createTable() {
    int product;
    for (int i = 0; i < NUM_ROWS; ++i) {
        const std::vector<int>& input = inputs[i];
        std::vector<int> inputCombination{1};
        for (const std::vector<int>& indexCombination : indexCombinations) {
            product = 1;
            for (int index : indexCombination)
                product *= input[index];
            inputCombination.push_back(product);
        }
        std::vector<int> results;
        for (int j = 0; j < r; ++j)
            results.push_back(resultsMatrix[j][i]);
        signTable.push_back(std::make_tuple(i+1, inputCombination, results, resultsMean[i]));
    }
}

void FactorialDesign::showSignTable() {
    int MAX_LABEL_LENGTH = 0;
    for (const std::string& label : columnLabels)
        MAX_LABEL_LENGTH = std::max(MAX_LABEL_LENGTH, static_cast<int>(label.length()));

    auto center = [MAX_LABEL_LENGTH](const std::string& s) {
        std::string padded = s;
        int maxLength = std::max(MAX_LABEL_LENGTH, static_cast<int>(s.length()));
        int numSpaces = maxLength - s.length();
        int leftSpaces = numSpaces / 2;
        int rightSpaces = numSpaces - leftSpaces;
        return std::string(leftSpaces, ' ') + padded + std::string(rightSpaces, ' ');
    };

    for (const std::string& label : columnLabels)
        std::cout << center(label) << " ";
    std::cout << std::endl;

    for (const auto& row : signTable) {
        std::cout << center(std::to_string(std::get<0>(row))) << " ";
        const std::vector<int>& inputCombination = std::get<1>(row);
        for (int input : inputCombination)
            std::cout << center(std::to_string(input)) << " ";
        const std::vector<int>& results = std::get<2>(row);
        for (int result : results)
            std::cout << center(std::to_string(result)) << " ";
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2) << std::get<3>(row);
        std::cout << center(ss.str()) << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void FactorialDesign::calculateStatistics() {
    effectSums.resize(NUM_ROWS);
    for (int j = 0; j < NUM_ROWS; ++j)
        for (int i = 0; i < NUM_ROWS; ++i)
            effectSums[j] += std::get<1>(signTable[i])[j] * resultsMean[i];

    for (int i = 0; i < NUM_ROWS; ++i)
        effectAverages.push_back(effectSums[i] / NUM_ROWS);

    squareSums.resize(NUM_ROWS - 1);
    for (int i = 1; i < NUM_ROWS; ++i)
        squareSums[i - 1] = NUM_ROWS * r * pow(effectAverages[i], 2);

    squareSumError = 0;
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < NUM_ROWS; ++j)
            squareSumError += pow(resultsMatrix[i][j] - resultsMean[j], 2);

    squareSumTotal = std::accumulate(squareSums.begin(), squareSums.end(), 0);
    squareSumTotal += squareSumError;

    factorVariations.resize(NUM_ROWS - 1);
    for (int i = 0; i < NUM_ROWS - 1; ++i)
        factorVariations[i] = round(10000.0 * squareSums[i] / squareSumTotal) / 100.0;

    errorVariation = round(10000.0 * squareSumError / squareSumTotal) / 100.0;
}

double FactorialDesign::linearRegression(const std::vector<int>& xs) {
    double y = effectAverages[0];
    double product;
    for (int i = 1; i < NUM_ROWS; ++i) {
        product = effectAverages[i];
        const std::vector<int>& indexCombination = indexCombinations[i-1];
        for (int index : indexCombination)
            product *= xs[index];
        y += product;
    }
    return y;
}

void FactorialDesign::showImpacts() {
    for (int i = 0; i < NUM_ROWS - 1; ++i)
        std::cout << "Impact of " << columnLabels[i+2] << ": " << factorVariations[i] << "%" << std::endl;
    std::cout << "Impact of experimental error: " << errorVariation << "%" << std::endl << std::endl;
}

void FactorialDesign::showLinearRegressionModel() {
    std::cout << "R^ = " << effectAverages[0];
    for (int i = 1; i < NUM_ROWS; ++i) {
        std::cout << " + " << effectAverages[i];
        const std::vector<int>& indexCombination = indexCombinations[i-1];
        for (int index : indexCombination)
            std::cout << "*x" << factorLabels[index];
    }
    std::cout << std::endl << std::endl;
}

void FactorialDesign::showResiduals() {
    std::cout << "Residuals:" << std::endl;
    double residual;
    for (int i = 0; i < NUM_ROWS; ++i) {
        std::vector<int> input = inputs[i];
        double modelResult = linearRegression(input);
        for (int j = 0; j < r; ++j) {
            residual = resultsMatrix[j][i] - modelResult;
            std::cout << " " << residual;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}

// Constructor
FactorialDesign::FactorialDesign() {
    std::cout << "======== Factorial Design of Experiments ========" << std::endl << std::endl;

    std::cout << "Enter the number of factors (k): ";
    std::cin >> k;
    std::cout << "Enter the number of replications (r): ";
    std::cin >> r;

    NUM_ROWS = 1 << k;
    std::vector<std::vector<int>> resultsMatrix(r, std::vector<int>(NUM_ROWS));

    bool generateRandom = false;
    std::string command;
    do {
        std::cout << std::endl;
        std::cout << "(1) Manually enter experimental results" << std::endl;
        std::cout << "(2) Generate random experimental results" << std::endl;
        std::cin >> command;
        if (command == "2") generateRandom = true;
    } while (command != "1" && command != "2");

    if (generateRandom) {
        srand(static_cast<unsigned>(time(0)));
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < NUM_ROWS; ++j)
                resultsMatrix[i][j] = rand() % 100;
    } else {
        std::cout << std::endl;
        std::cout << "Enter " << r << " sets of " << NUM_ROWS <<
        " experimental results, separated by spaces:" << std::endl;
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < NUM_ROWS; ++j)
                std::cin >> resultsMatrix[i][j];
    }
    std::cout << std::endl;

    generateIndexCombinations();
    createColumnLabels();
    generateInputs();
    calculateResultsMean();
    createTable();
    calculateStatistics();

    std::cout << "Factorial design created successfully." << std::endl << std::endl;
    do {
        std::cout << "(1) Show sign table" << std::endl;
        std::cout << "(2) Show impact of factors, interactions and experimental error" << std::endl;
        std::cout << "(3) Show linear regression model" << std::endl;
        std::cout << "(4) Show residuals of linear regression model" << std::endl;
        std::cout << "(5) Exit" << std::endl;
        std::cin >> command;
        std::cout << std::endl;
        if (command == "1") {
            bool show = true;
            if (k > 4) {
                std::cout << "Sign table is too big. (" << NUM_ROWS << " rows, " <<
                1 + NUM_ROWS << " columns)" << std::endl;
                std::cout << "Are you sure you want to show it? (y/n)" << std::endl;
                do {
                    std::cin >> command;
                } while (command != "y" && command != "n");
                std::cout << std::endl;
                if (command == "n") show = false;
            }
            if (show) showSignTable();
        }
        else if (command == "2") showImpacts();
        else if (command == "3") showLinearRegressionModel();
        else if (command == "4") showResiduals();
    } while (command != "5");

    std::cout << "Factorial design of experiments ended." << std::endl << std::endl;
}
