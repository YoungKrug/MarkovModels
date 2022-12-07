#pragma once
#include <string>
#include <unordered_map>

#include "Matrix.h"

class TransitionCalculations
{
public:
    Matrix::Matrix<std::string>GetTransitionTable(Matrix::Matrix<std::string>, Matrix::Matrix<std::string>);
    std::unordered_map<std::string, double> _transverseValues;
};
