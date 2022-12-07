#include "TransitionCalculations.h"

#include <unordered_map>

Matrix::Matrix<std::string> TransitionCalculations::GetTransitionTable(Matrix::Matrix<std::string> firstMarkovChain,
                                                                       Matrix::Matrix<std::string> otherMarkovChain)
{
    //Tables are the same size so
    Matrix::Matrix<std::string> newTable = firstMarkovChain;
    std::string delim = "->";
    for(int i = 1; i < newTable.matrix.size(); i++)
    {
        for(int j = 1; j < newTable.matrix[0].size(); j++)
        {
         
            double val1 = std::stod(firstMarkovChain.matrix[i][j]);
            double val2 = std::stod(otherMarkovChain.matrix[i][j]);
            std::string codonKey = "";
            codonKey = (newTable.matrix[i][0]) + (delim) + (newTable.matrix[0][j]);
            double constant = 0.0;
            if(val1 == 0.0 && val2 == 0.0)
            {
                //if NORF is 0, we its undefined isnt it?
                // if ORF is 0 its a domain issue

                val1 += 0.00001;
                val2 += 0.0001;
            }
            else if (val1 == 0.0 || val2 == 0.0)
            {
                if(val1 != 0.0)
                {
                    constant = val1/2.0;
                    val2 += constant;
                }
                else
                {
                    constant = val2/2.0;
                    val1 += constant;
                }
                
            }
            double newVale = std::log2(val1/val2);
            newTable.matrix[i][j] = std::to_string(newVale);
            _transverseValues.insert({codonKey, newVale});
        }
    }
    return newTable;
}
