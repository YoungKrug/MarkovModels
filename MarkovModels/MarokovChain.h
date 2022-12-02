#pragma once
#include <string>
#include <unordered_map>

#include "Matrix.h"
#include "SequenceData.h"

//How often it goes from one codon to another codon
//Create a bias table..
class MarokovChain
{
public:
    void ReadFiles(std::vector<std::string>);
    void ConstructMarkovChain();
    void ChooseListToAddToo(bool, std::string);
    void DisplayData() const;
    void CreateMatrix();
    MarokovChain(){}
private:
     enum States
     {
         ORF,
         NORF
     };
    States _state;
    std::unordered_map<std::string, SequenceData> _sequenceInformation;
    Matrix::Matrix<std::string> _transitionMatrix;
    std::unordered_map<std::string, int> _codonAmountsORF;
    std::unordered_map<std::string, int> _codonAmountsNORF;
    std::unordered_map<std::string, int> _codonToNumber;
    std::unordered_map<int, std::string> _numberToCodon;
    std::unordered_map<std::string,int> _codonToCodonAmounts;
    int _numberOfCodons = 0;
    std::string mapDelimiter = "->";
};
