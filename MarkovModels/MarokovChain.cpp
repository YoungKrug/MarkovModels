#include "MarokovChain.h"

#include <regex>

#include "fstream"

void MarokovChain::ReadFiles(std::vector<std::string> paths)
{
   
    std::string line;
    std::regex clockWiseOrNot("Counterclockwise"); // matches delimiters or consecutive non-delimiters
    int i = 0;
    for (auto path : paths)
    {
        std::ifstream files;
        files.open(path);
        while (!files.is_open())
        {
            std::cout << "Could not find file, Re enter file" << std::endl;
            std::cin >> path;
            files.open(path);
        }
        std::string seq;
        std::string name;
        while (std::getline(files, line))
        {
            if (line.find('>') != std::string::npos)
            {
                if(!seq.empty())
                {
                    std::smatch sm;
                    SequenceData temp;
                    if(i<= 0)
                        temp.isORF = true;
                    temp.sequenceType = SequenceData::ClockWise;
                    temp.name = name;
                    temp.Sequence = seq;
                    if(std::regex_search(name, sm, clockWiseOrNot))
                        temp.sequenceType = SequenceData::CounterClockWise;
                    _sequenceInformation.insert({name, temp});
                }
                name = line;
                seq.erase(seq.begin(), seq.end());
                continue;
            }
            seq += line;
        }
        
            
        SequenceData data;
        if(i<= 0)
            data.isORF = true;
        data.name = name;
        data.Sequence = seq;
        data.sequenceType = SequenceData::ClockWise;
        std::smatch sm;
        data.Sequence = seq;
        if(std::regex_search(name, sm, clockWiseOrNot))
            data.sequenceType = SequenceData::CounterClockWise;
        _sequenceInformation.insert({name, data});
        i++;
    }
}

void MarokovChain::ConstructMarkovChain()
{
    for(auto sequence: _sequenceInformation)
    {
        std::string seq = sequence.second.Sequence;
        std::string previousCodon = "";
        //std::cout <<match_results.size();
        for(int i = 0; i < static_cast<int>(seq.size()); i+=3)
        {
            std::string codon = seq.substr(i, 3);
            if(codon.size() < 3)
                continue;
            ChooseListToAddToo(sequence.second.isORF, codon);
            if(previousCodon == "")
            {
                previousCodon = codon;
                continue;
            }
            std::string codonToCodonKey;
            codonToCodonKey.append(std::to_string(_codonToNumber[previousCodon]));
            codonToCodonKey.append(mapDelimiter).append(std::to_string(_codonToNumber[codon]));
            if(_codonToCodonAmounts.find(codonToCodonKey) == _codonToCodonAmounts.end())
            {
                _codonToCodonAmounts.insert({codonToCodonKey, 1});
            }
            else
                _codonToCodonAmounts[codonToCodonKey]++;
            previousCodon = codon;
        }
    }
    CreateMatrix();
   // DisplayData();
}
void MarokovChain::DisplayData() const
{
    std::cout << "ORF Table" << std::endl;
    for(auto seqInfo : _codonAmountsORF)
    {
        std::cout << seqInfo.first  << " : " << seqInfo.second << std::endl;
    }
    std::cout << "\nNORF Table" << std::endl;
    for(auto seqInfo : _codonAmountsNORF)
    {
        std::cout << seqInfo.first  << " : " << seqInfo.second << std::endl;
    }
}

void MarokovChain::CreateMatrix()
{
    double probOfCodonInORF = 0.0;
    double probOfCodonInNORF = 0.0;
    int length = static_cast<int>(_codonToNumber.size());
    Matrix::Matrix<std::string> initMatrix(length + 1, length + 1);
    _transitionMatrix = initMatrix;
    for(int j = 1; j < length; j++)
    {
        _transitionMatrix.matrix[j][0] = _numberToCodon[j - 1];
    }
    for(int i = 1; i < length; i++)
    {
        _transitionMatrix.matrix[0][i] = _numberToCodon[i - 1];
    }
    for(int j = 1; j < length; j++)
    {
        bool restart = true;
        double totalAmount = 0;
        for(int i = 1; i < length; i++)
        {
            std::string key = std::to_string(j) + mapDelimiter + std::to_string(i);
            if(restart)
            {
                _transitionMatrix.matrix[j][i] = std::to_string(_codonToCodonAmounts[key]);
                totalAmount +=_codonToCodonAmounts[key];
                if(i + 1 >= length && restart)
                {
                    i = 0;
                    restart = false;
                }
            }
            else
            {
                double num = std::stoi(_transitionMatrix.matrix[j][i]);
                _transitionMatrix.matrix[j][i] = std::to_string(num / totalAmount);
            }
            
        }
    }
    _transitionMatrix.PrintToOutputFile("TrainingData/OutputFile.txt");
}

void MarokovChain::ChooseListToAddToo(bool val, const std::string codon)
{
    if(val)
    {
        if(_codonAmountsORF.find(codon) == _codonAmountsORF.end())
        {
            _codonAmountsORF.insert({codon, 1});
        }
        else
            _codonAmountsORF[codon]++;
    }
    else
    {
        if(_codonAmountsNORF.find(codon) == _codonAmountsNORF.end())
        {
            _codonAmountsNORF.insert({codon, 1});
        }
        else
            _codonAmountsNORF[codon]++;
    }
    if(_codonToNumber.find(codon) == _codonToNumber.end())
    {
        _codonToNumber.insert({codon, _numberOfCodons});
        _numberToCodon.insert({_numberOfCodons, codon});
        _numberOfCodons++;
    }
}
