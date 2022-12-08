#include "MarokovChain.h"

#include <map>
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
        if(sequence.second.sequenceType == SequenceData::CounterClockWise)
            seq = ReverseTranscription(seq);
        //std::cout <<match_results.size();
        int length = static_cast<int>(seq.size()) / 3;
        for(int i = 0; i < length; i++)
        {
            if(seq[i*3 + 3] == '\0')
            {
                continue;
            }
            std::string codon = seq.substr(i * 3, 3);
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
}
void MarokovChain::DisplayData() const
{
    std::string path;
    std::cout << "Please enter a output file for your distribution table.\n";
    std::cin >> path;
    std::ofstream output;
    output.open(path);
    if(!output.is_open())
    {
        std::cout << "Invalid output file, please enter a new one";
        std::cin >> path;
        output.open(path);
    }
    for(auto val : _sequenceScores)
    {
        std::string seq = val.first;
        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.cend());
        output << seq<< "," << val.second << std::endl;
    }
}

void MarokovChain::CreateMatrix()
{
    double probOfCodonInORF = 0.0;
    double probOfCodonInNORF = 0.0;
    int length = static_cast<int>(_codonToNumber.size());
    Matrix::Matrix<std::string> initMatrix(length + 1, length + 1);
    _transitionMatrix = initMatrix;
    for(int j = 1; j < length + 1; j++)
    {
        _transitionMatrix.matrix[j][0] = _numberToCodon[j - 1];
    }
    for(int i = 1; i < length + 1; i++)
    {
        _transitionMatrix.matrix[0][i] = _numberToCodon[i - 1];
    }
    for(int j = 1; j < length + 1; j++)
    {
        bool restart = true;
        double totalAmount = 0;
        for(int i = 1; i < length + 1; i++)
        {
            std::string key = std::to_string(j - 1) + mapDelimiter + std::to_string(i - 1);
            if(restart)
            {
                int val = _codonToCodonAmounts[key];
                _transitionMatrix.matrix[j][i] = std::to_string(val);
                totalAmount +=_codonToCodonAmounts[key];
                if(i >= length && restart)
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
}

void MarokovChain::GenerateScores(const TransitionCalculations scoreMatrix)
{
    for(auto seq : _sequenceInformation)
    {
        std::string sequence = seq.second.Sequence;
        std::string previousCodon = "";
        double score = 0.0;
        int length = static_cast<int>(sequence.size()) / 3;
        
        for(int i = 0; i < length; i++)
        {
            if(sequence[i*3 + 3] == '\0')
            {
                continue;
            }
            std::string codon = sequence.substr(i * 3, 3);
            if(previousCodon == "")
            {
                previousCodon = codon;
                continue;
            }
            std::string codonToCodonKey;
            codonToCodonKey = (previousCodon) + (mapDelimiter) + (codon);
            double probability;
            probability = scoreMatrix._transverseValues.at(codonToCodonKey);
            score += probability;
            previousCodon = codon;
        }
        double doubLength = static_cast<double>(length);
        double finalScore = score / doubLength;
        _sequenceScores.emplace_back(std::pair<std::string, double>(seq.second.name, finalScore));
    }
    DisplayData();
}

std::string MarokovChain::ReverseTranscription(std::string str)
{
    std::string newString ="";
    std::unordered_map<char, char> complements;
    complements.insert({'G', 'C'});
    complements.insert({'C', 'G'});
    complements.insert({'A', 'T'});
    complements.insert({'T', 'A'});
    for(int i = static_cast<int>(str.size()); i > 0; i--)
    {
        if(str[i] == '\0')
            continue;
        char val = complements[str[i]];
        newString += val;
        
    }
    return newString;
}

Matrix::Matrix<std::string> MarokovChain::GetTransitionTable()
{
    return _transitionMatrix;
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
