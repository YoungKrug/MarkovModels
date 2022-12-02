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
        //Find the window of exons
        //WE can use Regex to find windows
        std::string seq = sequence.second.Sequence;
        std::regex startCodons("ATG");
        std::smatch match_results;
        if(std::regex_search(seq,match_results, startCodons))
        {
            //std::cout <<match_results.size();
            for(int i = static_cast<int>(match_results.position()); i < static_cast<int>(seq.size()); i+=3)
            {
               
                std::string codon = seq.substr(i, 3);
                if(codon.size() < 3)
                    continue;
                ChooseListToAddToo(sequence.second.isORF, codon);
                if(codon._Equal("TGA") ||
              codon._Equal("TAA") ||
              codon._Equal("TAG"))
                    break; // Break the window, no more reading
            }
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
    for (auto codons: _codonAmountsORF)
    {
        probOfCodonInORF += codons.second;
    }
    for (auto codons: _codonAmountsNORF)
    {
        probOfCodonInNORF += codons.second;
    }
    std::cout << "Probability in ORF\n";
    for (auto codons: _codonAmountsORF)
    {
        std::cout << codons.first << " Has a probability of: " << static_cast<double>(codons.second) / probOfCodonInORF << " Amount: " << codons.second << std::endl;
    }
    std::cout << "\n\nProbability in NORF\n";
    for (auto codons: _codonAmountsNORF)
    {
        std::cout << codons.first << " Has a probability of: " << static_cast<double>(codons.second) / probOfCodonInNORF << " Amount: " << codons.second << std::endl;
    }
    
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
}
