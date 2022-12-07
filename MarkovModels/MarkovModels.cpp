
#include "MarokovChain.h"
#include "TransitionCalculations.h"

int main(int argc, char* argv[])
{
    MarokovChain chain;
    MarokovChain chain2;
    std::vector<std::string> paths1;
    std::vector<std::string> paths2;
    // // std::string path;
    // // std::cout <<"Please enter the path to the ORF data\n";
    // // std::cin >> path;
    // // paths.emplace_back(path);
    // // std::cout <<"Please enter the path to the NORF data\n";
    // // std::cin >> path;
    // paths.emplace_back(path);
    paths1.emplace_back("TrainingData/ORF.txt");
    paths2.emplace_back("TrainingData/NORF.txt");
    chain.ReadFiles(paths1);
    chain2.ReadFiles(paths2);
    chain.ConstructMarkovChain();
    chain2.ConstructMarkovChain();
    TransitionCalculations calculations;
    calculations.GetTransitionTable(chain.GetTransitionTable(), chain2.GetTransitionTable());
    chain.GenerateScores(calculations);
    chain2.GenerateScores(calculations);
    return 0;
}
