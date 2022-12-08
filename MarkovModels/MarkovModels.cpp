
#include "MarokovChain.h"
#include "TransitionCalculations.h"

int main(int argc, char* argv[])
{
    MarokovChain chain("ORF");
    MarokovChain chain2("NORF");
    std::vector<std::string> paths1;
    std::vector<std::string> paths2;
    std::string path;
    std::cout <<"Please enter the path to the ORF data\n";
    std::cin >> path;
    paths1.emplace_back(path);
    std::cout <<"Please enter the path to the NORF data\n";
    std::cin >> path;
    paths2.emplace_back(path);
    chain.ReadFiles(paths1);
    chain2.ReadFiles(paths2);
    chain.ConstructMarkovChain();
    chain2.ConstructMarkovChain();
    TransitionCalculations calculations;
    calculations.GetTransitionTable(chain.GetTransitionTable(), chain2.GetTransitionTable());
    chain.GenerateScores(calculations);
    chain2.GenerateScores(calculations);
    chain.DisplayData(false);
    chain2.DisplayData(false);
    std::vector<std::string> userPaths;
    std::cout << "Enter a fast A file to predict ORF and NORF sequences:\n ";
    std::cin >> path;
    userPaths.emplace_back(path);
    MarokovChain userChain("User");
    userChain.ReadFiles(userPaths);
    userChain.ConstructMarkovChain();
    userChain.GenerateScores(calculations);
    userChain.DisplayData(true);
    return 0;
}
