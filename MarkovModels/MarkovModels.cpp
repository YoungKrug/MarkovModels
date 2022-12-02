
#include "MarokovChain.h"

int main(int argc, char* argv[])
{
    MarokovChain chain;
    std::vector<std::string> paths;
    // // std::string path;
    // // std::cout <<"Please enter the path to the ORF data\n";
    // // std::cin >> path;
    // // paths.emplace_back(path);
    // // std::cout <<"Please enter the path to the NORF data\n";
    // // std::cin >> path;
    // paths.emplace_back(path);
    paths.emplace_back("TrainingData/ORF.txt");
    paths.emplace_back("TrainingData/NORF.txt");
    chain.ReadFiles(paths);
    chain.ConstructMarkovChain();
    return 0;
}
