#pragma once
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

namespace Matrix
{
    template <class T>
    class Matrix
    {
    public:
        std::vector<std::vector<T>> matrix;
        Matrix<T>(int ySize, int xSize)
        {
            std::vector <std::vector<T>> sMatrix(ySize, std::vector<T>(xSize));
            matrix = sMatrix;
        }
        Matrix<T>(){}
        void DisplayMatrix()
        {
            matrix[0][0] = " ";
          for(int j = 0; j < matrix.size(); j++)
          {
              for(int i = 0; i < matrix[0].size(); i++)
              {
                  std::cout << matrix[i][j] << " ";
              }
              std::cout << std::endl;
          }
            
        }
        void PrintToOutputFile(std::string path)
        {
            std::ofstream output;
            output.open(path);
            if(!output.is_open())
            {
                std::cout << "Invalid output file, please enter a new one";
                std::cin >> path;
                output.open(path);
            }
            for(int j = 0; j < matrix.size(); j++)
            {
                for(int i = 0; i < matrix[0].size(); i++)
                {
                    output << matrix[i][j] << "\t";
                }
                output << std::endl;
            }
            
        }
    };
}
