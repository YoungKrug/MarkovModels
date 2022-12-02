#pragma once
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
    };
}
