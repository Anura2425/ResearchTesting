#include <iostream>
#include <string>
#include <math.h>
#include "../MatrixClass/matrix2d.h"

using namespace std;

template <class T>
void PrintMatrix(matrix2D<T> matrix){
    int nRows = matrix.getNumRows();
    int nCols = matrix.getNumCols();
    for (int row = 0; row<nRows; ++row){
        for (int col = 0; col<nCols; ++col){
            cout << matrix.getElement(row, col) << " ";
        }
        cout << endl;
    }
}

int main(){


    double xs[12] {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

    matrix2D<double> testMatrix(3, 4, xs);

    PrintMatrix(testMatrix);

    return 0;
}