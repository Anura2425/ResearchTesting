#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

template <class T>
class matrix2D{
    public:
        matrix2D();
        matrix2D(int rows, int cols);
        matrix2D(int rows, int cols, const T *inputData);
        matrix2D(const matrix2D<T>& inputMatrix); 
        matrix2D(int rows, int cols, const std::vector<T> *inputData);

        ~matrix2D();

        // config methods
        bool resize(int numRows, int numCols);
        void setToIdentity();

        // element access methods
        T getElement(int row, int col);
        bool setElement(int row, int col, T elementValue);
        int getNumRows();
        int getNumCols();

        // manipulation methods
        bool inverse();

        // overload operators for matrix performing operations
        bool operator== (const matrix2D<T>& rhs);
        bool Compare (const matrix2D<T>& matrix1, double tolerance);

        template <class U> friend matrix2D<U> operator+ (const matrix2D<U>& lhs, const matrix2D<U>& rhs); // matrix + matrix
        template <class U> friend matrix2D<U> operator+ (const U& lhs, const matrix2D<U>& rhs);           // scalar + matrix
        template <class U> friend matrix2D<U> operator+ (const matrix2D<U>& lhs, const U& rhs);           // matrix + scalar

        template <class U> friend matrix2D<U> operator- (const matrix2D<U>& lhs, const matrix2D<U>& rhs); // matrix - matrix
        template <class U> friend matrix2D<U> operator- (const U& lhs, const matrix2D<U>& rhs);           // scalar - matrix
        template <class U> friend matrix2D<U> operator- (const matrix2D<U>& lhs, const U& rhs);           // matrix - scalar

        template <class U> friend matrix2D<U> operator* (const matrix2D<U>& lhs, const matrix2D<U>& rhs); // matrix * matrix
        template <class U> friend matrix2D<U> operator* (const U& lhs, const matrix2D<U>& rhs);           // scalar * matrix
        template <class U> friend matrix2D<U> operator* (const matrix2D<U>& lhs, const U& rhs);           // matrix * scalar

        bool Seperate(matrix2D<T> *matrix1, matrix2D<T> *matrix2, int colNum);
    //private:
    public: //pub for testing, otherwise should be private
        int sub2Ind(int row, int col);
        bool isSquare();
        bool closeEnough(T f1, T f2);
        void swapRow(int i, int j);
        void multAdd(int i, int j, T multFactor);
        void multRow(int i, T multFactor);
        bool Join(const matrix2D<T>& matrix2);
        int findRowWithMaxElement(int colNumber, int startingRow);
        void printMatrix();

    private: 
        T *m_matrixData;
        int m_rows, m_cols, m_elements;
};

//****************************//
//   CONSTRUCTOR FUNCTIONS    //
//****************************//

// default constructor
template <class T>
matrix2D<T>::matrix2D(){
    m_rows = 1;
    m_cols = 1;
    m_elements = 1;
    m_matrixData = new T[m_elements];
    m_matrixData[0] = 0.0;
}

// construct empty matrix (zeroed out)
template <class T>
matrix2D<T>::matrix2D(int rows, int cols){
    m_rows = rows;
    m_cols = cols;
    m_elements = m_rows * m_cols;
    m_matrixData = new T[m_elements];
    for (int i = 0; i < m_elements; ++i){
        m_matrixData[i] = 0.0;
    }
}

// construct from linear array of data
template <class T>
matrix2D<T>::matrix2D(int rows, int cols, const T *inputData){
    m_rows = rows;
    m_cols = cols;
    m_elements = m_rows * m_cols;
    m_matrixData = new T[m_elements];
    for (int i = 0; i < m_elements; ++i){
        m_matrixData[i] = inputData[i];
    }
}

// copy constructor
template <class T>
matrix2D<T>::matrix2D(const matrix2D<T>& inputMatrix){
    m_rows = inputMatrix.m_rows;
    m_cols = inputMatrix.m_cols;
    m_elements = inputMatrix.m_elements;
    m_matrixData = new T[m_elements];
    for (int i = 0; i < m_elements; ++i){
        m_matrixData[i] = inputMatrix.m_matrixData[i];
    }
}

// Construct from a vector
template <class T>
matrix2D<T>::matrix2D(int rows, int cols, const std::vector<T> *inputData){
    m_rows = rows;
    m_cols = cols;
    m_elements = m_rows * m_cols;
    m_matrixData = new T[m_elements];
    for (int i=0; i<m_elements; ++i){
        m_matrixData[i] = inputData->at(i);
    }
}

// destructor
template <class T>
matrix2D<T>::~matrix2D(){
    if (m_matrixData != nullptr){
        delete[] m_matrixData;
    }
}

//******************************//
//   CONFIGURATION FUNCTIONS    //
//******************************//

template <class T>
bool matrix2D<T>::resize(int numRows, int numCols){
    m_rows = numRows;
    m_cols = numCols;
    m_elements = (m_rows * m_cols);
    delete[] m_matrixData;
    m_matrixData = new T[m_elements];
    if (m_matrixData != nullptr){
        for (int i = 0; i < m_elements; ++i){
            m_matrixData[i] = 0.0;
        }
        return true;
    } else {
        return false;
    }
}

template <class T>
void matrix2D<T>::setToIdentity(){
    if (!isSquare())
        throw std::invalid_argument("Cannot form an indenity matrix from a non-square matrix");
    
    for (int row=0; row<m_rows; ++row){
        for (int col=0; col<m_cols; ++col){
            if (col == row){
                m_matrixData[sub2Ind(row, col)] = 1.0;
            } else {
                m_matrixData[sub2Ind(row, col)] = 0.0;
            }
        }
    }
}

//******************************//
//      ELEMENT FUNCTIONS       //
//******************************//

template <class T>
T matrix2D<T>::getElement(int row, int col){
    int linearIndex = sub2Ind(row, col);
    if (linearIndex >= 0){
        return m_matrixData[linearIndex];
    } else{
        return 0.0;
    }
}

template <class T>
bool matrix2D<T>::setElement(int row, int col, T elementValue){
    int linearIndex = sub2Ind(row, col);
    if (linearIndex >= 0){
        m_matrixData[linearIndex] = elementValue;
        return true; 
    } else {
        return false;
    }
}

template <class T>
int matrix2D<T>::getNumRows(){
    return m_rows;
}

template <class T>
int matrix2D<T>::getNumCols(){
    return m_cols;
}

template <class T>
bool matrix2D<T>::Compare(const matrix2D<T>& matrix1, double tolerance){
    int numRows1 = matrix1.m_rows;
    int numCols1 = matrix1.m_cols;
    if ((numRows1 != m_rows) || (numCols1 != m_cols))
        return false;
    
    double cumulativeSum = 0.0;
    for (int i=0; i<m_elements; ++i){
        T element1 = matrix1.m_matrixData[i];
        T element2 = m_matrixData[i];
        cumulativeSum += ((element1 - element2) * (element1 - element2));
    }
    double finalValue = sqrt(cumulativeSum / ((numRows1 * numCols1)-1));
    if (finalValue < tolerance){
        return true;
    } else {
        return false;
    }
}

//***************************//
//     INVERSE FUNCTION      //
//***************************//

template <class T>
bool matrix2D<T>::inverse(){
    if(!isSquare()){
        throw std::invalid_argument("Invalid: cannot compute the inverse of a non-square matrix");
    }

    matrix2D<T> identityMatrix(m_rows, m_cols);
    identityMatrix.setToIdentity();

    int originalNumCols = m_cols;
    Join(identityMatrix);

    // begin main process
    int cRow, cCol;
    int maxCount = 100;
    int count = 0;
    bool completeFlag = false;
    while ((!completeFlag) && (count < maxCount)){
        for (int diagIndex = 0; diagIndex<m_rows; ++diagIndex){
            // loop over the diagonal of the matrix to check if the elements are 1's
            cRow = diagIndex;
            cCol = diagIndex;

            // find the indext of the max elem in curr column
            int maxIndex = findRowWithMaxElement(cCol, cRow);

            // if this isnt the curr row, then swap
            if (maxIndex!=cRow){
                swapRow(cRow, maxIndex);
            }

            if (m_matrixData[sub2Ind(cRow, cCol)] != 1.0){
                T multFactor = 1.0 / m_matrixData[sub2Ind(cRow, cCol)];
                multRow(cRow, multFactor);
            }

            for (int rowIndex = cRow + 1; rowIndex<m_rows; ++rowIndex){
                if (!closeEnough(m_matrixData[sub2Ind(rowIndex, cCol)], 0.0)) {
                    int rowOneIndex = cCol;

                    T currentElementValue = m_matrixData[sub2Ind(rowIndex, cCol)];

                    T rowOneValue = m_matrixData[sub2Ind(rowOneIndex, cCol)];

                    if (!closeEnough(rowOneValue, 0.0)){
                        T correctionFactor = -(currentElementValue / rowOneValue);
                        multAdd(rowIndex, rowOneIndex, correctionFactor);
                    }
                }
            }
            for (int colIndex = cCol+1; colIndex<originalNumCols; ++colIndex){
                if (!closeEnough(m_matrixData[sub2Ind(cRow, colIndex)], 0.0)){
                    int rowOneIndex = colIndex;
                    
                    T currentElementValue = m_matrixData[sub2Ind(cRow, colIndex)];

                    T rowOneValue = m_matrixData[sub2Ind(rowOneIndex, colIndex)];

                    if (!closeEnough(rowOneValue, 0.0)){
                        T correctionFactor = -(currentElementValue / rowOneValue);

                        multAdd(cRow, rowOneIndex, correctionFactor);
                    }
                }
            }
        }
        matrix2D<T> leftHalf;
        matrix2D<T> rightHalf;
        this->Seperate(&leftHalf, &rightHalf, originalNumCols);

        if (leftHalf == identityMatrix){
            completeFlag = true;

            m_cols = originalNumCols;
            m_elements = m_rows * m_cols;
            delete[] m_matrixData;
            m_matrixData = new T[m_elements];
            for (int i = 0; i<m_elements; ++i){
                m_matrixData[i] = rightHalf.m_matrixData[i];
            }
        }
        count++;
    }
    return completeFlag;
}

//*********************************//
//       OVERLOADED OPERATORS      //
//*********************************//

// matrix + matrix
template <class T>
matrix2D<T> operator+ (const matrix2D<T>& lhs, const matrix2D<T>& rhs){
    int numRows = lhs.m_rows;
    int numCols = lhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i<numElements; ++i){
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar + matrix
template <class T>
matrix2D<T> operator+ (const T& lhs, const matrix2D<T>& rhs){
    int numRows = rhs.m_rows;
    int numCols = rhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = lhs + rhs.m_matrixData[i];
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix + scalar
template <class T>
matrix2D<T> operator+ (const matrix2D<T>& lhs, const T& rhs){
    int numRows = lhs.m_rows;
    int numCols = lhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = lhs.m_matrixData[i] + rhs;
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix - matrix
template <class T>
matrix2D<T> operator- (const matrix2D<T>& lhs, const matrix2D<T>& rhs){
    int numRows = lhs.m_rows;
    int numCols = lhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i<numElements; ++i){
        tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar - matrix
template <class T>
matrix2D<T> operator- (const T& lhs, const matrix2D<T>& rhs){
    int numRows = rhs.m_rows;
    int numCols = rhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = lhs - rhs.m_matrixData[i];
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix - scalar
template <class T>
matrix2D<T> operator- (const matrix2D<T>& lhs, const T& rhs){
    int numRows = lhs.m_rows;
    int numCols = lhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = lhs.m_matrixData[i] - rhs;
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}


// matrix * matrix
template <class T>
matrix2D<T> operator* (const matrix2D<T>& lhs, const matrix2D<T>& rhs){
    int r_numRows = rhs.m_rows;
    int r_numCols = rhs.m_cols;
    int l_numRows = lhs.m_rows;
    int l_numCols = lhs.m_cols;
    
    if (l_numCols == r_numRows){
        // The output is the same size as the RHS in matrix multiplication
        T *tempResult = new T[lhs.m_rows * rhs.m_cols];
        for (int lhsRow=0; lhsRow<l_numRows; lhsRow++){
            for (int rhsCol=0; rhsCol<r_numCols; rhsCol++){
                T elementResult = 0.0;
                for (int lhsCol=0; lhsCol<l_numCols; lhsCol++){
                    int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
                    int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;
                    // perform calculation on the elements at the linear indices
                    elementResult+=(lhs.m_matrixData[lhsLinearIndex] * rhs.m_matrixData[rhsLinearIndex]);
                }
                // store result
                int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
                tempResult[resultLinearIndex] = elementResult;
            }
        }
        matrix2D<T> result(l_numRows, r_numCols, tempResult);
        delete[] tempResult;
        return result;
    } else {
        matrix2D<T> result(1, 1);
        return result;
    }
}

// scalar * matrix
template <class T>
matrix2D<T> operator* (const T& lhs, const matrix2D<T>& rhs){
    int numRows = rhs.m_rows;
    int numCols = rhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = lhs * rhs.m_matrixData[i];
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix * scalar
template <class T>
matrix2D<T> operator* (const matrix2D<T>& lhs, const T& rhs){
    int numRows = lhs.m_rows;
    int numCols = lhs.m_cols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; ++i){
        tempResult[i] = lhs.m_matrixData[i] * rhs;
    }

    matrix2D<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix == matrix
template <class T>
bool matrix2D<T>::operator== (const matrix2D<T>& rhs){
    if ((this->m_rows != rhs.m_rows) && (this->m_cols != rhs.m_cols)){
        return false;
    }

    bool flag = true;
    for (int i = 0; i<this->m_elements; ++i){
        if (!closeEnough(this->m_matrixData[i], rhs.m_matrixData[i])){
            flag = false;
        }
    }
    return flag;
}

template <class T>
bool matrix2D<T>::Seperate(matrix2D<T> *matrix1, matrix2D<T> *matrix2, int colNum){
    int numRows = m_rows;
    int numCols1 = colNum;
    int numCols2 = m_cols - colNum;

    // resize matrices
    matrix1->resize(numRows, numCols1);
    matrix2->resize(numRows, numCols1);

    for (int row = 0; row < m_rows; ++row){
        for (int col = 0; col < m_cols; ++col){
            if (col < colNum){
                matrix1->setElement(row, col, this->getElement(row,col));
            } else {
                matrix2->setElement(row, col-colNum, this->getElement(row, col));
            }
        }
    }
    return true;
}

template <class T>
bool matrix2D<T>::Join(const matrix2D<T>& matrix2){
    int numRows1 = m_rows;
    int numRows2 = matrix2.m_rows;
    int numCols1 = m_cols;
    int numCols2 = matrix2.m_cols;

    if (numRows1 != numRows2){
        throw std::invalid_argument("Invalid: cannot join matrices with different numbers of rows");
    }

    // allocate mem
    T* newMatrixData = new T[numRows1*(numCols1+numCols2)];

    // copy the matrices into the new one
    int linearIndex, resultLinearIndex;
    for (int i=0; i<numRows1; ++i){
        for (int j=0; j<(numCols1+numCols2); ++j){
            resultLinearIndex = (i*(numCols1+numCols2)) + j;
            
            if (j<numCols1){
                linearIndex = (i * numCols1) + j;
                newMatrixData[resultLinearIndex] = m_matrixData[linearIndex];
            } else {
                linearIndex = (i * numCols2) + (j - numCols1);
                newMatrixData[resultLinearIndex] = matrix2.m_matrixData[linearIndex];
            } 
        }
    }

    m_cols = numCols1+numCols2;
    m_elements = m_rows * m_cols;
    delete[] m_matrixData;
    m_matrixData = new T[m_elements];
    for (int i=0; i<m_elements; ++i){
        m_matrixData[i] = newMatrixData[i];
    }

    delete[] newMatrixData;
    return true;
}

//******************************//
//       PRIVATE FUNCTIONS      //
//******************************//

template <class T>
int matrix2D<T>::sub2Ind(int row, int col){
    if ((row < m_rows) && (row >= 0) && (col < m_cols) && (col >= 0)){
        return (row * m_cols) + col;
    } else {
        return -1;
    }
}

template <class T>
bool matrix2D<T>::isSquare(){
    if (m_cols == m_rows){
        return true;
    } else {
        return false;
    }
}

template <class T>
void matrix2D<T>::swapRow(int i, int j){
    T *tempRow = new T[m_cols];
    for (int k = 0; k < m_cols; ++k){
        tempRow[k] = m_matrixData[sub2Ind(i,k)];
    }

    for (int k = 0; k < m_cols; ++k){
        m_matrixData[sub2Ind(i,k)] = m_matrixData[sub2Ind(j,k)];
    }

    for (int k = 0; k < m_cols; ++k){
        m_matrixData[sub2Ind(j,k)] = tempRow[k];
    }

    delete[] tempRow;
}

template <class T>
void matrix2D<T>::multAdd(int i, int j, T multFactor){
    for (int k=0; k<m_cols; ++k){
        m_matrixData[sub2Ind(i,k)] += (m_matrixData[sub2Ind(j,k)] * multFactor);
    }
}

template <class T>
int matrix2D<T>::findRowWithMaxElement(int colNumber, int startingRow){
    T tempValue = m_matrixData[sub2Ind(startingRow, colNumber)];
    int rowIndex = startingRow;
    for (int k=startingRow+1; k<m_rows; ++k){
        if (fabs(m_matrixData[sub2Ind(k, colNumber)]) > fabs(tempValue)){
            rowIndex = k;
            tempValue = m_matrixData[sub2Ind(k, colNumber)];
        }
    }
    return rowIndex;
}

template <class T>
void matrix2D<T>::multRow(int i, T multFactor){
    for (int k=0; k<m_cols; ++k){
        m_matrixData[sub2Ind(i,k)] *= multFactor;
    }
}

template <class T>
void matrix2D<T>::printMatrix(){
    int nRows = this->getNumRows();
    int nCols = this->getNumCols();
    for (int row = 0; row<nRows; ++row){
        for (int col = 0; col<nCols; ++col){
            std::cout << std::fixed << std::setprecision(3) << this->getElement(row, col) << "  ";
        }
        std::cout << std::endl;
    }
}

template <class T>
bool matrix2D<T>::closeEnough(T f1, T f2){
    return fabs(f1 - f2) < 1e-9;
}

#endif