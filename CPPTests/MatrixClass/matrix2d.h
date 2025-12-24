#ifndef MATRIX2D_H
#define MATRIX2D_H

template <class T>
class matrix2D{
    public:
        matrix2D();
        matrix2D(int rows, int cols);
        matrix2D(int rows, int cols, const T *inputData);
        matrix2D(const matrix2D<T>& inputMatrix); 

        ~matrix2D();


        // config methods
        bool resize(int numRows, int numCols);

        // element access methods
        T getElement(int row, int col);
        bool setElement(int row, int col, T elementValue);
        int getNumRows();
        int getNumCols();

        // overload operators for matrix performing operations
        bool operator== (const matrix2D<T>& rhs);

        template <class U> friend matrix2D<U> operator+ (const matrix2D<U>& lhs, const matrix2D<U>& rhs); // matrix + matrix
        template <class U> friend matrix2D<U> operator+ (const U& lhs, const matrix2D<U>& rhs);           // scalar + matrix
        template <class U> friend matrix2D<U> operator+ (const matrix2D<U>& lhs, const U& rhs);           // matrix + scalar

        template <class U> friend matrix2D<U> operator- (const matrix2D<U>& lhs, const matrix2D<U>& rhs); // matrix - matrix
        template <class U> friend matrix2D<U> operator- (const U& lhs, const matrix2D<U>& rhs);           // scalar - matrix
        template <class U> friend matrix2D<U> operator- (const matrix2D<U>& lhs, const U& rhs);           // matrix - scalar

        template <class U> friend matrix2D<U> operator* (const matrix2D<U>& lhs, const matrix2D<U>& rhs); // matrix * matrix
        template <class U> friend matrix2D<U> operator* (const U& lhs, const matrix2D<U>& rhs);           // scalar * matrix
        template <class U> friend matrix2D<U> operator* (const matrix2D<U>& lhs, const U& rhs);           // matrix * scalar

    private:
        int sub2Ind(int row, int col);

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
        if (this->m_matrixData[i] != rhs.m_matrixData[i]){
            flag = false;
        }
    }
    return flag;
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

#endif