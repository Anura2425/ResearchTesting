# GOAL: find the recursive relation for generating the inverse of any n x n square matrices
import numpy as np

# Manually assign each value of the matrix
# matrix = np.array(np.asmatrix('1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25')) # interesting pattern that when the matrix has a linear progression it has no inverse (det is always 0?)

# Random matrix (change dimensions in size variable)
matrix = np.random.randint(10, size=(4, 4))
print("Initial Matrix (A):")
print(matrix)



def inverse(matrix):
    '''
    Title:          inverse()
    Params:         matrix - 2D array
    Output:         inverse - 2D array
    Description:    Finds the inverse of an input array
    Requirements:   Must be a square matrix, determinant cant = 0
    '''
    det = determinant(matrix)
    if det == 0:
        print("Determinant is 0, no valid inverse")
        return None
    
    adj = adjugate(matrix)
    inv = adj/det
    return inv


def determinant(matrix):
    '''
    Title:          determinant()
    Description:    Recursively computes determinant using Laplace expansion
    '''
    n = len(matrix)
    if n == 1:
        return matrix[0][0]

    if n == 2:
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0] # ad-bc
    
    # Laplace expansion (note: extremely inefficient with an n of higher value, on the order of O(n!))
    det = 0
    for j in range(n):
        submatrix = np.delete(np.delete(matrix, 0, axis=0), j, axis=1) # splice out the items in the same row and column as the cofactor
        cofactor = ((-1)**j) * matrix[0][j] * determinant(submatrix)
        det+=cofactor
        
    return det


def adjugate(matrix):
    '''
    Title:          adjugate()
    Description:    Computes adj(A) = transpose of cofactor matrix
    '''
    n = len(matrix)
    cofactor_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            sub = np.delete(np.delete(matrix, i, axis=0), j, axis = 1)
            cofactor_matrix[i][j] = ((-1) ** (i + j)) * determinant(sub)
            
    return cofactor_matrix.T # transpose of the cofactor matrix 

inv = inverse(matrix)
print("\nInverse of Matrix (A**-1)")
print(inv)

# Check correctness using A * A^{-1} = I
if inv is not None:
    print("\nCheck A * A^{-1} = I:")
    identity_check = matrix @ inv
    print(np.round(identity_check).astype(int)) # @ matrix multiplication op