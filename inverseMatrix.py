# GOAL: find the recursive relation for generating the inverse of any n x n square matrices


'''
Title:          Inverse
Params:         matrix - 2D array
Output:         inverse - 2D array
Description:    Finds the inverse of an input array
Requirements:   Must be a square matrix, 

'''

import numpy as np

matrix = np.random.randint(10, size=(3, 3))
print("Initial Matrix (A):")
print(matrix)

def inverse(matrix):
    if(determinant == 0):
        print("Determinant is 0, no valid inverse")
        return
    determinant(matrix)
    adjugate(matrix)

def determinant(matrix):
    if (matrix):

    pass

def adjugate(matrix):
    pass