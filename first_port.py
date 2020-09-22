"""
Main module

This small program can calculate the solution system of position functions
for a simple system of two spring coupled masses. Later version will take
bigger than 2x2 matrices for the mass and spring constants.

Reduced row echelon form function (rref) adapted from Wikipedia pseudocode
"""

from math import sqrt, cos, sin, pow
from numpy import arange, array, csingle, divide, identity, linalg, matmul, multiply, subtract
# import matplotlib?
# TODO: refactor variable names to conform to PEP 8


def rref(matrix: array):
    lead = 0
    row_count, column_count = matrix.shape
    for r in range(row_count):
        if column_count <= lead:
            return matrix
        i = r
        while matrix[i, lead] == 0:  # search for leading integer
            i += 1
            if row_count == i:
                i = r
                lead += 1
                if column_count == lead:
                    return matrix
        if i != r:
            matrix[[i, r]] = matrix[[r, i]]  # swap rows i and r
        matrix[r] = divide(matrix[r], matrix[r, lead])
        for y in range(row_count):
            if y != r:
                matrix[y] = subtract(matrix[y], multiply(matrix[r], matrix[y, lead]))
        lead += 1


# Variable definitions
mass_one = 10.0  # later user defined mass
mass_two = 20.0  # masses
mass_matrix = array([[mass_one, 0], [0, mass_two]])  # Mass matrix, mass_one, mass_two 2x2
spring_one = 100.0  # spring coefficients?
spring_two = 150.0
spring_three = 125.0
# obtained by drawing diagram of system
spring_matrix = array([[-(spring_one + spring_two), spring_two], [spring_two, -(spring_two + spring_three)]])
A = matmul(linalg.inv(mass_matrix), spring_matrix)
# inverted mass_matrix = [[.1,0],[0,.5]]
# Find eigenvalues. Write separate function for this eventually
a = A[0, 0]  # The following 4 will be input args
b = A[0, 1]
c = A[1, 0]
f = A[1, 1]
# Function will return two values (scalars):
e_val_1_A = (((a + f) + sqrt(pow((-a - f), 2) - (4 * (a * f - b * c)))) * .5)  # Eigenvalue 1 should be -7.3691
e_val_2_A = (((a + f) - sqrt(pow((-a - f), 2) - (4 * (a * f - b * c)))) * .5)  # Eigenvalue 2 should be -31.381
om1 = sqrt(abs(e_val_1_A))  # omega om1 = 2.7146
om2 = sqrt(abs(e_val_2_A))  # is this needed? om2 = 5.6019
# v =  # Find eigenvector of e-val 1
time_matrix = arange(0, 11, 1)  # example. Possible user input
B1 = rref(A - (e_val_1_A * identity(2)))  # B1 = [[1, -0.8508],[0,0]]
B2 = rref(A - (e_val_2_A * identity(2)))  # B2 = [[1, 2.3508],[0,0]]
if B1[0, 1] == 0 or B2[0, 1] == 0:
    print('Something is wrong')
# change to endif if running on octave
# ev1 = [B1(1, 2)*B1(2, 2);B1(2, 2)]; %eigenvectors
# ev2 = [B2(1, 2)*B2(2, 2);B2(2, 2)];
ev1 = array([B1[0, 1], 1])  # ev1 = [-0.8508,1]
ev2 = array([B2[0, 1], 1])  # ev2 = [2.3508,1]
c1 = 1
c3 = c1
c2 = (0+1j)
c4 = c2
x = [ev1 * (c1 * cos(om1 * t) + c2 * sin(om1 * t)) + ev2 * (c3 * cos(om2 * t) + c4 * sin(om2 * t)) for t in time_matrix]
# x, ev1, t, and ev2 are arrays/vectors
# c1, om1, c2, c3, om2, and c4 are floats/scalars
# (14) c2, c4 have i in them
print(x)