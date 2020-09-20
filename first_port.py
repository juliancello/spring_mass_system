"""
Main module

This small program can calculate the solution system of position functions
for a simple system of two spring coupled masses. Later version will take
bigger than 2x2 matrices for the mass and spring constants.

Reduced row echelon function (rref) adapted from Wikipedia pseudocode
"""

from math import sqrt, cos, sin, pow
from numpy import arange, array, divide, identity, linalg, multiply, subtract
# import matplotlib?
# TODO: refactor variable names to conform to PEP 8

def rref(matrix: array):
    lead = 0
    row_count, column_count = matrix.shape
    for r in range(row_count):
        if column_count <= lead:
            return matrix
        i = r
        while matrix[i, lead] == 0: # search for leading integer
            i += 1
            if row_count == i:
                i = r
                lead += 1
                if column_count == lead:
                    return matrix
        if i != r:
            matrix[[i, r]] = matrix[[r, i]] #swap rows i and r
        matrix[r] = divide(matrix[r], matrix[r,lead])
        for j in range(row_count):
            if j != r:
                matrix[j] = subtract(matrix[j], multiply(matrix[r], matrix[j, lead]))
        lead += 1


#Variable definitions
m1=10.0 #later user defined
m2=20.0 # masses
m=array([[m1,0],[0,m2]]) #Mass matrix, m1, m2 2x2
k1=100.0 # spring coefficients?
k2=150.0
k3=125.0
k=array([[-(k1+k2),k2],[k2,-(k2+k3)]]) #obtained by drawing diagram of system
iden=identity(2)
#n=[m,iden]
#rref(n)
A=linalg.inv(m)*k

#Find eigenvalues. Write separate function for this eventually
a=A[0,0] #The following 4 will be input args
b=A[0,1]
c=A[1,0]
f=A[1,1]
#Function will return two values:
e_val_1_A=(((a+f)+sqrt(pow((-a-f),2)-(4*(a*f-b*c))))*.5) #Eigenvalue 1
e_val_2_A=(((a+f)-sqrt(pow((-a-f),2)-(4*(a*f-b*c))))*.5) #Eigenvalue 2
om1=sqrt(abs(e_val_1_A)) #omega
om2=sqrt(abs(e_val_2_A)) #is this needed?
#v= #Find eigenvector of e-val 1
t=arange(0,10,1) #example. Possible user input
B1=rref(A-(e_val_1_A*iden))
B2=rref(A-(e_val_2_A*iden))
if B1(1,2)==0 or B2(1,2)==0:
  print('Something is wrong')
#change to endif if running on octave
#ev1=[B1(1,2)*B1(2,2)B1(2,2)] #eigenvectors
#ev2=[B2(1,2)*B2(2,2)B2(2,2)]
ev1=[B1(1,2)1]
ev2=[B2(1,2)1]
c1=1
c2=c1
c3=c1
c4=c1 #I do not know about multivariable calc and only a bit about partial derivatives
x=ev1*(c1*cos(om1*t)+c2*sin(om1*t))+ ev2*(c3*cos(om2*t)+c4*sin(om2*t))  #(14) c2, c4 have i in them

