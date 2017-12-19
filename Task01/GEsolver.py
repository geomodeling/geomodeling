import numpy as np

##########################################
# Gaussian elimination solver
##########################################
def GEsolver(A, b):
    n = len(A)
    # Forward elimination
    for index in range(n-1):
        for row in range(index+1, n):
            multiplier = A[row,index]/A[index,index]
            A[row, index:] = A[row, index:] - multiplier*A[index, index:]
            b[row] = b[row] - multiplier*b[index]
    
    # Back Substitution
    x = np.zeros(n)
    for index in range(n-1, -1, -1):
        x[index] = (b[index] - np.dot(A[index,index+1:],x[index+1:]))/A[index,index]
    return x

##########################################
# Test data
##########################################
A = np.array([
        [1,-1,1,-1],
        [1,0,0,0],
        [1,1,1,1],
        [1,2,4,8]
])
b = np.array([14,4,2,2])

##########################################
# Show result
##########################################
print("result from lib function:")
print(np.linalg.solve(A, b))
print("result from custom function:")
print(GEsolver(A,b))
