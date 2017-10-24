#Gaussian elimination solver
def GEsolver(A, b):
    np =  len(A)
    for index in range(np-1):
        #Eliminate
        for row in range(index+1, np):
            multiplier = A[row,index]/A[index,index]
            A[row, index:] = A[row, index:] - multiplier*A[index, index:]
            b[row] = b[row] - multiplier*b[index]
    
    # Back Substitution
    x = np.zeros(np)
    for index in range(np-1, -1, -1):
        x[index] = (b[index] - np.dot(A[index,index+1:],x[index+1:]))/A[index,index]
    return x

A = np.array([[1.,-1.,1.,-1.],[1.,0.,0.,0.],[1.,1.,1.,1.],[1.,2.,4.,8.]])
b =  np.array([[14.],[4.],[2.],[2.]])
print GEsolver(A,b)
