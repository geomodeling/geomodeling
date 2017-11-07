import numpy as np

def TDMAsolver(A, res):
    a = A.diagonal(-1)
    b = A.diagonal(0)
    c = A.diagonal(1)

    nf = len(res) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, res)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

##########################################
# Test data
##########################################
A = np.array([
    [10., 2, 0, 0], 
    [3, 10, 4, 0], 
    [0, 1, 7, 5], 
    [0, 0, 3, 4]
])
b = np.array([3.,4,5,6])
##########################################
print("result from lib function:")
print(np.linalg.solve(A, b))
print("result from custom function:")
print(TDMAsolver(A, b))
##########################################
