import numpy as np
from LUP_decomposition import LUP
# from LUP_decomposition2 import LUP2 # The modified LUP function if you have it.

def computeDet(A):
    n = 
    L, U, P, success = LUP(A)
    if ... :
        print("This matrix appears to be singular!")
        det = 0.0
    else:
        det = 1.0
        for i in range(n):
            ...
    return det

# Small test case:
n = 8
A = np.random.rand(n,n)
det1 = computeDet(A)
det2 = np.linalg.det(A)
print("Det computed by me is %e and by NumPy is %e - should match up to the sign!" % (det1,det2))

def computeDet2(A):
    n = ...
    L, U, P, success, parity = LUP2(A)
    if ... :
        print("This matrix appears to be singular!")
        det = 0.0
    else:
        det = ...
        for i in range(n):
            ...
    return det

det1 = computeDet2(A)
print("Det computed by me is %e and by NumPy is %e - should match!" % (det1,det2))

