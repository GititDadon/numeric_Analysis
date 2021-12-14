import numpy as np
from termcolor import colored
# Name: Gitit Dadon ID:212280911
def find_det(A):
    """"" Returns The Determinant Of A Matrix"""
    S=np.array(A)
    return np.linalg.det(S)

def copyM(A, size):
    newMat = []
    for i in range(size):
        rowList = []
        for j in range(size):
            rowList.append(A[i][j])
        newMat.append(rowList)
    return newMat

def Is_Dominant(A, n):
    mat = copyM(A, n)
    for i in range(0, n):
        sum = 0
        for j in range(0, n):
            sum = sum + abs(mat[i][j])
        sum = sum - abs(mat[i][i])
        if abs(A[i][i]) < sum:
            return False
    return True

def initialize_D(A, n):
    newmat = copyM(A, n)
    for i in range(n):
        val = abs(newmat[i][i])
        for j in range(n):
            (newmat[i][j]) = 0
        (newmat[i][i]) = val
    return newmat
def calcM1(A, B):
    """"" Returns The Result Matrix Of 2 Matrix Calculation"""

    matrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                matrix[i][j] += A[i][k] * B[k][j]
    return matrix
def calcM2(A, b):
    """"" Calculates 3x1 matrix"""
    newmat = [[0, ], [0, ], [0, ]]
    for i in range(len(A)):
        for j in range(len(b[0])):
            for k in range(len(A)):
                newmat[i][j] += (A[i][k] * b[k][j])

    return newmat
def inverse_Matrix(A):
    """"" Calculating An Inverse Of 4X4 Matrix According To Gauss Elimination Method Using Elementary Matrices.
    """
    if find_det(A)!=0:
        I=np.identity(len(A))
        I_copy=I
        A_copy=A
        indicies=list(range(len(A))) # indicies=[0,1,2,3]
        for i in range(len(A)): # any element in range 0-3
            div=1.0/A_copy[i][i]
            for j in range(len(A)):
                A_copy[i][j]*=div
                I_copy[i][j]*=div
            for s in indicies[0:i]+indicies[i+1:]:
                #from index 0 untill i (not including i) and from index i+1 untill the end
                val=A_copy[s][i]
                for x in range(len(A)):
                    A_copy[s][x]-=val*A_copy[i][x]
                    I_copy[s][x]-=val*I_copy[i][x]
    return I_copy

def LU_decomp(A, n):
    L = [[0, 0, 0],[0, 0, 0],[0, 0, 0]]
    U = [[0, 0, 0],[0, 0, 0],[0, 0, 0]]
    for i in range(len(A)):
        for j in range(len(A)):
            if i != j:
                if i > j:
                    L[i][j] = A[i][j]
                else:
                    U[i][j] = A[i][j]

    return L, U
def negativeM(A):
    for i in range(len(A)):
        for j in range(len(A[i])):
            A[i][j] = -A[i][j]
    return A
def sumMat(A, B):
    for i in range(len(A)):
        for j in range(len(A[0])):
            A[i][j] += B[i][j]
    return A
def LU_soloution(mat, n):
    L, U = LU_decomp(mat, n)
    D = initialize_D(mat, n)
    D1 = inverse_Matrix(D)
    return L, U, D, D1
def subMat(A, B):
    sub = 0
    for i in range(len(A)):
        for j in range(1):
            sub += A[i][j] - B[i][j]
    return sub

def Gauss_Seidel_Method(A, n, B):
    """"" Gets  NxN Matrix and returns Calculation Approximation According To Gauss Seidel Method"""
    L,U=LU_decomp(A,n)
    D=[[4,0,0],[0,10,0],[0,0,5]]
    # Diagnol
    Xr = [[0], [0], [0]]
    Xr1 = [[0], [0], [0]]
    sumLD=sumMat(L,D)  #(L+D)
    H=inverse_Matrix(sumLD) #(L+D)^-1
    m=negativeM(H) # -(L+D)^-1
    G=calcM1(m,U) #-(L+D)^-1 * U
    Xr1 = sumMat(calcM2(G, Xr), calcM1(H, B))# According to Formula
    mat=[]
    Xr2=[]
    while (abs(subMat(Xr1, Xr)) > 0.001):
        Xr=Xr1
        Xr1 = sumMat(calcM2(G, Xr), calcM2(H, B))
        for i in Xr1:
            z=[np.absolute(i) for i in Xr1]
        for x,j in enumerate(z):
            print(j)
            if x%3==2:
                print("\t")
    return

def Yakobi_Method(mat, n, B):
    Xr = [[0], [0], [0]]
    Xr1 = [[0], [0], [0]]
    L, U, D, D1 = LU_soloution(mat, n)
    H = D1
    G = negativeM(calcM1(D1, sumMat(L, U)))
    Xr1 = sumMat(calcM2(G, Xr), calcM2(H, B))
    while abs(subMat(Xr1, Xr)) > 0.001:
        Xr = Xr1
        Xr1 = sumMat(calcM2(G, Xr), calcM2(H, B))
        for i in Xr1:
            z = [np.absolute(i) for i in Xr1]
        for x, j in enumerate(z):
            print(j)
            if x % 3 == 2:
                print("\t")
    return 


def main():
    """"" Activating All Functions together"""

    A = [[4, 2, 0],[2, 10, 4], [0, 4, 5]] # Eti's Example From Lesson.
    B = [[2],[6],[5]]
    n = len(A)
    D = initialize_D(A, n)
    if not Is_Dominant(A, n):
        print("No Dominant Diagnol Was Found.")
    else:
        print(colored("Yakobi_Method:","red",attrs=['bold']))
        Yakobi_Method(A, n, B)
        print(colored("Gauss_Seidel_Method:","red",attrs=['bold']))
        print(Gauss_Seidel_Method(A, n, B))


main()
