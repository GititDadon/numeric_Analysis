import numpy as np
def find_det(A):
    """"" Returns The Determinant Of A Matrix"""
    S=np.array(A)
    return np.linalg.det(S)
def inverse_Matrix(A):
    """"" Calculating An Inverse Of 4X4 Matrix According To Gauss Elimination Method Using Elementary Matrices.
    """
# Matrix multiplication 4X4#
def calc_Matrix(A, B):
	result = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
	for i in range(len(A)):
		# iterate  columns
		for j in range(len(B[0])):
			# iterate  rows
			for k in range(len(A)):
				result[i][j] += A[i][k] * B[k][j]

	return result


# Matrix multiplication 4X1#
def Multiply_Matrix(A, B):
	result = [[0, ], [0, ], [0, ], [0, ]]
	for i in range(len(A)):
		# iterate through columns of Y
		for j in range(len(B[0])):
			# iterate through rows of Y
			for k in range(len(A)):
				result[i][j] += (A[i][k] * B[k][j])

	return result


# Inversion matrix calculation#
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

# Decomposition into L and U matrices#
def LU_dec():
	X = [[1,5, 1, 0],
		 [1, 2, 0, 0],
		 [0, 1, 1, 3],
		 [0, 1, 0, 4]]
	I=np.identity(4)
	#Intialziation
	L = calc_Matrix(I, I)
	for j in range(len(X)):
		for i in range(len(X)):
			if i > j:
				I=np.identity(4)
				if X[i][j] != 0:
					I[i][j] = -X[i][j] / X[j][j]
					X = calc_Matrix(I, X)
					I[i][j] = -I[i][j]
					L = calc_Matrix(L, I)
	return L, X


# Matrix solution using the LU method#
def LU_Soloution(S):
	L, U = LU_dec()
	L1 = inverse_Matrix(L)
	U1 = inverse_Matrix(U)
	z = calc_Matrix(L1, S)
	m = calc_Matrix(U1, z)
	print("The Soloution Of LU Linear Equations System Is :")
	for i in m:
		print(i)
A=[[1,0,4,-6],[2,5,0,3],[-1,2,3,5],[2,1,-2,3]]
print("The Inverse Of The Matrix Is:",inverse_Matrix(A))
print(LU_Soloution(A))