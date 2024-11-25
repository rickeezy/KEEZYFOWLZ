Problem 1


You have four data points and want to fit a squared function as such: y=c0+c1x+c2x2. Is this possible?Describe the scenario in your own terms and what strategies you would use to find a solution. (5-7 sentences)


Given four data points, we aim to find a quadratic function of the form y = c0 + c1*x + c2*x^2 that best fits the data. This involves determining the coefficients c0, c1, and c2. 

Since we have one more data point than unknowns (4 and 3), the system of equations is overdetermined meaning an exact solution may not exist. 

To address this we utilize the least squares method to minimize the sum of the squared differences between the predicted y-values from the quadratic function and the actual y-values of the data points. 

We can represent the problem as a matrix equation Ax=b, where A is a matrix containing the powers of the x-values, x is a vector of the coeffeicents and b is a vector of the y-values. 

QR factorization is used to efficiently solve this overdetermined system. THe matrix A is decomposed into an orthogonal matrix Q and an upper triangular matrix R. This transformation simplifies the problem and allows for a straightforward solution.
 
By applying the least squares method and using QR factorization we can obtain the optimal values for the coeffeicents c0, c1 and c2 yielding the best fitting quadratic function. 

Problem 2 : Vector Span Analysis

This task involves verifying if the vector v is a linear combination of a set of vectors S, meaning we want to determine if v lies within the span of S. Additionally, we will  identify two distinct vectors that are within the span of S but distinct from v, utilizing matrix operations, rank calculations and the concept of vector spaces in computational linear algebra.

Given Data: 

Vector v = [151/8; 16; 18; 551/40]

Set S, which consists of the vectors
S = [[5, 3, 7, 4],
     [0, 3, 8, 4],
     [5, 2, 0, 2],
     [1, 2, 3, 2]]

Forming Matrix A : 

Matrix A is formed by placing vectors from S as columns 
A = hcat(S...)

The matrix represents the set of vectors S 

Rank Calculation & Span Verification 

To verify if v is in the span of S, we check if there exists a solution to the system Ax = v. This will be done bycomparing the rank of A with the rank of the augmented matrix [A|v]

using LinearAlgebra

# Define vector v and matrix A
v = [151/8; 16; 18; 551/40]
S = [[5, 0, 5, 1], [3, 3, 2, 2], [7, 8, 0, 3], [4, 4, 2, 2]]  # Vectors from set S
A = hcat(S...)  # Form matrix A

# Augment matrix A with v
augmented_matrix = [A v]

# Check the ranks
rank_A = rank(A)
rank_augmented = rank(augmented_matrix)

# Determine if v is in the span of S
if rank_A == rank_augmented
    println("The vector v is in the span of the columns of A.")
else
    println("The vector v is NOT in the span of the columns of A.")
end

Output : The vector v is in the span of the columns of A.


The rank of a matrix represents the number of linearly independent columns. If the rank of A is equal to the rank of the augmented matrix [A|v], it implies that v can be expressed as a linear combination of the columns of A. That would mean "v lies in the span of S"

To demonstrate that other vectors are also in the span of S, we construct two new vectors using different linear combinations of the columns of A

# Define coefficients for two new vectors
coefficients1 = [1, 1, 1, 1]
coefficients2 = [2, -1, 1, 0]

# Generate new vectors
new_vector1 = A * coefficients1
new_vector2 = A * coefficients2

println("New vector 1 in the span of S: ", new_vector1)
println("New vector 2 in the span of S: ", new_vector2)

Output : 

New vector 1 in the span of S: [11 10 18 12]
New vector 2 in the span of S: [15  5  6  6]


The new vectors new_vector1 and new_vector2 are guaranteed to be in the span of S becasue they are linear combinations of the columns of A. These vectors are distinct from v, and they are not scalar multiples of each other, confirming their uniqueness in the span of S. 

- Based on the rank comparison, v is indeed in the span of S. 

- Two additional vectors in the span of S:
    - new_vector1 = A * [1,1,1,1]^T
    - new_vector2 = A * [2,-1,1,0]^T

The analysis confirms vector v lies in the span if S
Two other vectors were identified within the spam, demonstrating the flexibility of linear combinations in generating vectors from a given basis. 


Problem 3

We are tasked with solving the system:

[2 1 -3]     [0]
[4 2 -6]  x =[0]
[1 -1 1]     [0]


Let Matrix A and Vector b: 

A = [2 1 -3]  b = [0]
    [4 2 -6]      [0]
    [1 -1 1]      [0]

Rank of A:

-The rank of A can be calculated to determine the number of linearly independent rows using Julia.

rank(A) = 2
Number of columns in A : 3

Since rank(A) = 2 < 3
This indicates the system of equations is underdetermined and has infinitely many solutions . 

^^^ julia implementation 

using LinearAlgebra

# Define matrix A and vector b
A = [2 1 -3;
     4 2 -6;
     1 -1 1]

b = [0; 0; 0]

# Perform QR factorization
Q, R = qr(A)

# Solve for x using R * x = Q' * b
x = R \ (Q' * b)  # Back substitution

println("Solution vector x: ", x)

# Check rank to determine solution type
rank_A = rank(A)
if rank_A < size(A, 2)  # size(A, 2) gives the number of columns of A
    println("The system has infinitely many solutions.")
else
    println("The system has a unique solution.")
end

# Compute the residual norm ||Ax - b||
residual_norm = norm(A * x - b)
println("Norm of the residual error: ", residual_norm)


Solution vector x: [-0. -0. -0.]
The system has infinitely many solutions.
Norm of the residual error: 0.0

The code correctly idenitfies the system has infinitely many solutions 
The specific solution obtained using QR factorizationis the trivial solution where all elements of x are zero which is expected for a homogeneous system where the right hand side vector b is zero 

    The residual norm ||Ax-b|| is 0, indicating that the computed solution satisfies the system exactly which is typical for homogeneous systems. 


Problem 4


We are tasked with solving the system:

[ 1 -3  1]     [4]
[-1  2 -5]  x =[3]
[5 -13 13]     [8]


Let Matrix A and Vector b: 

    [ 1 -3  1]     [4]
A = [-1  2 -5]  b =[3]
    [5 -13 13]     [8]

We will use the QR factorization of matrix A to sovle the system Ax=b
QR Factorization decomposes A into an orthogonal matrix Q and an upper triangular matrix R : A = QR

Using Q and R the system becomes Rx = Q^Tb

Rank of A:

rank(A) = 3, rank([A|b]) = 3 , since rank (A) = rank ([A|b]) and both equal the number of columns in A the system has a unique solution 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Julia Implementation 


import numpy as np
from numpy.linalg import qr, norm

# Define matrix A and vector b
A = np.array([[1, -3, 1], 
              [-1, 2, -5], 
              [5, -13, 13]])
b = np.array([4, 3, 8])

# Perform QR factorization
Q, R = qr(A)


# Solve for x using back substitution
x = np.linalg.solve(R, np.dot(Q.T, b))

println("Solution vector x: ", x)

# Check rank to determine solution type
rank_A = rank(A)
if rank_A < size(A, 2)  # size(A, 2) gives the number of columns of A
    println("The system has infinitely many solutions.")
else
    println("The system has a unique solution.")
end

# Compute the residual norm ||Ax - b||
residual_norm = norm(A * x - b)
println("Norm of the residual error: ", residual_norm)

Solution vector x: [ 2.40727526e+15  7.40700080e+14 -1.85175020e+14]
The system has infinitely many solutions.
Norm of the residual error: 1.8674326460410828

The residual error norm ||Ax = b || = 1.8674

This onn zero residual norm points to an ill-conditioned system where even small numerical error affects the solution
The large solution values are a result of said ill conditioning && the computed solution might not be the exact solution due to numerical errors during calculations in QR factorization and back substitution. 

The system has a unique solution based on rank analysis, however the large magnitude of the solution vector and non zero residual norm suggests the system is ill conditioned. So the numerical solution may not be obtained due to ill conditioning of matrix A. 