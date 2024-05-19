The power method provides an intuitive way of finding the dominant eigenvalue of a
square matrix A of size n by n. In this method, one relies on the fact that any arbitrary vector
b0 of size n by 1, upon repeated multiplication with the matrix A, tends to align the vector
in the direction of the dominant eigenvalue. Thus, starting with an arbitrary guess vector
b0 one calculates Ab0,A^2b0,A^3b0,... with the expectation that the sequence will converge to
the dominant eigenvector. To approximate the eigenvalue at any step, there is a scalar value
calculated using the matrix A and the vector x obtained at ith power iteration given by <x,Ax>/
<x,x> also known as the Rayleigh quotient for matrix A and vector x.
