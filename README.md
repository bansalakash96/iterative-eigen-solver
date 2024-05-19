# iterative-eigen-solver
This repository contains iterative techniques for solving eigenvalue problems, with a focus on large sparse matrices. Various methods are employed, including power iterations and Krylov subspace methods such as the Lanczos algorithm.
Power iterations: The power method provides an intuitive way of finding the dominant eigenvalue of a
square matrix A of size n by n. In this method, one relies on the fact that any arbitrary vector
b0 of size n by 1, upon repeated multiplication with the matrix A, tends to align the vector
in the direction of the dominant eigenvalue.
Lanczos algorithm: Lanczos algorithm is an effective tool for constructing an approximate tridiagonalization of a
symmetric matrix. The basic procedure creates a set of vectors and makes them orthogonal
to create an orthonormal basis of the Krylov subspace. While the original formulation failed
in finite precision, subsequent modification where one reorthogonalizes the set of vectors
is quite stable. In its current form, it is known for its efficiency in computing a subset
of eigenvalues and eigenvectors for large sparse symmetric matrices via this tridiagonal
representation.
