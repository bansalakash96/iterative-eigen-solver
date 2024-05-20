//Dated: 16 May, 2024(Thursday)
#include "conventional_lanczos_functions.h"


int main(){

int n=503;

/* Diagonal Matrix Defined */
double array[n] ;

for(int i=0;i< n-2;i++){
	array[i] = i;
	}

	array[n-2] = 550;
	array[n-1] = 600;

// std::cout << " array is " << std::endl;

// for(int i=0;i<n;i++){
// 	std::cout << "array[" << i << "] is " << array[i] << std::endl;
// }

Matrix a = Matrix(n,n);
for(int i=0;i<n;i++){
	a(i,i) = array[i];
}

// std::cout << "a is " << std::endl;
// std::cout << a << std::endl;
// std::cout << std::endl;

int Num_of_vectors_in_updated_B = 1;//Total Size of orthogonal matrix to be passed to Gram-Schmidt

/**********************************************************************************************************************************/

ComplexMatrix B(n,Num_of_vectors_in_updated_B);
ComplexMatrix B_conventional(n,Num_of_vectors_in_updated_B);
ComplexMatrix eigen_vectors_after_swap(n,Num_of_vectors_in_updated_B);
ComplexMatrix eigen_values(Num_of_vectors_in_updated_B,1);
ComplexMatrix eigen_vectors_after_swap_conventional(n,Num_of_vectors_in_updated_B);
ComplexMatrix eigen_values_conventional(n,Num_of_vectors_in_updated_B);
ComplexMatrix M_conventional_using_kahan_dots(Num_of_vectors_in_updated_B,Num_of_vectors_in_updated_B);


/*******************************************/

//Computing Random vector d0
ComplexMatrix d0(n,1);
for(int i=0;i<n;i++){
	// d0(i,0) = 1.0;
	d0(i,0) = uniform_random(-1,1);
}
/*******************************************/

/*********************************************************************************/
//Computing vector b0
ComplexMatrix b0(n,1);
b0 = d0;
double norm_b0 = NormComplex_kahan(b0,n).real();
b0 = b0/norm_b0;
for(int i=0;i<n;i++){
	B_conventional(i,0) = b0(i,0);
}
M_conventional_using_kahan_dots(0,0) = DotComplex_kahan(column_vector(B_conventional,n,0),a*column_vector(B_conventional,n,0),n);
eigen_compute(B_conventional,M_conventional_using_kahan_dots,eigen_values_conventional,eigen_vectors_after_swap_conventional,
	          n,Num_of_vectors_in_updated_B);
/*********************************************************************************************************************************/

ComplexMatrix b1_conventional(n,1);

b1_conventional = a*b0 - DotComplex_kahan(a*b0,b0,n)*b0;
double norm_b1_conventional = NormComplex_kahan(b1_conventional,n).real();
b1_conventional = b1_conventional/norm_b1_conventional;

/*******************************************************************************************************/	
/* MAKING THE COLUMNS OF MATRIX B as 2 */
Num_of_vectors_in_updated_B +=1;
store_vector_in_matrix_B(B_conventional,b1_conventional,n,Num_of_vectors_in_updated_B);

store_previous_elements_in_M(M_conventional_using_kahan_dots,n,Num_of_vectors_in_updated_B);
store_last_row_colum_in_M_tridiagonal(a,M_conventional_using_kahan_dots,B_conventional,n,Num_of_vectors_in_updated_B);
eigen_compute(B_conventional,M_conventional_using_kahan_dots,eigen_values_conventional,eigen_vectors_after_swap_conventional,
	          n,Num_of_vectors_in_updated_B);
std::cout << eigen_values_conventional << std::endl;
write_eigen_values( Num_of_vectors_in_updated_B, eigen_values_conventional);
// std::cout << eigen_vectors_after_swap_conventional << std::endl;
/*******************************************************************************************************/		

/*******************************************************************************************************/
ComplexMatrix b2_conventional(n,1);

b2_conventional = a*b1_conventional - DotComplex_kahan(a*b1_conventional,b1_conventional,n)*b1_conventional
								    - DotComplex_kahan(a*b1_conventional,b0,n)*b0;
double norm_b2_conventional = NormComplex_kahan(b2_conventional,n).real();
b2_conventional = b2_conventional/norm_b2_conventional;

/*******************************************************************************************************/	
/* MAKING THE COLUMNS OF MATRIX B as 3 */
Num_of_vectors_in_updated_B +=1;
store_vector_in_matrix_B(B_conventional,b2_conventional,n,Num_of_vectors_in_updated_B);
store_previous_elements_in_M(M_conventional_using_kahan_dots,n,Num_of_vectors_in_updated_B);
store_last_row_colum_in_M_tridiagonal(a,M_conventional_using_kahan_dots,B_conventional,n,Num_of_vectors_in_updated_B);
eigen_compute(B_conventional,M_conventional_using_kahan_dots,eigen_values_conventional,eigen_vectors_after_swap_conventional,
	          n,Num_of_vectors_in_updated_B);
write_eigen_values( Num_of_vectors_in_updated_B, eigen_values_conventional);
std::cout << eigen_values_conventional << std::endl;
// std::cout << eigen_vectors_after_swap_conventional << std::endl;
/*******************************************************************************************************/		


ComplexMatrix b_k_minus_1_conventional(n,1);
ComplexMatrix b_k_minus_2_conventional(n,1);
ComplexMatrix b_k_conventional(n,1);
double alpha_conventional;
double beta_conventional;
double norm_b_k_conventional;

for(int k=0;k<100;k++){

b_k_minus_1_conventional = column_vector(B_conventional,n,Num_of_vectors_in_updated_B-1);
b_k_minus_2_conventional = column_vector(B_conventional,n,Num_of_vectors_in_updated_B-2);

alpha_conventional = DotComplex_kahan(a*b_k_minus_1_conventional,b_k_minus_1_conventional,n).real();
beta_conventional  = DotComplex_kahan(a*b_k_minus_1_conventional,b_k_minus_2_conventional,n).real();

b_k_conventional = a*b_k_minus_1_conventional-alpha_conventional*b_k_minus_1_conventional-beta_conventional*b_k_minus_2_conventional;

norm_b_k_conventional = NormComplex_kahan(b_k_conventional,n).real();

std::cout << "Norm of b_k_conventional is " << norm_b_k_conventional << std::endl;
std::cout << std::endl;

b_k_conventional = b_k_conventional/norm_b_k_conventional;

Num_of_vectors_in_updated_B +=1;
store_vector_in_matrix_B(B_conventional,b_k_conventional,n,Num_of_vectors_in_updated_B);
store_previous_elements_in_M(M_conventional_using_kahan_dots,n,Num_of_vectors_in_updated_B);
store_last_row_colum_in_M_tridiagonal(a,M_conventional_using_kahan_dots,B_conventional,n,Num_of_vectors_in_updated_B);
eigen_compute(B_conventional,M_conventional_using_kahan_dots,eigen_values_conventional,eigen_vectors_after_swap_conventional,
	          n,Num_of_vectors_in_updated_B);
write_eigen_values( Num_of_vectors_in_updated_B, eigen_values_conventional);

// std::cout << eigen_values_conventional << std::endl;

}

	return 0;
}