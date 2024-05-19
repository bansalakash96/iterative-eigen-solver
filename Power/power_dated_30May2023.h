#include <octave/oct.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <chrono>


/*******************************************************************************************************************************************************/
/* Routine for finding random numbers uniformly distributed between (a,b) */
std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count()) ;
const long int r_max= 4294967295;   

double uniform_random(double a,double b)
{
    double res;
  

    
    res=((double)rng()/(r_max))*(b-a)+a;
 
    return res;
  
}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
std::complex <double> NormComplex(ComplexMatrix x,int n){ 

	std::complex <double> ans=0.0;//Important to initialise with 0

	for(int i=0;i<n;i++){

		ans += x(i,0)*x(i,0); 
	}

	return std::sqrt(ans);
}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
std::complex <double> DotComplex(ComplexMatrix x,ComplexMatrix y,int n){

	std::complex <double> ans=0.0;//replace by complex

	for(int i=0;i<n;i++){

		ans += x(i,0)*y(i,0); 
	}

	return ans;
}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
/* This routine is used to extract a column vector from a Matrix
/* Column_Num will start from 0 to n-1 */
ComplexMatrix column_vector(ComplexMatrix a,int n,int Column_Num){

	ComplexMatrix result(n,1);

	for(int i=0;i<n;i++){

		result(i,0) = a(i,Column_Num);
	}

	return result;

}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
/* This routine is used to extract a column vector from a Matrix
/* Column_Num will start from 0 to n-1 */
Matrix column_vector_real(ComplexMatrix a,int n,int Column_Num){

	Matrix result(n,1);

	for(int i=0;i<n;i++){

		result(i,0) = a(i,Column_Num).real();
	}

	return result;

}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
/* In this routine, we divide by n i.e. the number of components of the vector also as Prof. told you */
double Vector_Error_Norm( ComplexMatrix x_previous, ComplexMatrix x_next, int n ) {

	ComplexMatrix Error_vector(n,1);

	for(int i=0; i<n ;i++){
		Error_vector(i,0) = x_next(i,0) - x_previous(i,0);
	}

	std::complex <double> ans=0.0;//Important to initialise with 0

	for(int i=0; i<n; i++){
		ans += Error_vector(i,0)*Error_vector(i,0); 
	}

	std::cout << " square of vector_error norm ::   is "  <<	ans << std::endl;
	std::cout << " vector_error norm ::   is "  <<	std::sqrt( ans.real() ) << std::endl;

	return std::sqrt( ans.real() / n );

}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
/* In this routine, we divide by n i.e. the number of components of the vector also as Prof. told you */
double absolute_difference_Vector_Error_Norm( ComplexMatrix x_previous, ComplexMatrix x_next, int n ) {

	Matrix Error_vector(n,1);

	for(int i=0; i<n ;i++){
		Error_vector(i,0) = fabs(x_next(i,0).real()) - fabs(x_previous(i,0).real());
	}

	double ans=0.0;//Important to initialise with 0

	for(int i=0; i<n; i++){
		ans += Error_vector(i,0)*Error_vector(i,0); 
	}

	// std::cout << " square of vector_error norm ::   is "  <<	ans << std::endl;
	// std::cout << " vector_error norm ::   is "  <<	std::sqrt( ans.real() ) << std::endl;

	return std::sqrt( ans / n );

}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
/* This routine is used for Matrix-Vector multiplication of a Symmetric-Toeplitz Matrix where we are storing matrix A in a vector form*/ 
ComplexMatrix toeplitz_symmetric_mat_vec_multiply( Matrix a , ComplexMatrix x, int n ){

	ComplexMatrix result(n,1);

	/* RULE 1 */

	for( int j=1;j<n+1;j++){
		for(int i=1;i<n-j+2;i++){

		result(j-1,0) += a(i-1,0)*x(i + j -1 - 1 , 0 );
		}
	}

	/* RULE 2 */

	for( int j=2; j < n+1 ; j++){
		for(int i=2; i < j+1 ; i++){
			result(j-1,0) += a(i-1,0)*x(j-i , 0);
		}
	}

	return result;
	}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
/* Basic Power Method */
ComplexMatrix power( Matrix a,ComplexMatrix x,std::complex <double> &lambda,int n,int loops ){

	ComplexMatrix y(n,1);

	for( int i=0;i<loops;i++ ){

		// y = toeplitz_symmetric_mat_vec_multiply(a,x,n);
		y = a*x;
		lambda = DotComplex(y,x,n)/DotComplex(x,x,n); 
		y = y/NormComplex(y,n);


	}
	return y;
}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
void power_tolerance( Matrix a, ComplexMatrix x, ComplexMatrix &eigen_values, ComplexMatrix &eigen_vectors, int n, int  m, int &iter, double toler ){


	ComplexMatrix temp(n,1);
	std::complex <double> lambda_previous;
	std::complex <double> temp_lambda;		
	double rel_error = 1000;
	double vector_error_norm = 1000;
	 std::complex <double> inner_product;

	 double Error_VECTOR;

	 std:: ofstream file;
	 char fileName[250];
	 // sprintf(fileName,"./power_method_files/example1/eigen_value_%d.csv", m+1 ) ;

	 sprintf(fileName,"./eigen_values_data/example3/eigen_value_%d.csv", m+1 ) ;
  	 
  	 file.open(fileName) ;

  	 	file << m+1 << " " << 0 << " " << 0 << " " << 0  << std::endl;
  	 	file << std::endl;

	// while (rel_error > toler  || Error_VECTOR > toler){

	while( rel_error > toler ){	
		
		// temp = toeplitz_symmetric_mat_vec_multiply(a,x,n);
		temp = a*x;
		temp_lambda = DotComplex(temp,x,n)/DotComplex(x,x,n); 
		temp = temp/NormComplex(temp,n);

		std::cout << "temp_lambda is " << temp_lambda << std::endl;
		rel_error = fabs( ( lambda_previous.real() - temp_lambda.real() )/ temp_lambda.real()  );
		vector_error_norm = Vector_Error_Norm(x , temp , n);

		inner_product = DotComplex(x , temp , n);

		Error_VECTOR = fabs(1 - fabs(inner_product)); 


		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << " inner_product :::: " << inner_product << std::endl;
		std::cout << " ERROR VECTOR is " << Error_VECTOR << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		// vector_error_norm = absolute_difference_Vector_Error_Norm(x , temp , n);
			

		// std::cout << " THE 2 VECTORS ARE : " << std::endl;
		// std::cout << "x is " << x << std::endl;
		// std::cout << "temp vector CHANGING signs is: " << temp << std::endl;

		// std::cout << " vector_error_norm is: " << vector_error_norm << std::endl;
		// std::cout << "rel_error is " << rel_error << std::endl;
		lambda_previous = temp_lambda;
		x = temp;
		++iter;

		file << iter << " " << lambda_previous.real() << " " << rel_error << " " << Error_VECTOR << std::endl;

	}

	file.close();

	std::cout << std::endl;
	std::cout << std::endl;
		std::cout << std::endl;
	std::cout << std::endl;
		std::cout << std::endl;
	std::cout << std::endl;
	std::cout << " Number of iterations is " << iter << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
		std::cout << std::endl;
	std::cout << std::endl;
		std::cout << std::endl;
	std::cout << std::endl;


	eigen_values(m,0) = temp_lambda;

	for(int i=0;i<n;i++){
		eigen_vectors(i,m) = temp(i,0);
	}

}
/*******************************************************************************************************************************************************/

/*******************************************************************************************************************************************************/
void general_serial_power( Matrix  &first, int n, ComplexMatrix &eigen_values, ComplexMatrix &eigen_vectors,  int m, double toler){

/*Declarations*/

int loopcount[m];

/*Initialising all loops to 0*/
for(int i=0;i<m;i++)
loopcount[i]=0;


ComplexMatrix x(n,1);
/*Initialising the vectors with starting columns of matrix itself*/
for(int i=0;i<n;i++)
x(i,0) = uniform_random(-1,1);

for(int i=0;i<m;i++){

power_tolerance( first, x, eigen_values, eigen_vectors, n, i, loopcount[i], toler );

first = first - (  (column_vector_real(eigen_vectors, n, i)*column_vector_real(eigen_vectors, n, i).transpose()) * eigen_values(i,0).real()  );

}


}
/*******************************************************************************************************************************************************/
