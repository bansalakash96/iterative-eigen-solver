#include <octave/oct.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <chrono>
#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

//Changes in this Header file than the version1 :
//Change 1 : We will multiply step2 matrix : eigen_vectors_after_swap's 0th column by our matrix a for faster convergence of step 2. 
//Change 2 : As Ganesh said and you also agree, so you should find norm of difference of vectors e_1 and e_0 rather than first taking 
// the absolute value and also divide by n before taking squareroot


// std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count()) ;

std::mt19937 rng(42);
const long int r_max= 4294967295;   

double uniform_random(double a,double b)
{
    double res;
  

    
    res=((double)rng()/(r_max))*(b-a)+a;
 
    return res;
  
}

/* REMEMBER : We are using routines of octave and INDEX starts from 0 HERE AS WELL */

// array       : It should be updated with the values in the absolute Descending Order when this array comes out of this function
// array_index : This should give you the indices of array where the elements were stored before sorting them.
void SORT_absolute_Descending(double arr[],int indices[], int size) {
    // Calculate absolute values and populate the indices array
    for (int i = 0; i < size; ++i) {
        indices[i] = i;
    }

    // Custom sorting based on absolute values (bubble sort)
    for (int i = 0; i < size - 1; ++i) {
        for (int j = 0; j < size - i - 1; ++j) {
            if (std::fabs(arr[indices[j]]) < std::fabs(arr[indices[j + 1]])) {
                // Swap indices if absolute values are in descending order
                int temp = indices[j];
                indices[j] = indices[j + 1];
                indices[j + 1] = temp;
            }
        }
    }

    // Update the original array with the sorted values
    double sortedArr[size];
    for (int i = 0; i < size; ++i) {
        sortedArr[i] = arr[indices[i]];
    }
    for (int i = 0; i < size; ++i) {
        arr[i] = sortedArr[i];
    }
}



std::complex <double> NormComplex(ComplexMatrix x,int n){ 

	std::complex <double> ans=0.0;//Important to initialise with 0

	for(int i=0;i<n;i++){

		ans += x(i,0)*x(i,0); 
	}

	return std::sqrt(ans);
}





 std::complex <double> DotComplex(ComplexMatrix x,ComplexMatrix y,int n){

	std::complex <double> ans;//replace by complex

	for(int i=0;i<n;i++){

		ans += x(i,0)*y(i,0); 
	}

	return ans;
}

std::complex <double> sum_kahan(ComplexMatrix x, int n){

	std:: complex <double> S, C, Y,T;

	S = x(0,0);
	C = 0;

	for(int j=1; j<n; j++){

		Y = x(j,0) - C;
		T = S + Y;
		C = (T-S) - Y;
		S = T;
	}

	return S;
}

 std::complex <double> DotComplex_kahan(ComplexMatrix x,ComplexMatrix y,int n){

	std::complex <double> ans;//replace by complex


	ComplexMatrix products_matrix(n,1);

	for(int i=0;i<n;i++){

		products_matrix(i,0) = x(i,0)*y(i,0);
	}

	ans = sum_kahan( products_matrix, n);

	return ans;
}




std::complex <double> NormComplex_kahan(ComplexMatrix x,int n){ 

	std::complex <double> ans=0.0;//Important to initialise with 0

	ComplexMatrix products_matrix(n,1);

	for(int i=0;i<n;i++){

		products_matrix(i,0) = x(i,0)*x(i,0); 
	}

	ans = sum_kahan( products_matrix, n);

	return std::sqrt(ans);
}



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

/* Basic Power Method */
ComplexMatrix power( Matrix a,ComplexMatrix x,std::complex <double> &lambda,int n,int loops ){



	ComplexMatrix y(n,1);

	for( int i=0;i<loops;i++ ){

		y = toeplitz_symmetric_mat_vec_multiply(a,x,n);

		lambda = DotComplex(y,x,n)/DotComplex(x,x,n); 
		y = y/NormComplex(y,n);


	}
	return y;
}


/* Column_Mum will start from 0 to n-1 */
ComplexMatrix column_vector(ComplexMatrix a,int n,int Column_Num){

	ComplexMatrix result(n,1);

	for(int i=0;i<n;i++){

		result(i,0) = a(i,Column_Num);
	}

	return result;

}

std::complex <double> Vector_Error_Norm( ComplexMatrix x_previous, ComplexMatrix x_next, int n ) {


	ComplexMatrix Error_vector(n,1);

	for(int i=0; i<n ;i++){

		Error_vector(i,0) = x_next(i,0) - x_previous(i,0);
	}

/*****************************************************************************************************************/
	// Dated :: 8 April, 2023 (Saturday)

	std::complex <double> ans=0.0;//Important to initialise with 0

	for(int i=0;i<n;i++){
		ans += Error_vector(i,0)*Error_vector(i,0); 
	}

	return std::sqrt( ans.real()/n );

/*****************************************************************************************************************/


}


ComplexMatrix workinggramschmidt(ComplexMatrix obj, int n, int n_start ){

int row=n;
int col=n_start;
int i,j,k;

	ComplexMatrix base(n,n_start);

 	ComplexMatrix r(row,col);
 	ComplexMatrix v(row,col); 	


    v=obj;

for(j=0;j<col;j++){

for(i=0;i<row;i++){


std::complex <double> temp = NormComplex( column_vector(v,n,j),n );
                      
                  base(i,j) = v(i,j)/temp;

// double temp=v(j).normcomplex();
// base(i,j)=v(i,j)/temp;

}

for(k=j+1;k<col;k++){

std::complex <double> temp2 = DotComplex(column_vector(v,n,k),column_vector(base,n,j),n);
    // complex temp2=dotcomplex(v(k),base(j));
for(i=0;i<row;i++){



v(i,k)=v(i,k)- base(i,j)*temp2;


}
}

}

return base;


}


//UPDATE :: 2 MAY,2023 :: This routine has some issue for sure.
//The last vector is loosing its orthogonality as we keep on adding vectors 
// BUT :: There is no loss in workinggramschmidt ROUTINE :: STRANGE :: But must be a reason.

// ONLY VALID WHEN YOU HAVE ALREADY MADE NORM FOR 0th Vector in your  matrix as UNITY.
// Writing a routine which will be making the last vector of matrix orthogonal to all the previous vectors and then making its norm as 1.
void single_vector_gramschmidt(ComplexMatrix &obj, int n, int n_start ){

	Complex dot_temp;

	for(int j=0;j<n_start-1;j++){

		dot_temp = DotComplex( column_vector( obj , n , j) , column_vector( obj , n , n_start - 1) , n );
		
		for(int i=0;i<n;i++){
			// obj(i, 1) = obj(i,1) - dot_temp*obj(i,0);
			obj(i, n_start-1) = obj(i,n_start-1) - obj(i,j)*dot_temp;
		}
	}

	Complex norm_temp = NormComplex( column_vector( obj , n , n_start-1) ,n);
	// std::cout << " CHECK NORM IS : " << norm_temp << std::endl;

	for(int i=0;i<n;i++){
		obj(i, n_start-1) = obj( i , n_start - 1 )/norm_temp;
	}
}







void function_sort_descending_absolute_magnitude( int n , int p , ComplexMatrix &eigen_values, ComplexMatrix &eigen_vectors, 
	                                             ComplexMatrix& eigen_vectors_after_swap ){

/* Finding real part of eigenvalues in order to use a function for absolute values of eigenvalues in the descending order*/
/************************************************************************************/
double eigen_real_Step1[p];
int array_index_Step1[p];
for(int i=0;i<p;i++){
	eigen_real_Step1[i] = eigen_values(i,0).real();
}

SORT_absolute_Descending( eigen_real_Step1 , array_index_Step1 , p);

for(int i=0;i<p;i++){
	eigen_values(i,0) = eigen_real_Step1[i];
}

/* Sorting eigenvectors in descending order of absolute magnitude of eigenvalues */
/******************************************************************************************/
	for(int j=0;j<p;j++){
			for(int i=0;i<n;i++){
				eigen_vectors_after_swap(i,j) = eigen_vectors(i,array_index_Step1[j]);
								}
							  }

}





ComplexMatrix toeplitz_symmetric_mat_mat_multiply( Matrix a , ComplexMatrix Q, int n , int p ){

	// Important : p is number of vectors in Q (Complex Matrix)
	ComplexMatrix result(n,p);
	
	ComplexMatrix temp(n,1);

	for(int m=0;m<p;m++){

			// column_vector(result , n , m)	= toeplitz_symmetric_mat_vec_multiply( a , column_vector(Q,n,m) , n );
	

		temp = toeplitz_symmetric_mat_vec_multiply( a , column_vector(Q,n,m) , n );
		
		for(int i=0;i<n;i++){
			result(i,m) = temp(i,0);
		}
						}

	return result;

}	

/********************************************************************************************************************************************/




void conventional_lanczos_main ( Matrix a ,  int n , int &Num_of_vectors_in_updated_B , ComplexMatrix &B , ComplexMatrix &eigen_vectors_after_swap ,
 						  			  ComplexMatrix &eigen_values, ComplexMatrix &alpha, ComplexMatrix &beta, ComplexMatrix &w, ComplexMatrix &v){

	std::cout << std::endl;
	std::cout << "Num_of_vectors_in_updated_B " << Num_of_vectors_in_updated_B << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

/*********************************************************************************************************************************************/

	ComplexMatrix temp_B = B;

	B = ComplexMatrix( n, Num_of_vectors_in_updated_B );

	for(int j=0;j<Num_of_vectors_in_updated_B-1;j++){
		for(int i=0;i<n;i++){
			B(i,j) = temp_B(i,j);
		}
	}

 
/*********************************************************************************************************************************************/
/*********************************************************************************************************************************************/

	/* Actual Lanczos implemented in Golub book */ 

	Complex temp;

	for(int i=0;i<n;i++){

		temp = w(i,0);

		w(i,0) = v(i,0)/beta(Num_of_vectors_in_updated_B-2,0);
		
		B(i,Num_of_vectors_in_updated_B-1) = w(i,0);

		v(i,0) = -beta(Num_of_vectors_in_updated_B-2,0)*temp;

		}

	v = v + a*w;	

	alpha(Num_of_vectors_in_updated_B-1,0) = DotComplex( w , v , n );
	v = v - alpha(Num_of_vectors_in_updated_B-1,0)*w;
	beta( Num_of_vectors_in_updated_B-1 , 0 ) = NormComplex( v , n );

	// std::cout << "alpha " << alpha << std::endl;
	std::cout << "beta " << beta(Num_of_vectors_in_updated_B-1,0).real() << std::endl;
	// std::cout << "B is " << B << std::endl;
	// std::cout << std::endl;

	// for(int i=0;i< Num_of_vectors_in_updated_B - 1 ;i++){

	// 	std::cout << DotComplex( column_vector(B,n,i) ,  column_vector(B,n,Num_of_vectors_in_updated_B -1 ) , n )  << "    " ;

	// 	}

	// 	std::cout << std::endl;
	// 	std::cout << std::endl;

/*********************************************************************************************************************************************/
/*********************************************************************************************************************************************/		
	ComplexMatrix M(Num_of_vectors_in_updated_B,Num_of_vectors_in_updated_B);

	for (int i=0;i<Num_of_vectors_in_updated_B;i++){

		M(i,i) = alpha(i,0);
	}

	for(int i=0;i<Num_of_vectors_in_updated_B-1;i++){
		M(i+1,i) = beta(i,0);
		M(i,i+1) = beta(i,0);
	}


	// std::cout << " Tridiagonal matrix M is " << std::endl;
	// std::cout << M << std::endl;	

/*********************************************************************************************************************************************/

	EIG eig2;
	eig2 = EIG(M);
	eigen_values = eig2.eigenvalues();
	ComplexMatrix eigen_vectors = eig2.right_eigenvectors();

	// std::cout << " eigen_vectors is : " << eigen_vectors << std::endl;
	// std::cout << std::endl;


	// std::cout <<  M*column_vector(eigen_vectors , Num_of_vectors_in_updated_B , 0 ) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << eigen_values(0,0)*column_vector(eigen_vectors , Num_of_vectors_in_updated_B , 0 ) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;

	/* Step 2e : Finding eigenvectors closer to actual eigenvectors of a: */
 	ComplexMatrix Eig_iter = B*eigen_vectors; //Very Important Step

 	eigen_vectors_after_swap = ComplexMatrix( n, Num_of_vectors_in_updated_B );//IMPORTANT : Actually, the size increases as add n_add vectors.
    
    function_sort_descending_absolute_magnitude( n ,  Num_of_vectors_in_updated_B , eigen_values, Eig_iter, eigen_vectors_after_swap );

    // std::cout << " Num_of_vectors_in_updated_B is : " << Num_of_vectors_in_updated_B << std::endl;
    // std::cout << "eigenvalues for A in LANCZOS " << std::endl;
    // std::cout << eigen_values << std::endl;

    // std::cout << " eigen_vectors_after_swap for A in LANCZOS " << std::endl;
    // std::cout << eigen_vectors_after_swap << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
/*********************************************************************************************************************************************/

	}


void toeplitz_lanczos_main ( Matrix a ,  int n , int &Num_of_vectors_in_updated_B , ComplexMatrix &B , ComplexMatrix &eigen_vectors_after_swap ,
 						  			  ComplexMatrix &eigen_values, ComplexMatrix &alpha, ComplexMatrix &beta, ComplexMatrix &w, ComplexMatrix &v){

	std::cout << std::endl;
	std::cout << "Num_of_vectors_in_updated_B " << Num_of_vectors_in_updated_B << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

/*********************************************************************************************************************************************/

	ComplexMatrix temp_B = B;

	B = ComplexMatrix( n, Num_of_vectors_in_updated_B );

	for(int j=0;j<Num_of_vectors_in_updated_B-1;j++){
		for(int i=0;i<n;i++){
			B(i,j) = temp_B(i,j);
		}
	}

 
/*********************************************************************************************************************************************/
/*********************************************************************************************************************************************/

	/* Actual Lanczos implemented in Golub book */ 

	Complex temp;

	for(int i=0;i<n;i++){

		temp = w(i,0);

		w(i,0) = v(i,0)/beta(Num_of_vectors_in_updated_B-2,0);
		
		B(i,Num_of_vectors_in_updated_B-1) = w(i,0);

		v(i,0) = -beta(Num_of_vectors_in_updated_B-2,0)*temp;

		}

	v = v + toeplitz_symmetric_mat_vec_multiply(a,w,n);	

	alpha(Num_of_vectors_in_updated_B-1,0) = DotComplex( w , v , n );
	v = v - alpha(Num_of_vectors_in_updated_B-1,0)*w;
	beta( Num_of_vectors_in_updated_B-1 , 0 ) = NormComplex( v , n );

	// std::cout << "alpha " << alpha << std::endl;
	std::cout << "beta " << beta << std::endl;
	// std::cout << "B is " << B << std::endl;
	// std::cout << std::endl;

	for(int i=0;i< Num_of_vectors_in_updated_B - 1 ;i++){

		std::cout << DotComplex( column_vector(B,n,i) ,  column_vector(B,n,Num_of_vectors_in_updated_B -1 ) , n )  << "    " ;

		}

		std::cout << std::endl;
		std::cout << std::endl;

/*********************************************************************************************************************************************/
/*********************************************************************************************************************************************/		
	ComplexMatrix M(Num_of_vectors_in_updated_B,Num_of_vectors_in_updated_B);

	for (int i=0;i<Num_of_vectors_in_updated_B;i++){

		M(i,i) = alpha(i,0);
	}

	for(int i=0;i<Num_of_vectors_in_updated_B-1;i++){
		M(i+1,i) = beta(i,0);
		M(i,i+1) = beta(i,0);
	}


	// std::cout << " Tridiagonal matrix M is " << std::endl;
	// std::cout << M << std::endl;	

/*********************************************************************************************************************************************/

	EIG eig2;
	eig2 = EIG(M);
	eigen_values = eig2.eigenvalues();
	ComplexMatrix eigen_vectors = eig2.right_eigenvectors();

	// std::cout << " eigen_vectors is : " << eigen_vectors << std::endl;
	// std::cout << std::endl;


	// std::cout <<  M*column_vector(eigen_vectors , Num_of_vectors_in_updated_B , 0 ) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;
	// std::cout << eigen_values(0,0)*column_vector(eigen_vectors , Num_of_vectors_in_updated_B , 0 ) << std::endl;
	// std::cout << std::endl;
	// std::cout << std::endl;

	/* Step 2e : Finding eigenvectors closer to actual eigenvectors of a: */
 	ComplexMatrix Eig_iter = B*eigen_vectors; //Very Important Step

 	eigen_vectors_after_swap = ComplexMatrix( n, Num_of_vectors_in_updated_B );//IMPORTANT : Actually, the size increases as add n_add vectors.
    
    function_sort_descending_absolute_magnitude( n ,  Num_of_vectors_in_updated_B , eigen_values, Eig_iter, eigen_vectors_after_swap );

    // std::cout << " Num_of_vectors_in_updated_B is : " << Num_of_vectors_in_updated_B << std::endl;
    std::cout << "eigenvalues for A in LANCZOS " << std::endl;
    std::cout << eigen_values << std::endl;

    // std::cout << " eigen_vectors_after_swap for A in LANCZOS " << std::endl;
    // std::cout << eigen_vectors_after_swap << std::endl;
    std::cout << std::endl;
    // std::cout << std::endl;
/*********************************************************************************************************************************************/

	}

		void store_vector_in_matrix_B(ComplexMatrix &B,ComplexMatrix &vector,int n, int Num_of_vectors_in_updated_B){

			std::cout << "Num_of_vectors_in_updated_B is " << Num_of_vectors_in_updated_B << std::endl;
	        std::cout << std::endl;
			ComplexMatrix temp_B = B;

			B = ComplexMatrix( n, Num_of_vectors_in_updated_B );

		for( int j=0;j < Num_of_vectors_in_updated_B-1; j++){
			for(int i=0;i<n;i++){

				B(i,j) = temp_B(i,j);
			}
		}

		for(int i=0;i<n;i++){
			B(i,Num_of_vectors_in_updated_B-1) = vector(i,0);
                            }
	}

/*****************************************************************************************************************************/
/* COMPUTING EIGENPAIRS USING GNU OCTAVE */
void eigen_compute(ComplexMatrix B, ComplexMatrix M, ComplexMatrix &eigen_values, ComplexMatrix &eigen_vectors_after_swap, int n,
	               int Num_of_vectors_in_updated_B){
EIG eig;
eig = EIG(M);

eigen_values = eig.eigenvalues();
ComplexMatrix eigen_vectors = eig.right_eigenvectors();

/* Step 2e : Finding eigenvectors closer to actual eigenvectors of a: */
ComplexMatrix Eig_iter = B*eigen_vectors; //Very Important Step

eigen_vectors_after_swap = ComplexMatrix( n, Num_of_vectors_in_updated_B );

//IMPORTANT : Actually, the size increases as add n_add vectors.
function_sort_descending_absolute_magnitude( n ,  Num_of_vectors_in_updated_B , eigen_values, Eig_iter
	, eigen_vectors_after_swap );

}
/*******************************************************************************************************************************/


void store_previous_elements_in_M(ComplexMatrix &M,int n, int Num_of_vectors_in_updated_B){

			ComplexMatrix temp_M = M;

			M = ComplexMatrix( Num_of_vectors_in_updated_B, Num_of_vectors_in_updated_B );

			for(int i=0;i<Num_of_vectors_in_updated_B-1;i++){
		for( int j=0;j < Num_of_vectors_in_updated_B-1; j++){

				M(i,j) = temp_M(i,j);
			}
		}

	}


void store_last_row_colum_in_M_tridiagonal(Matrix a,ComplexMatrix &M,ComplexMatrix &B,int n,int Num_of_vectors_in_updated_B){

M(Num_of_vectors_in_updated_B-1,Num_of_vectors_in_updated_B-1) = 
 DotComplex_kahan(column_vector(B,n,Num_of_vectors_in_updated_B-1),a*column_vector(B,n,Num_of_vectors_in_updated_B-1),n);

M(Num_of_vectors_in_updated_B-2,Num_of_vectors_in_updated_B-1) = 
 DotComplex_kahan(column_vector(B,n,Num_of_vectors_in_updated_B-2),a*column_vector(B,n,Num_of_vectors_in_updated_B-1),n);

M(Num_of_vectors_in_updated_B-1,Num_of_vectors_in_updated_B-2) = 
 DotComplex_kahan(column_vector(B,n,Num_of_vectors_in_updated_B-1),a*column_vector(B,n,Num_of_vectors_in_updated_B-2),n);

}

void Updated_store_last_row_colum_in_M_tridiagonal
(ComplexMatrix abi,ComplexMatrix abi_minus_1, ComplexMatrix &M,ComplexMatrix &B,int n,int Num_of_vectors_in_updated_B){

M(Num_of_vectors_in_updated_B-1,Num_of_vectors_in_updated_B-1) = 
 DotComplex_kahan(column_vector(B,n,Num_of_vectors_in_updated_B-1), abi ,n);

M(Num_of_vectors_in_updated_B-2,Num_of_vectors_in_updated_B-1) = 
 DotComplex_kahan(column_vector(B,n,Num_of_vectors_in_updated_B-1),abi_minus_1,n);

M(Num_of_vectors_in_updated_B-1,Num_of_vectors_in_updated_B-2) = 
	M(Num_of_vectors_in_updated_B-2,Num_of_vectors_in_updated_B-1);

}


void write_eigen_values(int Num_of_vectors_in_updated_B,ComplexMatrix eigen_values) {

std::ofstream file;
char fileName[250];

sprintf( fileName, "./eigen_values_conventional/%d.csv", Num_of_vectors_in_updated_B ) ;

file.open(fileName) ;
file << Num_of_vectors_in_updated_B   ;

for(int i=0;i<Num_of_vectors_in_updated_B;i++){
	file << " " << eigen_values(i,0).real()  ;
	}

file << std::endl;
file.close();


}

void compute_eigen_vectors(Matrix a_original,ComplexMatrix &eigen_vectors_original,int n){

ComplexMatrix eigen_values(n,1);

EIG eig;
eig = EIG(a_original);

eigen_values = eig.eigenvalues();
eigen_vectors_original = eig.right_eigenvectors();

}


void matrix_a_original(int n, Matrix &a_original){

for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		a_original(i,j) = 0.0;
	}
}

/* Diagonal Matrix Defined */
double array_original[n] ;

for(int i=0;i<n-2;i++){
	array_original[i] = i*1;
	}
	array_original[n-2] = 550;
	array_original[n-1] = 600;

for(int i=0;i<n;i++){
	a_original(i,i) = array_original[i];
}

}

double compute_2_vectors_error_norm( ComplexMatrix vector_1, ComplexMatrix vector_2, int n){

	std::complex <double> ans=0.0;//Important to initialise with 0
	ComplexMatrix products_matrix(n,1);

	// ComplexMatrix x = vector_1 - vector_2;


ComplexMatrix x(n/2,1);


for(int i=0;i<n/2;i++){
	x(i,0) = vector_1(i,0).real() - vector_2(i,0).real();
}


// for(int i=0;i<n/2;i++){
// 	x(i,0) = fabs(vector_1(i,0).real()) - fabs(vector_2(i,0).real());
// }


	for(int i=0;i<n;i++){
		products_matrix(i,0) = x(i,0)*x(i,0); 
	}

	ans = sum_kahan( products_matrix, n);

	return std::sqrt(ans.real());
}

// only valid for EVEN Num_of_vectors_in_updated_B
void compute_vector_error_norm( int Num_of_vectors_in_updated_B,int n,
	                      ComplexMatrix eigen_vectors_original,
	                      ComplexMatrix eigen_vectors_after_swap ){


Matrix norm_error_vector(Num_of_vectors_in_updated_B/2,1);
double error_norm;


for(int iter=0; iter<Num_of_vectors_in_updated_B/2 ;iter++){

// std::cout << "(" << 2*iter << ", " << n-iter << ") " << std::endl;

error_norm = compute_2_vectors_error_norm( column_vector(eigen_vectors_after_swap,n,2*iter), 
	                                       column_vector(eigen_vectors_original,n,n-iter), n);

norm_error_vector(iter,0) = error_norm;

}

std::cout << "matrix for error norm is " << std::endl;
std::cout << norm_error_vector << std::endl;
std::cout << std::endl;
}



