# include "block_lanczos.h"

int main(){

int n=503;
int p=2;//Number of random vectors in the starting block 
int number_of_blocks_in_Xs=1;
int num_of_vectors_in_Xs = p*number_of_blocks_in_Xs;
ComplexMatrix Xs;
ComplexMatrix Ms;
std:: ofstream file;
char fileName[250];

ComplexMatrix eigen_values;
ComplexMatrix eigen_vectors_after_swap;

/**********************************************************************************************************************************/
/* Diagonal Matrix Defined */
double array[n] ;

for(int i=0;i< n-2;i++){
	array[i] = i*1;
}

array[n-1] = 550;
array[n-2] = 600;

// array[n-1] = 1050;
// array[n-2] = 1100;


Matrix a = Matrix(n,n);

for(int i=0;i<n;i++){
	a(i,i) = array[i];
}

std::cout << "Diagonal Matrix A of size  (" << n << "x" << n << ") is : " << std::endl;
// std::cout << a << std::endl;
// std::cout << std::endl;
/**********************************************************************************************************************************/



/*****************************************************************************************************************************************/

ComplexMatrix X1(n,p);

X1 = intialize_random_vectors(n,p);

// std::cout << " InitiaL BLOCK of Random vectors is: " << std::endl;
// std::cout << X1 << std::endl;
// std::cout << std::endl;


/* QR Test2 */
// for(int i=0;i<p;i++){

// std::cout << "Norm of column number " << i << " for matrix q_test is " << NormComplex(column_vector(X1,n,i),n)  << std::endl;

// }
// std::cout << std::endl;

/* QR Test3 */
// std::cout << "All Dotproducts: " << std::endl;

// for(int i=0;i< n -1;i++){
// 	for(int j=i+1;j< p ;j++){

// std::cout << DotComplex( column_vector(X1,n,i) ,  column_vector(X1,n,j) , n )  << "    " ;
// 	}
// }
// std::cout << std::endl;


ComplexMatrix M1(p,p);
ComplexMatrix Z2(n,p);
ComplexMatrix X2;
ComplexMatrix R2; 

// std::cout << "transpose of X1 is " << std::endl;
// std::cout << transpose(X1,n,p) << std::endl;
// std::cout << std::endl;

M1 = transpose(X1,n,p)*a*X1;

/* Computing Z2 */
Z2 = a*X1 - X1*M1;

// std::cout << " Z2 is "  << std::endl;
// std::cout << Z2 << std::endl;
// std::cout << std::endl;

octave::math::qr<ComplexMatrix> qrObject( Z2 , octave::math::qr<ComplexMatrix>::economy);
X2 = qrObject.Q();
R2 = qrObject.R();

// std::cout << "  BLOCK X2 is: " << std::endl;
// std::cout << X2 << std::endl;
// std::cout << std::endl;

// std::cout << " Upper triangular matrix R2 is " << std::endl;
// std::cout << R2 << std::endl;
// std::cout << std::endl;

/* QR Test2 */
// for(int i=0;i<p;i++){

// std::cout << "Norm of column number " << i << " for matrix q_test is " << NormComplex(column_vector(X2,n,i),n)  << std::endl;

// }
// std::cout << std::endl;

/* QR Test3 */
// std::cout << "All Dotproducts: " << std::endl;

// for(int i=0;i< n -1;i++){
// 	for(int j=i+1;j< p ;j++){

// std::cout << DotComplex( column_vector(X2,n,i) ,  column_vector(X2,n,j) , n )  << "    " ;
// 	}
// }
// std::cout << std::endl;


/* Computed X1 and X2 blocks */
/**********************************************************************************************************************************************/
Ortho_check_2_blocks(X1, X2, n, p);

Xs = ComplexMatrix(X1);


for(int j=0;j<p;j++){
	for(int i=0;i<n;i++){
	
		Xs(i,j) = X1(i,j);
	
	}
}
// std::cout << " Initial Xs is " << std::endl;
// std::cout << Xs << std::endl;
// std::cout << std::endl;

number_of_blocks_in_Xs  +=1;
num_of_vectors_in_Xs += p;
Concatenate(Xs,X2,n,p,number_of_blocks_in_Xs, num_of_vectors_in_Xs);

// std::cout << "  Xs after X1 and X2 is " << std::endl;
// std::cout << Xs << std::endl;
// std::cout << std::endl;


Ms = transpose(Xs, n, num_of_vectors_in_Xs)*a*Xs;

std::cout << std::endl;
std::cout << " Block tridiagonal matrix Ms is " << std::endl;
std::cout << Ms << std::endl;

/*********************************************************************************************************************************************/

	EIG eig2;
	eig2 = EIG(Ms);
	eigen_values = eig2.eigenvalues();
	ComplexMatrix eigen_vectors = eig2.right_eigenvectors();

 	ComplexMatrix Eig_iter = Xs*eigen_vectors; //Very Important Step

 	eigen_vectors_after_swap = ComplexMatrix( n, num_of_vectors_in_Xs );//IMPORTANT : Actually, the size increases as add n_add vectors.
    
    function_sort_descending_absolute_magnitude( n ,  num_of_vectors_in_Xs , eigen_values, Eig_iter, eigen_vectors_after_swap );

    std::cout << "eigenvalues for A in Block LANCZOS " << std::endl;
    std::cout << eigen_values << std::endl;

    // std::cout << " eigen_vectors_after_swap for A in LANCZOS " << std::endl;
    // std::cout << eigen_vectors_after_swap << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
/*********************************************************************************************************************************************/
sprintf( fileName, "./block_eigen_values_lanczos_main/update_test1_M_PENTADIAGONAL/%d.csv", num_of_vectors_in_Xs ) ;


file.open(fileName) ;

file << num_of_vectors_in_Xs   ;

for(int i=0;i<num_of_vectors_in_Xs;i++){
	file << " " << eigen_values(i,0).real()  ;
	}

file << std::endl;
file.close();




/*********************************************************************************************************************************************/

ComplexMatrix Xi;
ComplexMatrix Ximinus1;
ComplexMatrix Mi;
ComplexMatrix Ri;
ComplexMatrix Ziplus1;

/* Let's start the loop to Compute blocks: X3,X4,...,Xs-1 */

Xi = X2;
Ximinus1 = X1;
Mi =  transpose(Xi,n,p)*a*Xi;
Ri = R2;

/* LOOPS starts */ 
for(int q=0;q<150;q++){
// for(int q=0;q<5;q++){

function_loops_block_lanczos( a, Xs, Ziplus1, Xi, Ximinus1, Mi, Ri, X2, X1, R2, n, p, number_of_blocks_in_Xs, num_of_vectors_in_Xs,
							 eigen_values, eigen_vectors, eigen_vectors_after_swap, Eig_iter );


sprintf( fileName, "./block_eigen_values_lanczos_main/update_test1_M_PENTADIAGONAL/%d.csv", num_of_vectors_in_Xs ) ;
file.open(fileName) ;
file << num_of_vectors_in_Xs   ;


std::cout << " num_of_vectors_in_Xs is  "<< num_of_vectors_in_Xs << std::endl;

for(int i=0;i<num_of_vectors_in_Xs;i++){
	file << " " << eigen_values(i,0).real()  ;
	}

file << std::endl;
file.close();

}

	return 0;

}




