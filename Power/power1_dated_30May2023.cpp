# include "power_dated_30May2023.h"

int main(){

int n=203;
int n_start = 3;
ComplexMatrix eigen_values(n_start,1);
ComplexMatrix eigen_vectors(n,n_start);

Matrix a = Matrix(n,n);
/*******************************************************************************************************************************/
/* SYMMETRIC TOEPLITZ MATRIX DEFINED */

// double values[n];
// for(int i=0;i<n;i++)
// {
//     values[i] = i + 1;

// }

// for(int i=0;i<n;i++){
// 	a(0,i) = values[i];	
// }

// for(int i=0;i<n;i++){
// 	for(int j=0;j<n;j++){
// 		if( (i-j) == (0))  
// 			a(i,j) = a(0,0);
// }}

// for(int k=1;k<n;k++){
// for(int i=1;i<n;i++){
// for(int j=1;j<n;j++){
//  if( (i-j)==(-k) )  
//  a(i,j) = a(0,k);
// }}}

// for(int k=1;k<n;k++){
// for(int i=0;i<n;i++){
// for(int j=0;j<n;j++){
//  if( (i-j) == (k) )  
//     a(i,j) = a(j,i);
// }}}

// std::cout << "Matrix A of size  (" << n << "x" << n << ") is : " << std::endl;
// std::cout << a << std::endl;

/*******************************************************************************************************************************/

/* Diagonal Matrix Defined */
double array[n] ;

for(int i=0;i< 201;i++){

	array[i] = i*0.0001;
	// array[i] = i+1;
	// array[i] = (i*i + 2*(i- 0.25)) / (i + 1.92);
		// array[i] = ( i*i*(-0.000896) + 0.000000196*(i - 0.25)) / (i + 1.92);	//TRY FOR OTHER LANCZOS
	
}

	array[201] = 2.5;
	array[202] = 3.0;

	// array[198]= 3;
	// array[199]=3.3;


for(int i=0;i<n;i++){
	std::cout << array[i] << std::endl;
}

// Matrix a = Matrix(n,n);



for(int i=0;i<n;i++){

	a(i,i) = array[i];

}

std::cout << "Diagonal Matrix A of size  (" << n << "x" << n << ") is : " << std::endl;
/**********************************************************************************************************************************/






// // int iter=0;
double toler = 1e-6;

// // power_tolerance(a, x, eigen_values, eigen_vectors, n, 0, iter, toler );



general_serial_power(a, n, eigen_values, eigen_vectors, n_start, toler);



std::cout << " eigen values are : " << std::endl;
std::cout << eigen_values << std::endl;

// std::cout << "eigen_vectors are : " << std::endl;
// std::cout << eigen_vectors << std::endl;








	return 0;
}