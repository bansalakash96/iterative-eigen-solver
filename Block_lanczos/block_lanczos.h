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

std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count()) ;
const long int r_max= 4294967295;   

double uniform_random(double a,double b)
{
    double res;
  

    
    res=((double)rng()/(r_max))*(b-a)+a;
 
    return res;
  
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


/* Column_Mum will start from 0 to n-1 */
ComplexMatrix column_vector(ComplexMatrix a,int n,int Column_Num){

        ComplexMatrix result(n,1);

        for(int i=0;i<n;i++){

                result(i,0) = a(i,Column_Num);
        }

        return result;

}

/* REMEMBER : We are using routines of octave and INDEX starts from 0 HERE AS WELL */

// array       : It should be updated with the values in the absolute Descending Order when this array comes out of this function
// array_index : This should give you the indices of array where the elements were stored before sorting them.
void SORT_absolute_Descending( double array[] , int array_index[] ,int n ){

        double abs_array[n];
        double sorted_array[n];

        for(int i=0;i<n;i++){

                        abs_array[i] = fabs(array[i]);

                   }
 
        for( int loop=0; loop<n; loop++ ){

                        double max_value = ( *( std::max_element( abs_array , abs_array + n ) ) );

                                        for( int i=0;i<n;i++ ){

                                                        if( abs_array[i] == max_value ){

                                                                sorted_array[loop] = array[i];
                                                                array_index[loop] = i;
           
                                                        abs_array[i] = 0;
                                                                                   }    
                                                             }

                                                             }

        for(int i=0;i<n;i++){

                array[i] = sorted_array[i];

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


ComplexMatrix transpose(ComplexMatrix a, int &n, int &p){

        ComplexMatrix result(p,n);


        for(int i=0;i<p;i++){
                for(int j=0;j<n;j++){
                        result(i,j) = a(j,i);
                }
        }

        return result;
        
        }


/* Initializing the starting block with p random vectors */
ComplexMatrix intialize_random_vectors(int n,int p){

ComplexMatrix result_random(n,p);

for(int i=0;i<n;i++){
        for(int j=0;j<p;j++){
                result_random(i,j) = uniform_random(-1,1);
        }
}

// std::cout << result_random << std::endl;

ComplexMatrix q_test;
ComplexMatrix r_test;

octave::math::qr<ComplexMatrix> qrObject(result_random, octave::math::qr<ComplexMatrix>::economy);
q_test = qrObject.Q();
r_test = qrObject.R();

// std::cout << r_test << std::endl;

return q_test;

}

void Ortho_check_2_blocks( ComplexMatrix block1, ComplexMatrix block2,int n, int p){

        for(int i=0;i<p;i++){

                for(int j=0;j< p ;j++){

                        std::cout << DotComplex( column_vector(block1,n,i) ,  column_vector(block2,n,j) , n )  << "    " ;
        
                                      }

                        std::cout << std::endl;

                            }
        }


ComplexMatrix Concatenate(ComplexMatrix &Xs, ComplexMatrix X2, int n, int p, int number_of_blocks_in_Xs, int num_of_vectors_in_Xs){

        

        ComplexMatrix temp_Xs = Xs;

        // std::cout << " num_of_vectors_in_Xs is " << num_of_vectors_in_Xs << std::endl;

        Xs = ComplexMatrix(n, num_of_vectors_in_Xs);

        for(int j=0;j< num_of_vectors_in_Xs - p;j++){
                for(int i=0;i<n;i++){
                        Xs(i,j) = temp_Xs(i,j);
                }
        }



        for(int j=num_of_vectors_in_Xs - p; j < num_of_vectors_in_Xs ; j++ ){

                for(int i=0;i<n;i++){

                        // std::cout << "(i,j) is (" <<  i << "," << j << ") "<< std::endl;

                        Xs(i,j) = X2(i,j - (num_of_vectors_in_Xs - p) );
                }
        }

        return Xs;

   }


void function_loops_block_lanczos( Matrix a, ComplexMatrix &Xs,  ComplexMatrix &Ziplus1, ComplexMatrix & Xi, 
                                   ComplexMatrix & Ximinus1,
                                   ComplexMatrix & Mi,
                                   ComplexMatrix &Ri, 
                                  ComplexMatrix &X2, ComplexMatrix &X1, ComplexMatrix &R2,
                                  int n, int p, int &number_of_blocks_in_Xs, int &num_of_vectors_in_Xs,
                                  ComplexMatrix &eigen_values, ComplexMatrix &eigen_vectors, ComplexMatrix &eigen_vectors_after_swap, 
                                  ComplexMatrix &Eig_iter ){

Ziplus1 = a*Xi - Xi*Mi - Ximinus1*transpose(Ri, p, p);

octave::math::qr<ComplexMatrix> qrObject( Ziplus1 , octave::math::qr<ComplexMatrix>::economy);
X2 = qrObject.Q();
R2 = qrObject.R();

X1 = Xi;


// std::cout << "  BLOCK X2 is: " << std::endl;
// std::cout << X2 << std::endl;
// std::cout << std::endl;

// std::cout << " Upper triangular matrix R2 is " << std::endl;
// std::cout << R2 << std::endl;
// std::cout << std::endl;

/****************************************************************************************************************************************/
/* QR Test2 */
// for(int i=0;i<p;i++){

// std::cout << "Norm of column number " << i << " for matrix q_test is " << NormComplex(column_vector(X2,n,i),n)  << std::endl;

// }
// std::cout << std::endl;

/* QR Test3 */
// std::cout << "All Dotproducts: " << std::endl;

// for(int i=0;i< n -1;i++){
//         for(int j=i+1;j< p ;j++){

// std::cout << DotComplex( column_vector(X2,n,i) ,  column_vector(X2,n,j) , n )  << "    " ;
//         }
// }
// std::cout << std::endl;
/****************************************************************************************************************************************/

/* Computed X1 and X2 blocks */
/**********************************************************************************************************************************************/
// Ortho_check_2_blocks(X1, X2, n, p);

number_of_blocks_in_Xs  +=1;
num_of_vectors_in_Xs += p;
Concatenate(Xs,X2,n,p,number_of_blocks_in_Xs, num_of_vectors_in_Xs);

// std::cout << "  Xs after concatenating X2 is " << std::endl;
// std::cout << Xs << std::endl;
// std::cout << std::endl;

// std::cout << "All Dotproducts: " << std::endl;

// for(int i=0;i< num_of_vectors_in_Xs -1;i++){
//         for(int j=i+1;j< num_of_vectors_in_Xs ;j++){

                // std::cout << "(i,j) is    ( "  << i << ", " << j << " )" << std::endl;


// std::cout << DotComplex( column_vector(Xs,n,i) ,  column_vector(Xs,n,j) , n )  << "    " ;
//         }
// }
// std::cout << std::endl;



   ComplexMatrix Ms_full = transpose(Xs, n, num_of_vectors_in_Xs)*a*Xs;
        
     ComplexMatrix Ms(num_of_vectors_in_Xs, num_of_vectors_in_Xs);


    for(int i=0;i<num_of_vectors_in_Xs;i++){

        Ms(i,i)  =  Ms_full(i,i);

    }

        for(int i=0;i<num_of_vectors_in_Xs-1;i++){

        Ms(i,i+1)=  Ms_full(i,i+1);
        Ms(i+1,i)=  Ms_full(i+1,i);
            } 


        for(int i=0;i<num_of_vectors_in_Xs-2;i++){

        Ms(i,i+2)=  Ms_full(i,i+2);
        Ms(i+2,i)=  Ms_full(i+2,i);
            } 




// std::cout << std::endl;
// std::cout << " Block tridiagonal matrix Ms is " << std::endl;
// std::cout << Ms << std::endl;

/*********************************************************************************************************************************************/


/*********************************************************************************************************************************************/

        EIG eig2;
        eig2 = EIG(Ms);
        eigen_values = eig2.eigenvalues();
        eigen_vectors = eig2.right_eigenvectors();

        Eig_iter = Xs*eigen_vectors; //Very Important Step

        eigen_vectors_after_swap = ComplexMatrix( n, num_of_vectors_in_Xs );//IMPORTANT : Actually, the size increases as add n_add vectors.
    
    function_sort_descending_absolute_magnitude( n ,  num_of_vectors_in_Xs , eigen_values, Eig_iter, eigen_vectors_after_swap );

    // std::cout << "eigenvalues for A in Block LANCZOS " << std::endl;
    // std::cout << eigen_values << std::endl;

    // std::cout << " eigen_vectors_after_swap for A in LANCZOS " << std::endl;
    // std::cout << eigen_vectors_after_swap << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
/*********************************************************************************************************************************************/



/*********************************************/
Xi = X2;
Ximinus1 = X1;
Mi =  transpose(Xi,n,p)*a*Xi;
Ri = R2;
// /*********************************************/




   }


