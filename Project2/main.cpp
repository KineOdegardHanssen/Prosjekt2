#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

// Functions, in the order of which they are used
void    create_vectors(int n, double h, vec &a, vec &b, vec &c, vec &potential);
//mat     set_matrix(int n, vec &a, vec &b, vec &c);
mat     set_matrix(int n, double h, vec &potential);
vec     create_potential(vec &rho);
void    jacobis_method(int n, mat &A, mat &R);
double  maxoffdiagonal(int n, int &k, int &l, mat &A);
void    rotate(int n, int &k, int &l, mat &A, mat &R, double &max_off_diag);

int main()
{
    int n_exponent, n;
    double h, rho_min, rho_max, n_step;

    // Creating for-loop to try the solving mechanism for different size matrices
    for(n_exponent=1; n_exponent < 2; n_exponent++){
        // vectors of length n-1: Set n in program to n=n-1

        n = pow(10,n_exponent);                     // Dimension of matrix
        n_step = 100;
        n = n_step;
        mat A;
        vec a, b, c;                                // Vectors giving the three diagnonals of the matrix A
        vec rho = zeros(n+2);                       // Vector saving the dimentionless radial positions
        vec potential = zeros(n+2);                 // Vector saving the potential
        vec i = linspace(0, n_step, n+2);           // Vector making calculating rho easier

        // Assigning values to vector rho
        rho_min = 0;
        rho_max = 10;
        h = (rho_max-rho_min)/n_step;
        rho = rho_min + i*h;

        // Calling functions to create matrix
        potential = create_potential(rho);          // Finding the potential at each point
        //create_vectors(n, h, a, b, c, potential);   // Creating vectors for the tridiagonal matrix system
        //A = set_matrix(n,a,b,c);                    // Creating matrix A
        A = set_matrix(n, h, potential);
        A.save("A.txt", raw_ascii);                 // Saving matrix to be able to have a look at it

        // Creating matrix R, which stores the eigenvectors of A
        mat R(n,n);
        R.zeros();

        // Using Jacobi's method to solve the eigenvalue problem
        jacobis_method(n, A, R);

    } // End for loop over n values

    return 0;

} // End of program

vec create_potential(vec &rho)
{
    int n = rho.size();
    vec potential(n);

    for(int j = 0; j < n; j++)
    {
        potential[j] = rho[j]*rho[j];
    }

    return potential;
}   // End create_potential-function

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// PLEASE DELETE THESE FUNCTIONS IF YOU DECIDE NOT TO USE THEM
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
void create_vectors(int n, double h, vec &a, vec &b, vec &c, vec &potential){

    // Resizing the vector to size needed for the current for-loop
    a.resize(n+2);              // Vector storing elements a_ij (j = i-1) for tridiagonal matrix A
    b.resize(n+2);              // Vector storing elements a_ii           for tridiagonal matrix A
    c.resize(n+2);              // Vector storing elements a_ij (j = i+1) for tridiagonal matrix A

    // Filling all elements of the vectors with zeros
    a.zeros();
    b.zeros();
    c.zeros();

    // Assign values to the elements of a, b, c
    for(int i=1; i<=n+1; i++){
        if(i == 1){
            b[i] = 2 + h*h*potential[i];
            c[i] = -1;
        }else if(i == n+1){
            a[i] = -1;
            b[i] = 2 + h*h*potential[i];
        }else{
            a[i] = -1;
            b[i] = 2 + h*h*potential[i];
            c[i] = -1;
        } // End if statements
    } // End vector value assignment
} // End the create_vectors-function

mat set_matrix(int n, vec &a, vec &b, vec &c){

    // Can merge this function with create_vector when the code is running

    mat A(n,n);                 // Create full matris from tridiagonal matrix vectors
    A.zeros();                  // Fill matrix with zeros

    // Assign values to elements of matrix
    for (int i = 0; i<n; i++){
        if (i == 0){
            A(i,i)   = b[i+1];
            A(i,i+1) = c[i+1];
        }else if (i==n-1){
            A(i,i)   = b[i+1];
            A(i,i-1) = a[i+1];
        }else{
            A(i,i)   = b[i+1];
            A(i,i+1) = c[i+1];
            A(i,i-1) = a[i+1];
        } // End if statements
    } // End creating matrix
    return A;
}// End of function
*/

mat set_matrix(int n, double h, vec &potential){
    mat A(n,n);
    A.zeros();

    for(int i = 0; i<n; i++){
        if (i == 0){
            A(i,i)   = 2 + h*h*potential[i+1];
            A(i,i+1) = -1;
        }else if (i==n-1){
            A(i,i)   = 2 + h*h*potential[i+1];
            A(i,i-1) = -1;
        }else{
            A(i,i-1) = -1;
            A(i,i)   = 2 + h*h*potential[i+1];
            A(i,i+1) = -1;
        }   // End the various if-statements
    }   // End for-loop

    return A;

}   // End set_matrix-function

void jacobis_method(int n, mat &A, mat &R){
    double epsilon = pow(10,-8);
    double max_iteration_number = n*n*n;
    double max_off_diag;
    int k, l;
    int iteration_number = 0;

    // Setting up the eigenvector matrix
    for (int i=0; i<n; i++){
        R(i,i) = 1.0;
    }

    max_off_diag = maxoffdiagonal(n, k, l, A);

    while (fabs(max_off_diag) > epsilon && iteration_number <= max_iteration_number){
        max_off_diag = maxoffdiagonal(n, k, l, A);
        rotate(n, k, l, A, R, max_off_diag);
        iteration_number++;
    }
}

void rotate(int n, int &k, int &l, mat &A, mat &R, double &max_off_diag){
    int i;
    double c, s, t, tau;
    double all, alk, akk, ail, aik, rik, ril;

    // Finding values of the trigonomical functions
    if (max_off_diag !=0){
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if(tau < 0){
            t = -1.0/(-tau+sqrt(1+tau*tau));
        }else{
            t = 1.0/(tau+sqrt(1+tau*tau));
        }   // Ending if-statement
        c = 1/(1+t*t);
        s = t*c;
    }else{
        c = 1.0;
        s = 0.0;
    }   // Ending if-statement

    all = A(l,l);
    akk = A(k,k);
    alk = A(k,l);

    // Changing the matrix elements with indices k and l
    A(l,l) = all*c*c + 2*alk*c*s + akk*s*s;
    A(k,k) = akk*c*c - 2*alk*c*s + all*s*s;
    A(k,l) = 0;
    A(l,k) = 0;

    for (i=0; i<n; i++){
        if (i !=k && i !=l ){
            aik = A(i,k);
            ail = A(i,l);
            A(i,k) = aik*c - ail*s;
            A(i,l) = ail*c + aik*s;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);
        }   // Ending if-statement

        // New eigenvectors
        rik = R(i,k);
        ril = R(i,l);
        R(i,k) = c*rik - s*ril;
        R(i,l) = c*ril + s*rik;
    }   // Ending for-loop

    return;

}   // Ending rotate-function

double maxoffdiagonal(int n, int &k, int &l, mat &A){
    // In his lecture notes, MHJ said that this could be done more elegantly
    // This is pretty much what he did :/
    double maxval = 0.0;
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if(fabs(A(i,j)) > maxval){
                maxval = fabs(A(i,j));
                l = i;
                k = l;
            } // End if-test
        } // End j-loop
    } // End i-loop
    return maxval;
} // End maxoffdiagonal-function


