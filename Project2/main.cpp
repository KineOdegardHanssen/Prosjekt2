#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

// Functions, in the order of which they are used
mat     set_matrix(int n, double h, vec &potential);
void    create_rho_and_potential(int n, double h, double rho_min, vec &rho, vec &potential);
void    jacobis_method(int n, mat &A, mat &R);
double  maxoffdiagonal(int n, int &k, int &l, mat &A);
void    rotate(int n, int &k, int &l, mat &A, mat &R, double &max_off_diag);

int main()
{
    double h, rho_min, rho_max;

    int n_step = 100;
    int n = n_step;
    mat A;
    vec rho = zeros(n+2);                       // Vector saving the dimentionless radial positions
    vec potential = zeros(n+2);                 // Vector saving the potential
    mat R = eye(n,n);                           // Creating matrix R, which stores the eigenvectors of A


    // Calling functions to create matrix
    rho_min = 0.;
    rho_max = 10;
    h = (rho_max-rho_min)/(n_step+1);
    create_rho_and_potential(n, h, rho_min, rho, potential);          // Finding the potential at each point
    A = set_matrix(n, h, potential);
    A.save("A.txt", raw_ascii);                 // Saving matrix to be able to have a look at it

    // Using Jacobi's method to solve the eigenvalue problem
    jacobis_method(n, A, R);

    return 0;

} // End of program

void create_rho_and_potential(int n, double h, double rho_min,  vec &rho, vec &potential){
    for(int i = 1; i < n+2; i++){
        rho[i] = rho_min + i*h;
        potential[i] = rho[i]*rho[i];
    }   // End for-loop
}   // End create_potential-function

mat set_matrix(int n, double h, vec &potential){
    mat A(n,n);
    A.zeros();
    double offDiagonalValue = -1.0/h/h;
    A.diag(-1).fill(offDiagonalValue);
    A.diag(1).fill(offDiagonalValue);

    for(int i = 0; i<n; i++){
        A(i,i)   = 2/h/h + potential[i+1];
    }   // End for-loop

    return A;

}   // End set_matrix-function

void jacobis_method(int n, mat &A, mat &R){
    int k, l;
    int iteration_number = 0;
    double max_iteration_number = n*n;    // We stop diagonalising the matrix after this to not kill something
    double max_off_diag;
    double epsilon = pow(10,-8);            // If the maximum off-diagonal value is smaller than this,
                                            // we consider the matrix diagonalised

    max_off_diag = maxoffdiagonal(n, k, l, A);

    while (fabs(max_off_diag) > epsilon && iteration_number <= max_iteration_number){
        max_off_diag = maxoffdiagonal(n, k, l, A);
        rotate(n, k, l, A, R, max_off_diag);
        iteration_number++;
    }   // End while-loop

    // Print out error in case matrix doesn't get diagonalised, and the solution went to hell
    if(fabs(max_off_diag) > epsilon && iteration_number == max_iteration_number){
        cout << "Matrix was not succesfully diagonalised." << endl;
    }   // End if-statement
}   // End jacobis_method-function

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


