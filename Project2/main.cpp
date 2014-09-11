#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

//vec tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n);
vec create_potential(int j, vec &rho);
mat set_matrix(int n, vec &a, vec &b, vec &c);
void create_vectors(int n, double h, vec &a, vec &b, vec &c, vec &f, vec &u_real, vec &potential);
void jacobis_method(int n, mat &A, mat &R);
void rotate(int n, int &k, int &l, mat &A, mat &R, double &max_off_diag);
double maxoffdiagonal(int n, int &k, int &l, mat &A);

int main()
{
    int n_exponent, n, j;
    double start, finish, operation_time, relative_error, h, rho_min, rho_max, n_step,
           start_arma, finish_arma, operation_time_arma, relative_error_arma;


    // Creating a file to write the time and error to
    ofstream Amatrix;
    Amatrix.open("A.dat");
    //Amatrix.setf(ios::scientific);         // Forcing scientific notation
    /*
    // Creating a file to write the time and error to for armadillo solving
    ofstream timeanderror_arma;
    timeanderror_arma.open("timeanderror_arma.dat");
    timeanderror_arma.setf(ios::scientific);         // Forcing scientific notation

    // Creating a file to write the found function values to
    ofstream uvalues;
    uvalues.open("uvalues_8nr4.dat");
    uvalues.setf(ios::scientific);              // Forcing scientific notation
    */

    // Creating for-loop to try the solving mechanism for different size matrices
    for(n_exponent=1; n_exponent < 2; n_exponent++){
        // vectors of length n-1: Set n in program to n=n-1

        n = pow(10,n_exponent);                 // Dimension of matrix
        n_step = 100;
        vec a, b, c, f, u_real, u_arma;         // Please see sub functions for explanation
        vec u = zeros(n+2);                     // Vector to store the solution of the problem
        vec u_short = zeros(n);
        vec u_real_short = zeros(n);
        vec u_arma_short= zeros(n);
        vec rho = zeros(n+2);
        vec i = linspace(0, n_step, n+2);
        vec potential = zeros(n+2);
        mat A;
        double max_off_diag;

        rho_min = 0;
        rho_max = 10;
        h = (rho_max-rho_min)/n_step;
        rho = rho_min*i*h;

        potential = create_potential(j, rho);

        create_vectors(n, h, a, b, c, f, u_real, potential);  // Creating vectors for the tridiagonal matrix system

        A = set_matrix(n,a,b,c);

        A.save("A.txt", raw_ascii); // For checking the matrix
        mat R(n,n);
        R.zeros();

        jacobis_method(n, A, R);


        //start = clock();                        // Starting clock
        /*u = tridiagonal_matrix(a, b, c, f, n);  // Solving the tridiagonal matrix problem
        finish = clock();                       // Ending clock
        operation_time = (finish - start)/(double) CLOCKS_PER_SEC;  // Calculating time in seconds*/

        //u_short = u.subvec(1,n);                // Defining u-vector with non-zero elements
        //u_real_short = u_real.subvec(1,n);      // Defining u_real-vector with non-zero elements
        //relative_error = log10(max(abs((u_real_short-u_short)/u_real_short)));  // Calculating relative error
        /*
        // Solving problem using LU-decomposition, and printing u-vector to file, if the matrix is small
        if(n_exponent < 4){
            start_arma = clock();
            armadillo_solve(n, a, b, c, f, u_arma);
            finish_arma = clock();
            operation_time_arma = (finish_arma - start_arma)/(double) CLOCKS_PER_SEC;

            u_arma_short = u_arma.subvec(0,n-1);  // Defining u_arma-vector with non-zero elements
            //relative_error_arma = log10(max(abs((u_real_short-u_arma_short)/u_real_short)));

            // Writing time and errors for armadillo solutions to file
            timeanderror_arma << setiosflags(ios::showpoint | ios:: uppercase);
            timeanderror_arma << setw(10) << setprecision(8) << n_exponent << " "
                              << setw(10) << setprecision(8) << log10(max(abs((u_real_short-u_arma_short)/u_real_short))) << " "
                              << setw(10) << setprecision(8) << operation_time_arma << endl;

            // Writing found function values to file
            uvalues << setiosflags(ios::showpoint | ios:: uppercase);
            uvalues << setw(10) << setprecision(8) << u << endl;
        }*/

        // Writing error and time to file

        /*
        if(n_exponent == 1){                    // Setting information about the columns in the first row
            timeanderror << setiosflags(ios::showpoint | ios:: uppercase);
            timeanderror << setw(10) << setprecision(8) << "exponent"
                         << setw(10) << setprecision(8) << "log10 rel.error"
                         << setw(10) << setprecision(8) << "time" << endl;
        }

        timeanderror << setiosflags(ios::showpoint | ios:: uppercase);
        timeanderror << setw(10) << setprecision(8) << n_exponent << " "
                     << setw(10) << setprecision(8) << log10(max(abs((u_real_short-u_short)/u_real_short))) << " "
                     << setw(10) << setprecision(8) << operation_time << endl;*/


    } // End for loop over n values

    // Closing files
    /*timeanderror.close();
    uvalues.close();*/

    return 0;

} // End of program

vec create_potential(int j, vec &rho)
{
    int n = rho.size();
    vec potential = vec(n);

    for(j = 0; j < n; j++)
    {
        potential[j] = rho[j]*rho[j];
    }

    return potential;
}

void create_vectors(int n, double h, vec &a, vec &b, vec &c, vec &f, vec &u_real, vec &potential){

    // Resizing the vector to size needed for the current for-loop
    a.resize(n+2);              // Vector storing elements a_ij (j = i-1) for tridiagonal matrix A
    b.resize(n+2);              // Vector storing elements a_ii           for tridiagonal matrix A
    c.resize(n+2);              // Vector storing elements a_ij (j = i+1) for tridiagonal matrix A
    f.resize(n+2);              // Vector storing elements in vector f in the problem Au = f
    u_real.resize(n+2);         // Vector storing the values of the actual solution to the problem

    // Filling all elements of the vectors with zeros
    a.zeros();
    b.zeros();
    c.zeros();
    f.zeros();
    u_real.zeros();

    // Assign values to the elements of a, b, c
    for(int i=1; i<=n+1; i++){

        double x = i*h;         // Variable of functions

        if(i == 1){
            b[i] = 2+ h*h*potential[i];
            c[i] = -1;
            f[i] = h*h*100*exp(-10*x);
            u_real[i] = 1-(1-exp(-10))*x - exp(-10*x);
        }else if(i == n+1){
            a[i] = -1;
            b[i] = 2 + h*h*potential[i];
            f[i] = h*h*100*exp(-10*x);
            u_real[i] = 1-(1-exp(-10))*x - exp(-10*x);
        }else{
            a[i] = -1;
            b[i] = 2 + h*h*potential[i];
            c[i] = -1;
            f[i] = h*h*100*exp(-10*x);
            u_real[i] = 1-(1-exp(-10))*x - exp(-10*x);
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

void jacobis_method(int n, mat &A, mat &R){
    double epsilon = pow(10,-8);
    double max_off_diag;
    int k,l;

    // Setting up the eigenvector matrix
    for (int i=0; i<n; i++){
        R(i,i) = 1.0;
    }

    max_off_diag = maxoffdiagonal(n, k, l, A);

    while (fabs(max_off_diag) > epsilon ){
        // Could add a limit on the number of iterations in while-statement
        max_off_diag = max_off_diag = maxoffdiagonal(n, k, l, A);
        rotate(n, k, l, A, R, max_off_diag);
    }


}

double maxoffdiagonal(int n, int &k, int &l, mat &A){
    // In his lecture notes, MHJ said that this could be done more elegantly
    // This is pretty much what he did :/
    double maxval = 0.0;
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if (fabs(A(i,j)) > maxval){
                maxval = fabs(A(i,j));
                l = i;
                k = l;
            } // End if-test
        } // End j-loop
    } // End i-loop
    return maxval;
} // End maxoffdiagonal-function

void rotate(int n, int &k, int &l, mat &A, mat &R, double &max_off_diag){
    double c, s, t, tau, i;
    double all, alk, akk, ail, aik, rik, ril;
    if (max_off_diag !=0){
        tau = (A(l,l)-A(k,k))/(2*A(k,l));
        if (tau < 0){
            t = -1.0/(-tau+sqrt(1+tau*tau));}
        else{
            t = 1.0/(tau+sqrt(1+tau*tau));}
        c = 1/(1+t*t);
        s =t*c;


        }
    else{
        c = 1.0;
        s = 0.0;
    }
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
        }

    // New eigenvectors
    rik = R(i,k);
    ril = R(i,l);
    R(i,k) = c*rik - s*ril;
    R(i,l) = c*ril + s*rik;
    }
    return;
}

