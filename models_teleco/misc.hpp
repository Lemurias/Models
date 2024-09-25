#pragma once

#ifndef _MISC_H_
#define _MISC_

#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <complex>
#include <random>
#include <iterator>
#include <algorithm>
#include <numeric>  // For std::accumulate

// Define a type alias for convenience
using Complex = std::complex<double>;
using DoubleVector = std::vector<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<std::vector<Complex>>;
using DoubleMatrix = std::vector<std::vector<double>>;


/************************************* .
* Modulation:
Options: 'SP_BPSK', 'DP_BPSK', 'SP_QPSK', 'DP_QPSK'
* ************************************/
typedef enum
{
	SP_BPSK,
	DP_BPSK,
	SP_QPSK,
	DP_QPSK,
	PAM_4,
	NRZ
}modulation_t;

/*************************************
*  Units
*************************************/
double const	THZ = 1e12;
double const	GHZ = 1e9;
double const	MHZ = 1e6;
double const	kHZ = 1e3;
double const	km = 1e3;
double const	nm = 1e-9;
double const	us = 1e-6;
double const	ns = 1e-9;
double const	mW = 1e-3;
double const	nA = 1e-9;
double const	ps = 1e-12;

/*************************************/

/*************************************
* 		   Constants
**************************************/

const double pi = 3.14159265358979323846;	// pi value
const double q	= 1.602E-19;				// electron charge
const double kb = 1.3806504E-23;			// physconst('Boltzmann') [J/K]
const double c  = 299792458;				// Speed of light (ITU)

/*************************************.
 * Viewer parameters
*************************************/
typedef struct
{
	bool    displayFigures;         // Show figures (=1 yes / =0 no)
	bool    displayConstelation;    // Show constelation (=1 yes / =0 no)
	bool    save;                   // Save? (=1 yes / =0 no)
}parameters_type;

//typedef struct
//{
//	double              bits_per_symbol;
//	double              baudrate;
//	modulation_t        modulation;
//	double              bitrate;
//	std::uint16_t		n_symbols;
//	std::uint32_t       n_bits;
//	std::uint32_t       n_samples_per_bit;
//	std::uint32_t	    n_samples;				// Total number of samples
//	std::uint16_t	    n_samples_per_sym;		// Samples per symbol
//}grl_prm;

void saveMatrixInFile(const std::vector<std::vector<double>>& matrix, const std::string& fileName);
void saveMatrixInFile(const std::vector<std::vector<Complex>>& matrix, const std::string& fileName);
void saveMatrixInFile(std::unique_ptr<std::unique_ptr<Complex[]>[]>& matrix, std::uint32_t rows, std::uint32_t cols, const std::string& fileName);
void saveMatrixInFile(std::unique_ptr<std::unique_ptr<double[]>[]>& matrix, std::uint32_t rows, std::uint32_t cols, const std::string& fileName);
void saveMatrixInFile(std::unique_ptr<double[]>& vector, std::uint32_t cols, const std::string& fileName);
void saveMatrixInFile(const std::vector<double> &vect, const std::string& fileName);

/* ======================================================================
				Matlab functions to C++ functions
   ======================================================================= */
   // The mod function in MATLAB returns the remainder after division of x by y.
uint32_t mod(uint32_t x, uint32_t y);

// Function to generate normally distributed random numbers in an vector
std::vector<double> randn_vec(int size, double mean, double stddev);

// Function to generate normally distributed random numbers in an m x n matrix
std::vector<std::vector<double>> randn(int m, int n, double mean, double stddev);

// Function to compute the cumulative sum of a 2D array along rows
std::vector<std::vector<double>> cumsumRows(const std::vector<std::vector<double>>& A);

// Function to compute the cumulative sum of a 2D array along columns
std::vector<std::vector<double>> cumsumCols(const std::vector<std::vector<double>>& A);

std::vector<double> cumsum(const std::vector<double>& X);

// Function to perform element-wise multiplication of a scalar with a 2D array
std::vector<std::vector<double>> elementwise_multiply(const std::vector<std::vector<double>>& A, double scalar);

// Function to create an m x n matrix filled with ones
std::vector<std::vector<double>> ones(int m, int n);

// Function to create an m x n matrix filled with zeros
std::vector<std::vector<double>> zeros(int m, int n);

// Function to perform element-wise multiplication of two 2D vectors
std::vector<std::vector<double>> elementwise_multiply(const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2);

// Function to compute the element-wise complex exponential of a vector
std::vector<Complex> complex_exponential(const std::vector<double>& vect);

// Function to compute the element-wise complex exponential of a 2D vector
std::vector<std::vector<Complex>> complex_exponential(const std::vector<std::vector<double>>& mat);

// Function to perform the final operation
std::vector<std::vector<Complex>> compute_E_laser(const std::vector<std::vector<Complex>>& E_l, const std::vector<std::vector<double>>& phi_tx, const std::vector<std::vector<double>>& t);

// Function to create a 2D rotation matrix
std::vector<std::vector<double>> create_rotation_matrix(double theta_pol);

std::vector<double> sign(std::vector<double>& matx);

std::vector<std::vector<double>> kron(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);

// Function to compute the bit reversal for the FFT algorithm 
//void bit_reversal(ComplexVector& a);

// Function to compute the FFT algorithm
void fft(ComplexVector& a);

// Function to compute the IFFT algorithm
void ifft(ComplexVector& a);

// Function to compute the FFT algorithm optimized
void fft_o(ComplexVector& a);

// Function to compute the IFFT algorithm optimized
void ifft_o(ComplexVector& a);

// Function to convert a double vector to complex vector
ComplexVector convertToComplex(const DoubleVector& realVector);

// Function to print a std::complex<double>
void printComplexVector(const ComplexVector& vec);

// Function to multiply double matrix
DoubleMatrix multiply_d_matrices(const DoubleMatrix& A, const DoubleMatrix& B);

// Function to multiply complex matrix
ComplexMatrix multiply_matrices(const ComplexMatrix& A, const ComplexMatrix& B);

// multiply a double matrix by a complex matrix
ComplexMatrix multiply_matrices(const DoubleMatrix& A, const ComplexMatrix& B);

// Function to create a 2xN matrix from a complex vector
ComplexMatrix create_matrix(const ComplexVector& vec);

// Function to create a 2xN matrix from a double vector
ComplexMatrix convert_to_complex_matrix(const DoubleMatrix& double_matrix);

// Function to compute the FFTShift function
std::vector<double> fftshift(const std::vector<double>& input);

// Function to compute the FFTShift function
void fft_shift(std::vector<Complex>& input);

// Function to compute the FFTShift function
void fft_shift(std::vector<double>& input);

// Optimized function to perform ifftshift on a 1D array
std::vector<double> ifftshift(const std::vector<double>& data);

// Optimized function to perform ifftshift on a 1D array
std::vector<Complex> ifftshift(const std::vector<Complex>& data);

// Optimized function to perform ifftshift on a 2D matrix
std::vector<std::vector<double>> ifftshift(const std::vector<std::vector<double>>& data);

// Function to perform element-wise division of a 2D matrix by a scalar
std::vector<std::vector<double>> elementWiseDivision(const std::vector<std::vector<double>>& matrix, double scalar);

double variance(const std::vector<double>& data, bool isSample = false);

#endif