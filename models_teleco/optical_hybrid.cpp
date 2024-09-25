#include "optical_hybrid.hpp"

void optical_hybrid(ComplexMatrix& esig, ComplexMatrix& elo, receiver_params_t& rx_params, ComplexMatrix& E1, ComplexMatrix& E2, ComplexMatrix& E3, ComplexMatrix& E4)
{
	size_t rows_A;
	size_t cols_A;
	size_t rows_B;
	size_t cols_B;
	
	ComplexMatrix esig_slow(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));
	ComplexMatrix esig_fast(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));;
	ComplexMatrix elo_slow(2, std::vector<Complex>(elo[0].size(), Complex(0.0, 0.0)));;
	ComplexMatrix elo_fast(2, std::vector<Complex>(elo[0].size(), Complex(0.0, 0.0)));;

	for (size_t j = 0; j < 2; j++) {
		for (size_t k = 0; k < esig[0].size(); k++) {
			esig_slow[0][k] = esig[0][k];		// first coupler input in slow axis pol
			esig_slow[1][k] = 0.0,				// second coupler input in slow axis pol

			esig_fast[0][k] = esig[1][k];		// first coupler input in fast axis pol
			esig_fast[1][k] = 0.0;				// second coupler input in fast axis pol
		}
	}

	for (size_t j = 0; j < 2; j++) {
		for (size_t k = 0; k < elo[0].size(); k++) {
			elo_slow[1][k] = elo[0][k];			// first coupler input in slow axis pol
			elo_slow[0][k] = 0.0;				// second coupler input in slow axis pol

			elo_fast[1][k] = elo[1][k];		// first coupler input in fast axis pol
			elo_fast[0][k] = 0.0;				// second coupler input in fast axis pol
		}
	}

	/************************************************************************************/
	
	double p = -1;								// p = -1 invert the conjugated phase
	double a_dB = rx_params.coup_excess_loss_dB;// Excess loss of couplers
	double a = pow(10, -a_dB / 20);

	/************************************************************************************/
	ComplexMatrix E_A_slow(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));
	ComplexMatrix E_A_fast(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));

	double aux = a / sqrt(2);

	ComplexMatrix aux_matrix(2, std::vector<Complex>(2));
	aux_matrix[0][0] = aux;
	aux_matrix[0][1] = aux * Complex(0, p);
	aux_matrix[1][0] = aux * Complex(0, p);
	aux_matrix[1][1] = aux;

	
	rows_A = aux_matrix.size();
	cols_A = aux_matrix[0].size();
	rows_B = esig_slow.size();
	cols_B = esig_slow[0].size();

	// check if it´s posible to multiply 
	if (cols_A != rows_B) {
		throw std::invalid_argument("Number of columns in A must be equal to the number of rows in B.");
	}

	// Perform matrix multiplication
	for (size_t j = 0; j < rows_A; j++) {
		for (size_t k = 0; k < cols_B; k++) {
			for (size_t l = 0; l < cols_A; l++) {
				/* First coupler for signal path */
				E_A_slow[j][k] += aux_matrix[j][l] * esig_slow[l][k];
				E_A_fast[j][k] += aux_matrix[j][l] * esig_fast[l][k];
			}
		}
	}

	ComplexMatrix E_C_slow(2, std::vector<Complex>(elo[0].size(), Complex(0.0, 0.0))); 
	ComplexMatrix E_C_fast(2, std::vector<Complex>(elo[0].size(), Complex(0.0, 0.0)));

	rows_B = elo_slow.size();
	cols_B = elo_slow[0].size();

	// check if it´s posible to multiply 
	if (cols_A != rows_B) {
		throw std::invalid_argument("Number of columns in A must be equal to the number of rows in B.");
	}

	// Perform matrix multiplication
	for (size_t j = 0; j < rows_A; ++j) {
		for (size_t k = 0; k < cols_B; ++k) {
			for (size_t l = 0; l < cols_A; ++l) {
				/* First coupler for LO path */
				E_C_slow[j][k] += aux_matrix[j][l] * elo_slow[l][k];
				E_C_fast[j][k] += aux_matrix[j][l] * elo_fast[l][k];
			}
		}
	}

	/* Second coupler for the upper path */
	ComplexMatrix E1_slow(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));
	ComplexMatrix E1_fast(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));

	// Perform matrix multiplication
	for (size_t j = 0; j < rows_A; ++j) {
		for (size_t k = 0; k < cols_B; ++k) {
				E1_slow[j][k] = aux_matrix[j][0] * E_A_slow[0][k] + aux_matrix[j][1] * E_C_slow[0][k];
				E1_fast[j][k] = aux_matrix[j][0] * E_A_fast[0][k] + aux_matrix[j][1] * E_C_fast[0][k];
		}
	}
		
	/* Second coupler for the lower path */
	double	phase = pi / 2;
	ComplexMatrix E3_slow(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));
	ComplexMatrix E3_fast(2, std::vector<Complex>(esig[0].size(), Complex(0.0, 0.0)));

	// Perform matrix multiplication
	for (size_t j = 0; j < rows_A; ++j) {
		for (size_t k = 0; k < cols_B; ++k) {
			E3_slow[j][k] = aux_matrix[j][0] * E_A_slow[1][k] + aux_matrix[j][1] * E_C_slow[1][k] * std::exp(Complex(0, phase)); 
			E3_fast[j][k] = aux_matrix[j][0] * E_A_fast[1][k] + aux_matrix[j][1] * E_C_fast[1][k] * std::exp(Complex(0, phase));
		}
	}

	/* Jones Vector of the hybrid's output */
	E1.resize(esig_slow.size());
	E2.resize(esig_slow.size());
	E3.resize(esig_slow.size());
	E4.resize(esig_slow.size());
	
	for (size_t k = 0; k < esig_slow.size(); k++) {
		E1[k].resize(esig_slow[0].size());
		E2[k].resize(esig_slow[0].size());
		E3[k].resize(esig_slow[0].size());
		E4[k].resize(esig_slow[0].size());
	}

	for (size_t k = 0; k < cols_B; ++k) {
		E1[0][k] = E1_slow[0][k];
		E1[1][k] = E1_fast[0][k];
		E2[0][k] = E1_slow[1][k];
		E2[1][k] = E1_fast[1][k];
		E3[0][k] = E3_slow[0][k];
		E3[1][k] = E3_fast[0][k];
		E4[0][k] = E3_slow[1][k];
		E4[1][k] = E3_fast[1][k];
	}
}