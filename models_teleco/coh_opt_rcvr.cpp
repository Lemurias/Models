#include "misc.hpp"
#include "setup.hpp"
#include "coh_opt_rcvr.hpp"
#include "optical_hybrid.hpp"
#include "balanced_detector.hpp"
#include <array>
#include <cmath> 

/*****************************************************************************
* Coherent_optical_receiver
* Inputs: -Received optical electric field(Ein_rx)
* Local oscilator optical electric field(Elo_)
* Output : -Optical signals at the output of the coherent receiver signal(x8)
*****************************************************************************/
void coherent_optical_receiver(ComplexMatrix& ein_rx, ComplexMatrix& elo_rx, receiver_params_t& rcv_prms, modulation_setup_t& main_prms, ComplexMatrix& E1_8)
{
	uint32_t i, j;

	/****************************************************************************/
	/*  Pol Beam Splitter */

	double k_pbs = 0.5;							// power coupling factor
	double aux1 = std::sqrt(1 - k_pbs);
	double aux2 = std::sqrt(k_pbs);

	
	ComplexMatrix Ein_rx_up(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));
	ComplexMatrix Ein_rx_down(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));

	for (i = 0; i < 2; i++) {
		for (j = 0; j < main_prms.n_samples; j++) {
			// Rx signal to uppper 90� hybrid
			Ein_rx_up[i][j] = aux1 * ein_rx[i][j];
			// Rx signal to uppper 90� hybrid
			Ein_rx_down[i][j] = aux2 * ein_rx[i][j];
		}
	}
	
	double theta_slow = 0 * pi / 180;			// supose slow axis is aligned with theta = 0[degree]
	double theta_fast = theta_slow + pi / 2;	// orthogonal polarization component
	
	// Linear polarizer matrix 
	std::array<std::array<double, 2>, 2> lpm_slow {0.0};
	lpm_slow[0][0] = std::pow(std::cos(theta_slow), 2);
	lpm_slow[0][1] = std::sin(theta_slow) * std::cos(theta_slow);
	lpm_slow[1][0] = std::sin(theta_slow) * std::cos(theta_slow);
	lpm_slow[1][1] = std::pow(std::sin(theta_slow), 2);

	// Linear polarizer matrix 
	std::array<std::array<double, 2>, 2> lpm_fast{ 0.0 };
	lpm_fast[0][0] = std::pow(std::cos(theta_fast), 2);
	lpm_fast[0][1] = std::sin(theta_fast) * std::cos(theta_fast);
	lpm_fast[1][0] = std::sin(theta_fast) * std::cos(theta_fast);
	lpm_fast[1][1] = std::pow(std::sin(theta_fast), 2);

	
	ComplexMatrix Ein_after_pbs_slow(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));
	ComplexMatrix Ein_after_pbs_fast(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));

	// Ein_after_pbs_slow = LPM_slow * Ein_rx_up;
	// Ein_after_pbs_fast = LPM_fast * Ein_rx_down;
	// Perform matrix multiplication
	for (i = 0; i < 2; ++i) {
		for (j = 0; j < main_prms.n_samples; ++j) {
			for (size_t k = 0; k < 2; ++k) {
				Ein_after_pbs_slow[i][j] += lpm_slow[i][k] * Ein_rx_up[k][j];
				Ein_after_pbs_fast[i][j] += lpm_fast[i][k] * Ein_rx_down[k][j];
			}
		}
	}
	
	double theta_pol = 90 * pi / 180;

	std::array<std::array<double, 2>, 2> rot_pol;
	rot_pol[0][0] = std::cos(theta_pol);
	rot_pol[0][1] = -std::sin(theta_pol);
	rot_pol[1][0] = std::sin(theta_pol);
	rot_pol[1][1] = std::cos(theta_pol);
	
	ComplexMatrix Esig_up(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));
	ComplexMatrix Esig_down(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));
	for (i = 0; i < 2; ++i) {
		for (j = 0; j < main_prms.n_samples; ++j) {
			Esig_up[i][j] = Ein_after_pbs_slow[i][j]; // Signal field to upper 90� hybrid
		}
	}

	// rotate the polarization aligned with fast axis
	for (i = 0; i < 2; ++i) {
		for (j = 0; j < main_prms.n_samples; ++j) {
			for (size_t k = 0; k < 2; ++k) {
				Esig_down[i][j] += rot_pol[i][k] * Ein_after_pbs_fast[k][j];
			}
		}
	}
	/****************************************************************************
	 *  LO Beam Splitter 
	****************************************************************************/
	
	double k_bs = 0.5;	// power coupling factor
	aux1 = std::sqrt(1 - k_bs);
	aux2 = std::sqrt(k_bs);
	
	ComplexMatrix Elo_up(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));
	ComplexMatrix Elo_down(2, std::vector<Complex>(main_prms.n_samples, Complex(0.0, 0.0)));

	for (i = 0; i < 2; i++) {
		for (j = 0; j < main_prms.n_samples; j++) {
			// LO field to upper 90� hybrid
			Elo_up[i][j] = aux1 * elo_rx[i][j];
			// LO field to lower 90� hybrid
			Elo_down[i][j] = aux2 * elo_rx[i][j];
		}
	}
	
	/****************************************************************************/
	
	ComplexMatrix E1, E2, E3, E4, E5, E6, E7, E8;
	optical_hybrid(Esig_up, Elo_up, rcv_prms, E1, E2, E3, E4);
	optical_hybrid(Esig_down, Elo_down, rcv_prms, E5, E6, E7, E8);

	E1_8.resize(16);
	for (i = 0; i < 16; i++) {
		E1_8[i].resize(main_prms.n_samples);
	}
	
	for (size_t i = 0; i < main_prms.n_samples; i++) {
		E1_8[0][i] = E1[0][i];
		E1_8[1][i] = E1[1][i];
		E1_8[2][i] = E2[0][i];
		E1_8[3][i] = E2[1][i];
		E1_8[4][i] = E3[0][i];
		E1_8[5][i] = E3[1][i];
		E1_8[6][i] = E4[0][i];
		E1_8[7][i] = E4[1][i];
		E1_8[8][i] = E5[0][i];
		E1_8[9][i] = E5[1][i];
		E1_8[10][i] = E6[0][i];
		E1_8[11][i] = E6[1][i];
		E1_8[12][i] = E7[0][i];
		E1_8[13][i] = E7[1][i];
		E1_8[14][i] = E8[0][i];
		E1_8[15][i] = E8[1][i];
	}
}
