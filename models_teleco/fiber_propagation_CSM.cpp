#include "fiber_propagation_CSM.hpp"

#include <random>
#include <cmath>

void fiber_propagation_CSM(ComplexMatrix& signal_in, channel_params_t& ch_prms, modulation_setup_t& mod_stp, laser_params_t& laser_tx)
{
	double alpha = ch_prms.alpha_dB / (10 * std::log10(std::exp(1)));			// [1/km]
	double D = ch_prms.D * ps / (nm * km);
	double beta2 = -(std::pow(laser_tx.Lambda, 2)) / (2 * pi * c) * D * 1000;	// [ps^2/km]
	
	/************************************************
	* ---- (only for run quick tests) ----
	* 
	************************************************/ 
	ch_prms.length = 10;
	ch_prms.dz = ch_prms.length / 2;
	alpha = 0;
	beta2 = 0 * 2.6E-24;
	ch_prms.beta3 = 0.0;
	ch_prms.pmd = 3.75 * 0.5;					// [ps]
	ch_prms.sopmd = 9.4308 * 1 / sqrt(10);		// [ps ^ 2]

	std::vector<double> w(mod_stp.size_of_f);

	double two_pi = 2 * pi;
	for (size_t i = 0; i < mod_stp.size_of_f; i++){
		w[i] = two_pi * mod_stp.f[i];			// Angular frequency axis
	}

	w = fftshift(w);							// Required fftshift (!!!)

	// Initialize output
	ComplexVector signalX_fft(signal_in[0].size());
	ComplexVector signalY_fft(signal_in[1].size());
	
	for (size_t i = 0; i < mod_stp.size_of_f; i++) {
		signalX_fft[i] = signal_in[0][i];
		signalY_fft[i] = signal_in[1][i];
	}
	
	fft_o(signalX_fft);
	fft_o(signalY_fft);

	// Initialize fiber transfer amtrix
	std::vector<double> M11(mod_stp.n_samples, 1.0);
	std::vector<double> M12(mod_stp.n_samples, 0.0);
	std::vector<double> M21(mod_stp.n_samples, 0.0);
	std::vector<double> M22(mod_stp.n_samples, 1.0);

	// Dispersion and attenuation effects
	ComplexVector D_Att_GVD(mod_stp.n_samples);

	for (size_t i = 0; i < mod_stp.size_of_f; i++) {
		D_Att_GVD[i] = -alpha / 2 * ch_prms.length - Complex(0.0, (beta2 / 2) * ch_prms.length * std::pow(w[i], 2)) - Complex(0.0, (ch_prms.beta3 / 6) * ch_prms.length * std::pow(w[i], 3));
	}
	
	// DGD and SOPMD of each propagation step
	double dgd_dz	= ch_prms.pmd * std::sqrt(ch_prms.dz) * ps;			// [s]
	double sodgd_dz = ch_prms.sopmd * std::sqrt(ch_prms.dz) * ps * ps;	// [s ^ 2]

	std::cout << "------------------------------------------" << '\n';
	std::cout << "Number of simulation fiber segments = " << ch_prms.length / ch_prms.dz << '\n';
	std::cout << "DGD per segment = " << dgd_dz * 1E12 << " [ps]" << '\n';
	std::cout << "SODGD per segment = " << sodgd_dz * 1E24 << " [ps^2]" << '\n';

	// Obtaining fiber transfer matrix
	size_t i_dz = static_cast<size_t>(ch_prms.dz);
	size_t i_L_minus_dz = static_cast<size_t>(ch_prms.length - ch_prms.dz);

	std::vector<Complex> Dx(mod_stp.n_samples);
	std::vector<Complex> Dy(mod_stp.n_samples);

	// Initialize random number generator
	std::random_device rd;							// Obtain a random number from hardware
	std::mt19937 gen(rd());							// Seed the generator
	std::uniform_real_distribution<> dis(0.0, 1.0); // Define the range [0, 1)

	std::vector<Complex> M_dz11(mod_stp.n_samples);
	std::vector<Complex> M_dz12(mod_stp.n_samples);
	std::vector<Complex> M_dz21(mod_stp.n_samples);
	std::vector<Complex> M_dz22(mod_stp.n_samples);

	for (size_t i = 0; i < i_L_minus_dz; i = i + i_dz) {
		// DGD and SOPMD

		for (size_t j = 0; j < mod_stp.n_samples; j++) {
			Dx[j] = Complex(0.0, 0.5 * dgd_dz * w[j] + 0.25 * sodgd_dz * w[j]);
			Dy[j] = Complex(0.0, -0.5 * dgd_dz * w[j] - 0.25 * sodgd_dz * w[j]);
		}

		// Generate a random unit vector for the PMD direction
		// Generate a random angle in radians
		double theta = 2 * pi * dis(gen);	// Random angle
		double phi	 = 2 * pi * dis(gen);	// Random phase

		// Segment matrix transfer function (dz)
		/**********************************
		* 
		* M_dz = [M_dz11    M_dz12]
		*        [M_dz21    M_dz22]
		* 
		* ********************************/
		for (size_t j = 0; j < mod_stp.n_samples; j++) {
			M_dz11[j] = std::exp(Dx[i]) * std::cos(theta) * std::exp(Complex(0.0, phi));
			M_dz12[j] = std::exp(Dx[i]) * std::sin(theta) * std::exp(Complex(0.0, phi));
			M_dz21[j] = - std::exp(Dy[i]) * std::sin(theta);
			M_dz11[j] = std::exp(Dy[i]) * std::cos(theta);
		}
	
		// New accumulated fiber transfer function
		// T = T * M_dz

	}







}
