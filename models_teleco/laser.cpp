#ifndef _LASER_H_
#define _LASER_

#include "setup.hpp"
#include "laser.hpp"
#include "instruments.hpp"

/**********************************************************************************
 * Returns the complex electric field emited by the laser
 *
 * NOTE: This function works for the Tx and Rx laser (in the later case, the input
 * parameter must be Param_Rx instead of Param_Tx.
 * *******************************************************************************/                                                                                
//void compute_laser_out(modulation_setup_t& mod_stp, laser_params_t& laser_prm, std::unique_ptr<std::unique_ptr<Complex[]>[]>& elaser_out, int N, int M) // 2D dynamic array )
void compute_laser_out(modulation_setup_t& mod_stp, laser_params_t& laser_prm, ComplexMatrix& elaser_out, int N, int M) // 2D dynamic array )
{
    // index variables
    uint32_t i;

    double R = mod_stp.baudrate * GHZ;
    //uint32_t n_Symbols = get_nSymbols();

    uint32_t n_samples_tot = mod_stp.n_symbols * mod_stp.n_samples_per_sym;
    double P_avg = pow(10, ( laser_prm.P_dBm / 10.0)) * mW;
    double delta_nu = laser_prm.LW;                                           // Linewidth in Hz

    double sigma = sqrt(2 * pi * delta_nu * mod_stp.dt);
    double theta_pol = 0 * pi / 180;                                            // Polarization rotation angle[rad]

    // Create a random number generator and normal distribution
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::normal_distribution<> d(0.0, 1.0); // Normal distribution with specified mean and stddev

    //const uint32_t nTotal = n_samples_tot;
    //std::array<double, nTotal> np;   
    std::unique_ptr<double[]> np;    // unique pointer for dinamic array
    np = std::make_unique<double[]>(n_samples_tot);  // create an array of n_samples values

    //size_t              size_of_np;  // size of variable t

    for (i = 0; i < n_samples_tot; ++i) {
        np[i] = d(gen) * sigma;         // Generate a random number from normal distribution
    }

    std::unique_ptr<double[]> phi;    // unique pointer for dinamic array
    phi = std::make_unique<double[]>(n_samples_tot);  // create an array of n_samples values
    phi[0] = np[0];

    for (i = 1; i < n_samples_tot; ++i) {
        phi[i] = phi[i - 1] + np[i];
    }

    /* ======================================================================
     *                    Relative Intensity Noise
     * RIN noise is modeled as complex AWGN noise described by a PSD Param_Tx.Laser_RIN_dBHz.
     * RIN_dBHz is defined as one sided PSD
     * The RIN electric field is included in the same polarization axis than that of the
     * laser, which is assumed to emit linearly polarized light.
     *
     ======================================================================= */

    // Convert RIN from dB/Hz to power variance
    double RIN_lin = pow(10, laser_prm.RIN_dBHz / 10.0);   // [1/Hz]  //* pow(P_avg, 2);    // [W^2/Hz]
    double RIN_var = RIN_lin * P_avg / (2 * mod_stp.dt);
    double noise_power = sqrt(RIN_var/2);                   // [W]

    ComplexMatrix noise_field(2, std::vector<Complex>(n_samples_tot));
    //std::vector<Complex> noise_field(n_samples_tot, 0.0);
    //std::unique_ptr<Complex[]> noise_field;    // unique pointer for dinamic array
    //noise_field = std::make_unique<Complex[]>(n_samples_tot); // Crear array de punteros a filas
    for (uint32_t i = 0; i < n_samples_tot; ++i) {
        noise_field[0][i] = Complex(d(gen), d(gen)) * noise_power;   // Generate a random number from normal distribution
        noise_field[0][i] = Complex(0.0,0.0);
    }

    // Laser output
    double aux = sqrt(P_avg);

    for (uint32_t i = 0; i < n_samples_tot; ++i) {
        elaser_out[0][i] = (aux + noise_field[0][i]) * (std::exp(Complex(0, -2 * pi * laser_prm.FO * mod_stp.t[i] + phi[i])));
        elaser_out[1][i] = 0.0;
    }

    return;

    // Rot_pol = [cos(theta_pol) -sin(theta_pol);sin(theta_pol) cos(theta_pol)]; % rotation matrix
    // std::vector<std::vector<double>> Rot_pol = create_rotation_matrix(theta_pol);
    // ComplexMatrix rot_pol_cx = convert_to_complex_matrix(Rot_pol);
    // e_laser = multiply_matrices(rot_pol_cx, e_laser);

    double Power_meas_dBm = Power_meter(noise_field);
    DoubleVector intensity(elaser_out[0].size());

    for (uint32_t i = 0; i < n_samples_tot; ++i) {
        intensity[i] = std::pow(std::abs(elaser_out[0][i]),2) + std::pow(std::abs(elaser_out[1][i]),2);   
    }

    double rin2 = 10 * std::log10(variance(intensity) * mod_stp.dt / std::pow(P_avg,2));
}

#endif
