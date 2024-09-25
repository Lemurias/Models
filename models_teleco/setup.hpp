#pragma once

#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <fstream>

#include "misc.hpp"

typedef struct
{
    double              bits_per_symbol;
    double              baudrate;
    modulation_t        modulation;
    double              bitrate;
    std::uint16_t		n_symbols;
    std::uint32_t       n_bits;
    std::uint32_t       n_samples_per_bit;
    std::uint32_t	    n_samples;				// Total number of samples
    std::uint16_t	    n_samples_per_sym;		// Samples per symbol

    double				dt;						// [s]
    //std::vector<double>	t;
    std::unique_ptr<double[]> t;    // unique pointer for dinamic array
    size_t              size_of_t;  // size of variable t
    double				df;
    //std::vector<double>	f;
    std::unique_ptr<double[]> f;    // unique pointer for dinamic array
    size_t              size_of_f;  // size of variable f
    double				siganlrollOffFactor;
}modulation_setup_t;

/*************************************.
* Transmitter's parameters
*************************************/
typedef struct
{
    // Laser
    double	P_dBm;				// Pavg of the laser [dBm]
    double	Lambda;				// wavelength
    double	LW;					// Line Width[Hz]
    double	FO;					// Frequency offset [Hz]
    double	RIN_dBHz;			// Relative Intensity Noise [dB/Hz] 
    double	Freq_var;			// Frequency fluctuations
}laser_params_t;

typedef struct
{
    // Pulse shaping(raise - cosine filter)
    bool	p_shaping_flag;			// = true enable pulse shaping
    double	rolloff;				// 
}pulse_shaping_t;

typedef struct
{
    // RF driver
    std::vector<std::vector<double>> Vpi_dc;// { 2, std::vector<double>(3, 0.0) };
    std::vector<std::vector<double>> Vpi_rf;// { 2, std::vector<double>(2, 0.0) };

    double	Gain_dB;			// Gain of the electrical amplifier
    std::vector<std::vector<double>> u0;// { 2, std::vector<double>(2, 0.0) };	// Peak voltaje of the electrical signal at the driver's output
    double	NF_dB;			// Amplifier noise figure in dB
    double	BW_3dB;			// Amplifier bandwidth
    bool	Vpi_f_flag;				// = true Enable frequency dependency of Vpi
}rf_driver_t;

typedef struct
{
    // Optical modulator (DP-IQ modulator)
    double	IL_dB_mod;				// DP-IQ modulator insertion loss [dB]
    std::vector<std::vector<double>> IQ_phase_error;// { 2, std::vector<double>(4, 0.1) };	// I / Q Quadrature Phase error[degree]
    double	ER_p_dB;				// Extinction ratio for parent MZI (Amplitude imbalance)
    double	ER_c_dB;				// Extinction ratio for child MZI  (Amplitude imbalance)
    std::vector<std::vector<double>> Udc;//{ 2, std::vector<double>(2,0.0) };	// [Udc_Ix Udc_Qx; Udc_Iy Udc_Qy] DC bias voltage refered to Vpi_DC
    std::vector<std::vector<double>> u_ph;// { 2, std::vector<double>(1, 0.0) };		// I / Q phase voltage refered to Vpi_dc
}optical_modulator_t;

/************************************************************************
* Receiver parameters
************************************************************************/
typedef struct
{
    // Laser
    laser_params_t laser_prm_rcv;
    // Coherent receiver
    double PD_Resp;             // Photodiodes responsivity[A / W]
    double PD_id;               // Photodiodes dark current
    bool   filter_on;           // = true enable the PD filter
    double PD_BW;               // Photodiodes bandwidth[Hz]
    double PD_T;                // PD noise temperature[K]
    double PD_RL;               // PD load impedance[Ohms]
    bool   DC_block_on;         // = true enable the PD DC coupled output

    double coup_excess_loss_dB; // Excess loss of hybrid's couplers
}receiver_params_t;
/************************************************************************/

/************************************************************************
* Channel parameters
************************************************************************/
typedef struct
{
    // fiber parameters
    double  length;             // Fiber length [km]
    double  alpha_dB;           // Fiber atenuation [dB/km]
    double  D;                  // GVD dispersion coef. [ps/(nm km)] - It is later converted to beta2
    // double S;                // Dispersion slope
    double  beta3;              // Third order dispersion coef. [ps^3/km]; 
    double  pmd;                // PMD coeficient  [ps/sqrt(km)]
    double  sopmd;              // Second-order PMD coeficient [ps^2/sqrt(km)]
    double  PDL_dB;             // Polarization-dependent loss [dB]
    double  dz;                 // Step size of Coarse Step Method [km]
}channel_params_t;
/************************************************************************/


void modulat_setup(modulation_t modul, double bit_rate, std::uint32_t nbits, std::uint32_t  samples_pr_bits, modulation_setup_t& modu_stp);
void laser_setup(double	pdbm, double wl, double lw, double	fo, double	rin_dbhz, double fvar, laser_params_t& laser_params);
void saveParametersTxInFile(modulation_setup_t& modul_stp, const std::string& fileName);
void saveParametersRxInFile(receiver_params_t& rx_stp, const std::string& fileName);
void receiver_setup(laser_params_t& laser_rcv, double pd_resp, double pd_id, bool fter_on, double pd_bw, double pd_t, double pd_rl, bool dc_block_on, double coup_excess, receiver_params_t& coh_rcv);