/*
    version 1.0
    Transcript from matlab to C++: Leonardo Ortiz
    Last update: 08/05/24
*/

#include <iostream>
#include <cstdint>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>
#include <complex>
#include <algorithm> // for std::max_element

#include "instruments.hpp"
#include "misc.hpp"
#include "setup.hpp"
#include "laser.hpp"
#include "digital_symbols_gen.hpp"
#include "pulse_shaping.hpp"
#include "rf_driver.hpp"
#include "optical_modulator.hpp"
#include "coh_opt_rcvr.hpp"
#include "optical_hybrid.hpp"
#include "balanced_detector.hpp"


modulation_setup_t  modulation_stp;
laser_params_t      laser_stp;
laser_params_t      laser_rcv;
pulse_shaping_t     pulse_shng;
rf_driver_t         rf_drv;
optical_modulator_t opt_modulator;
receiver_params_t   rcv_params;
channel_params_t    channel_params;

// Funci�n para crear una matriz de NxM
std::unique_ptr<std::unique_ptr<Complex[]>[]> createMatrix(int N, int M);
//std::unique_ptr<std::unique_ptr<double[]>[]> symbols;
std::vector<std::vector<double>> symbols;
std::vector<std::vector<Complex>> ef_laser;
std::vector<std::vector<std::complex<double>>> Eout_Tx;


int main()
{
    uint32_t i, j;

    modulat_setup(DP_QPSK, 960, pow(2, 15), pow(2, 1), modulation_stp);

    saveParametersTxInFile(modulation_stp, "Modulation Setup.txt");

    laser_setup(13.0, 1550.0, 15*kHZ, 0.0, -145.0, 0.0, laser_stp);

    // Create a matrix
    uint32_t field_cols = modulation_stp.n_symbols * modulation_stp.n_samples_per_sym;
    //auto EF_laser_out = createMatrix(2, field_cols);

    ComplexMatrix EF_laser_out(2, std::vector<Complex>(field_cols, Complex(0.0, 0.0)));
    compute_laser_out(modulation_stp, laser_stp, EF_laser_out, 2, field_cols);

    double Power_avg_Laser_dBm = Power_meter(EF_laser_out);     //% Test point
    double Power_avg_Laser_dBm_new = Power_meter_new(EF_laser_out);     //% Test point

    std::cout << "Power out laser: " << Power_avg_Laser_dBm << " dBm\n";
    
    /***********************************************************************************/
    saveMatrixInFile(EF_laser_out, "Tx_Elaser.txt");
    /***********************************************************************************/


    /****************************************************************
      ****************** IQ Data Generation **************************
      *
      * ---------------- - Digital stage----------------------------
      *  1) Generation of digital symbols
      * *************************************************************/
    std::uint32_t symbols_rows;
    std::uint32_t symbols_cols; 
    digital_symbols_gen(modulation_stp, symbols, symbols_rows, symbols_cols);
   
    std::string text;

    switch (modulation_stp.modulation)
    {
    case SP_BPSK:
        text = "SP_BPSK_" + std::to_string(modulation_stp.bitrate) + "Gbps.txt";
        break;
    case DP_BPSK:
        text = "DP_BPSK_" + std::to_string(modulation_stp.bitrate) + "Gbps.txt";
        break;
    case SP_QPSK:
        text = "SP_QPSK_" + std::to_string(modulation_stp.bitrate) + "Gbps.txt";
        break;
    case DP_QPSK:
        text = "DP_QPSK_" + std::to_string(modulation_stp.bitrate) + "Gbps.txt";
        break;
    default:
        break;

    }

    /******************************************************************/
    saveMatrixInFile(symbols, text);
    /******************************************************************/

    /****************************************************************
     ****************** Pulse shaping **************************
     * 
     * 2) Pulse shaping and generation of analog electrical signals
     ***************************************************************/
    std::vector<std::vector<double>> dac_out(4, std::vector<double>(modulation_stp.n_samples, 0.0));

    pulse_shng.p_shaping_flag = true;     // enable pulse shaping
    pulse_shng.rolloff = 0.5;
    
    if(pulse_shng.p_shaping_flag == true)
    {
        dac_out = pulse_shaping(symbols, modulation_stp, pulse_shng);
    }

    // save file for matlab inspection/verification EYE diagram     
    saveMatrixInFile(dac_out, "dac_out.txt");

    /****************************************************************
     ----------------------- RF Driver -----------------------------
     ***************************************************************/
    rf_drv.Vpi_dc.resize(2);
    for (i = 0; i < 2; ++i)
        rf_drv.Vpi_dc[i].resize(static_cast<double>(3));

    rf_drv.Vpi_rf.resize(2);
    for (i = 0; i < 2; ++i)
        rf_drv.Vpi_rf[i].resize(static_cast<double>(2));

    rf_drv.u0.resize(2);
    for (i = 0; i < 2; ++i)
        rf_drv.u0[i].resize(static_cast<double>(2));

     // RF driver
    rf_drv.Vpi_dc[0][0] = 3.34;   // Vpi_dc_Ix
    rf_drv.Vpi_dc[0][1] = 3.34;   // Vpi_dc_Qx
    rf_drv.Vpi_dc[0][2] = 3.34;   // Vph_x
    rf_drv.Vpi_dc[1][0] = 3.34;   // Vpi_dc_Iy
    rf_drv.Vpi_dc[1][1] = 3.34;   // Vpi_dc_Qy
    rf_drv.Vpi_dc[1][2] = 3.34;   // Vph_y

    // Constant RF Vpi
    rf_drv.Vpi_rf[0][0] = 2.6;   // Vpi_rf_Ix
    rf_drv.Vpi_rf[0][1] = 2.6;   // Vpi_rf_Qx
    rf_drv.Vpi_rf[1][0] = 2.6;   // Vpi_rf_Iy
    rf_drv.Vpi_rf[1][1] = 2.6;   // Vpi_rf_Qy

    rf_drv.Gain_dB   = 20;  // Gain of the electrical amplifier

    // Peak voltaje of the electrical signal at the driver's output
    for (i = 0; i < rf_drv.Vpi_rf.size(); i++) {
        for (j = 0; j < rf_drv.Vpi_rf[0].size(); j++) {
            rf_drv.u0[i][j] = 0.95 * rf_drv.Vpi_rf[i][j];
        }
    }

    rf_drv.NF_dB = 3.0;    // Amplifier noise figure in dB       
    rf_drv.BW_3dB = 1.75 * modulation_stp.baudrate * GHZ; // Amplifier bandwidth
    DoubleMatrix ampl_signal = rf_driver(dac_out, rf_drv);

    // save file for matlab inspection/verification EYE diagram     
    saveMatrixInFile(ampl_signal, "ampl_signal.txt");

    /****************************************************************
      Pre distortion filter for the frequency dependence of the Vpi
     ***************************************************************/

    rf_drv.Vpi_f_flag = false;    // = true Enable frequency dependency of Vpi

    // if(Param_Tx.Vpi_f_flag == true)
    //     [ampl_signal] = RF_signal_normalized(ampl_signal, Param_Tx, Param_main, Units);


    /****************************************************************
     ----------------- Optical Modulator ----------------------------
     ***************************************************************/

     // Optical modulator (DP-IQ modulator)
    opt_modulator.IL_dB_mod = 0;    // DP-IQ modulator insertion loss [dB]
    opt_modulator.ER_p_dB = 22;     // Extinction ratio for parent MZI(Amplitude imbalance)
    opt_modulator.ER_c_dB = 20;     // Extinction ratio for child MZI(Amplitude imbalance)

    // I/Q Quadrature Phase error [degree]
    opt_modulator.IQ_phase_error.resize(2);
    for (i = 0; i < 2; ++i)
    opt_modulator.IQ_phase_error[i].resize(static_cast<double>(4));

    for (i = 0; i < opt_modulator.IQ_phase_error.size(); i++) {
        for (j = 0; j < opt_modulator.IQ_phase_error[0].size(); j++) {
            opt_modulator.IQ_phase_error[i][j] = 0.1;
        }
    }

    opt_modulator.IQ_phase_error = elementwise_multiply(opt_modulator.IQ_phase_error, 0.0);
    /*
    opt_modulator.IQ_phase_error[0][0] = 0.01;
    opt_modulator.IQ_phase_error[0][1] = 0.01;
    opt_modulator.IQ_phase_error[0][2] = 0.01;
    opt_modulator.IQ_phase_error[0][3] = 0.01;
    opt_modulator.IQ_phase_error[1][0] = 0.01;
    opt_modulator.IQ_phase_error[1][1] = 0.01;
    opt_modulator.IQ_phase_error[1][2] = 0.01;
    opt_modulator.IQ_phase_error[1][3] = 0.01;
    */

    opt_modulator.Udc.resize(2);
    for (i = 0; i < 2; ++i)
        opt_modulator.Udc[i].resize(static_cast<double>(2));

    // [Udc_Ix Udc_Qx; Udc_Iy Udc_Qy] DC bias voltage refered to Vpi_DC
    opt_modulator.Udc[0][0] = rf_drv.Vpi_dc[0][0];
    opt_modulator.Udc[0][1] = rf_drv.Vpi_dc[0][1];
    opt_modulator.Udc[1][0] = rf_drv.Vpi_dc[1][0];
    opt_modulator.Udc[1][1] = rf_drv.Vpi_dc[1][1];

    opt_modulator.u_ph.resize(2);
    for (i = 0; i < 2; ++i)
        opt_modulator.u_ph[i].resize(static_cast<double>(1));

       
    // I/Q phase voltage refered to Vpi_dc
    opt_modulator.u_ph[0][0] = 0.5 * rf_drv.Vpi_dc[0][2];
    opt_modulator.u_ph[1][0] = 0.5 * rf_drv.Vpi_dc[1][2];

    DP_I_Q_Modulator(EF_laser_out, ampl_signal, opt_modulator, rf_drv, Eout_Tx);
    
    double Power_avg_mod_out_dBm = Power_meter(Eout_Tx);

    std::cout << "Output power Tx: " << Power_avg_mod_out_dBm << '\n';


    saveMatrixInFile(Eout_Tx, "E_mod_out.txt");
    /********************************************************************/

    /*********************************************************************
    *******************  CHANNEL  ****************************************
    *********************************************************************/
    channel_params.length = 10.0;                       // Fiber length[km]
    channel_params.alpha_dB = 0.2;                      // Fiber atenuation [dB/km]
    channel_params.D = 17.5;                            // GVD dispersion coef. [ps/(nm km)] - It is later converted to beta2
    channel_params.beta3 = 0.0;                         // Third order dispersion coef. [ps^3/km]; 
    channel_params.pmd = 0.1;                           // PMD coeficient  [ps/sqrt(km)]
    channel_params.sopmd = 0.0;                         // Second-order PMD coeficient [ps^2/sqrt(km)]
    channel_params.PDL_dB = 0.0;                        // Polarization-dependent loss [dB]
    channel_params.dz = 0.01 * channel_params.length;   // Step size of Coarse Step Method [km]


    /*[Eout_channel, T, M] = fiber_propagation_CSM(Eout_tx, Param_fiber, Param_Tx, Param_main, Units);
    PMDcoef_estimator;

    Power_avg_fiber_out_dBm = Power_meter(Eout_channel, Units, 'output of fiber propagation');% Test point
        Optical_spectrum_analyzer(Eout_channel, 0, Param_main, Units, 'Output of fiber prop. - Optical spectrum');*/





    /*********************************************************************
    * RECEIVER
    *********************************************************************/

    /*********************************************************************
    * Local Oscilator Laser 
    *********************************************************************/
    laser_setup(8.0, 1550.0, 15 * kHZ, 0.0, -145.0, 0.0, laser_rcv);

    receiver_setup(laser_rcv, 1, 10E-9, true, 0.75 * modulation_stp.baudrate * GHZ, 290, 50, true, 0.0, rcv_params);
    
    saveParametersRxInFile(rcv_params, "Receiver_setup.txt");
    
    // Create a matrix
    uint32_t rcv_field_cols = modulation_stp.n_symbols * modulation_stp.n_samples_per_sym;
    //auto EF_rcv_laser_out = createMatrix(2, rcv_field_cols);
    ComplexMatrix EF_rcv_laser_out(2, std::vector<Complex>(rcv_field_cols, Complex(0.0, 0.0)));
    compute_laser_out(modulation_stp, laser_rcv, EF_rcv_laser_out, 2, rcv_field_cols);

    saveMatrixInFile(EF_rcv_laser_out, "Elo_.txt");

    //Power_avg_Laser_dBm = Power_meter(EF_rcv_laser_out, 2, rcv_field_cols);     //% Test point
    Power_avg_Laser_dBm = Power_meter(EF_rcv_laser_out);     //% Test point

    /*********************************************************************
    * Coherent Receiver
    *********************************************************************/
    ComplexMatrix E1_8;
    
    // coherent_optical_receiver()
    //E1_8 = Coherent_optical_receiver(Ein_rx, Elo_, Param_Rx);
    coherent_optical_receiver(Eout_Tx, EF_rcv_laser_out, rcv_params, modulation_stp, E1_8);

    saveMatrixInFile(E1_8, "E1_8.txt");

    /*********************************************************************
    * Photodetectors %
    *********************************************************************/
    uint32_t n_cols = E1_8[0].size();
	std::vector<double> I_X(n_cols);
	std::vector<double> Q_X(n_cols);
	std::vector<double> I_Y(n_cols);
	std::vector<double> Q_Y(n_cols);

	balanced_detector(E1_8[0].data(), E1_8[1].data(), E1_8[2].data(), E1_8[3].data(), n_cols, rcv_params, modulation_stp, I_X);
    balanced_detector(E1_8[4].data(), E1_8[5].data(), E1_8[6].data(), E1_8[7].data(), n_cols, rcv_params, modulation_stp, Q_X);
    balanced_detector(E1_8[8].data(), E1_8[9].data(), E1_8[10].data(), E1_8[11].data(), n_cols, rcv_params, modulation_stp, I_Y);
    balanced_detector(E1_8[12].data(), E1_8[13].data(), E1_8[14].data(), E1_8[15].data(), n_cols, rcv_params, modulation_stp, Q_Y);
    
	if(rcv_params.DC_block_on == true)
	{
		double mean_ix = 0.0;
		double mean_qx = 0.0;
		double mean_iy = 0.0;
		double mean_qy = 0.0;
		for(size_t i = 0; i < n_cols; i++){
			mean_ix += I_X[i];
			mean_qx += Q_X[i];
			mean_iy += I_Y[i];
			mean_qy += Q_Y[i];
		}
		mean_ix /= n_cols;
		mean_qx /= n_cols;
		mean_iy /= n_cols;
		mean_qy /= n_cols;

		for(size_t i = 0; i < n_cols; i++){
			I_X[i] -=  mean_ix;
			Q_X[i] -=  mean_qx;
			I_Y[i] -=  mean_iy;
			Q_Y[i] -=  mean_qy;
		}
	}

    saveMatrixInFile(I_X, "I_X_output_receptor.txt");
    saveMatrixInFile(Q_X, "Q_X_output_receptor.txt");
    saveMatrixInFile(I_Y, "I_Y_output_receptor.txt");
    saveMatrixInFile(Q_Y, "Q_Y_output_receptor.txt");

	// if Param_Rx.TIA_on == 1
     
	// 	Nn          = length(I_X);
	// 	kb          = physconst('Boltzmann');
	// 	G_tia       = 10^(Param_Rx.TIA_G_dB/20);
	// 	NF_tia      = 10^(Param_Rx.TIA_NF_dB/10);
	// 	Te          = (NF_tia-1)*290; % Equivalent noise temeperature
	// 	B           = Param_ADC.Antialiasing_filter_f3dB;
		
	// 	P_noise_out = (G_tia^2)*Te*kb*B; % Noise power added by the TIA
		
	// 	I_X = I_X*G_tia + sqrt(P_noise_out*50).*randn(1,Nn);  
	// 	Q_X = Q_X*G_tia + sqrt(P_noise_out*50).*randn(1,Nn);
	// 	I_Y = I_Y*G_tia + sqrt(P_noise_out*50).*randn(1,Nn);
	// 	Q_Y = Q_Y*G_tia + sqrt(P_noise_out*50).*randn(1,Nn);
	//  else
     
 	// end

    

    return 0;
}


// Funci�n para crear una matriz de NxM
std::unique_ptr<std::unique_ptr<Complex[]>[]> createMatrix(int N, int M) {
    // Asignaci�n din�mica de memoria para la matriz
    std::unique_ptr<std::unique_ptr<Complex[]>[]> matrix(new std::unique_ptr<Complex[]>[N]);
    for (int i = 0; i < N; ++i) {
        matrix[i] = std::unique_ptr<Complex[]>(new Complex[M]);
    }
    return matrix;
}

