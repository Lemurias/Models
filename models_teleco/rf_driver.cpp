#include "rf_driver.hpp"


/* ======================================================================
                          RF driver
 *
 * @brief   compute rf_driver function.
 * @param   matrix double of signals
 * @retval  matrix <double>
 *
 ======================================================================= */
std::vector<std::vector<double>> rf_driver(DoubleMatrix& sgnls, rf_driver_t& rf_drv)
{
    uint32_t i, j;

    std::vector<std::vector<double>> u0(4, std::vector<double>(1));
    std::vector<std::vector<double>> rf_result(4, std::vector<double>(sgnls[0].size()));

    u0[0][0] = rf_drv.u0[0][0];
    u0[1][0] = rf_drv.u0[0][1];
    u0[2][0] = rf_drv.u0[1][0];
    u0[3][0] = rf_drv.u0[1][1];

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < sgnls[0].size(); ++j) {
            rf_result[i][j] = u0[i][0] * sgnls[i][j];
        }
    }

    return(rf_result);
}

/* ======================================================================
                          RF signal normalized
 *
 * @brief   RF signal normalized.
 * @param   double matrix of amplifier signal
 * @param   cTransmitter_params struct
 * @param   cParam_signals struct
 * @retval  RF signal normalized
 *
 ======================================================================= */
void RF_signal_normalized(DoubleMatrix& ampl_sgnl, modulation_setup_t& modul_stp)
{
    uint32_t rows_symbols = ampl_sgnl.size();
    uint16_t i;

    std::vector<double> Vpi_f_meas{ 25 };
    std::vector<double> f_meas{ 25 };

    Vpi_f_meas[0] = 0.0;
    Vpi_f_meas[1] = 3.34;
    Vpi_f_meas[2] = 3.34;
    Vpi_f_meas[3] = 3.4;
    Vpi_f_meas[4] = 3.45;
    Vpi_f_meas[5] = 3.5;
    Vpi_f_meas[6] = 3.55;
    Vpi_f_meas[7] = 3.6;
    Vpi_f_meas[8] = 3.7;
    Vpi_f_meas[9] = 3.73;
    Vpi_f_meas[10] = 3.8;
    Vpi_f_meas[11] = 3.85;
    Vpi_f_meas[12] = 3.87;
    Vpi_f_meas[13] = 3.96;
    Vpi_f_meas[14] = 4.0;
    Vpi_f_meas[15] = 4.12;
    Vpi_f_meas[16] = 4.20;
    Vpi_f_meas[17] = 4.25;
    Vpi_f_meas[18] = 4.30;
    Vpi_f_meas[19] = 4.40;
    Vpi_f_meas[20] = 4.45;
    Vpi_f_meas[21] = 4.505;
    Vpi_f_meas[22] = 4.61;
    Vpi_f_meas[23] = 4.70;
    Vpi_f_meas[24] = 4.77;

    f_meas[0] = 0.0;
    for (i = 1; i < 25; i++)
    {
        f_meas[i] = (f_meas[i - 1] + 5.0) * GHZ;
    }
    f_meas[0] = 0.02 * GHZ;

    /*

    Vpi_f_meas_bi = [fliplr(Vpi_f_meas(1, 1:end)) Vpi_f_meas];% Create two sided spectrum
        f_bi = [fliplr(-f_meas(1, 1:end)) f_meas];

    V_pi_f_fit = interp1(f_bi, Vpi_f_meas_bi, f, 'linear', 'extrap');% Interpolate acording to the simulation frequency axis

        V_pi_f = repmat(V_pi_f_fit, r, 1);% This line could be replaced for the Vpi meas for each MZM

        symbols_F = fftshift(fft(symbols, [], 2).*Param_main.dt, 2);

    H_f = [Vpi_rf(1, 1). / V_pi_f(1, :); Vpi_rf(1, 2). / V_pi_f(2, :); ...
        Vpi_rf(2, 1). / V_pi_f(3, :); Vpi_rf(2, 2). / V_pi_f(4, :)];

    signal_normalized = real(ifft(ifftshift(symbols_F.*H_f, 2), [], 2). / Param_main.dt);

    EYE_diagram(Param_main, signal_normalized(1, :));
    */
}