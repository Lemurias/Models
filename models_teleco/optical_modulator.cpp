#include "optical_modulator.hpp"
#include "misc.hpp"
#include "setup.hpp"
#include "instruments.hpp"

/* ======================================================================
                          DP_I_Q_Modulator
 *
 * @brief   compute DP_I_Q_Modulator.
 * @param   complex matrix of Laser's Electric field
 * @param   double matrix of amplifier signal
 * @param   cTransmitter_params struct
 * @retval  complex matrix of E_out
 *
 ======================================================================= */
//void DP_I_Q_Modulator(std::vector<std::vector<Complex>>& elaser, std::vector<std::vector<double>>& ampl_signal, optical_modulator_t& opt_mod, rf_driver_t& rf_drv, std::vector<std::vector<Complex>>& e_out)
void DP_I_Q_Modulator(ComplexMatrix& elaser, std::vector<std::vector<double>>& ampl_signal, optical_modulator_t& opt_mod, rf_driver_t& rf_drv, std::vector<std::vector<Complex>>& e_out)
{
    double aux{ 0.0 };
    uint32_t i, j;

    /*   Elaser = Elaser(1, :);

       Ix = symbols(1, :);
       Qx = symbols(2, :);
       Iy = symbols(3, :);
       Qy = symbols(4, :);

       IL_dB_mod = Param_Tx.IL_dB_mod;

       Vpi_dc = Param_Tx.Vpi_dc;
       Vpi_rf = Param_Tx.Vpi_rf;

   */
    double er_p = pow(10, (opt_mod.ER_p_dB / 10));
    double er_c = pow(10, (opt_mod.ER_c_dB / 10));

    double k_p = 0.5 * (1 - 1 / sqrt(er_p));        // power combining ratio(or splitting) for the Y branch at the parent MZI
    double k_c = 0.5 * (1 - 1 / sqrt(er_c));        // power combining ratio(or splitting) for the Y branch at the child MZI

    std::vector<std::vector<double>> theta_err = elementwise_multiply(opt_mod.IQ_phase_error, (pi / 180));   // [th_Ix_u th_Ix_d th_Qx_u th_Qx_d; th_Iy_u th_Iy_d th_Qy_u th_Qy_d]

    uint32_t N_tot = ampl_signal[0].size();               //length(Ix);

    /* ------------------ Modulated phases of MZM ---------------------- */

    aux = opt_mod.Udc[0][0] / rf_drv.Vpi_dc[0][0];

    // Copy the first row from the matrix

    std::vector<double> U_ix1(N_tot, 0.0);
    std::vector<double> U_ix2(N_tot, 0.0);

    for (j = 0; j < N_tot; ++j)
        U_ix1[j] = (ampl_signal[0][j] / rf_drv.Vpi_rf[0][0] + aux) / 2;

    // U_ix2 = -U_ix1;
    for (i = 0; i < U_ix1.size(); ++i)
        U_ix2[i] = -U_ix1[i];

    // Copy the first row from the matrix

    aux = opt_mod.Udc[0][1] / rf_drv.Vpi_dc[0][1];

    std::vector<double> U_qx1(N_tot, 0.0);
    std::vector<double> U_qx2(N_tot, 0.0);

    for (j = 0; j < N_tot; ++j)
        U_qx1[j] = (ampl_signal[1][j] / rf_drv.Vpi_rf[0][1] + aux) / 2;

    // U_qx2 = -U_qx1;
    for (i = 0; i < U_qx1.size(); ++i)
        U_qx2[i] = -U_qx1[i];

    // Copy the second row from the matrix
    aux = opt_mod.Udc[1][0] / rf_drv.Vpi_dc[1][0];

    std::vector<double> U_iy1(N_tot, 0.0);
    std::vector<double> U_iy2(N_tot, 0.0);

    for (j = 0; j < N_tot; ++j)
        U_iy1[j] = (ampl_signal[2][j] / rf_drv.Vpi_rf[1][0] + aux) / 2;

    // U_iy2 = -U_iy1;
    for (i = 0; i < U_iy1.size(); ++i)
        U_iy2[i] = -U_iy1[i];

    // Copy the third row from the matrix
    aux = opt_mod.Udc[1][1] / rf_drv.Vpi_dc[1][1];

    std::vector<double> U_qy1(N_tot, 0.0);
    std::vector<double> U_qy2(N_tot, 0.0);

    for (j = 0; j < N_tot; ++j)
        U_qy1[j] = (ampl_signal[3][j] / rf_drv.Vpi_rf[1][1] + aux) / 2;

    // U_qy2 = -U_qy1;
    for (i = 0; i < U_qy1.size(); ++i)
        U_qy2[i] = -U_qy1[i];


    std::vector<double> Pd_p{ sqrt(k_p), sqrt(1 - k_p) };   // vector of the power divider(parent MZI)
    std::vector<double> Pc_p{ sqrt(k_p), sqrt(1 - k_p) };   // vector of the power combiner(parent MZI)

    ComplexMatrix E_a(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));

    // First power divider (parent)
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < N_tot; ++j) {
            E_a[i][j] = Pd_p[i] * elaser[0][j];
        }
    }

    ComplexMatrix E_Ix_Qx_in(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));
    ComplexMatrix E_Iy_Qy_in(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));

    // Second power divider(parent)
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < N_tot; ++j) {
            E_Ix_Qx_in[i][j] = Pd_p[i] * E_a[0][j];
        }
    }

    // Second power divider(parent)
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < N_tot; ++j) {
            E_Iy_Qy_in[i][j] = Pd_p[i] * E_a[1][j];
        }
    }

    /*********************************************************************************************/

    std::vector<double> aux_vector(N_tot, 0.0);
    std::vector<Complex> aux_complex_vector1(N_tot);
    std::vector<Complex> aux_complex_vector2(N_tot);
    std::vector<Complex> EI_x(N_tot);
    std::vector<Complex> EQ_x(N_tot);
    std::vector<Complex> EI_y(N_tot);
    std::vector<Complex> EQ_y(N_tot);

    double aux_double1;
    std::complex<double> aux_complex1;
    std::complex<double> aux_complex2;

    // exp(1i * pi * u_ph(1, 1) / Vpi_dc(1, 3))
    aux_complex2 = std::exp(Complex(0, pi * opt_mod.u_ph[0][0] / rf_drv.Vpi_dc[0][2]));

    for (j = 0; j < N_tot; ++j)
    {
        // (theta_err(1,1) + pi.*U_ix1)
        aux_double1 = U_ix1[j] * pi + theta_err[0][0];
        // k_c.*exp(1i*(aux_vector))
        aux_complex1 = std::exp(Complex(0.0, aux_double1)) * k_c;
        // theta_err(1,2) + pi.*U_ix2
        aux_double1 = U_ix2[j] * pi + theta_err[0][1];
        // (1-k_c).*exp(1i*(aux_vector))
        aux_complex1 += std::exp(Complex(0, aux_double1)) * (1 - k_c);

        EI_x[j] = aux_complex1 * E_Ix_Qx_in[0][j];
    }

    for (j = 0; j < N_tot; ++j)
    {
        // theta_err(1,3) + pi.*U_qx1
        aux_double1 = U_qx1[j] * pi + theta_err[0][2];
        // k_c.*exp(1i*(theta_err(1,3) + pi.*U_qx1))
        aux_complex1 = std::exp(Complex(0.0, aux_double1)) * k_c;
        // theta_err(1,4) + pi.*U_qx2
        aux_double1 = U_qx2[j] * pi + theta_err[0][3];
        // (1-k_c).*exp(1i*(theta_err(1,4) + pi.*U_qx2))
        aux_complex1 += std::exp(Complex(0, aux_double1)) * (1 - k_c);

        EQ_x[j] = aux_complex1 * E_Ix_Qx_in[1][j] * aux_complex2;
    }

    for (j = 0; j < N_tot; ++j)
    {
        // theta_err(2,1) + pi.*U_iy1
        aux_double1 = U_iy1[j] * pi + theta_err[1][0];
        // k_c.*exp(1i*(theta_err(2,1) + pi.*U_iy1))
        aux_complex1 = std::exp(Complex(0.0, aux_double1)) * k_c;
        // theta_err(2,2) + pi.*U_iy2
        aux_double1 = U_iy2[j] * pi + theta_err[1][1];
        // (1-k_c).*exp(1i*(theta_err(2,2) + pi.*U_iy2))
        aux_complex1 += std::exp(Complex(0, aux_double1)) * (1 - k_c);

        EI_y[j] = aux_complex1 * E_Iy_Qy_in[0][j];
    }


    // exp(1i * pi * u_ph(2, 1) / Vpi_dc(2, 3));
    aux_complex2 = std::exp(Complex(0, pi * opt_mod.u_ph[1][0] / rf_drv.Vpi_dc[1][2]));

    for (j = 0; j < N_tot; ++j)
    {
        // theta_err(2,3) + pi.*U_qy1
        aux_double1 = U_qy1[j] * pi + theta_err[1][2];
        // k_c.*exp(1i*(theta_err(2,3) + pi.*U_qy1))
        aux_complex1 = std::exp(Complex(0.0, aux_double1)) * k_c;
        // theta_err(2,4) + pi.*U_qy2
        aux_double1 = U_qy2[j] * pi + theta_err[1][3];
        // (1-k_c).*exp(1i*(theta_err(2,4) + pi.*U_qy2))
        aux_complex1 += std::exp(Complex(0, aux_double1)) * (1 - k_c);

        EQ_y[j] = aux_complex1 * E_Iy_Qy_in[1][j] * aux_complex2;
    }

    /*********************************************************************************************/

    // E_out_x = Pc_p * [EI_x; EQ_x];
    // First power combiner(parent)
    std::vector<std::vector<Complex>> E_out_x(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));
    for (i = 0; i < N_tot; i++)
    {
        E_out_x[0][i] = Pc_p[0] * EI_x[i] + Pc_p[1] * EQ_x[i];
    }

    // E_out_y = Pc_p * [EI_y; EQ_y]; 
    // First power combiner(parent)
    std::vector<std::vector<Complex>> E_out_y(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));
    for (i = 0; i < N_tot; i++)
    {
        E_out_y[0][i] = Pc_p[0] * EI_y[i] + Pc_p[1] * EQ_y[i];
    }

    double theta_pol = 90 * pi / 180;   // Polarization rotation angle[rad]

    // Rot_pol = [cos(theta_pol) -sin(theta_pol);sin(theta_pol) cos(theta_pol)];    // rotation matrix
    std::vector<std::vector<double>> Rot_pol = create_rotation_matrix(theta_pol);

    std::vector<std::vector<Complex>> E_out_y_rot(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));
    //std::vector<std::vector<Complex>> E_out(2, std::vector<Complex>(N_tot, Complex(0.0, 0.0)));

    E_out_y_rot = multiply_matrices(Rot_pol, E_out_y);  // Pol rotation

    //E_out = Pc_p(1, 1).*E_out_x + Pc_p(1, 2).*E_out_y_rot; 
    // Second power combiner(parent)
    aux_double1 = -opt_mod.IL_dB_mod / (20 * log10(exp(1)));

    e_out.resize(2, std::vector<Complex>(N_tot));
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < N_tot; j++)
        {
            e_out[i][j] = (Pc_p[0] * E_out_x[i][j] + Pc_p[1] * E_out_y_rot[i][j]) * exp(aux_double1);
        }
    }

//    double Power_avg_mod_out_dBm = Power_meter(e_out);
}
