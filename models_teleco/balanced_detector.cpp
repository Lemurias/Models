#include "balanced_detector.hpp"
#include "detector_filter.hpp"


//void balanced_detector(DoubleMatrix& ein1, DoubleMatrix& ein2, receiver_params_t& rx_prms, modulation_setup_t& mod_stp, std::vector<double>& iout)
//void balanced_detector(Complex* block_1, Complex* block_2, size_t nCols,receiver_params_t& rx_prms, modulation_setup_t& mod_stp, std::vector<double>& iout)
void balanced_detector(Complex* row1_in1, Complex* row2_in1, Complex* row1_in2, Complex* row2_in2, size_t nCols, receiver_params_t& rx_prms, modulation_setup_t& mod_stp, std::vector<double>& iout)
{
    //double resp = rx_prms.PD_Resp;
    //double id   = rx_prms.PD_id;
    //double delta_f = rx_prms.PD_BW;
    //double T = rx_prms.PD_T;
    //double Rload = rx_prms.PD_RL;

    std::vector<double> Is1(nCols, 0.0);
    std::vector<double> Is2(nCols, 0.0);

    double mean_1 = 0.0;
    double mean_2 = 0.0;
    //double auxiliar = 0.0;

    for (size_t i = 0; i < nCols; i++)
    {
        Is1[i] = std::pow(std::abs(row1_in1[i]), 2) + std::pow(std::abs(row2_in1[i]), 2);
        mean_1 += std::pow(Is1[i], 2);

        Is2[i] = std::pow(std::abs(row1_in2[i]), 2) + std::pow(std::abs(row2_in2[i]), 2);
        mean_2 += std::pow(Is2[i], 2);

        /*  if (i >= (nCols - 10))
              auxiliar = mean_1 + mean_2;*/
    }

    double Is1_rms = std::sqrt(mean_1 / nCols);
    double Is2_rms = std::sqrt(mean_1 / nCols);

    /*************************************************************************************/
    /*************************** Shot Noise **********************************************/

    // shot noise variance of PD1 [A^2]
    double var_shot1 = 2 * q * (Is1_rms + rx_prms.PD_id) * rx_prms.PD_BW;

    // shot noise variance of PD2 [A^2]
    double var_shot2 = 2 * q * (Is2_rms + rx_prms.PD_id) * rx_prms.PD_BW;

    /*************************************************************************************/
    /*************************** Thermal noise  ******************************************/

    double var_ther = 4 * kb * rx_prms.PD_T * rx_prms.PD_BW / rx_prms.PD_RL;

    /*************************************************************************************/
    /*************************** Total noise  ********************************************/

    double var_total_pd1 = var_shot1 + var_ther;   // total noise variance of PD1 [A^2]
    double var_total_pd2 = var_shot2 + var_ther;   // total noise variance of PD2 [A^2]


    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::normal_distribution<> d(0.0, 1.0); // Normal distribution with specified mean and stddev

    std::vector<double> noise_curr_pd1(nCols);
    std::vector<double> noise_curr_pd2(nCols);

    for (int i = 0; i < nCols; ++i) {
        noise_curr_pd1[i] = std::sqrt(var_total_pd1) * d(gen); // noise current of PD1 [A]
        noise_curr_pd2[i] = std::sqrt(var_total_pd2) * d(gen); // noise current of PD2 [A]
    }

    /*************************************************************************************/

    std::vector<double> I1(nCols, 0.0);
    std::vector<double> I2(nCols, 0.0);
    std::vector<double> Iout(nCols, 0.0);

    for (size_t i = 0; i < nCols; i++)
    {
        I1[i] = Is1[i] + noise_curr_pd1[i];
        I2[i] = Is2[i] + noise_curr_pd2[i];
        Iout[i] = I1[i] - I2[i];
    }
    //rx_prms.filter_on = false;
    if (rx_prms.filter_on == true)
    {
        detector_filter(Iout, rx_prms, mod_stp, iout);
    }

    /*for (size_t i = 0; i < nCols; i++)
    {
        iout[i] = Iout[i];
    }*/

}