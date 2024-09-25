#include "misc.hpp"
#include "pulse_shaping.hpp"
#include "digital_symbols_gen.hpp"
#include "setup.hpp"

/* ======================================================================
                          pulse_shaping
 *
 * @brief   compute the pulse_shaping.
 * @param   symbols matrix
 * @param   cTransmitter_params struct
 * @param   cParam_signals struc
 * @retval  pulse shaping matrix
 *
 ======================================================================= */
std::vector<std::vector<double>> pulse_shaping(std::unique_ptr<std::unique_ptr<double[]>[]>& symbs, modulation_setup_t& modul_stp, pulse_shaping_t& puls_sh)
{
    uint32_t i, j, index;

    // ch1, ch2, ch3, ch4 = symbols vector ( one sample vector / symbol)
    uint16_t n_symbols = modul_stp.n_symbols;
    uint16_t n_sampl_per_sym = modul_stp.n_samples_per_sym;

    uint32_t vector_size = modul_stp.n_symbols * modul_stp.n_samples_per_sym;

    std::vector<double> ch1_upsampled(vector_size, 0.0);
    std::vector<double> ch2_upsampled(vector_size, 0.0);
    std::vector<double> ch3_upsampled(vector_size, 0.0);
    std::vector<double> ch4_upsampled(vector_size, 0.0);

 /*   symbs = std::make_unique<std::unique_ptr<double[]>[]>(4);
    for (int i = 0; i < 4; ++i) {
        symbs[i] = std::make_unique<double[]>(modul_stp.n_symbols);
    }*/

    for (i = 0; i < n_symbols; i++)
    {
        index = i * n_sampl_per_sym;
        ch1_upsampled[index] = symbs[0][i];
        ch2_upsampled[index] = symbs[1][i];
        ch3_upsampled[index] = symbs[2][i];
        ch4_upsampled[index] = symbs[3][i];

    }

    //std::unique_ptr<std::unique_ptr<double[]>[]> signal = std::make_unique<std::unique_ptr<double[]>[]>(4);
    //for (int i = 0; i < 4; ++i) {
    //    signal[i] = std::make_unique<double[]>(matrix_cols_size);
    //}
    /*std::vector<double> ch1_upsampled(vector_size, 0.0);
    std::vector<double> ch2_upsampled(vector_size, 0.0);
    std::vector<double> ch3_upsampled(vector_size, 0.0);
    std::vector<double> ch4_upsampled(vector_size, 0.0);*/

    //for (i = 0; i < 4; i++) {
    //    for (j = 0; j < matrix_cols_size; j++) {
    //        signal[i][j] = 0.0;
    //    }
    //}

    //for (i = 0; i < n_symbols; i++)
    //{
    //    j = i * n_sampl_per_sym;
    //    signal[0][j] = symbs[0][i];
    //    signal[1][j] = symbs[1][i];
    //    signal[2][j] = symbs[2][i];
    //    signal[3][j] = symbs[3][i];
    //}

    /****************************************************************************************
     * Lo siguiente es porque necesito el vector de frecuencia para
     * generar el filtro coseno elevado-- - si a la funcion
     * Electrical_signal_gen.m se le pasa como argumento el vector de
     * frecuencias de la simulacion(lo que llamo "f_vect"), lo siguiente se podría borrar.
    ****************************************************************************************/

    std::vector<double> RC_transferfunction = raised_cosine_filter(modul_stp, 1 / (modul_stp.baudrate * GHZ), puls_sh.rolloff);

    fft_shift(RC_transferfunction);
    ComplexVector cx_shiftRC_transferFunction = convertToComplex(RC_transferfunction);

    std::vector<std::vector<double>> signal(4, std::vector<double>(modul_stp.n_samples, 0.0));

    /***************** Ix ********************************/
    ComplexVector cx_ch1_upsampled = convertToComplex(ch1_upsampled);
    fft_o(cx_ch1_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        cx_ch1_upsampled[i] *= cx_shiftRC_transferFunction[i];

    ifft_o(cx_ch1_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        signal[0][i] = cx_ch1_upsampled[i].real() * modul_stp.n_samples_per_sym;

    /***************** Qx ********************************/
    ComplexVector cx_ch2_upsampled = convertToComplex(ch2_upsampled);
    fft_o(cx_ch2_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        cx_ch2_upsampled[i] *= cx_shiftRC_transferFunction[i];

    ifft_o(cx_ch2_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        signal[1][i] = cx_ch2_upsampled[i].real() * modul_stp.n_samples_per_sym;

    /***************** Iy ********************************/
    ComplexVector cx_ch3_upsampled = convertToComplex(ch3_upsampled);
    fft_o(cx_ch3_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        cx_ch3_upsampled[i] *= cx_shiftRC_transferFunction[i];

    ifft_o(cx_ch3_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        signal[2][i] = cx_ch3_upsampled[i].real() * modul_stp.n_samples_per_sym;

    /***************** Qy ********************************/
    ComplexVector cx_ch4_upsampled = convertToComplex(ch4_upsampled);
    fft_o(cx_ch4_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        cx_ch4_upsampled[i] *= cx_shiftRC_transferFunction[i];

    ifft_o(cx_ch4_upsampled);
    for (i = 0; i < modul_stp.n_samples; i++)
        signal[3][i] = cx_ch4_upsampled[i].real() * modul_stp.n_samples_per_sym;

    /*****************************************************/
    return(signal);
}

/* ======================================================================
                          raised_cosine_filter
 *
 * @brief   compute raised cosine filter.
 * @param   cParam_signals struc
 * @param   double Tsymb
 * @param   double beta
 * @retval  raised cosine filter  vector<double>
 *
 ======================================================================= */
std::vector<double> raised_cosine_filter(modulation_setup_t& modul_stp, double Tsymb, double beta)
{
    /**********************************************************
     * f: Frequency vector
     * Tsymb : Symbol time
     * beta : Roll - off factor
     *
     * Initialize the frequency response vector
     * H = zeros(size(f));
     *
     * Compute the frequency response for each frequency point
     *
     * *******************************************************/
    uint32_t i;

    // Initialize the frequency response vector
    std::vector<double> H(modul_stp.n_samples, 0.0);

    // Compute the frequency response for each frequency point
    double aux = (1 - beta) / (2 * Tsymb);

    for (i = 0; i < modul_stp.n_samples; i++)
    {
        if (fabs(modul_stp.f[i]) <= aux)
            H[i] = Tsymb;
        else
            if (fabs(modul_stp.f[i]) <= ((1 + beta) / (2 * Tsymb)))
            {
                double aux1 = fabs(modul_stp.f[i]);
                double aux2 = pi * Tsymb / beta;
                H[i] = (Tsymb / 2) * (1 + std::cos(aux2 * (aux1 - aux)));
            }
            else
                H[i] = 0;
    }

    double Hmax = *std::max_element(H.begin(), H.end());

    if (Hmax != 0)
    {   // Avoid division by zero
        // Divide each element form vector by Hmax
        for (double& val : H) {
            val /= Hmax;
        }
    }

    return (H);

}