#include "instruments.hpp"
#include "misc.hpp"

/************************************************************************
 * Function that measures the power of the input signal Ein
 * Ein is the field defined in two orthogonal polarizations (Jones Vector)
***************************************************************************/
double Power_meter(std::vector<std::vector<std::complex<double>>>& Ein)
{
    uint32_t n_rows = static_cast<uint32_t>(Ein.size());           // Number of rows in Ein
    uint32_t n_cols = static_cast<uint32_t>(Ein[0].size());        // Number of columns in Ein
    uint32_t  n_tot = n_cols;               // The number of columns in each row (assuming N_rows is 2)

    double sum_power = 0.0;

    //double Power_meas_W = 1 / N_tot.*sum((abs(Ein(1, :))). ^ 2 + (abs(Ein(2, :))). ^ 2);

    // Sum the magnitude squared of elements from both rows
    for (uint32_t i = 0; i < n_rows; ++i) {
        for (uint32_t j = 0; j < n_cols; ++j) {
            sum_power += pow(abs(Ein[i][j]), 2);
        }
    }

    // Compute the power measurement
    double Power_meas_W = sum_power / n_tot;
    double Power_meas_dBm = 10 * log10(Power_meas_W / mW);

    return (Power_meas_dBm);
}

//double Power_meter(std::unique_ptr<std::unique_ptr<Complex[]>[]> e_in, std::uint32_t rows, std::uint32_t cols)
double Power_meter(std::unique_ptr<std::unique_ptr<Complex[]>[]>& e_in, std::uint32_t rows, std::uint32_t cols)
{
    //uint32_t n_rows = static_cast<uint32_t>(Ein.size());           // Number of rows in Ein
    //uint32_t n_cols = static_cast<uint32_t>(Ein[0].size());        // Number of columns in Ein
    //uint32_t  n_tot = n_cols;               // The number of columns in each row (assuming N_rows is 2)

    double sum_power = 0.0;

    //double Power_meas_W = 1 / N_tot.*sum((abs(Ein(1, :))). ^ 2 + (abs(Ein(2, :))). ^ 2);
    
    
    // Sum the magnitude squared of elements from both rows
    for (uint32_t i = 0; i < rows; ++i) {
        for (uint32_t j = 0; j < cols; ++j) {
            sum_power += pow(abs(e_in[i][j]), 2);
        }
    }

    // Compute the power measurement
    double Power_meas_W = sum_power / cols;
    double Power_meas_dBm = 10 * log10(Power_meas_W / mW);

    return (Power_meas_dBm);
}

double Power_meter_new(std::vector<std::vector<std::complex<double>>>& Ein)
{
    uint32_t n_tot = static_cast<uint32_t>(Ein[0].size());        // The number of columns in each row (assuming N_rows is 2)
    
    double power_meas_W_pol_X = 0.0;
    double power_meas_W_pol_Y = 0.0;
    double power_meas_W = 0.0;

    // Sum the magnitude squared of elements from both rows
    for (uint32_t i = 0; i < n_tot; ++i) {
        power_meas_W_pol_X += pow(abs(Ein[0][i]), 2);
        power_meas_W_pol_Y += pow(abs(Ein[1][i]), 2);
    }
   
    // Compute the power measurement
    power_meas_W_pol_X /= n_tot;
    power_meas_W_pol_Y /= n_tot;
    power_meas_W = power_meas_W_pol_X + power_meas_W_pol_Y;

    double power_meas_dBm_pol_X = 10 * std::log10(power_meas_W_pol_X / mW);
    double power_meas_dBm_pol_Y = 10 * std::log10(power_meas_W_pol_Y / mW);
    double power_meas_dBm = 10 * std::log10(power_meas_W / mW);

    return (power_meas_dBm);
}