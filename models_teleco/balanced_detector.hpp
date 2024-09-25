#include "misc.hpp"
#include "setup.hpp"


//void balanced_detector(Complex* block_1, Complex* block_2, size_t nCols, receiver_params_t& rx_prms, modulation_setup_t& mod_stp, std::vector<double>& iout);
void balanced_detector(Complex* row1_in1, Complex* row2_in1, Complex* row1_in2, Complex* row2_in2, size_t nCols, receiver_params_t& rx_prms, modulation_setup_t& mod_stp, std::vector<double>& iout);