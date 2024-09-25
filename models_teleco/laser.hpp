
#include "misc.hpp"

//void compute_laser_out(modulation_setup_t& mod_stp, laser_params_t& laser_prm, std::unique_ptr<std::unique_ptr<Complex[]>[]> elaser_out);
//void compute_laser_out(modulation_setup_t& mod_stp, laser_params_t& laser_prm, std::unique_ptr<std::unique_ptr<Complex[]>[]>& elaser_out, int N, int M); // 2D dynamic array )
void compute_laser_out(modulation_setup_t& mod_stp, laser_params_t& laser_prm, ComplexMatrix& elaser_out, int N, int M);