#pragma once

#ifndef _OPTICAL_MODULATOR_H_
#define _OPTICAL_MODULATOR_

#include "setup.hpp"
#include "misc.hpp"


//void DP_I_Q_Modulator(std::vector<std::vector<Complex>>& elaser, std::vector<std::vector<double>>& ampl_signal, optical_modulator_t& opt_mod, rf_driver_t& rf_drv, std::vector<std::vector<Complex>>& e_out);
void DP_I_Q_Modulator(ComplexMatrix& elaser, std::vector<std::vector<double>>& ampl_signal, optical_modulator_t& opt_mod, rf_driver_t& rf_drv, std::vector<std::vector<Complex>>& e_out);
#endif