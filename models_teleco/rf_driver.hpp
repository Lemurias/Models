#pragma once

#ifndef _RF_DRIVER_H_
#define _RF_DRIVER_

#include "setup.hpp"
#include "misc.hpp"

std::vector<std::vector<double>> rf_driver(DoubleMatrix& sgnls, rf_driver_t& rf_drv);
void RF_signal_normalized(DoubleMatrix& ampl_sgnl, modulation_setup_t& modul_stp);

#endif