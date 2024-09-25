#pragma once

#ifndef _PULSE_SHAPING_H_
#define _PULSE_SHAPING_

#include "setup.hpp"
#include "misc.hpp"

//std::vector<std::vector<double>> pulse_shaping(std::unique_ptr<std::unique_ptr<double[]>[]>& symbs, modulation_setup_t& modul_stp, pulse_shaping_t& puls_sh);
std::vector<std::vector<double>> pulse_shaping(std::vector<std::vector<double>>& symbs, modulation_setup_t& modul_stp, pulse_shaping_t& puls_sh);
std::vector<double> raised_cosine_filter(modulation_setup_t& modul_stp, double Tsymb, double beta);


#endif