#ifndef _DIGITAL_SYMBOL_GEN_H_
#define _DIGITAL_SYMBOL_GEN_

#include "setup.hpp"

//void digital_symbols_gen(modulation_setup_t& mod_stp, std::unique_ptr<std::unique_ptr<double[]>[]>& symbs, std::uint32_t& n, std::uint32_t& m);
void digital_symbols_gen(modulation_setup_t& mod_stp, std::vector<std::vector<double>>& symbs, std::uint32_t& n, std::uint32_t& m);

#endif