#ifndef INSTRUMENTS_HPP
#define INSTRUMENTS_HPP

#include<iostream>
#include<vector>
#include<complex>

#include "misc.hpp"

// [Power_avg_tx_dBm] = Power_meter(Eout_tx, Units); % Test point
double Power_meter(std::vector<std::vector<std::complex<double>>>& Eout_tx);

//double Power_meter(std::unique_ptr<std::unique_ptr<Complex[]>[]> e_in, std::uint32_t rows, std::uint32_t cols);
double Power_meter(std::unique_ptr<std::unique_ptr<Complex[]>[]>& e_in, std::uint32_t rows, std::uint32_t cols);
//double Power_meter(std::unique_ptr<Complex[]>& e_in, std::uint32_t rows, std::uint32_t cols);

double Power_meter_new(std::vector<std::vector<std::complex<double>>>& Ein);

#endif