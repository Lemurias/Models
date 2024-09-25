#include <memory>

#include "misc.hpp"
#include "digital_symbols_gen.hpp"
#include "setup.hpp"


//void digital_symbols_gen(modulation_setup_t& mod_stp, std::unique_ptr<std::unique_ptr<double[]>[]>& symbs, std::uint32_t& n, std::uint32_t& m)
void digital_symbols_gen(modulation_setup_t& mod_stp, std::vector<std::vector<double>>& symbs, std::uint32_t& n, std::uint32_t& m)
{
    uint32_t i, k;
    uint32_t j;

    // Create a random number generator and normal distribution
    std::random_device rd;  // Obtain a random number from hardware
    std::mt19937 gen(rd()); // Seed the generator
    std::normal_distribution<> d(0.0, 1.0); // Normal distribution with specified mean and stddev

    //std::unique_ptr<double[]> data_sgn;    // unique pointer for dinamic array
    //data_sgn = std::make_unique<double[]>(mod_stp.n_bits);  // create an array of n_samples values
    std::vector<double> data_sgn(mod_stp.n_bits);  // create an array of n_samples values

    for (i = 0; i < mod_stp.n_bits; ++i) {
        data_sgn[i] = (d(gen) >= 0) ? 1 : -1; // Apply the sign function
    }
    
    n = 4;
    m = static_cast<uint32_t>(mod_stp.n_symbols);

    /*symbs = std::make_unique<std::unique_ptr<double[]>[]>(4);
    for (int i = 0; i < 4; ++i) {
        symbs[i] = std::make_unique<double[]>(mod_stp.n_symbols);
    }*/
    symbs.resize(4);
    for (uint32_t i = 0; i < 4; i++) {
        symbs[i].resize(mod_stp.n_symbols);
    }

    switch (mod_stp.modulation)
    {
        case SP_BPSK:
            /*
            ch1 = data;
            ch2 = zeros(1, N_symbols);
            ch3 = zeros(1, N_symbols);
            ch4 = zeros(1, N_symbols);
            */
            for (i = 0; i < mod_stp.n_bits; i++)
                symbs[0][i] = data_sgn[i];

            for (j = 1; j < n; j++)
            {
                for (i = 0; i < mod_stp.n_symbols; i++)
                    symbs[j][i] = 0.0;
            }
            break;
        case DP_BPSK:
            /*
            ch1 = data(1:2 : end);
            ch2 = zeros(1, N_symbols);
            ch3 = data(2:2 : end);
            ch4 = zeros(1, N_symbols);
            */
            k = 0;
            for (i = 0; i < mod_stp.n_bits; i = i + 2)
            {
                symbs[0][i / 2] = data_sgn[i];
                symbs[2][(i+1) / 2] = data_sgn[i+1];
            }
            for (j = 1; j < n; j = j + 2)
            {
                for (i = 0; i < mod_stp.n_symbols; i++)
                    symbs[j][i] = 0.0;
            }
            break;
        case SP_QPSK:
            /*
            ch1 = data(1:2 : end);
            ch2 = data(2:2 : end);
            ch3 = zeros(1, N_symbols);
            ch4 = zeros(1, N_symbols);
            */
            k = 0;
            for (i = 0; i < mod_stp.n_bits; i = i + 2)
            {
                symbs[0][i / 2] = data_sgn[i];
                symbs[1][(i+1) / 2] = data_sgn[i+1];
            }
            for (j = 2; j < n; j++)
            {
                for (i = 0; i < mod_stp.n_symbols; i++)
                    symbs[j][i] = 0.0;
            }
            break;
        case DP_QPSK:
            /*
            ch1 = data(1:4 : end);
            ch2 = data(2:4 : end);
            ch3 = data(3:4 : end);
            ch4 = data(4:4 : end);
            */
            for (i = 0; i < mod_stp.n_bits; i = i + 4)
            {
                symbs[0][i / 4] = data_sgn[i];
                symbs[1][i / 4] = data_sgn[i + 1];
                symbs[2][i / 4] = data_sgn[i + 2];
                symbs[3][i / 4] = data_sgn[i + 3];
            }
            break;
        case PAM_4:

            break;
        case NRZ:

            break;
        default:

            break;
    }

    n = 4;
    m = mod_stp.n_symbols;
}