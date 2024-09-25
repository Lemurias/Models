#include "setup.hpp"


void modulat_setup(modulation_t modul, double bit_rate, std::uint32_t nbits, std::uint32_t  samples_pr_bits, modulation_setup_t& modu_stp)
{
    std::int32_t i = { 0 };
    std::int32_t j = { 0 };
    std::int32_t half_n_samples = { 0 };

    switch (modul)
    {
    case SP_BPSK:
        modu_stp.bits_per_symbol = 1;
        break;
    case DP_BPSK:
        modu_stp.bits_per_symbol = 2;
        break;
    case SP_QPSK:
        modu_stp.bits_per_symbol = 2;
        break;
    case DP_QPSK:
        modu_stp.bits_per_symbol = 4;
        break;
    case PAM_4:
        modu_stp.bits_per_symbol = 2;
        break;
    case NRZ:
        modu_stp.bits_per_symbol = 1;
        break;
    default:
        modu_stp.bits_per_symbol = 1;
        break;
    }

    modu_stp.modulation = modul;
    modu_stp.bitrate = bit_rate;
    modu_stp.n_bits = nbits;
    modu_stp.n_samples_per_bit = samples_pr_bits;
    modu_stp.baudrate = modu_stp.bitrate / modu_stp.bits_per_symbol;   // Baud rate (in Gbaud)

    modu_stp.n_symbols = static_cast<uint16_t>((modu_stp.n_bits - mod(modu_stp.n_bits, modu_stp.bits_per_symbol)) / modu_stp.bits_per_symbol);

    modu_stp.n_samples = modu_stp.n_bits * modu_stp.n_samples_per_bit; // Total number of samples	
    modu_stp.n_samples_per_sym = static_cast<uint16_t>(modu_stp.n_samples_per_bit * modu_stp.bits_per_symbol); // Samples per symbol

    modu_stp.dt = 1 / (modu_stp.n_samples_per_sym * (modu_stp.bitrate * GHZ / modu_stp.bits_per_symbol));   // in [s]

    modu_stp.size_of_t = modu_stp.n_samples;
    modu_stp.t = std::make_unique<double[]>(modu_stp.size_of_t);  // create an array of n_samples values

    for (i = 0; i < modu_stp.n_samples; ++i) {
        modu_stp.t[i] = i * modu_stp .dt;
    }
    
    modu_stp.df = 1 / (modu_stp.n_samples * modu_stp.dt);                  // (Bitrate*Units.GHz)/N_bits;

    half_n_samples = modu_stp.n_samples / 2;

    modu_stp.size_of_f = modu_stp.n_samples;
    modu_stp.f = std::make_unique<double[]>(modu_stp.size_of_f);  // create an array of n_samples values

    for (i = 0; i < modu_stp.n_samples; ++i) {
        modu_stp.f[i] = static_cast<double>(i - half_n_samples) * modu_stp.df;
    }
}


void laser_setup(double	pdbm, double wl, double lw, double	fo, double	rin_dbhz, double fvar, laser_params_t& laser_params)
{
    laser_params.P_dBm       = pdbm;     // Pavg of the laser [dBm]
    laser_params.Lambda      = wl;		// wavelength
    laser_params.LW          = lw;		// Line Width[Hz]
    laser_params.FO          = fo;		// Frequency offset [Hz]
    laser_params.RIN_dBHz    = rin_dbhz;	// Relative Intensity Noise [dB/Hz] 
    laser_params.Freq_var    = fvar;		// Frequency fluctuations
}

void saveParametersTxInFile(modulation_setup_t& modul_stp, const std::string& fileName)
{
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    file << "Parameters transmitter:" << "\n";
    file << "Bitrate: " << modul_stp.bitrate << "\n";
    file << "Baudrate: " << modul_stp.baudrate << "\n";
    file << "Bits per symbol: " << modul_stp.bits_per_symbol << "\n";
    switch (modul_stp.modulation)
    {
    case SP_BPSK:
        file << "Modulation: SP_BPSK" << "\n";
        break;
    case DP_BPSK:
        file << "Modulation: DP_BPSK" << "\n";
        break;
    case SP_QPSK:
        file << "Modulation: SP_QPSK" << "\n";
        break;
    case DP_QPSK:
        file << "Modulation: DP_QPSK" << "\n";
        break;
    default:
        break;
    }
    file << "N bits: " << modul_stp.n_bits << "\n";
    file << "N samples: " << modul_stp.n_samples << "\n";
    file << "N samples per bits: " << modul_stp.n_samples_per_bit << "\n";
    file << "N samples per sym: " << modul_stp.n_samples_per_sym << "\n";
    file << "N symbols: " << modul_stp.n_symbols << "\n";
    //file << "Signal Roll Off Factor: " << modul_stp.siganlrollOffFactor << "\n"; 
    file << "df: " << modul_stp.df << "\n";
    file << "dt: " << modul_stp.dt << "\n";
    file << "\n\r";        // New line after each row

    uint32_t cols_f = static_cast<uint32_t>(modul_stp.size_of_f);
    file << "f[i] ; t[i] \n";
    for (uint32_t i = 0; i < cols_f; ++i)
    {
        file << modul_stp.f[i] << ';' << modul_stp.t[i] << '\n';
    }

    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;
}

void saveParametersRxInFile(receiver_params_t& rx_stp, const std::string& fileName)
{
    std::ofstream file(fileName);

    if (!file) {
        std::cerr << "The file could not be opened for writing.." << std::endl;
        return;
    }

    file << "Parameters transmitter:" << "\n";
    file << "Laser Rcv Power: " << rx_stp.laser_prm_rcv.P_dBm << "dBm" << "\n";
    file << "Laser Rcv Lambda: " << rx_stp.laser_prm_rcv.Lambda << "nm" << "\n";
    file << "Laser Rcv LW: " << rx_stp.laser_prm_rcv.LW << "Hz" << "\n";
    file << "Laser Rcv FO: " << rx_stp.laser_prm_rcv.FO << "Hz" << "\n";
    file << "Laser Rcv RIN: " << rx_stp.laser_prm_rcv.RIN_dBHz << "dBHz" << "\n";
    //file << "Laser Rcv: " << rx_stp.laser_prm_rcv.Freq_var << "dBm" << "\n";
    file << "Laser Rcv BW: " << rx_stp.PD_BW << "Hz" << "\n";
    file << "Laser Rcv ID: " << rx_stp.PD_id << "A" << "\n";
    file << "Laser Rcv Resp: " << rx_stp.PD_Resp << "A/W" << "\n";
    file << "Laser Rcv RL: " << rx_stp.PD_RL << "ohms" << "\n";
    file << "Laser Rcv T: " << rx_stp.PD_T << "ºK" << "\n";
    file << "Laser Rcv coup excess loss: " << rx_stp.coup_excess_loss_dB << "dB" << "\n";
    file << "Laser Rcv DC Block ON: " << rx_stp.DC_block_on << "\n";
    file << "Laser Filter ON: " << rx_stp.filter_on << "\n";
    file << "\n\r";        // New line after each row

    file.close();
    std::cout << "Matrix saved in " << fileName << std::endl;
}


void receiver_setup(laser_params_t& laser_rcv, double pd_resp, double pd_id, bool fter_on, double pd_bw, double pd_t, double pd_rl, bool dc_block_on, double coup_excess, receiver_params_t& coh_rcv)
{
    coh_rcv.laser_prm_rcv.P_dBm = laser_rcv.P_dBm;
    coh_rcv.laser_prm_rcv.Lambda = laser_rcv.Lambda;
    coh_rcv.laser_prm_rcv.LW = laser_rcv.LW;
    coh_rcv.laser_prm_rcv.FO = laser_rcv.FO;
    coh_rcv.laser_prm_rcv.RIN_dBHz = laser_rcv.RIN_dBHz;
    coh_rcv.laser_prm_rcv.Freq_var = laser_rcv.Freq_var;

    coh_rcv.PD_Resp = pd_resp;
    coh_rcv.PD_id = pd_id;
    coh_rcv.filter_on = fter_on;
    coh_rcv.PD_BW = pd_bw;
    coh_rcv.PD_T = pd_t;
    coh_rcv.PD_RL = pd_rl;
    coh_rcv.DC_block_on = dc_block_on;
    coh_rcv.coup_excess_loss_dB = coup_excess;
}