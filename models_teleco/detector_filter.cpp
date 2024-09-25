#include "detector_filter.hpp"

void detector_filter(std::vector<double>& signal_in, receiver_params_t& rx_prms, modulation_setup_t& mod_stp, std::vector<double>& signal_filtered)
{
    ComplexVector H(signal_in.size());
    ComplexVector sig_in(signal_in.size());
    ComplexVector tf_signal_in(signal_in.size());
    ComplexVector tf_signal_out(signal_in.size());

    /*
    % w0 = 2*pi*f0;
    % w  = 2*pi.*f;
    % s  = 1j*w;

    % filter_order = 5;
    % [b,a] = besself(filter_order,w0);   % Bessel analog filter design
    % H     = polyval(b,s)./p
    */

    // Low pass filter for the PD
    for(size_t i = 0; i < signal_in.size(); i++){
        H[i] = 1.0/Complex(1, mod_stp.f[i]/rx_prms.PD_BW);
    }

    for(size_t i = 0; i < signal_in.size(); i++){
        tf_signal_in[i] = Complex(signal_in[i],0.0);
    }

    fft_o(tf_signal_in);
    fft_shift(tf_signal_in);

    for(size_t i = 0; i < signal_in.size(); i++){
        tf_signal_out[i] = H[i] * tf_signal_in[i];
    }

    tf_signal_out = ifftshift(tf_signal_out);
    ifft_o(tf_signal_out);

    for(size_t i = 0; i < signal_in.size(); i++){
        signal_filtered[i] = tf_signal_out[i].real();
    }
     



}