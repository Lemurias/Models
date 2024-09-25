#include "misc.hpp"
#include "setup.hpp"

/*************************************************************************************************************
* Fiber_propagation_CMS
* 
* The script simulates the propagation of the optical signal using the
* coarse - step method(CSM).
* Input: Optical electric field as a 2xN matrix(1st row : X polarization, 2nd row : Y polarization)
* Outputs : -signal_out = Optical electric field after fiber propagation
* -T = Total fiber matrix transfer function -> struct formed by 4 elements
* -M = Polarization fiber matrix transfer function(PMD, SOPMD, RSOP)
* 
* ***********************************************************************************************************/

/************************************************************************************************************
* Fiber parameters, included in Param_fiber.
* 
* L : Fiber length(km)
* alpha : Attenuation coefficient(dB / km)
* beta2 : Group velocity dispersion coefficient(ps ^ 2 / km)
* beta3 : Third - order dispersion coefficient(ps ^ 3 / km)
* pmd : PMD parameter, polarization mode dispersion(ps / sqrt(km))
* sopmd : Second order PMD parameter(ps ^ 2 / sqrt(km))
* PDL_dB : polarization dependent losses[dB]
* dz : Integration step along the fiber length(km)->must be equal or
* larger than the fiber correlation length
* 
* **********************************************************************************************************/

void fiber_propagation_CSM(ComplexMatrix& signal_in, channel_params_t& ch_prms, modulation_setup_t& mod_stp);