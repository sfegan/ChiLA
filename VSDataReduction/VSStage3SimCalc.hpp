//-*-mode:c++; mode:font-lock;-*-

/*! \file VSStage3SimCalc.hpp

  Helper class for performing various calculations on simulations.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSStage3SimCalc.hpp,v 3.1 2009/05/13 05:48:44 matthew Exp $

*/

#ifndef VSSTAGE3SIMCALC_HPP
#define VSSTAGE3SIMCALC_HPP

#include <VSResultsSimData.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSSimEnergyWeightCalc.hpp>

namespace VERITAS
{
  class VSStage3SimCalc
  {
  public:
    VSStage3SimCalc() { }
    ~VSStage3SimCalc() { }


    void analyze(VSStage3SimDatum& data);
    void analyze(std::vector<VSStage3SimArrayTableDatum*>& tables);

    void calc(const VSLimitedErrorsHist<double, double>& h,
	      const VSLimitedErrorsHist<double, double>& htot,
	      const VSLimitedErrorsHist<double, double>& harea,
	      VSLimitedErrorsHist<double,double>& effarea,
	      VSLimitedErrorsHist<double,double>& diffrate);

    double fitThreshold(const VSLimitedErrorsHist<double, double>& h);

    void calcRate(double& rate, double& rate_err,
		  const VSLimitedErrorsHist<double, double>& effarea,
		  VSSpectrumFn* spectrum);
    void calcKernel(const VSSimple2DHist<double, double>& h,
		    VSNSpace& egy_kernel);

    void calcEnergyBias(const VSSimple2DHist<double, double>& h,
			VSLimitedErrorsHist<double, double>& egy_bias_hist,
			VSLimitedErrorsHist<double, double>& egy_rms_hist,
			VSLimitedErrorsHist<double, double>& 
			egy_log10_bias_hist,
			VSLimitedErrorsHist<double, double>& 
			egy_log10_rms_hist);


    void calcEffarea(VSStage3SimArrayTableDatum& d);
  };
}

#endif // VSSTAGE3SIMCALC_HPP
