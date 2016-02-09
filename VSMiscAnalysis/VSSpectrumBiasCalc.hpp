//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSpectrumBiasCalc.hpp

  Class for evaluating performance of different spectral
  reconstruction methods.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSSpectrumBiasCalc.hpp,v 1.1 2010/01/19 22:44:07 matthew Exp $

*/

#ifndef VSSPECTRUMBIASCALC_HPP
#define VSSPECTRUMBIASCALC_HPP

#include <VSOptions.hpp>
#include <RandomNumbers.hpp>
#include <VSSpectrumCalc.hpp>
#include <VSSourceInjector.hpp>

namespace VERITAS
{
  class VSHistRNG
  {
  public:
    VSHistRNG(const VSLimitedErrorsHist<double,double>& hist);

    double rnd(RandomNumbers* rng) const;

    void fill(RandomNumbers* rng,
	      VSLimitedErrorsHist<double,double>& h, 
	      double mu) const;

  private:
    std::vector<std::pair<double,double> > m_icdf;
  };


  class VSSpectrumBiasCalc 
  {
  public:
    struct Options
    {
      Options(): niter(1) { }
      
      unsigned niter;
    };


    class Data
    {
    public:

      VSLimitedErrorsHist<double,double> hist;
      VSSpectrumData                     sp_data;
      VSSimple2DHist<double,double>      emc_erec_hist;
      VSLimitedErrorsHist<double,double> bias_hist;
      VSLimitedErrorsHist<double,double> rbias_hist;

      void save(VSOctaveH5Writer* writer) const;

    };

    VSSpectrumBiasCalc(RandomNumbers* rng,
		       const Options& opt = s_default_options);

    void analyze(const std::string& s3_file);

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");


    void save(VSOctaveH5Writer* writer) const;

  private:

    RandomNumbers*                    m_rng;

    VSSpectrumCalc*                   m_spectrum_calc;

    std::vector<Data>                  m_data;
    std::vector<VSLimitedErrorsHist<double,double> > m_bias_hist;
    std::vector<VSSimpleStat2<double,double> > m_bias_stat;
    
    VSLimitedErrorsHist<double,double> m_mean_bias_hist;
    VSLimitedErrorsHist<double,double> m_mean_bias_dev_hist;
    VSLimitedErrorsHist<double,double> m_mean_dev_hist;
    VSLimitedErrorsHist<double,double> m_sim_flux_hist;
    VSLimitedErrorsHist<double,double> m_sim_dfde_hist;
    VSLimitedErrorsHist<double,double> m_sim_edfde_hist;
    VSLimitedErrorsHist<double,double> m_sim_e2dfde_hist;
    VSLimitedErrorsHist<double,double> m_ntot_hist;
    VSLimitedErrorsHist<double,double> m_ncov_hist;
    VSLimitedErrorsHist<double,double> m_coverage_hist;
    

    Options                           m_options;

    static Options                    s_default_options;

  };
}


#endif // VSSPECTRUMBIASCALC_HPP
