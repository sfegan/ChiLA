//-*-mode:c++; mode:font-lock;-*-

/*! \file VSIntegralAnalysisData.hpp

  Class to store stage3 integral analysis results.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.7 $
  \date       07/29/2007

  $Id: VSIntegralAnalysisData.hpp,v 3.7 2010/10/20 03:46:18 matthew Exp $

*/

#ifndef VSINTEGRALANALYSISDATA_HPP
#define VSINTEGRALANALYSISDATA_HPP

#include <Astro.h>
#include <VSOctaveH5Reader.hpp>
#include <VSSimpleGraph.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSSimple2DHist.hpp>
#include <VSGenSpace.hpp>
#include <VSAAlgebra.hpp>
#include <VSScaledParameterCalc.hpp>
#include <VSEnergyCalcLT.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSMergedCalibrationData.hpp>
#include <VSDiagnosticsData.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSIntegralAnalysisDatum
  // ==========================================================================
  class VSIntegralAnalysisDatum 
  {
  public:

    VSIntegralAnalysisDatum();    
    ~VSIntegralAnalysisDatum() { }

    std::string       method;

    double            src_dec_J2000;
    double            src_ra_J2000;

    VSAAlgebra::Vec2D src_xy;

    double      elaptime_min;
    double      livetime_min;

    // Counts in the aperture region ------------------------------------------
    unsigned    on_counts;
    double      on_counts_err;
    double      on_rate;
    double      on_rate_err;

    // Counts in the off region ------------------------------------------
    unsigned    off_counts;
    double      off_counts_err;
    double      off_rate;
    double      off_rate_err;

    // Number of counts in the background model inside the aperture region ----
    double      bkgnd;
    double      bkgnd_err;
    double      bkgnd_rate;
    double      bkgnd_rate_err;

    // Background per unit solid angle at the source position -----------------
    double      bkgnd_density;
    double      bkgnd_density_err;
    double      bkgnd_density_rate;
    double      bkgnd_density_rate_err;

    // Total excess counts ----------------------------------------------------
    double      excess;
    double      excess_err;
    double      excess_rate;
    double      excess_rate_err;

    // Integral flux ----------------------------------------------------------
    double      flux;
    double      flux_err;
    double      flux100;
    double      flux100_err;
    double      flux316;
    double      flux316_err;
    double      flux1000;
    double      flux1000_err;
    double      flux_eth;
    double      flux_eth_err;

    // Differential flux at various energies ----------------------------------
    double      dfde;
    double      dfde_err;
    double      dfde100;
    double      dfde100_err;
    double      dfde316;
    double      dfde316_err;
    double      dfde1000;
    double      dfde1000_err;
    double      dfde_eth;
    double      dfde_eth_err;

    double      egy_threshold;

    double      alpha;
    double      significance;
    double      sigma_sqrthr;
    double      excess_rate_ul95;
    double      excess_rate_ul99;

    // Flux upper limits ------------------------------------------------------
    double      flux_ul95;
    double      flux_ul99;
    double      flux100_ul95;
    double      flux100_ul99;
    double      flux316_ul95;
    double      flux316_ul99;
    double      flux1000_ul95;
    double      flux1000_ul99;
    double      flux_eth_ul95;
    double      flux_eth_ul99;

    double      dfde_ul95;
    double      dfde_ul99;
    double      dfde100_ul95;
    double      dfde100_ul99;
    double      dfde316_ul95;
    double      dfde316_ul99;
    double      dfde1000_ul95;
    double      dfde1000_ul99;
    double      dfde_eth_ul95;
    double      dfde_eth_ul99;

    double      e2dfde_ul95;
    double      e2dfde_ul99;
    double      e2dfde100_ul95;
    double      e2dfde100_ul99;
    double      e2dfde316_ul95;
    double      e2dfde316_ul99;
    double      e2dfde1000_ul95;
    double      e2dfde1000_ul99;
    double      e2dfde_eth_ul95;
    double      e2dfde_eth_ul99;

    VSLimitedErrorsHist<double, double>  th2_bkgnd_counts_model_hist;
    VSLimitedErrorsHist<double, double>  th2_source_counts_model_hist;
    VSLimitedErrorsHist<double, double>  th2_source_model_hist;
    VSLimitedErrorsHist<double, double>  th2_counts_model_hist;

    // Effective area ---------------------------------------------------------
    VSNSpace    effarea;
    VSNSpace    psf;
    VSNSpace    drde;
    VSNSpace    edrde;

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,method);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,src_dec_J2000);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,src_ra_J2000);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,elaptime_min);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,livetime_min);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,on_counts);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,on_counts_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,on_rate);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,on_rate_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,off_counts);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,off_counts_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,off_rate);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,off_rate_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_rate);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_rate_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_density);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_density_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_density_rate);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,bkgnd_density_rate_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,excess);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,excess_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,excess_rate);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,excess_rate_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux100);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux100_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux316);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux316_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux1000);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux1000_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_eth);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_eth_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde100);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde100_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde316);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde316_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde1000);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde1000_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_eth);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_eth_err);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,egy_threshold);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,alpha);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,significance);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,sigma_sqrthr);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,excess_rate_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,excess_rate_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux100_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux100_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux316_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux316_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux1000_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux1000_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_eth_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,flux_eth_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde100_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde100_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde316_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde316_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde1000_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde1000_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_eth_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,dfde_eth_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde100_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde100_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde316_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde316_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde1000_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde1000_ul99);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde_eth_ul95);
      H5_ADDMEMBER(c,VSIntegralAnalysisDatum,e2dfde_eth_ul99);
    }
  };

  // ==========================================================================
  // VSIntegralAnalysisRunDatum
  // ==========================================================================
  class VSIntegralAnalysisRunDatum : public VSIntegralAnalysisDatum
  {
  public:
    VSIntegralAnalysisRunDatum(unsigned run_number = 0);
    virtual ~VSIntegralAnalysisRunDatum();

    unsigned    run_number;
    std::string run_start_time_string;
    std::string run_stop_time_string;

    double      run_start_time_mjd;
    double      run_stop_time_mjd;
    double      obs_dec_J2000;
    double      obs_ra_J2000;

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,run_number);
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,run_start_time_string);
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,run_stop_time_string);
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,run_start_time_mjd);
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,run_stop_time_mjd);
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,obs_dec_J2000);
      H5_ADDMEMBER(c,VSIntegralAnalysisRunDatum,obs_ra_J2000);  
    }
  };

  // ==========================================================================
  // VSIntegralAnalysisData
  // ==========================================================================
  class VSIntegralAnalysisData 
  {
  public:
    VSIntegralAnalysisData();
    virtual ~VSIntegralAnalysisData();


    // Various maps -----------------------------------------------------------
    // sky -- derotated with origin at source position
    // fov -- derotated with origin at camera center
    // cam -- non-derotated with origin at camera center

    VSSimple2DHist<double, double>       sky_on_counts_hist;
    VSSimple2DHist<double, double>       sky_off_counts_hist;
    VSSimple2DHist<double, double>       sky_bkgnd_hist;
    VSSimple2DHist<double, double>       sky_bkgnd_rate_hist;
    VSSimple2DHist<double, double>       sky_bkgnd_density_rate_hist;
    VSSimple2DHist<double, double>       sky_excess_hist;
    VSSimple2DHist<double, double>       sky_excess_rate_hist;
    VSSimple2DHist<double, double>       sky_source_hist;
    VSSimple2DHist<double, double>       sky_flux_hist;
    VSSimple2DHist<double, double>       sky_flux_ul95_hist;
    VSSimple2DHist<double, double>       sky_alpha_hist;
    VSSimple2DHist<double, double>       sky_domega_hist;
    VSSimple2DHist<double, double>       sky_significance_hist;
    VSSimple2DHist<double, double>       sky_livetime_hist;
    VSSimple2DHist<double, double>       sky_acceptance_hist;

    VSLimitedErrorsHist<double, double>  significance_hist;
    VSLimitedErrorsHist<double, double>  significance_excluded_hist;

    VSIntegralAnalysisDatum              src_data;
    VSIntegralAnalysisDatum              maxsig_data;

    double                               livetime_min;

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    void load(const std::vector<VSIntegralAnalysisDatum*>& results);
      
    static void _compose(VSOctaveH5CompositeDefinition& c)
    {

    }

    std::vector<unsigned>             m_nchan;
  };

}  // namespace VERITAS

#endif // VSINTEGRALANALYSISDATA_HPP
