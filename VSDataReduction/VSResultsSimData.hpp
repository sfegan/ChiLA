//-*-mode:c++; mode:font-lock;-*-

/*! \file VSStage3SimData.hpp

  Class to store results of the stage 3 analysis

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.10 $
  \date       08/16/2008

  $Id: VSResultsSimData.hpp,v 3.10 2009/11/12 23:58:41 matthew Exp $

*/

#ifndef VSSTAGE3SIMDATA_HPP
#define VSSTAGE3SIMDATA_HPP

#include <VSSimpleGraph.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSSimple2DHist.hpp>
#include <VSNSpace.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSEventData.hpp>
#include <VSSimulationData.hpp>
#include <VSAnalysisStage3Data.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSStage3SimScopeDatum
  // ==========================================================================
  class VSStage3SimScopeDatum
  {
  public:
    VSStage3SimScopeDatum();
    ~VSStage3SimScopeDatum();
    
    // Histograms of Focal Plane Parameters -----------------------------------
    VSLimitedErrorsHist<double, double> log10_N_hist;
    VSLimitedErrorsHist<double, double> fp_width_hist;
    VSLimitedErrorsHist<double, double> fp_length_hist;
    VSLimitedErrorsHist<double, double> fp_dist_hist;
    VSSimple2DHist<double, double>      fp_centroid_hist;

    // Energy Bias Histograms -------------------------------------------------
    VSSimple2DHist<double, double>      core_log10_egyerr_hist;
    VSSimple2DHist<double, double>      log10_N_log10_egyerr_hist;
    VSSimple2DHist<double, double>      fp_dist_log10_egyerr_hist;
    VSSimple2DHist<double, double>      fp_disp_log10_egyerr_hist;
    VSSimple2DHist<double, double>      nimage_log10_egyerr_hist;

    void merge(VSStage3SimScopeDatum* sim);
    void accumulate(const VSArraySimulationDatum& sim_data,
		    const VSEventScopeDatum& scope_data,
		    unsigned event_code,
		    double weight = 1.0);

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;
  };

  // ==========================================================================
  // VSStage3SimArrayDatumBase
  // ==========================================================================
  class VSStage3SimArrayDatumBase
  {
  public:
    VSStage3SimArrayDatumBase(const std::vector<unsigned>& nchan = 
			       std::vector<unsigned>());
    ~VSStage3SimArrayDatumBase();
    
    // Event Counts -----------------------------------------------------------
    double                              ntotal;
    double                              ntriggered;
    double                              nreconstructed;
    double                              ncuts_selected;
    double                              nselected;

    // TeV resolution ---------------------------------------------------------
    double                              th68;
    double                              th68_err;
    double                              th90;
    double                              th90_err;
    double                              th95;
    double                              th95_err;
    double                              thsq68;
    double                              thsq68_err;
    double                              thsq90;
    double                              thsq90_err;
    double                              thsq95;
    double                              thsq95_err;

    // Theta Histograms -------------------------------------------------------
    VSLimitedErrorsHist<double, double> thetasq_hist;
    VSSimple2DHist<double, double>      thetasq_core_hist;
    VSSimple2DHist<double, double>      thetasq_coremc_hist;

    // Energy Bias histograms -------------------------------------------------
    VSSimple2DHist<double, double>      core_log10_egyerr_hist;
    VSSimple2DHist<double, double>      coremc_log10_egyerr_hist;
    VSSimple2DHist<double, double>      log10_N2_log10_egyerr_hist;
    VSSimple2DHist<double, double>      thetasq_log10_egyerr_hist;

    // Direction Histograms ---------------------------------------------------
    VSSimple2DHist<double, double>      fov_total_hist;
    VSSimple2DHist<double, double>      fov_triggered_hist;
    VSSimple2DHist<double, double>      fov_reconstructed_hist;
    VSSimple2DHist<double, double>      fov_cuts_selected_hist;
    VSSimple2DHist<double, double>      fov_selected_hist;

    // Core Position Histograms -----------------------------------------------
    VSSimple2DHist<double, double>      core_triggered_hist;
    VSSimple2DHist<double, double>      core_reconstructed_hist;
    VSLimitedErrorsHist<double, double> core_reconstructed_err_x_hist;
    VSLimitedErrorsHist<double, double> core_reconstructed_err_y_hist;
    VSLimitedErrorsHist<double, double> core_reconstructed_err_r_hist;
    VSSimple2DHist<double, double>      core_cuts_selected_hist;
    VSLimitedErrorsHist<double, double> core_cuts_selected_err_x_hist;
    VSLimitedErrorsHist<double, double> core_cuts_selected_err_y_hist;
    VSLimitedErrorsHist<double, double> core_cuts_selected_err_r_hist;
    VSSimple2DHist<double, double>      core_selected_hist;
    VSLimitedErrorsHist<double, double> core_selected_err_x_hist;
    VSLimitedErrorsHist<double, double> core_selected_err_y_hist;
    VSLimitedErrorsHist<double, double> core_selected_err_r_hist;

    // Temporary data ---------------------------------------------------------
    std::vector< double >               m_thetasq;

    std::vector<VSStage3SimScopeDatum*>      scope;

    void setNChan(const std::vector<unsigned>& nchan);

    void merge(VSStage3SimArrayDatumBase* sim);
    void finalize();
    void accumulate(const VSArraySimulationDatum& sim_data,
		    const VSEventArrayDatum& event_data,
		    unsigned event_code,
		    double weight = 1.0);

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,ntotal);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,ntriggered);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,nreconstructed);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,ncuts_selected);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,nselected);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,th68);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,th68_err);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,th90);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,th90_err);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,th95);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,th95_err);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,thsq68);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,thsq68_err);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,thsq90);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,thsq90_err);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,thsq95);
      H5_ADDMEMBER(c,VSStage3SimArrayDatumBase,thsq95_err);    
    }
  };

  // ==========================================================================
  // VSStage3SimArrayTableDatum
  // ==========================================================================
  class VSStage3SimArrayTableDatum : public VSStage3SimArrayDatumBase
  {
  public:
    VSStage3SimArrayTableDatum(const std::vector<unsigned>& nchan = 
				std::vector<unsigned>());
    ~VSStage3SimArrayTableDatum();

    unsigned                            table_index;
    double                              egy_tev;

    double                              sampling_area;
    double                              effarea_triggered;
    double                              effarea_triggered_err;
    double                              effarea_reconstructed;
    double                              effarea_reconstructed_err;
    double                              effarea_cuts_selected;
    double                              effarea_cuts_selected_err;
    double                              effarea_selected;
    double                              effarea_selected_err;

    // Energy Reconstruction --------------------------------------------------
    VSLimitedErrorsHist<double, double> egy_counts_hist; 

    void merge(VSStage3SimArrayTableDatum* sim);
    void finalize();
    void accumulate(const VSArraySimulationDatum& sim_data,
		    const VSEventArrayDatum& event_data,
		    unsigned event_code,
		    double weight = 1.0);

    void setEnergyBinSize(double egy_binsize, double egy_lo, double egy_hi);
    
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,table_index);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,egy_tev);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,sampling_area);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_triggered);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_triggered_err);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_reconstructed);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_reconstructed_err);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_cuts_selected);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_cuts_selected_err);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_selected);
      H5_ADDMEMBER(c,VSStage3SimArrayTableDatum,effarea_selected_err);
    }
  };

  // ==========================================================================
  // VSStage3SimDatum
  // ==========================================================================
  class VSStage3SimDatum : public VSStage3SimArrayDatumBase
  {
  public:
    VSStage3SimDatum(const std::vector<unsigned>& nchan = 
		      std::vector<unsigned>());
    virtual ~VSStage3SimDatum();
    
    double                              src_ra_deg;
    double                              src_dec_deg;
    double                              obs_ra_deg;
    double                              obs_dec_deg;

    double                              gamma_rate_triggered;
    double                              gamma_rate_triggered_err;
    double                              gamma_rate_reconstructed;
    double                              gamma_rate_reconstructed_err;
    double                              gamma_rate_cuts_selected;
    double                              gamma_rate_cuts_selected_err;
    double                              gamma_rate_selected;
    double                              gamma_rate_selected_err;
    double                              crab_rate_triggered;
    double                              crab_rate_triggered_err;
    double                              crab_rate_reconstructed;
    double                              crab_rate_reconstructed_err;
    double                              crab_rate_cuts_selected;
    double                              crab_rate_cuts_selected_err;
    double                              crab_rate_selected;
    double                              crab_rate_selected_err;

    double                              proton_rate_triggered;
    double                              proton_rate_triggered_err;
    double                              proton_rate_reconstructed;
    double                              proton_rate_reconstructed_err;
    double                              proton_rate_cuts_selected; 
    double                              proton_rate_cuts_selected_err; 
    double                              proton_rate_selected;
    double                              proton_rate_selected_err;

    double                              egy_threshold_triggered;
    double                              egy_threshold_reconstructed;
    double                              egy_threshold_cuts_selected;
    double                              egy_threshold_selected;

    double                              log10_egy_binsize;
    double                              log10_egy_min;
    double                              log10_egy_max;

    // Event Counts -----------------------------------------------------------
    VSSimpleHist<unsigned>              table_id_hist;
    VSLimitedErrorsHist<double,double > egymc_count_hist;
    VSLimitedErrorsHist<double,double > egymc_total_hist;
    VSLimitedErrorsHist<double,double > egymc_triggered_hist;
    VSLimitedErrorsHist<double,double > egymc_reconstructed_hist;
    VSLimitedErrorsHist<double,double > egymc_cuts_selected_hist;
    VSLimitedErrorsHist<double,double > egymc_selected_hist;
    VSLimitedErrorsHist<double,double > egymc_fluence_hist;
    VSLimitedErrorsHist<double,double > egymc_sampling_area_hist;

    // Effective Area ---------------------------------------------------------
    VSLimitedErrorsHist<double, double> effarea_triggered_hist;
    VSLimitedErrorsHist<double, double> effarea_reconstructed_hist;
    VSLimitedErrorsHist<double, double> effarea_cuts_selected_hist;
    VSLimitedErrorsHist<double, double> effarea_selected_hist;
    
    // Differential Rate ------------------------------------------------------
    VSLimitedErrorsHist<double, double> diffrate_triggered_hist;
    VSLimitedErrorsHist<double, double> diffrate_reconstructed_hist;
    VSLimitedErrorsHist<double, double> diffrate_cuts_selected_hist;
    VSLimitedErrorsHist<double, double> diffrate_selected_hist;

    // Trigger and Selection efficiencies -------------------------------------
    VSSimple2DHist<double, double>      egymc_core_total_hist;
    VSSimple2DHist<double, double>      egymc_core_triggered_hist;
    VSSimple2DHist<double, double>      egymc_core_triggered_norm_hist;
    VSSimple2DHist<double, double>      egymc_core_triggered_eff_hist;
    VSSimple2DHist<double, double>      egymc_core_reconstructed_hist;
    VSSimple2DHist<double, double>      egymc_core_reconstructed_norm_hist;
    VSSimple2DHist<double, double>      egymc_core_reconstructed_eff_hist;
    VSSimple2DHist<double, double>      egymc_core_selected_hist;
    VSSimple2DHist<double, double>      egymc_core_selected_norm_hist;
    VSSimple2DHist<double, double>      egymc_core_selected_eff_hist;

    // TeV resolution ---------------------------------------------------------
    VSSimple2DHist<double, double>      egymc_thetasq_hist;
    VSSimple2DHist<double, double>      egymc_thetasq_norm_hist;
    VSSimple2DHist<double, double>      egymc_theta_hist;
    VSLimitedErrorsHist<double, double> egymc_thsq68_hist;
    VSLimitedErrorsHist<double, double> egymc_thsq90_hist;
    VSLimitedErrorsHist<double, double> egymc_thsq95_hist;

    // Energy Reconstruction --------------------------------------------------
    VSSimple2DHist<double, double>      egymc_egy_hist;

    VSNSpace                            egy_kernel;
    VSNSpace                            effarea;

    VSSimple2DHist<double, double>      egymc_log10_egyerr_hist;

    VSLimitedErrorsHist<double, double> egymc_bias_hist;
    VSLimitedErrorsHist<double, double> egymc_rms_hist;
    VSLimitedErrorsHist<double, double> egymc_log10_bias_hist;
    VSLimitedErrorsHist<double, double> egymc_log10_rms_hist;


    //    std::vector<VSStage3SimArrayTableDatum*> table;

    void merge(VSStage3SimDatum* sim);
    void finalize();
    void accumulate(const VSArraySimulationDatum& sim_data,
		    const VSEventArrayDatum& event_data,
		    unsigned event_code,
		    double weight = 1.0);

    void setEnergyBinSize(double egy_binsize, double egy_lo, double egy_hi);

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSStage3SimDatum,src_ra_deg);
      H5_ADDMEMBER(c,VSStage3SimDatum,src_dec_deg);
      H5_ADDMEMBER(c,VSStage3SimDatum,obs_ra_deg);
      H5_ADDMEMBER(c,VSStage3SimDatum,obs_dec_deg);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_triggered);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_triggered_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_reconstructed);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_reconstructed_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_cuts_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_cuts_selected_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,gamma_rate_selected_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_triggered);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_triggered_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_reconstructed);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_reconstructed_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_cuts_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_cuts_selected_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,crab_rate_selected_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_triggered);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_triggered_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_reconstructed);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_reconstructed_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_cuts_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_cuts_selected_err);
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_selected);      
      H5_ADDMEMBER(c,VSStage3SimDatum,proton_rate_selected_err);      
      H5_ADDMEMBER(c,VSStage3SimDatum,egy_threshold_triggered);
      H5_ADDMEMBER(c,VSStage3SimDatum,egy_threshold_reconstructed);
      H5_ADDMEMBER(c,VSStage3SimDatum,egy_threshold_cuts_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,egy_threshold_selected);
      H5_ADDMEMBER(c,VSStage3SimDatum,log10_egy_binsize);
      H5_ADDMEMBER(c,VSStage3SimDatum,log10_egy_min);
      H5_ADDMEMBER(c,VSStage3SimDatum,log10_egy_max);
    }

  private:
    // Copy Constructor -------------------------------------------------------
    VSStage3SimDatum(const VSStage3SimDatum& o);

    // Assignment Operator ----------------------------------------------------
    VSStage3SimDatum& operator= (const VSStage3SimDatum& o);
  };
}

#endif // VSSTAGE3SIMDATA_HPP
