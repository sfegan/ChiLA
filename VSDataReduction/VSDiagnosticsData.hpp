//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDiagnosticsData.hpp
  Diagnostics data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/09/2007

  $Id: VSDiagnosticsData.hpp,v 3.22 2009/11/05 22:52:26 matthew Exp $

*/

#ifndef VSDIAGNOSTICSDATA_HPP
#define VSDIAGNOSTICSDATA_HPP

#include<vector>

#include<VSOctaveIO.hpp>
#include<VSSimpleHist.hpp>

namespace VERITAS
{

  class VSPartialScopeDiagnosticsBase
  {
  public:
    VSPartialScopeDiagnosticsBase(unsigned nchan=0);
    virtual ~VSPartialScopeDiagnosticsBase();

    void clear();
    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;

    unsigned peak_hist_iel(unsigned ichan, unsigned ipeak) const
    {
      // Try to keep down the number of cache misses the event loop.
      // The simple "return ichan*256+ipeak;" has horrible performance
      // since each channel will lie far from the last in memory whereas
      // ipeak is usual approximately the same from channel to channel
      // so storing it as "return ichan+m_nchan*ipeak;" is much better.
      // The solution here where the "peaks" are stored in chunks of
      // four elements seems even better (thanks valgrind!)
      return (ichan+(ipeak>>2)*m_nchan)<<2 | (ipeak&0x03);
    }

    // DATA -------------------------------------------------------------------

    std::vector<unsigned> channel_in_image;
    std::vector<unsigned> channel_is_triggered;
    std::vector<unsigned> channel_is_triggered_in_largest_region;
    std::vector<unsigned> channel_is_triggered_isolated;
    std::vector<unsigned> channel_is_triggered_in_sub_threshold_event;
    std::vector<unsigned> channel_is_triggered_no_l2;
    std::vector<unsigned> channel_is_triggered_in_largest_region_no_l2;
    std::vector<unsigned> channel_is_triggered_isolated_no_l2;
    std::vector<unsigned> channel_is_triggered_in_sub_threshold_event_no_l2;
    std::vector<unsigned> channel_is_hit;
    std::vector<unsigned> channel_is_lo_gain;
    std::vector<unsigned> channel_is_hi_gain;
    std::vector<unsigned> channel_is_lo_gain_no_l2;
    std::vector<unsigned> channel_is_hi_gain_no_l2;
    std::vector<VSSimpleHist<double> > channel_raw_tij;
    std::vector<VSSimpleHist<double> > channel_image_tij;
    std::vector<VSLimitedHist<double> > channel_raw_log10_sij;
    std::vector<VSLimitedHist<double> > channel_image_log10_sij;
    std::vector<VSLimitedHist<double> > channel_logain_log10_sij;
    std::vector<unsigned> channel_peak_hist_trig;
    std::vector<unsigned> channel_peak_hist_notr;
    std::vector<unsigned> channel_peak_hist_lowg;

  protected:
    unsigned m_nchan;

  private:
    VSPartialScopeDiagnosticsBase(const VSPartialScopeDiagnosticsBase&);
    VSPartialScopeDiagnosticsBase& 
    operator=(const VSPartialScopeDiagnosticsBase&);
  };

  class VSScopeDiagnosticsBase
  {
  public:
    VSScopeDiagnosticsBase(unsigned nchan=0, unsigned nsample=0, 
			   unsigned nscope=0);
    virtual ~VSScopeDiagnosticsBase();

    void clear();
    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;

    // DATA -------------------------------------------------------------------

    // SCOPE PREFIX -----------------------------------------------------------
    unsigned scope_trigger;
    unsigned scope_sent_l3;
    unsigned scope_has_event;
    unsigned scope_event_processed;
    unsigned scope_image;
    unsigned scope_has_muon;
    VSSimpleHist<int32_t> scope_trigger_ev_num_hist;
    VSSimpleHist<int32_t> scope_sent_l3_ev_num_hist;
    VSSimpleHist<int32_t> scope_has_event_ev_num_hist;
    std::vector<VSLimitedHist<double> > scope_tdc_diff;
    VSLimitedHist<double> scope_l3_ticks_veto_vdaq_hist;

    // L2 PREFIX --------------------------------------------------------------
    VSLimitedHist<double> l2_rate_hist;
    VSLimitedHist<double> l2_high_res_rate_hist;

    // GPS PREFIX -------------------------------------------------------------
    VSLimitedHist<int32_t> gps_scope_event_dt;
    std::vector<double> gps_scope_diff_mean_ev_hist;
    std::vector<double> gps_scope_diff_dev_ev_hist;
    std::vector<double> gps_scope_diff_rms_ev_hist;

    // CHANNEL PREFIX ---------------------------------------------------------
    std::vector<unsigned> channel_is_raw_max1;
    std::vector<unsigned> channel_is_sig_max1;
    std::vector<unsigned> channel_is_raw_top3;
    std::vector<unsigned> channel_is_sig_top3;
    std::vector<unsigned> channel_is_raw_max1_no_l2;
    std::vector<unsigned> channel_is_sig_max1_no_l2;
    std::vector<unsigned> channel_is_raw_top3_no_l2;
    std::vector<unsigned> channel_is_sig_top3_no_l2;
    double* channel_mean_lotrace_mtx;
    double* channel_mean_hitrace_mtx;
    double* channel_mean_lotrace_mtx_no_l2;
    double* channel_mean_hitrace_mtx_no_l2;
    double* channel_pedestal_covariance_mtx;
    std::vector<double> channel_cfd_threshold;

    double channeMeanLotrace(unsigned ichan, unsigned isamp) const
    {
      return channel_mean_lotrace_mtx[isamp*m_nchan + ichan];
    }

    double channeMeanHitrace(unsigned ichan, unsigned isamp) const
    {
      return channel_mean_hitrace_mtx[isamp*m_nchan + ichan];
    }

    double channeMeanLotraceNoL2(unsigned ichan, unsigned isamp) const
    {
      return channel_mean_lotrace_mtx_no_l2[isamp*m_nchan + ichan];
    }

    double channeMeanHitraceNoL2(unsigned ichan, unsigned isamp) const
    {
      return channel_mean_hitrace_mtx_no_l2[isamp*m_nchan + ichan];
    }

    double channelPedestalCovariance(unsigned ichan, unsigned isamp) const
    {
      return channel_pedestal_covariance_mtx[isamp*m_nchan + ichan];
    }
    
    // CAMERA PREFIX ----------------------------------------------------------
    VSSimpleHist<int32_t> camera_ntrigger;
    VSSimpleHist<int32_t> camera_ntrigger_largest_region;
    VSSimpleHist<int32_t> camera_nimage_minus_ntrigger_largest_region;
    VSSimpleHist<int32_t> camera_ntrigger_no_l2;
    VSSimpleHist<int32_t> camera_ntrigger_largest_region_no_l2;
    VSSimpleHist<int32_t> camera_nimage_minus_ntrigger_largest_region_no_l2;
    VSLimitedHist<double> camera_log10_charge;
    VSLimitedHist<double> camera_log10_charge_no_l2;
    VSLimitedHist<double> camera_log10_charge_logain;
    VSLimitedHist<double> camera_log10_max1;
    VSLimitedHist<double> camera_log10_max1_no_l2;
    VSLimitedHist<double> camera_log10_sij_raw;
    VSLimitedHist<double> camera_log10_sij_image;
    VSLimitedHist<double> camera_log10_sij_logain;
    VSSimpleHist<int16_t> camera_nimage;
    VSSimpleHist<int16_t> camera_nimage_no_l2;

    VSSimpleHist<double> camera_width;
    VSSimpleHist<double> camera_length;
    VSSimpleHist<double> camera_psi;
    VSSimpleHist<double> camera_xc;
    VSSimpleHist<double> camera_yc;
    VSSimpleHist<double> intrinsic_width;
    VSSimpleHist<double> intrinsic_length;

    // SCOPE PREFIX -----------------------------------------------------------
    std::vector<double> scope_pos_t;
    std::vector<double> scope_pos_az;
    std::vector<double> scope_pos_zn;
    std::vector<double> scope_pos_ra;
    std::vector<double> scope_pos_dec;
    std::vector<double> scope_pos_ra_dev;
    std::vector<double> scope_pos_dec_dev;
    std::vector<double> scope_pos_l;
    std::vector<double> scope_pos_b;
    std::vector<double> scope_pos_l_dev;
    std::vector<double> scope_pos_b_dev;
    std::vector<double> scope_moon_separation;

    // MEDIAN PREFIX ----------------------------------------------------------
    std::vector<double> median_l1_rate_time;
    std::vector<double> median_l1_rate;
    std::vector<double> median_current_time;
    std::vector<double> median_current;
    std::vector<double> median_pedvar_time;
    std::vector<double> median_pedvar;

  protected:
    unsigned m_nchan;
    unsigned m_nsample;

  private:
    VSScopeDiagnosticsBase(const VSScopeDiagnosticsBase&);
    VSScopeDiagnosticsBase& operator=(const VSScopeDiagnosticsBase&);
  };

  class VSScopeDiagnosticsData: 
    public VSPartialScopeDiagnosticsBase, public VSScopeDiagnosticsBase
  {
  public:
    VSScopeDiagnosticsData(unsigned nchan=0, unsigned nsample=0, 
			   unsigned nscope=0);
    virtual ~VSScopeDiagnosticsData();

    void clear();
    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;
  };

  class VSArrayDiagnosticsBase
  {
  public:
    VSArrayDiagnosticsBase(unsigned nscope);
    virtual ~VSArrayDiagnosticsBase();

    void clear();
    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;

    // DATA -------------------------------------------------------------------

    // PACKETS AND EVENTS -----------------------------------------------------
    unsigned packets_found;
    unsigned events_found;
    unsigned events_processed;
    unsigned events_reconstructed;
    unsigned events_written;
    unsigned events_ped;
    unsigned events_l2;
    unsigned events_failed_software_trigger;

    // GPS PREFIX -------------------------------------------------------------
    uint64_t gps_ticks_elapsed;
    double gps_elaptime_sec;
    double gps_livetime_sec;
    VSLimitedHist<int32_t> gps_l3_event_dt;
    std::vector<double> gps_l3_diff_mean_ev_hist;
    std::vector<double> gps_l3_diff_dev_ev_hist;
    std::vector<double> gps_l3_diff_rms_ev_hist;

    // L3 PREFIX --------------------------------------------------------------
    unsigned has_l3;
    VSSimpleHist<int32_t> has_l3_ev_num_hist;
    VSLimitedHist<double> l3_dt;
    VSLimitedHist<double> l3_log10_dt;
    VSLimitedHist<double> l3_rate_hist;
    std::vector<unsigned> l3_events_missing;
    uint64_t l3_ticks_elapsed;
    uint64_t l3_ticks_veto_vdaq;
    uint64_t l3_ticks_veto_lev3;
    uint64_t l3_ticks_veto_both;
    VSLimitedHist<double> l3_ticks_elapsed_hist;
    VSLimitedHist<double> l3_ticks_veto_vdaq_hist;
    VSLimitedHist<double> l3_ticks_veto_lev3_hist;
    VSLimitedHist<double> l3_ticks_veto_both_hist;
    double l3_elaptime_sec;
    double l3_livetime_sec;

    // SCOPE PREFIX -----------------------------------------------------------
    VSSimpleHist<int8_t> scope_nimage;
    VSSimpleHist<int8_t> scope_nquality;
    VSSimpleHist<int8_t> scope_ntrigger;
    VSSimpleHist<int8_t> scope_nhas_event;
    VSSimpleHist<int8_t> scope_nsent_l3;

    // MOON PREFIX
    double moon_el;
    double moon_az;
    double moon_ra;
    double moon_dec;
    std::string moon_ra_string;
    std::string moon_dec_string;
    double moon_phase;
    double moon_angle;
    double moon_dphase_dt;

    // FIR PREFIX -------------------------------------------------------------

    struct FirDatum
    {
      double time;
      double ambient;
      double sky;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,FirDatum,time);
	H5_ADDMEMBER(c,FirDatum,ambient);
	H5_ADDMEMBER(c,FirDatum,sky);
      };
    };

    typedef std::vector<FirDatum> FirDataSet;

    std::vector<FirDataSet> fir;

  protected:
    unsigned m_nmask;

  private:
    VSArrayDiagnosticsBase(const VSArrayDiagnosticsBase&);
    VSArrayDiagnosticsBase& operator=(VSArrayDiagnosticsBase&);
  };

  class VSArrayDiagnosticsData: public VSArrayDiagnosticsBase
  {
  public:
    VSArrayDiagnosticsData(unsigned nscope=0);
    VSArrayDiagnosticsData(const std::vector<bool>& config_mask,
			   const std::vector<unsigned>& nchan, 
			   const std::vector<unsigned>& nsample);
    virtual ~VSArrayDiagnosticsData();

    void clear();
    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;

    // DATA -------------------------------------------------------------------

    std::vector<VSScopeDiagnosticsData*> scope;

  private:
    VSArrayDiagnosticsData(const VSArrayDiagnosticsData&);
    VSArrayDiagnosticsData operator=(const VSArrayDiagnosticsData&);
  };

}

#endif // defined VSDIAGNOSTICSDATA_HPP
