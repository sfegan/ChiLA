//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFLaserCalc.cpp

  Calculate timing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFLaserCalc.hpp,v 3.4 2007/06/22 22:29:01 matthew Exp $

*/

#ifndef VBFLASERCALC_HPP
#define VBFLASERCALC_HPP

#include<vector>

#include "VSSimpleStat.hpp"
#include "VSSimpleHist.hpp"
#include "VSTimingCalc.hpp"
#include "VSSimpleVBF.hpp"
#include "VSAnalysisStage1.hpp"
#include "VSLaserData.hpp"

// ============================================================================
// VBFLaserCalc
// ============================================================================

namespace VERITAS
{

  class VBFLaserCalc: public VSSimpleVBFVisitor
  {
  public:
    struct ChanData
    {
      ChanData(double hist_res_time, double hist_res_signal):
	time_stat(), signal_stat(), signal_corr_stat(), nflash(),
	global_uncalibrated(), global_excluded(), global_ped(), 
	global_crate(), ev_hit(), ev_time(), ev_signal(), 
	time_hist(hist_res_time), 
	signal_hist(hist_res_signal), signal_corr_hist(hist_res_signal) { }
      
      // Statistics counters --------------------------------------------------
      VSSimpleStat2<double>   time_stat;         // Arrival time
      VSSimpleStat2<double>   signal_stat;       // Amplitude
      VSSimpleStat2<double>   signal_corr_stat;  // Amplitude corrected
                                                 // by average of all pixels
      
      // Data -----------------------------------------------------------------
      unsigned                nflash;            // Number flashes recorded

      // Global quantities ----------------------------------------------------
      bool                    global_uncalibrated;
      bool                    global_excluded;
      double                  global_ped;
      unsigned                global_crate;
      
      // Temporary event level quantities -------------------------------------
      bool                    ev_hit;     
      double                  ev_time;
      double                  ev_signal;
      
      // Histograms -----------------------------------------------------------
      VSSimpleHist<double>    time_hist;
      VSSimpleHist<double>    signal_hist;
      VSSimpleHist<double>    signal_corr_hist;
    };
    
    struct ScopeData
    {
      ScopeData(double _hist_res_time, double _hist_res_signal,
		const std::vector<unsigned>& _boards_per_crate,
		const std::vector<unsigned>& _l2_pulse_channel,
		unsigned nchan=0):
	chan(nchan,ChanData(_hist_res_time,_hist_res_signal)),
	mean_nflash_stat(), mean_signal_stat(), 
	mean_crate_time_stat(),
	max_nsample(0), 
	hist_res_time(_hist_res_time), hist_res_signal(_hist_res_signal),
	boards_per_crate(_boards_per_crate), 
	l2_pulse_channel(_l2_pulse_channel), ev_nchan_flash(),
	ev_nchan_logain(), ev_nchan_higain(),
	signal_hist(_hist_res_signal),nchan_flash_hist(1),
	nchan_logain_hist(1),nchan_higain_hist(1)
      { }
      std::vector<ChanData>               chan;

      // Statistics counters --------------------------------------------------
      VSSimpleStat2<double>               mean_nflash_stat;
      VSSimpleStat2<double>               mean_signal_stat;  
      std::vector<VSSimpleStat2<double> > mean_crate_time_stat;

      // Data -----------------------------------------------------------------
      unsigned                            max_nsample;
      double                              hist_res_time;
      double                              hist_res_signal;
      std::vector<unsigned>               boards_per_crate;
      std::vector<unsigned>               l2_pulse_channel;

      // Temporary event level quantities -------------------------------------
      unsigned                            ev_nchan_flash;
      unsigned                            ev_nchan_logain;
      unsigned                            ev_nchan_higain;

      // Histograms -----------------------------------------------------------
      VSSimpleHist<double>    signal_hist;
      VSSimpleHist<unsigned>  nchan_flash_hist;
      VSSimpleHist<unsigned>  nchan_logain_hist;
      VSSimpleHist<unsigned>  nchan_higain_hist;

    };

    typedef std::vector<ScopeData*> ArrayData;

    VBFLaserCalc(VSTimingCalc* calc,
		 double threshold_dc, unsigned threshold_nchan,
		 unsigned threshold_nscope, unsigned threshold_nflash,
		 double ped_suppress_lo, double ped_suppress_hi,
		 double singlepe_dev,
		 const std::vector<std::vector<unsigned> >& boards_per_crate,
		 const std::vector<std::vector<unsigned> >& l2_pulse_channel,
		 const VSAnalysisStage1Data& stage1_data,
		 bool fill_hist = false, bool include_lo_gain = false);

    
    virtual ~VBFLaserCalc();
    
    virtual void visitArrayTrigger(bool& veto_array_event, void* user_data,
				   uint32_t             event_num,
				   const VEventType&    event_type,
				   uint32_t             trigger_mask,
				   uint32_t             flags,
				   const VSTime&        raw_time,
				   uint32_t             at_flags,
				   uint32_t             config_mask,
				   uint32_t             num_telescopes,
				   uint32_t             num_trigger_telescopes,
				   uint32_t             run_number,
				   const uint32_t*      ten_mhz_clocks,
				   const uint32_t*      cal_count,
				   const uint32_t*      ped_count,
				   const VArrayTrigger* trigger);

    virtual void visitScopeEvent(bool& veto_scope_event, void* user_data,
				 uint32_t               event_num,
				 uint32_t               telescope_num, 
				 const VEventType&      event_type,
				 uint32_t               trigger_mask,
				 uint32_t               flags,
				 const VSTime&          raw_time,
				 uint32_t               num_samples,
				 uint32_t               num_channels_saved,
				 uint32_t               num_channels_total,
				 uint32_t               num_clock_trigger,
				 const VEvent*          event);
    
    virtual void visitChannel(bool& veto_channel, void* user_data,
			      uint32_t                  channel_num, 
			      bool                      hit, 
			      bool                      trigger);
    
    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);
    
    virtual void leaveScopeEvent(bool veto_scope_event, void* user_data);

    virtual void leaveArrayEvent(bool veto_array_event, void* user_data);
    
    void getData(VSLaserData& data) const;

    ArrayData               m_data;
    
  private:
    VBFLaserCalc(VBFLaserCalc&);
    VBFLaserCalc& operator= (const VBFLaserCalc&);

    VSTimingCalc*           m_calc;
    unsigned                m_scope_id;
    unsigned                m_nflash_scope;
    unsigned                m_nscope;

    // Settings ---------------------------------------------------------------
    double                  m_threshold_dc;
    unsigned                m_threshold_nchan;
    unsigned                m_threshold_nscope;
    unsigned                m_threshold_nflash;
    double                  m_singlepe_dev;
    bool                    m_fill_hist;
    bool                    m_include_lo_gain;

    const VSAnalysisStage1Data& m_stage1_data;
  };
  

}

#endif // VBFLASERCALC_HPP
