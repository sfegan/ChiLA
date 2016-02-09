//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFAnalysisStage2Diagnostics.cpp

  Stage 2 analysis (diagnostics calculation class)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 3.20 $
  \date       05/10/2006

  $Id: VBFAnalysisStage2Diagnostics.cpp,v 3.20 2010/09/19 18:26:22 matthew Exp $

*/

#include <Astro.h>
#include <VSDataReductionConstants.hpp>
#include <VSABracketedMonotonicRootFinder.hpp>

#include <VSAnalysisStage2.hpp>

using namespace VERITAS;
using namespace SEphem;

// ----------------------------------------------------------------------------
//  ___  _                       _   _
// |   \(_)__ _ __ _ _ _  ___ __| |_(_)__ ___
// | |) | / _` / _` | ' \/ _ (_-<  _| / _(_-<
// |___/|_\__,_\__, |_||_\___/__/\__|_\__/__/
//             |___/
// ----------------------------------------------------------------------------

template<typename HIST> static inline void 
spreadTheTicks(double t1, double t0, int64_t ticks, HIST& hist)
{
  // Spread the 10MHz ticks over the bins that it encompasses... say
  // for example if there is a 3 minute gap in the data then don't
  // just assign all of the deadtime since the last event to one bin
  // in the histogram.

  int bin1 = hist.bin(t1);
  int bin0 = hist.bin(t0);

  if(bin1==bin0)hist.accumulate(t1,ticks);
  else
    {
      while((bin0<bin1)&&(ticks))
	{
	  double tbin0 = hist.val(bin0);
	  double tbin1 = hist.val(bin0+1);
	  double frac = (tbin1-t0)/(t1-t0);
	  int64_t bin_ticks = int64_t(round(double(ticks)*frac));

	  if(bin_ticks > ticks)bin_ticks=ticks;
	  hist.accumulate(0.5*(tbin0+tbin1),bin_ticks); 
	  ticks -= bin_ticks;
	  t0=tbin1;
	  bin0++;
	}
      if(ticks)hist.accumulate(t1,ticks);
    }
}

VBFAnalysisStage2::PartialScopeDiagnostics::
PartialScopeDiagnostics(unsigned nchan, unsigned nsample):
  VSPartialScopeDiagnosticsBase(nchan),
  channel_mean_lo_trace(new uint32_t[nchan * nsample]),
  channel_mean_hi_trace(new uint32_t[nchan * nsample]),
  channel_mean_lo_trace_no_l2(new uint32_t[nchan * nsample]),
  channel_mean_hi_trace_no_l2(new uint32_t[nchan * nsample]),
  channel_pedestal_covariance(new uint32_t[nchan * nsample]),
  channel_pedestal_mean(new uint32_t[nchan * nsample]),
  npedestal(new uint32_t[nchan]),
  m_nsample(nsample)
{
  for(unsigned i=0;i<nchan*nsample;i++)channel_mean_lo_trace[i]=0;
  for(unsigned i=0;i<nchan*nsample;i++)channel_mean_hi_trace[i]=0;
  for(unsigned i=0;i<nchan*nsample;i++)channel_mean_lo_trace_no_l2[i]=0;
  for(unsigned i=0;i<nchan*nsample;i++)channel_mean_hi_trace_no_l2[i]=0;
  for(unsigned i=0;i<nchan*nsample;i++)channel_pedestal_covariance[i]=0;
  for(unsigned i=0;i<nchan*nsample;i++)channel_pedestal_mean[i]=0;
  for(unsigned i=0;i<nchan;i++)npedestal[i]=0;
}

VBFAnalysisStage2::PartialScopeDiagnostics::~PartialScopeDiagnostics()
{
  delete[] channel_mean_lo_trace;
  delete[] channel_mean_hi_trace;
  delete[] channel_mean_lo_trace_no_l2;
  delete[] channel_mean_hi_trace_no_l2;
  delete[] channel_pedestal_covariance;
  delete[] channel_pedestal_mean;
  delete[] npedestal;
}

VBFAnalysisStage2::PartialScopeDiagnostics& 
VBFAnalysisStage2::PartialScopeDiagnostics::
operator+=(const VBFAnalysisStage2::PartialScopeDiagnostics& o)
{
  for(unsigned ichan=0;ichan<m_nchan;ichan++)
    {
      channel_in_image[ichan]
	+= o.channel_in_image[ichan];
      channel_is_triggered[ichan]
	+= o.channel_is_triggered[ichan];
      channel_is_triggered_in_largest_region[ichan]
	+= o.channel_is_triggered_in_largest_region[ichan];
      channel_is_triggered_isolated[ichan]
	+= o.channel_is_triggered_isolated[ichan];
      channel_is_triggered_in_sub_threshold_event[ichan]
	+= o.channel_is_triggered_in_sub_threshold_event[ichan];
      channel_is_triggered_no_l2[ichan]
	+= o.channel_is_triggered_no_l2[ichan];
      channel_is_triggered_in_largest_region_no_l2[ichan]
	+= o.channel_is_triggered_in_largest_region_no_l2[ichan];
      channel_is_triggered_isolated_no_l2[ichan]
	+= o.channel_is_triggered_isolated_no_l2[ichan];
      channel_is_triggered_in_sub_threshold_event_no_l2[ichan]
	+= o.channel_is_triggered_in_sub_threshold_event_no_l2[ichan];
      channel_is_hit[ichan]
	+= o.channel_is_hit[ichan];
      channel_is_lo_gain[ichan]
	+= o.channel_is_lo_gain[ichan];
      channel_is_hi_gain[ichan]
	+= o.channel_is_hi_gain[ichan];
      channel_is_lo_gain_no_l2[ichan]
	+= o.channel_is_lo_gain_no_l2[ichan];
      channel_is_hi_gain_no_l2[ichan]
	+= o.channel_is_hi_gain_no_l2[ichan];

      channel_raw_tij[ichan] += o.channel_raw_tij[ichan];
      channel_image_tij[ichan] += o.channel_image_tij[ichan];
      channel_raw_log10_sij[ichan] += o.channel_raw_log10_sij[ichan];
      channel_image_log10_sij[ichan] += o.channel_image_log10_sij[ichan];
      channel_logain_log10_sij[ichan] += o.channel_logain_log10_sij[ichan];

      for(unsigned ipeak=0;ipeak<256;ipeak++)
	{
	  channel_peak_hist_trig[peak_hist_iel(ichan,ipeak)]
	    += o.channel_peak_hist_trig[peak_hist_iel(ichan,ipeak)];
	  channel_peak_hist_notr[peak_hist_iel(ichan,ipeak)]
	    += o.channel_peak_hist_notr[peak_hist_iel(ichan,ipeak)];
	  channel_peak_hist_lowg[peak_hist_iel(ichan,ipeak)]
	    += o.channel_peak_hist_lowg[peak_hist_iel(ichan,ipeak)];
	}

      for(unsigned isample=0;isample<m_nsample;isample++)
	{
	  channel_mean_lo_trace[ichan*m_nsample + isample]
	    += o.channel_mean_lo_trace[ichan*m_nsample + isample];
	  channel_mean_hi_trace[ichan*m_nsample + isample]
	    += o.channel_mean_hi_trace[ichan*m_nsample + isample];
	  channel_mean_lo_trace_no_l2[ichan*m_nsample + isample]
	    += o.channel_mean_lo_trace_no_l2[ichan*m_nsample + isample];
	  channel_mean_hi_trace_no_l2[ichan*m_nsample + isample]
	    += o.channel_mean_hi_trace_no_l2[ichan*m_nsample + isample];
	  channel_pedestal_covariance[ichan*m_nsample + isample]
	    += o.channel_pedestal_covariance[ichan*m_nsample + isample];
	  channel_pedestal_mean[ichan*m_nsample + isample]
	    += o.channel_pedestal_mean[ichan*m_nsample + isample];
	}

      npedestal[ichan] += o.npedestal[ichan];
    }

  return *this;
}

VBFAnalysisStage2::PartialArrayDiagnostics::
PartialArrayDiagnostics(const std::vector<bool>& config_mask,
			const std::vector<unsigned>& nchan, 
			const std::vector<unsigned>& nsample):
  scopes(nchan.size())
{
  unsigned nscope = nchan.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned ns = 0;
      if(iscope < nsample.size())ns=nsample[iscope];
      if(((iscope<config_mask.size())&&(config_mask[iscope]))||(nchan[iscope]))
	scopes[iscope]=
	  new PartialScopeDiagnostics(nchan[iscope], ns);
    }
}

VBFAnalysisStage2::PartialArrayDiagnostics::~PartialArrayDiagnostics()
{
  unsigned nscope = scopes.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)delete scopes[iscope];
}

VBFAnalysisStage2::ScopeDiagnostics::
ScopeDiagnostics(unsigned nchan, unsigned nsample, unsigned nscope):
  PartialScopeDiagnostics(nchan, nsample),
  VSScopeDiagnosticsBase(nchan, nsample, nscope),
  scope_zero_set(false), 
  scope_az_zero(), 
  scope_ra_zero(), 
  scope_l_zero(), 
  scope_zn(), 
  scope_az(), 
  scope_ra(), 
  scope_dec(), 
  scope_l(), 
  scope_b(),
  scope_gps_diff_hist()
{
  // nothing to see here
}

VBFAnalysisStage2::ScopeDiagnostics::~ScopeDiagnostics()
{
  // nothing to see here
}

void VBFAnalysisStage2::ScopeDiagnostics::
integrateEventIED(IED* ied, ScopeIED* sied)
{
  const unsigned nchan = VSPartialScopeDiagnosticsBase::m_nchan;

  int32_t event_num = ied->event_num;

  if(sied->has_l3)
    {
      if(sied->l3_l2_received)
	{
	  scope_trigger++;
	  scope_trigger_ev_num_hist.accumulate(event_num);
	}

      if(sied->counts_l2>=0)
	{
	  l2_rate_hist.accumulate(ied->event_time_hist, 
				  unsigned(sied->counts_l2));
	  if((ied->has_event_time)
	     &&(ied->event_time_sec<3600)&&(ied->event_time_sec>=0))
	    spreadTheTicks(ied->event_time_sec, ied->last_event_time_sec,
			   sied->counts_l2, l2_high_res_rate_hist);
	}
      else l2_rate_hist.accumulate(-2, 1); //unsigned(sied->counts_l2));
    }

  if(ied->has_l3)
    {
      if(sied->l3_sent)
	{
	  scope_sent_l3++;
	  scope_sent_l3_ev_num_hist.accumulate(event_num);
	}
    }

  if(sied->has_vdaq)
    {
      scope_has_event++;
      scope_has_event_ev_num_hist.accumulate(event_num);

      if(ied->has_event_time)
	{
	  int32_t dt = sied->gps_dt;
	  if(sied->gps_dt > int64_t(gps_scope_event_dt.hiLimit()))
	    dt = gps_scope_event_dt.hiLimit()+1;
	  else if(sied->gps_dt < int64_t(gps_scope_event_dt.loLimit()))
	    dt = gps_scope_event_dt.loLimit()-1;
	  gps_scope_event_dt.accumulate(dt);

	  unsigned iebin = ied->event_num/1000;
	  if(iebin<100000)
	    {
	      static unsigned mgd_warn = 0;
	      if(iebin >= scope_gps_diff_hist.size())
		scope_gps_diff_hist.resize(iebin+1);
	      scope_gps_diff_hist[iebin].accumulate(sied->gps_dt);
	      if((fabs(sied->gps_dt)>1e9)&&(mgd_warn++<100))
		std::cerr << "MASSIVE GPS DIFF: " << ied->event_num 
			  << ' ' << ied->event_time.toString() 
			  << " T" << sied->iscope+1
			  << " = " << sied->gps_dt << " ns\n";
	    }
	}

      if(sied->processed_scope_event)
	{
	  scope_event_processed++; // For use by scope veto

	  if(sied->l3_l2_received)
	    {
	      if(sied->sij_max1_raw_j != nchan)
		{
		  channel_is_raw_max1[sied->sij_max1_raw_j]++;
		  channel_is_raw_top3[sied->sij_max1_raw_j]++;
		  if(sied->sij_max2_raw_j != nchan)
		    {
		      channel_is_raw_top3[sied->sij_max2_raw_j]++;
		      if(sied->sij_max3_raw_j != nchan)
		      channel_is_raw_top3[sied->sij_max3_raw_j]++;
		    }
		}
	    }
	  else
	    {
	      if(sied->sij_max1_raw_j != nchan)
		{
		  channel_is_raw_max1_no_l2[sied->sij_max1_raw_j]++;
		  channel_is_raw_top3_no_l2[sied->sij_max1_raw_j]++;
		  if(sied->sij_max2_raw_j != nchan)
		    {
		      channel_is_raw_top3_no_l2[sied->sij_max2_raw_j]++;
		      if(sied->sij_max3_raw_j != nchan)
			channel_is_raw_top3_no_l2[sied->sij_max3_raw_j]++;
		    }
		}
	    }

	  if(sied->total_signal>0)
	    {
	      if(sied->nchan_logain)
		{
		  camera_log10_charge_logain.
		    accumulate(log10(sied->total_signal));
		}

	      if(sied->l3_l2_received)
		{
		  camera_log10_charge.accumulate(log10(sied->total_signal));
		}
	      else
		{
		  camera_log10_charge_no_l2.
		    accumulate(log10(sied->total_signal));
		}
	    }

	  int32_t blah = 
	    int32_t(sied->nchan_image)-int32_t(sied->largest_region_nchan);
	  if(sied->l3_l2_received)
	    {
	      camera_nimage_minus_ntrigger_largest_region.
		accumulate(blah);
	      camera_ntrigger_largest_region.
		accumulate(sied->largest_region_nchan);
	      camera_ntrigger.
		accumulate(sied->nchan_trigger);
	      camera_nimage.accumulate(sied->nchan_image);
	    }
	  else
	    {
	      camera_nimage_minus_ntrigger_largest_region_no_l2.
		accumulate(blah);
	      camera_ntrigger_largest_region_no_l2.
		accumulate(sied->largest_region_nchan);
	      camera_ntrigger_no_l2.
		accumulate(sied->nchan_trigger);
	      camera_nimage_no_l2.accumulate(sied->nchan_image);
	    }

	  if((sied->has_azzn)
	     &&(ied->has_event_time)
	     &&(ied->event_time_sec>=0)&&(ied->event_time_sec<3600))
	    {
	      unsigned ibin = unsigned(floor(ied->event_time_sec/10));

	      if(scope_zn.size() <= ibin)
		{
		  scope_zn.resize(ibin+1);
		  scope_az.resize(ibin+1);
		  scope_ra.resize(ibin+1);
		  scope_dec.resize(ibin+1);
		  scope_l.resize(ibin+1);
		  scope_b.resize(ibin+1);
		}
	      
	      if(!scope_zero_set)
		{
		  scope_az_zero   = sied->az_rad;
		  scope_ra_zero   = sied->ra_rad;
		  scope_l_zero    = sied->l_rad;
		  scope_zero_set  = true;
		}

	      Angle _az(sied->az_rad-scope_az_zero);
	      scope_zn[ibin].accumulate(sied->zn_rad);
	      scope_az[ibin].accumulate(_az.radPM());

	      Angle _ra(sied->ra_rad-scope_ra_zero);
	      scope_dec[ibin].accumulate(sied->dec_rad);
	      scope_ra[ibin].accumulate(_ra.radPM());

	      Angle _l(sied->l_rad-scope_l_zero);
	      scope_b[ibin].accumulate(sied->b_rad);
	      scope_l[ibin].accumulate(_l.radPM());
	    }

	  if(sied->has_image)
	    {
	      scope_image++;

	      if(sied->l3_l2_received)
		{
		  if(sied->sij_max1_sig_j != nchan)
		    {
		      channel_is_sig_max1[sied->sij_max1_sig_j]++;
		      channel_is_sig_top3[sied->sij_max1_sig_j]++;

		      if(sied->sij_max1_sig > 0)
			camera_log10_max1.
			  accumulate(log10(sied->sij_max1_sig));

		      if(sied->sij_max2_sig_j != nchan)
			{
			  channel_is_sig_top3[sied->sij_max2_sig_j]++;
			  if(sied->sij_max3_sig_j != nchan)
			    channel_is_sig_top3[sied->sij_max3_sig_j]++;
			}
		    }
		}
	      else
		{
		  if(sied->sij_max1_sig_j != nchan)
		    {
		      channel_is_sig_max1_no_l2[sied->sij_max1_sig_j]++;
		      channel_is_sig_top3_no_l2[sied->sij_max1_sig_j]++;

		      if(sied->sij_max1_sig > 0)
			camera_log10_max1_no_l2.
			  accumulate(log10(sied->sij_max1_sig));

		      if(sied->sij_max2_sig_j != nchan)
			{
			  channel_is_sig_top3_no_l2[sied->sij_max2_sig_j]++;
			  if(sied->sij_max3_sig_j != nchan)
			    channel_is_sig_top3_no_l2[sied->sij_max3_sig_j]++;
			}
		    }
		}

	      if(sied->has_moments && sied->used_in_reconstruction)
		{
		  camera_xc.accumulate(sied->fp_xc);
		  camera_yc.accumulate(sied->fp_yc);
		  camera_length.accumulate(sied->fp_length);
		  camera_width.accumulate(sied->fp_width);
		  if(sied->fp_length>0)
		    camera_psi.accumulate(sied->fp_psi);
		  intrinsic_length.accumulate(sied->intrinsic_length);
		  intrinsic_width.accumulate(sied->intrinsic_width);
		}
	    }
	}
    }

  if(sied->has_muon)
    {
      scope_has_muon++;
    }
}

void VBFAnalysisStage2::ScopeDiagnostics::
integratePartial(const std::list<PartialScopeDiagnostics*>& sdiags)
{
  // --------------------------------------------------------------------------
  // Integrate counts etc
  // --------------------------------------------------------------------------

  for(std::list<PartialScopeDiagnostics*>::const_iterator 
	isdiag = sdiags.begin(); isdiag != sdiags.end(); isdiag++)
    *this += **isdiag;
}
      
void VBFAnalysisStage2::ScopeDiagnostics::
finalize(const ScopeMergedCal* scal,
	 const VSRunInfoData& run_info,
	 const VSMiscellaneousDBData::Scope* db_scope_data,
	 const std::vector<VATime>* nsb_time,
	 const std::vector<double>* nsb_median,
	 const SEphem::SphericalCoords& earth_position)
{
  const unsigned nchan = VERITAS::VSScopeDiagnosticsBase::m_nchan;
  const unsigned nsample = VERITAS::VSScopeDiagnosticsBase::m_nsample;

  for(unsigned ichan=0;ichan<nchan;ichan++)
    {
      VERITAS::VSScopeDiagnosticsBase::camera_log10_sij_raw += 
	channel_raw_log10_sij[ichan];
      VERITAS::VSScopeDiagnosticsBase::camera_log10_sij_image += 
	channel_image_log10_sij[ichan];
      VERITAS::VSScopeDiagnosticsBase::camera_log10_sij_logain += 
	channel_logain_log10_sij[ichan];
    }

  // --------------------------------------------------------------------------
  // Integrate GPS difference histograms
  // --------------------------------------------------------------------------

  {
    unsigned nbin = scope_gps_diff_hist.size();
    gps_scope_diff_mean_ev_hist.resize(nbin);
    gps_scope_diff_dev_ev_hist.resize(nbin);
    gps_scope_diff_rms_ev_hist.resize(nbin);
    for(unsigned ibin=0;ibin<nbin;ibin++)
    {
      gps_scope_diff_mean_ev_hist[ibin] = scope_gps_diff_hist[ibin].mean();
      gps_scope_diff_dev_ev_hist[ibin]  = scope_gps_diff_hist[ibin].dev();
      gps_scope_diff_rms_ev_hist[ibin]  = 
	sqrt(scope_gps_diff_hist[ibin].chi2(0));
    }
  }

  // --------------------------------------------------------------------------
  // Integrate mean LO and HI gain pulse shapes
  // --------------------------------------------------------------------------

  for(unsigned ichan=0;ichan<nchan;ichan++)
    channel_is_hit[ichan] = 
      channel_is_lo_gain[ichan] + channel_is_hi_gain[ichan]
      + channel_is_lo_gain_no_l2[ichan] + channel_is_hi_gain_no_l2[ichan];

  channel_mean_lotrace_mtx = new double[nchan*nsample];
  channel_mean_hitrace_mtx = new double[nchan*nsample];
  channel_mean_lotrace_mtx_no_l2 = new double[nchan*nsample];
  channel_mean_hitrace_mtx_no_l2 = new double[nchan*nsample];
  channel_pedestal_covariance_mtx = new double[nchan*nsample];
  
  for(unsigned ichan=0;ichan<nchan;ichan++)
    for(unsigned isample=0;isample<nsample;isample++)
      {
	channel_mean_lotrace_mtx[isample*nchan+ichan] = 
	  double(channel_mean_lo_trace[ichan*nsample+isample])
	  /double(channel_is_lo_gain[ichan]);
	channel_mean_hitrace_mtx[isample*nchan+ichan] = 
	  double(channel_mean_hi_trace[ichan*nsample+isample])
	  /double(channel_is_hi_gain[ichan]);
	channel_mean_lotrace_mtx_no_l2[isample*nchan+ichan] = 
	  double(channel_mean_lo_trace_no_l2[ichan*nsample+isample])
	  /double(channel_is_lo_gain_no_l2[ichan]);
	channel_mean_hitrace_mtx_no_l2[isample*nchan+ichan] = 
	  double(channel_mean_hi_trace_no_l2[ichan*nsample+isample])
	  /double(channel_is_hi_gain_no_l2[ichan]);

	if(npedestal[ichan])
	  {
	    int64_t s2 = channel_pedestal_covariance[ichan*nsample+isample];

	    int64_t s1 = 0;
	    for(unsigned jsample=0;jsample<nsample-isample;jsample++)
	      {
		unsigned ibase = ichan*nsample+jsample;
		s1 += ( uint64_t(channel_pedestal_mean[ibase])
			* uint64_t(channel_pedestal_mean[ibase + isample]) );
	      }

	    s1 /= npedestal[ichan];

	    double N = npedestal[ichan]*(nsample-isample);

	    channel_pedestal_covariance_mtx[isample*nchan+ichan] 
	      = double(s2-s1)/N;
	  }
	else
	  {
	    channel_pedestal_covariance_mtx[isample*nchan+ichan] = 0;
	  }
      }

  // --------------------------------------------------------------------------
  // Calculate CFD trigger thresholds
  // --------------------------------------------------------------------------

  const double threshold_fraction = 0.5;
  const double cfd_mv_per_fadc_dc = 2000.0 / 255.0 * 4.24 / 4.85;

  if((scal)&&(scal->channel.size()>=nchan))
    for(unsigned ichan=0;ichan<nchan;ichan++)
      {
	unsigned last_dc = 255;
	double last_fraction = 1.0;
	for(unsigned idc=0;idc<256;idc++)
	  {
	    unsigned ntrig = 
	      channel_peak_hist_trig[peak_hist_iel(ichan,255-idc)];
	    unsigned ntot = 
	      ntrig + channel_peak_hist_notr[peak_hist_iel(ichan,255-idc)];
	    double fraction = 1.0;
	    if(ntot)fraction = double(ntrig)/double(ntot);
	    if(fraction < threshold_fraction)break;
	    last_dc = 255-idc;
	    last_fraction = fraction;
	  }

	double fraction = 0.0;
	if(last_dc)
	  {
	    unsigned ntrig = 
	      channel_peak_hist_trig[peak_hist_iel(ichan,255-last_dc)];
	    unsigned ntot = 
	      ntrig + channel_peak_hist_notr[peak_hist_iel(ichan,255-last_dc)];
	    if(ntot)fraction = double(ntrig)/double(ntot);
	  }

	double threshold_dc = 
	  double(last_dc) - 1.0
	  + (threshold_fraction-fraction)/(last_fraction-fraction) 
	  - scal->channel[ichan].ped_hi;

	channel_cfd_threshold[ichan] = cfd_mv_per_fadc_dc * threshold_dc;
      }

  // --------------------------------------------------------------------------
  // Median Current
  // --------------------------------------------------------------------------

  if((db_scope_data)&&(!db_scope_data->hv_status.empty()))
    {
      unsigned nmeas = db_scope_data->hv_status.size();
      median_current_time.resize(nmeas);
      median_current.resize(nmeas);
      for(unsigned imeas=0;imeas<nmeas;imeas++)
	{
	  unsigned nchan = db_scope_data->hv_status[imeas].chan.size();
	  std::vector<double> current(nchan);
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    current[ichan] = 
	      db_scope_data->hv_status[imeas].chan[ichan].current;
	  median_current_time[imeas] = 
	    double(db_scope_data->hv_status[imeas].timestamp
		   - run_info.lo_event_time)/60.0/1e9;
	  median_current[imeas]      = median(current);
	}
    }

  // --------------------------------------------------------------------------
  // Median pedestal variance
  // --------------------------------------------------------------------------

  if((nsb_time)&&(!nsb_time->empty())
     &&(nsb_median)&&(nsb_median->size() <= nsb_time->size()))
    {
      unsigned nmeas = nsb_time->size();
      median_pedvar_time.resize(nmeas);
      median_pedvar.resize(nmeas);
      for(unsigned imeas=0;imeas<nmeas;imeas++)
	{
	  median_pedvar_time[imeas] =
	    double((*nsb_time)[imeas] - run_info.lo_event_time)/60.0/1e9;
	  median_pedvar[imeas] = (*nsb_median)[imeas];
	}
    }

  // --------------------------------------------------------------------------
  // Median L1 Rate
  // --------------------------------------------------------------------------

  if((db_scope_data)&&(!db_scope_data->l1_rate.empty()))
    {
      unsigned nmeas = db_scope_data->l1_rate.size();
      median_l1_rate_time.resize(nmeas);
      median_l1_rate.resize(nmeas);
      for(unsigned imeas=0;imeas<nmeas;imeas++)
	{
	  median_l1_rate_time[imeas] = 
	    double(db_scope_data->l1_rate[imeas].timestamp
		   - run_info.lo_event_time)/60.0/1e9;
	  median_l1_rate[imeas]      = 
	    median(db_scope_data->l1_rate[imeas].rate);
	}
    }

  // --------------------------------------------------------------------------
  // L3 Veto From Telescope
  // --------------------------------------------------------------------------

  if((db_scope_data)&&(!db_scope_data->l3_scope.empty()))
    {
      unsigned nmeas = db_scope_data->l3_scope.size();
      uint32_t last_scaler=0;
      double last_time = 0;
      for(unsigned imeas=0;imeas<nmeas;imeas++)
	{
	  uint32_t cur_scaler =
	    db_scope_data->l3_scope[imeas].scaler_vdacq_busy;
	  double cur_time =
	    double(db_scope_data->l3_scope[imeas].timestamp
		   - run_info.lo_event_time)/60.0/1e9;
	  uint32_t ticks = 0;
	  if(cur_scaler >= last_scaler)ticks = cur_scaler - last_scaler;
	  else ticks = cur_scaler + (0xFFFFFFFF-last_scaler) + 1;
	  spreadTheTicks(cur_time, last_time, ticks,
			 scope_l3_ticks_veto_vdaq_hist);
	  last_time = cur_time;
	  last_scaler = cur_scaler;
	}
    }

  // --------------------------------------------------------------------------
  // Add max1 to top3
  // --------------------------------------------------------------------------

  for(unsigned ichan=0;ichan<nchan;ichan++)
    {
      channel_is_raw_top3[ichan] += channel_is_raw_max1[ichan];
      channel_is_sig_top3[ichan] += channel_is_sig_max1[ichan];
      channel_is_raw_top3_no_l2[ichan] += channel_is_raw_max1_no_l2[ichan];
      channel_is_sig_top3_no_l2[ichan] += channel_is_sig_max1_no_l2[ichan];
    }

  // --------------------------------------------------------------------------
  // Pointing
  // --------------------------------------------------------------------------

  scope_pos_t.resize(scope_zn.size());
  scope_pos_az.resize(scope_zn.size());
  scope_pos_zn.resize(scope_zn.size());
  scope_pos_ra.resize(scope_zn.size());
  scope_pos_dec.resize(scope_zn.size());
  scope_pos_ra_dev.resize(scope_zn.size());
  scope_pos_dec_dev.resize(scope_zn.size());
  scope_pos_l.resize(scope_zn.size());
  scope_pos_b.resize(scope_zn.size());
  scope_pos_l_dev.resize(scope_zn.size());
  scope_pos_b_dev.resize(scope_zn.size());
  scope_moon_separation.resize(scope_zn.size());

  bool az_wrap = false;
  bool ra_wrap = false;
  bool l_wrap = false;
  unsigned nbin=scope_zn.size();
  for(unsigned ibin=0;ibin<nbin;ibin++)
    {
      scope_pos_t[ibin]   = (double(ibin)+0.5)*10.0;
      scope_pos_az[ibin]  = Angle::toDeg(scope_az[ibin].mean()+scope_az_zero);
      scope_pos_zn[ibin]  = Angle::toDeg(scope_zn[ibin].mean());
      scope_pos_ra[ibin]  = Angle::toDeg(scope_ra[ibin].mean()+scope_ra_zero);
      scope_pos_dec[ibin] = Angle::toDeg(scope_dec[ibin].mean());
      scope_pos_ra_dev[ibin]  = Angle::toDeg(scope_ra[ibin].dev());
      scope_pos_dec_dev[ibin] = Angle::toDeg(scope_dec[ibin].dev());
      scope_pos_l[ibin]   = Angle::toDeg(scope_l[ibin].mean()+scope_l_zero);
      scope_pos_b[ibin]   = Angle::toDeg(scope_b[ibin].mean());
      scope_pos_l_dev[ibin] = Angle::toDeg(scope_l[ibin].dev());
      scope_pos_b_dev[ibin] = Angle::toDeg(scope_b[ibin].dev());
      if(scope_pos_az[ibin] >= 360.0)az_wrap = true;
      if(scope_pos_ra[ibin] >= 360.0)ra_wrap = true;
      if(scope_pos_l[ibin]  >= 360.0)l_wrap  = true;

      if(scope_ra[ibin].count() > 0)
	{
	  VSTime t = run_info.lo_event_time + 
	    int64_t(scope_pos_t[ibin])*INT64_C(1000000000);
	  double mjd = t.getMJDDbl();

	  SphericalCoords moon_radec;
	  Astro::moonRaDecApparent(mjd, moon_radec);

	  SphericalCoords scope_radec;
	  scope_radec.setLatLongDeg(scope_pos_dec[ibin],scope_pos_ra[ibin]);

	  scope_moon_separation[ibin] = 
	    scope_radec.separation(moon_radec).deg();
	}
      else
	{
	  scope_moon_separation[ibin] = -1;
	}
	
    }
  
  if(az_wrap)for(unsigned ibin=0;ibin<nbin;ibin++)scope_pos_az[ibin] -= 360.0;
  if(ra_wrap)for(unsigned ibin=0;ibin<nbin;ibin++)scope_pos_ra[ibin] -= 360.0;
  if(l_wrap)for(unsigned ibin=0;ibin<nbin;ibin++)scope_pos_l[ibin]   -= 360.0;
}

SphericalCoords VBFAnalysisStage2::ArrayDiagnostics::
getMeanRaDec() const
{
  VSSimpleStat2<double> all_ra;
  VSSimpleStat2<double> all_dec;

  bool all_zero_found = false;
  double all_ra_zero = 0;
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(scope[iscope] == 0)continue;

      double scope_ra_zero = scope[iscope]->scope_ra_zero;
      if(!all_zero_found)
	{
	  all_ra_zero = scope_ra_zero;
	  all_zero_found = true;
	}

      VSSimpleStat2<double> scope_ra;
      VSSimpleStat2<double> scope_dec;
      for(unsigned ibin=0;
	  ibin<scope[iscope]->scope_zn.size();ibin++)
	{
	  scope_ra   += scope[iscope]->scope_ra[ibin];
	  scope_dec  += scope[iscope]->scope_dec[ibin];
	}

      all_ra.accumulate(scope_ra.mean() + scope_ra_zero - all_ra_zero);
      all_dec.accumulate(scope_dec.mean());
    }

  double ra = all_ra.mean() + all_ra_zero;
  double dec = all_dec.mean();
  if(ra<0)ra+=2.0*M_PI;
  SphericalCoords radec;
  radec.setLatLongRad(dec,ra);

  return radec;
}

void VBFAnalysisStage2::ScopeDiagnostics::
save(VSOctaveH5WriterStruct* s, bool slow_diagnostics) const
{
  VSPartialScopeDiagnosticsBase::save(s, slow_diagnostics);
  VSScopeDiagnosticsBase::save(s, slow_diagnostics);
}

VBFAnalysisStage2::ArrayDiagnostics::
ArrayDiagnostics(const std::vector<bool>& config_mask,
		 const std::vector<unsigned>& nchan, 
		 const std::vector<unsigned>& nsample,
		 const std::pair<double,double>& limited_dt_range):
  VSArrayDiagnosticsBase(nchan.size()), m_limited_dt_range(limited_dt_range),
  m_limited_dt(), l3_gps_diff_hist(),
  scope()
{
  unsigned nscope = nchan.size();
  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned ns = 0;
      if(iscope < nsample.size())ns=nsample[iscope];
      if(((iscope<config_mask.size())&&(config_mask[iscope]))
	 ||(nchan[iscope]))
	scope[iscope]=
	  new ScopeDiagnostics(nchan[iscope], ns, nscope);
    }
}

VBFAnalysisStage2::ArrayDiagnostics::
~ArrayDiagnostics()
{
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)delete scope[iscope];
}

void VBFAnalysisStage2::ArrayDiagnostics::integrateEventIED(IED* ied)
{
  events_found++;

  if(ied->has_event_time)
    {
      double dt = ied->event_time_sec-ied->last_event_time_sec;
      if((dt>0)&&(dt<3600))
	{
	  l3_log10_dt.accumulate(log10(dt));
	  if(dt<10)l3_dt.accumulate(dt);
	  else l3_dt.accumulate(-1);
	}

      if((dt>m_limited_dt_range.first)&&(dt<m_limited_dt_range.second))
	m_limited_dt.accumulate(dt);

      l3_rate_hist.accumulate(ied->event_time_hist);
    }

  if(ied->event_type == ET_PED)
    events_ped++;
  else if(ied->event_type == ET_L2)
    events_l2++;

#warning TEMPORARY ASSERT
  //  vsassert((ied->event_sent_trigger_mask<m_nmask)
  //	 &&(ied->event_has_vdaq_mask<m_nmask));
  if((ied->event_sent_trigger_mask<m_nmask)
     &&(ied->event_has_vdaq_mask<m_nmask))
    {
//       events_mask_mtx[ied->event_sent_trigger_mask*m_nmask
// 		      + ied->event_has_vdaq_mask]++;
    }
  else
    {
      std::cout << "Bad event mask: " << ied->event_num
		<< ' ' << ied->event_sent_trigger_mask
		<< ' ' << ied->event_has_vdaq_mask
		<< ' ' << m_nmask << '\n';
    }

  if(ied->event_failed_software_trigger)
    events_failed_software_trigger++;

  if(ied->processed_array_event)
    {
      events_processed++;
      if(ied->reconstruction_successful)events_reconstructed++;
      if(ied->write_event)events_written++;  
    }

  if((ied->event_num - ied->last_event_num > 1)&&
     (ied->event_num - ied->last_event_num < 10000))
    for(unsigned imissing=ied->last_event_num+1;
	imissing<ied->event_num; imissing++)
      l3_events_missing.push_back(imissing);  
      
  if(ied->has_l3)
    {
      has_l3++;
      has_l3_ev_num_hist.accumulate(ied->event_num);
      if(ied->ticks_elap >= 0)l3_ticks_elapsed   += ied->ticks_elap;
      if(ied->ticks_vdaq >= 0)l3_ticks_veto_vdaq += ied->ticks_vdaq;
      if(ied->ticks_lev3 >= 0)l3_ticks_veto_lev3 += ied->ticks_lev3;
      if(ied->ticks_both >= 0)l3_ticks_veto_both += ied->ticks_both;

      if(ied->has_event_time)
	{
	  int32_t dt = ied->gps_dt;
	  if(ied->gps_dt > int64_t(gps_l3_event_dt.hiLimit()))
	    dt = gps_l3_event_dt.hiLimit()+1;
	  else if(ied->gps_dt < int64_t(gps_l3_event_dt.loLimit()))
	    dt = gps_l3_event_dt.loLimit()-1;
	  gps_l3_event_dt.accumulate(dt);

	  unsigned iebin = ied->event_num/1000;
	  if(iebin<100000)
	    {
	      if(iebin >= l3_gps_diff_hist.size())
		l3_gps_diff_hist.resize(iebin+1);
	      l3_gps_diff_hist[iebin].accumulate(ied->gps_dt);
	    }

	  spreadTheTicks(ied->event_time_hist, ied->last_event_time_hist, 
			 ied->ticks_elap, l3_ticks_elapsed_hist);
	  spreadTheTicks(ied->event_time_hist, ied->last_event_time_hist,
			 ied->ticks_vdaq, l3_ticks_veto_vdaq_hist);
	  spreadTheTicks(ied->event_time_hist, ied->last_event_time_hist,
			 ied->ticks_lev3, l3_ticks_veto_lev3_hist);
	  spreadTheTicks(ied->event_time_hist, ied->last_event_time_hist,
			 ied->ticks_both, l3_ticks_veto_both_hist);
	}
    }

  if(ied->processed_array_event)
    {
      scope_nimage.accumulate(ied->nscope_image);
      scope_nquality.accumulate(ied->nscope_quality);
    }

  if(ied->has_l3)
    {
      scope_ntrigger.accumulate(ied->nscope_trigger);
      scope_nsent_l3.accumulate(ied->nscope_sent_l3);
    }

  scope_nhas_event.accumulate(ied->nscope_has_event);

  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope])
      {
	scope[iscope]->integrateEventIED(ied,ied->scope[iscope]);
	if(ied->scope[iscope]->l3_l2_received)
	  for(unsigned jscope=0;jscope<iscope;jscope++)
	    if((scope[jscope])&&(ied->scope[jscope]->l3_l2_received))
	      {
		int32_t idiff = ( int32_t(ied->scope[iscope]->l3_tdc)
				  - int32_t(ied->scope[jscope]->l3_tdc) );
		double diff = double(idiff)*TDC_CONV;
		scope[iscope]->scope_tdc_diff[jscope].accumulate(diff);
		scope[jscope]->scope_tdc_diff[iscope].accumulate(-diff);
	      }
      }
}

void VBFAnalysisStage2::ArrayDiagnostics::
integratePartial(const std::list<PartialArrayDiagnostics*>& diags)
{
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope])
      {
	std::list<PartialScopeDiagnostics*> sdiags;
	for(std::list<PartialArrayDiagnostics*>::const_iterator 
	      idiag=diags.begin(); idiag!=diags.end(); idiag++)
	  sdiags.push_back((*idiag)->partialScopeDiagnostics(iscope));
	scope[iscope]->
	  integratePartial(sdiags);
      }
}

class TauOverDeltaFunction
{
public:
  TauOverDeltaFunction(double lhs): m_lhs(lhs) { }
  double operator() (double x) { return m_lhs - x + 1.0/(exp(1.0/x)-1); }
private:
  double m_lhs;
};

void VBFAnalysisStage2::ArrayDiagnostics::
finalize(const ArrayMergedCal* cal,
	 const VSRunInfoData& run_info,
	 const VSMiscellaneousDBData* db_data,
	 const VBFTimeDepNSB::Data* nsb,
	 const SEphem::SphericalCoords& earth_position)
{
  // --------------------------------------------------------------------------
  // FIR data
  // --------------------------------------------------------------------------

  if((db_data)&&(!db_data->fir.empty()))
    {
      unsigned nfir = db_data->fir.size();
      fir.resize(nfir);
      for(unsigned ifir=0;ifir<nfir;ifir++)
	if(!db_data->fir[ifir].empty())
	  {
	    unsigned nmeas = db_data->fir[ifir].size();
	    fir[ifir].resize(nmeas);
	    for(unsigned imeas=0;imeas<nmeas;imeas++)
	      {
		fir[ifir][imeas].time    =
		  double(db_data->fir[ifir][imeas].timestamp
			 - run_info.lo_event_time)/60.0/1e9;
		fir[ifir][imeas].ambient = db_data->fir[ifir][imeas].ambient;
		fir[ifir][imeas].sky     = db_data->fir[ifir][imeas].sky;
	      }
	  }
    }

  // --------------------------------------------------------------------------
  // Integrate GPS difference histograms
  // --------------------------------------------------------------------------

  l3_elaptime_sec = double(l3_ticks_elapsed)*1e-7;
  l3_livetime_sec = double(l3_ticks_elapsed-l3_ticks_veto_both)*1e-7;

  gps_ticks_elapsed = 
    (run_info.hi_event_time-run_info.lo_event_time)/UINT64_C(100);
  gps_elaptime_sec = double(gps_ticks_elapsed)*1e-7;
  gps_livetime_sec = 0;

  unsigned nbin = l3_gps_diff_hist.size();
  gps_l3_diff_mean_ev_hist.resize(nbin);
  gps_l3_diff_dev_ev_hist.resize(nbin);
  gps_l3_diff_rms_ev_hist.resize(nbin);
  for(unsigned ibin=0;ibin<nbin;ibin++)
    {
      gps_l3_diff_mean_ev_hist[ibin] = l3_gps_diff_hist[ibin].mean();
      gps_l3_diff_dev_ev_hist[ibin]  = l3_gps_diff_hist[ibin].dev();
      gps_l3_diff_rms_ev_hist[ibin]  = 
	sqrt(l3_gps_diff_hist[ibin].chi2(0));
    }
  
  // --------------------------------------------------------------------------
  // Calculate DT-based livetime
  // ------------------------------------------------------------------------

  double dt_mean = m_limited_dt.mean();
  double delta   = m_limited_dt_range.second - m_limited_dt_range.first;
  double lhs     = (dt_mean-m_limited_dt_range.first)/delta;
  if(lhs < 0.5)
    {
      try
	{
	  TauOverDeltaFunction f(lhs);
	  double tod_hi =  VSAMath::findMonotonicHiBracket(f, 0, lhs);
	  double tod = 
	    VSAMath::findBracketedMonotonicRoot(f, 0, tod_hi, lhs*0.00001);
	  gps_livetime_sec = double(events_found)*tod*delta;
	}
      catch(VSAMath::FunctionNotMonotonic)
	{
	  vstream << '\n'
		  << "WARN: Delta-T livetime calculation caught exception\n";
	}
    }

  // ------------------------------------------------------------------------
  // Finalize per-telescope diagnostics
  // ------------------------------------------------------------------------

  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope])
      {
	const VSMiscellaneousDBData::Scope* db_scope_data = 0;
	if((db_data)&&(iscope<db_data->scope.size())
	   &&(db_data->scope[iscope].has_scope))
	  db_scope_data=&db_data->scope[iscope];

	const std::vector<VATime>* nsb_time = 0;
	const std::vector<double>* nsb_median = 0;
	if((nsb)&&(iscope<nsb->scope.size()))
	  nsb_time = &nsb->slice_time, 
	    nsb_median = &nsb->scope[iscope].slice_median_dev;

	scope[iscope]->finalize(&cal->scope[iscope],
				run_info, db_scope_data, nsb_time, nsb_median,
				earth_position);
      }

  // ------------------------------------------------------------------------
  // Moon data
  // ------------------------------------------------------------------------

  VSTime mid_time = run_info.lo_event_time + 
    (run_info.hi_event_time-run_info.lo_event_time)/2;
  double mjd = mid_time.getMJDDbl();
  double mjd_1min = mjd + 1.0/1440.0;

  SphericalCoords moon_radec;
  Astro::moonRaDecApparent(mjd, moon_radec);

  Angle lmst = Astro::mjdToLMST(mjd, earth_position.longitudeRad());
  SphericalCoords moon_azel = moon_radec;
  Astro::raDecToAzEl(lmst, earth_position, moon_azel);

  moon_el         = moon_azel.latitudeDeg();
  moon_az         = moon_azel.longitudeDeg();
  moon_ra         = moon_radec.latitudeDeg();
  moon_dec        = moon_radec.latitudeDeg();
  moon_ra_string  = moon_radec.longitude().hmsString(1);
  moon_dec_string = moon_radec.latitude().dmsString(1);
  moon_phase      = Astro::moonPhase(mjd);
  moon_angle      = Astro::moonAngle(mjd).deg();
  moon_dphase_dt  = (Astro::moonPhase(mjd_1min) - moon_phase)/60;
}

void VBFAnalysisStage2::ArrayDiagnostics::
save(VSOctaveH5WriterStruct* s, bool slow_diagnostics) const
{
  VSArrayDiagnosticsBase::save(s, slow_diagnostics);
  unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* c = s->writeCellVector("t",nscope);
  vsassert(c);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope])
      {
	VSOctaveH5WriterStruct* ss = c->writeStruct(iscope);
	vsassert(ss);
	scope[iscope]->save(ss, slow_diagnostics);
	delete ss;
      }
  delete c;
}
