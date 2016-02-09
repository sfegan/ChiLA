//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDiagnosticsData.cpp
  Diagnostics data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/09/2007

  $Id: VSDiagnosticsData.cpp,v 3.29 2009/11/05 22:52:26 matthew Exp $

*/

#include <vsassert>

#include <VSDiagnosticsData.hpp>
#include <VSDataReductionConstants.hpp>

using namespace VERITAS;

static void writeXY(VSOctaveH5WriterStruct* io, const std::string& name,
		    const std::vector<double>& x, 
		    const std::vector<double>& y)
{
  VSOctaveH5WriterStruct* s = io->writeStruct(name);
  vsassert(s);
  s->writeVector("x",x);
  s->writeVector("y",y);
  delete s;
}

static void readXY(VSOctaveH5ReaderStruct* io, const std::string& name,
		   std::vector<double>& x, 
		   std::vector<double>& y)
{
  VSOctaveH5ReaderStruct* s = io->readStruct(name);
  vsassert(s);
  s->readVector("x",x);
  s->readVector("y",y);
  delete s;
}

// ============================================================================
//
// PARTIAL SCOPE DIAGNOSTICS BASE
//
// ============================================================================

VSPartialScopeDiagnosticsBase::
VSPartialScopeDiagnosticsBase(unsigned nchan):
  channel_in_image(nchan),
  channel_is_triggered(nchan),
  channel_is_triggered_in_largest_region(nchan),
  channel_is_triggered_isolated(nchan),
  channel_is_triggered_in_sub_threshold_event(nchan),
  channel_is_triggered_no_l2(nchan),
  channel_is_triggered_in_largest_region_no_l2(nchan),
  channel_is_triggered_isolated_no_l2(nchan),
  channel_is_triggered_in_sub_threshold_event_no_l2(nchan),
  channel_is_hit(nchan),
  channel_is_lo_gain(nchan),
  channel_is_hi_gain(nchan),
  channel_is_lo_gain_no_l2(nchan),
  channel_is_hi_gain_no_l2(nchan),
  channel_raw_tij(nchan,0.5),
  channel_image_tij(nchan,0.5),
  channel_raw_log10_sij(nchan,VSLimitedHist<double>(0.05,-1,10)),
  channel_image_log10_sij(nchan,VSLimitedHist<double>(0.05,-1,10)),
  channel_logain_log10_sij(nchan,VSLimitedHist<double>(0.05,-1,10)),
  channel_peak_hist_trig(nchan * 256),
  channel_peak_hist_notr(nchan * 256),
  channel_peak_hist_lowg(nchan * 256),
  m_nchan(nchan)
{
  // nothing to see here
}

VSPartialScopeDiagnosticsBase::~VSPartialScopeDiagnosticsBase()
{
  // nothing to see here
}

void VSPartialScopeDiagnosticsBase::clear()
{
  m_nchan                                          = 0;

  // CHANNEL PREFIX FROM PARTIAL ----------------------------------------------
  channel_in_image                                 .clear();
  channel_is_triggered                             .clear();
  channel_is_triggered_in_largest_region           .clear();
  channel_is_triggered_isolated                    .clear();
  channel_is_triggered_in_sub_threshold_event      .clear();
  channel_is_triggered_no_l2                       .clear();
  channel_is_triggered_in_largest_region_no_l2     .clear();
  channel_is_triggered_isolated_no_l2              .clear();
  channel_is_triggered_in_sub_threshold_event_no_l2.clear();
  channel_is_hit                                   .clear();
  channel_is_lo_gain                               .clear();
  channel_is_hi_gain                               .clear();
  channel_is_lo_gain_no_l2                         .clear();
  channel_is_hi_gain_no_l2                         .clear();
  channel_raw_tij                                  .clear();
  channel_image_tij                                .clear();
  channel_raw_log10_sij                            .clear();
  channel_image_log10_sij                          .clear();
  channel_logain_log10_sij                         .clear();
  channel_peak_hist_trig                           .clear();
  channel_peak_hist_notr                           .clear();
  channel_peak_hist_lowg                           .clear();
}

void VSPartialScopeDiagnosticsBase::load(VSOctaveH5ReaderStruct* s)
{
  clear();
  VSOctaveH5ReaderCellVector* cell;

  // CHANNEL PREFIX -----------------------------------------------------------
  s->readVector("channel_in_image",channel_in_image);
  s->readVector("channel_is_triggered",channel_is_triggered);
  s->readVector("channel_is_triggered_in_largest_region",
		channel_is_triggered_in_largest_region);
  s->readVector("channel_is_triggered_isolated",
		channel_is_triggered_isolated);
  s->readVector("channel_is_triggered_in_sub_threshold_event",
		channel_is_triggered_in_sub_threshold_event);
  s->readVector("channel_is_triggered_no_l2",channel_is_triggered_no_l2);
  s->readVector("channel_is_triggered_in_largest_region_no_l2",
		channel_is_triggered_in_largest_region_no_l2);
  s->readVector("channel_is_triggered_isolated_no_l2",
		channel_is_triggered_isolated_no_l2);
  s->readVector("channel_is_triggered_in_sub_threshold_event_no_l2",
		channel_is_triggered_in_sub_threshold_event_no_l2);
  s->readVector("channel_is_hit",channel_is_hit);
  s->readVector("channel_is_lo_gain",channel_is_lo_gain);
  s->readVector("channel_is_hi_gain",channel_is_hi_gain);
  s->readVector("channel_is_lo_gain_no_l2",channel_is_lo_gain_no_l2);
  s->readVector("channel_is_hi_gain_no_l2",channel_is_hi_gain_no_l2);

  if(s->isCellVector("channel_tij"))
    {
      cell = s->readCellVector("channel_tij");
      vsassert(cell);
      m_nchan = cell->dimensions();
      channel_image_tij.resize(m_nchan,0.5);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5ReaderStruct* cs = cell->readStruct(ichan);
	  channel_image_tij[ichan].load(cs);
	  delete cs;
	}
      delete cell;
    }

  if(s->isCellVector("channel_tij_raw"))
    {
      cell = s->readCellVector("channel_tij_raw");
      vsassert(cell);
      vsassert(m_nchan == cell->dimensions());
      channel_raw_tij.resize(m_nchan,0.5);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5ReaderStruct* cs = cell->readStruct(ichan);
	  channel_raw_tij[ichan].load(cs);
	  delete cs;
	}
      delete cell;
    }

  if(s->isCellVector("channel_log10_sij_raw"))
    {
      cell = s->readCellVector("channel_log10_sij_raw");
      vsassert(cell);
      vsassert(m_nchan == cell->dimensions());
      channel_raw_log10_sij.
	resize(m_nchan,VSLimitedHist<double>(0.05,-1.,10.));
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5ReaderStruct* cs = cell->readStruct(ichan);
	  channel_raw_log10_sij[ichan].load(cs);
	  delete cs;
	}
      delete cell;
    }

  if(s->isCellVector("channel_log10_sij_image"))
    {
      cell = s->readCellVector("channel_log10_sij_image");
      vsassert(cell);
      vsassert(m_nchan == cell->dimensions());
      channel_image_log10_sij.
	resize(m_nchan,VSLimitedHist<double>(0.05,-1.,10.));
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5ReaderStruct* cs = cell->readStruct(ichan);
	  channel_image_log10_sij[ichan].load(cs);
	  delete cs;
	}
      delete cell;
    }

  if(s->isCellVector("channel_log10_sij_logain"))
    {
      cell = s->readCellVector("channel_log10_sij_logain");
      vsassert(cell);
      vsassert(m_nchan == cell->dimensions());
      channel_logain_log10_sij.
	resize(m_nchan,VSLimitedHist<double>(0.05,-1.,10.));
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5ReaderStruct* cs = cell->readStruct(ichan);
	  channel_logain_log10_sij[ichan].load(cs);
	  delete cs;
	}
      delete cell;
    }

  if(s->isMatrix("channel_peak_hist_triggered"))
    {
      unsigned npeak;
      unsigned* peak_trig;
      s->readMatrix("channel_peak_hist_triggered",npeak,m_nchan,peak_trig);
      if(npeak>256)npeak=256;
      channel_peak_hist_trig.resize(m_nchan*npeak);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	for(unsigned ipeak=0;ipeak<npeak;ipeak++)
	  channel_peak_hist_trig[peak_hist_iel(ichan,ipeak)] =
	    peak_trig[ipeak*m_nchan+ichan];
      delete[] peak_trig;
    }

  if(s->isMatrix("channel_peak_hist_notrigger"))
    {
      unsigned npeak;
      unsigned* peak_notr;
      s->readMatrix("channel_peak_hist_notrigger",npeak,m_nchan,peak_notr);
      if(npeak>256)npeak=256;
      channel_peak_hist_notr.resize(m_nchan*npeak);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	for(unsigned ipeak=0;ipeak<npeak;ipeak++)
	  channel_peak_hist_notr[peak_hist_iel(ichan,ipeak)] =
	    peak_notr[ipeak*m_nchan+ichan];
      delete[] peak_notr;
    }

  if(s->isMatrix("channel_peak_hist_lowgain"))
    {
      unsigned npeak;
      unsigned* peak_lowg;
      s->readMatrix("channel_peak_hist_lowgain",npeak,m_nchan,peak_lowg);
      if(npeak>256)npeak=256;
      channel_peak_hist_lowg.resize(m_nchan*npeak);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	for(unsigned ipeak=0;ipeak<npeak;ipeak++)
	  channel_peak_hist_lowg[peak_hist_iel(ichan,ipeak)] =
	    peak_lowg[ipeak*m_nchan+ichan];
      delete[] peak_lowg;
    }
}

void VSPartialScopeDiagnosticsBase::
save(VSOctaveH5WriterStruct* s, bool slow_diagnostics) const
{
  VSOctaveH5WriterCellVector* cell;

  // CHANNEL PREFIX -----------------------------------------------------------
  s->writeVector("channel_in_image",channel_in_image);
  s->writeVector("channel_is_triggered",channel_is_triggered);

  s->writeVector("channel_is_triggered_no_l2",channel_is_triggered_no_l2);
  s->writeVector("channel_is_hit",channel_is_hit);
  s->writeVector("channel_is_lo_gain",channel_is_lo_gain);
  s->writeVector("channel_is_hi_gain",channel_is_hi_gain);
  s->writeVector("channel_is_lo_gain_no_l2",channel_is_lo_gain_no_l2);
  s->writeVector("channel_is_hi_gain_no_l2",channel_is_hi_gain_no_l2);

  if(slow_diagnostics)
    {
      s->writeVector("channel_is_triggered_in_largest_region",
		     channel_is_triggered_in_largest_region);
      s->writeVector("channel_is_triggered_isolated",
		     channel_is_triggered_isolated);
      s->writeVector("channel_is_triggered_in_sub_threshold_event",
		     channel_is_triggered_in_sub_threshold_event);
      s->writeVector("channel_is_triggered_in_largest_region_no_l2",
		     channel_is_triggered_in_largest_region_no_l2);
      s->writeVector("channel_is_triggered_isolated_no_l2",
		     channel_is_triggered_isolated_no_l2);
      s->writeVector("channel_is_triggered_in_sub_threshold_event_no_l2",
		     channel_is_triggered_in_sub_threshold_event_no_l2);

      cell = s->writeCellVector("channel_tij",m_nchan);
      vsassert(cell);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5WriterStruct* cs = cell->writeStruct(ichan);
	  channel_image_tij[ichan].save(cs);
	  delete cs;
	}
      delete cell;

      cell = s->writeCellVector("channel_tij_raw",m_nchan);
      vsassert(cell);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5WriterStruct* cs = cell->writeStruct(ichan);
	  channel_raw_tij[ichan].save(cs);
	  delete cs;
	}
      delete cell;

      cell = s->writeCellVector("channel_log10_sij_raw",m_nchan);
      vsassert(cell);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5WriterStruct* cs = cell->writeStruct(ichan);
	  channel_raw_log10_sij[ichan].save(cs);
	  delete cs;
	}
      delete cell;

      cell = s->writeCellVector("channel_log10_sij_image",m_nchan);
      vsassert(cell);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5WriterStruct* cs = cell->writeStruct(ichan);
	  channel_image_log10_sij[ichan].save(cs);
	  delete cs;
	}
      delete cell;

      cell = s->writeCellVector("channel_log10_sij_logain",m_nchan);
      vsassert(cell);
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	{
	  VSOctaveH5WriterStruct* cs = cell->writeStruct(ichan);
	  channel_logain_log10_sij[ichan].save(cs);
	  delete cs;
	}
      delete cell;

      const unsigned npeak = 256;
      unsigned* peak_trig = new unsigned[m_nchan*npeak];
      unsigned* peak_notr = new unsigned[m_nchan*npeak];
      unsigned* peak_lowg = new unsigned[m_nchan*npeak];
      for(unsigned ichan=0;ichan<m_nchan;ichan++)
	for(unsigned ipeak=0;ipeak<npeak;ipeak++)
	  {
	    peak_trig[ipeak*m_nchan+ichan] = 
	      channel_peak_hist_trig[peak_hist_iel(ichan,ipeak)];
	    peak_notr[ipeak*m_nchan+ichan] = 
	      channel_peak_hist_notr[peak_hist_iel(ichan,ipeak)];
	    peak_lowg[ipeak*m_nchan+ichan] = 
	      channel_peak_hist_lowg[peak_hist_iel(ichan,ipeak)];
	  }
      s->writeMatrix("channel_peak_hist_triggered",npeak,m_nchan,peak_trig);
      s->writeMatrix("channel_peak_hist_notrigger",npeak,m_nchan,peak_notr);
      s->writeMatrix("channel_peak_hist_lowgain",npeak,m_nchan,peak_lowg);
      delete[] peak_trig;
      delete[] peak_notr;
      delete[] peak_lowg;
    }
}

// ============================================================================
// 
// SCOPE DIAGNOSTICS BASE
// 
// ============================================================================

VSScopeDiagnosticsBase::
VSScopeDiagnosticsBase(unsigned nchan, unsigned nsample, unsigned nscope):
  // SCOPE PREFIX -------------------------------------------------------------
  scope_trigger(), 
  scope_sent_l3(), 
  scope_has_event(), 
  scope_event_processed(), 
  scope_image(), 
  scope_has_muon(),
  scope_trigger_ev_num_hist(1000), 
  scope_sent_l3_ev_num_hist(1000),
  scope_has_event_ev_num_hist(1000),
  scope_tdc_diff(nscope,VSLimitedHist<double>(TDC_RES,-1000.0,1000.0)),
  scope_l3_ticks_veto_vdaq_hist(HIST_PERIOD,-5,120),

  // L2 PREFIX ----------------------------------------------------------------
  l2_rate_hist(HIST_PERIOD,-5,120),
  l2_high_res_rate_hist(0.005,-5,120),

  // GPS PREFIX ---------------------------------------------------------------
  gps_scope_event_dt(100 /* nano-seconds */, -10000, 10000),
  gps_scope_diff_mean_ev_hist(),
  gps_scope_diff_dev_ev_hist(),
  gps_scope_diff_rms_ev_hist(),

  // CHANNEL PREFIX -----------------------------------------------------------
  channel_is_raw_max1(nchan),
  channel_is_sig_max1(nchan),
  channel_is_raw_top3(nchan),
  channel_is_sig_top3(nchan),
  channel_is_raw_max1_no_l2(nchan),
  channel_is_sig_max1_no_l2(nchan),
  channel_is_raw_top3_no_l2(nchan),
  channel_is_sig_top3_no_l2(nchan),
  channel_mean_lotrace_mtx(),
  channel_mean_hitrace_mtx(),
  channel_mean_lotrace_mtx_no_l2(),
  channel_mean_hitrace_mtx_no_l2(),
  channel_pedestal_covariance_mtx(),
  channel_cfd_threshold(nchan),

  // CAMERA PREFIX ------------------------------------------------------------
  camera_ntrigger(1),
  camera_ntrigger_largest_region(1),
  camera_nimage_minus_ntrigger_largest_region(1),
  camera_ntrigger_no_l2(1), 
  camera_ntrigger_largest_region_no_l2(1),
  camera_nimage_minus_ntrigger_largest_region_no_l2(1),
  camera_log10_charge(0.1,-1.0,10.0), 
  camera_log10_charge_no_l2(0.1,-1.0,10.0), 
  camera_log10_charge_logain(0.1,-1.0,10.0), 
  camera_log10_max1(0.02,-1.0,10.0), 
  camera_log10_max1_no_l2(0.02,-1.0,10.0), 
  camera_log10_sij_raw(0.02,-1.0,10.0), 
  camera_log10_sij_image(0.02,-1.0,10.0), 
  camera_log10_sij_logain(0.02,-1.0,10.0), 
  camera_nimage(1), 
  camera_nimage_no_l2(1),
  camera_width(0.01),
  camera_length(0.01), 
  camera_psi(5), 
  camera_xc(0.02), 
  camera_yc(0.02),
  intrinsic_width(0.01),
  intrinsic_length(0.01), 

  // SCOPE PREFIX -------------------------------------------------------------
  scope_pos_t(),
  scope_pos_az(),
  scope_pos_zn(),
  scope_pos_ra(),
  scope_pos_dec(),
  scope_pos_ra_dev(),
  scope_pos_dec_dev(),
  scope_pos_l(),
  scope_pos_b(),
  scope_pos_l_dev(),
  scope_pos_b_dev(),
  scope_moon_separation(),

  // MEDIAN PREFIX ------------------------------------------------------------
  median_l1_rate_time(), 
  median_l1_rate(), 
  median_current_time(), 
  median_current(), 
  median_pedvar_time(), 
  median_pedvar(),

  m_nchan(nchan),
  m_nsample(nsample)
{
  // nothing to see here
}

VSScopeDiagnosticsBase::~VSScopeDiagnosticsBase()
{
  delete[] channel_mean_lotrace_mtx;
  delete[] channel_mean_hitrace_mtx;
  delete[] channel_mean_lotrace_mtx_no_l2;
  delete[] channel_mean_hitrace_mtx_no_l2;
  delete[] channel_pedestal_covariance_mtx;
}

void VSScopeDiagnosticsBase::clear()
{
  m_nchan                                          = 0;
  m_nsample                                        = 0;

  // SCOPE PREFIX -------------------------------------------------------------
  scope_trigger                                    = 0;
  scope_sent_l3                                    = 0;
  scope_has_event                                  = 0;
  scope_event_processed                            = 0;
  scope_image                                      = 0;
  scope_has_muon                                   = 0;
  scope_trigger_ev_num_hist                        .clear();
  scope_sent_l3_ev_num_hist                        .clear();
  scope_has_event_ev_num_hist                      .clear();
  scope_tdc_diff                                   .clear();
  scope_l3_ticks_veto_vdaq_hist                    .clear();

  // L2 PREFIX ----------------------------------------------------------------
  l2_rate_hist                                     .clear();
  l2_high_res_rate_hist                            .clear();

  // GPS PREFIX ---------------------------------------------------------------
  gps_scope_event_dt                               .clear();
  gps_scope_diff_mean_ev_hist                      .clear();
  gps_scope_diff_dev_ev_hist                       .clear();
  gps_scope_diff_rms_ev_hist                       .clear();

  // CHANNEL PREFIX -----------------------------------------------------------
  channel_is_raw_max1                              .clear();
  channel_is_sig_max1                              .clear();
  channel_is_raw_max1_no_l2                        .clear();
  channel_is_sig_max1_no_l2                        .clear();
  channel_is_raw_top3                              .clear();
  channel_is_sig_top3                              .clear();
  channel_is_raw_top3_no_l2                        .clear();
  channel_is_sig_top3_no_l2                        .clear();
  delete[] channel_mean_hitrace_mtx;
  channel_mean_hitrace_mtx                         = 0;
  delete[] channel_mean_lotrace_mtx;
  channel_mean_lotrace_mtx                         = 0;
  delete[] channel_mean_hitrace_mtx_no_l2;
  channel_mean_hitrace_mtx_no_l2                   = 0;
  delete[] channel_mean_lotrace_mtx_no_l2;
  channel_mean_lotrace_mtx_no_l2                   = 0;
  delete[] channel_pedestal_covariance_mtx;
  channel_pedestal_covariance_mtx                  = 0;
  channel_cfd_threshold                            .clear();

  // CAMERA PREFIX ------------------------------------------------------------
  camera_ntrigger                                  .clear();
  camera_ntrigger_largest_region                   .clear();
  camera_nimage_minus_ntrigger_largest_region      .clear();
  camera_ntrigger_no_l2                            .clear();
  camera_ntrigger_largest_region_no_l2             .clear();
  camera_nimage_minus_ntrigger_largest_region_no_l2.clear();
  camera_log10_charge                              .clear();
  camera_log10_charge_no_l2                        .clear();
  camera_log10_charge_logain                       .clear();
  camera_log10_max1                                .clear();
  camera_log10_max1_no_l2                          .clear();
  camera_log10_sij_raw                             .clear();
  camera_log10_sij_image                           .clear();
  camera_log10_sij_logain                          .clear();
  camera_nimage                                    .clear();
  camera_nimage                                    .clear();
  camera_xc                                        .clear();
  camera_yc                                        .clear();
  camera_length                                    .clear();
  camera_width                                     .clear();
  camera_psi                                       .clear();
  intrinsic_length                                 .clear();
  intrinsic_width                                  .clear();

  // SCOPE PREFIX -------------------------------------------------------------
  scope_pos_t                                      .clear();
  scope_pos_az                                     .clear();
  scope_pos_zn                                     .clear();
  scope_pos_ra                                     .clear();
  scope_pos_dec                                    .clear();
  scope_pos_ra_dev                                 .clear();
  scope_pos_dec_dev                                .clear();
  scope_pos_l                                      .clear();
  scope_pos_b                                      .clear();
  scope_pos_l_dev                                  .clear();
  scope_pos_b_dev                                  .clear();
  scope_moon_separation                            .clear();

  // MEDIAN PREFIX ------------------------------------------------------------
  median_l1_rate_time                              .clear();
  median_current_time                              .clear();
  median_current_time                              .clear();
}

void VSScopeDiagnosticsBase::load(VSOctaveH5ReaderStruct* s)
{
  clear();
  unsigned nscope;
  VSOctaveH5ReaderCellVector* cell;

  // SCOPE PREFIX -------------------------------------------------------------
  s->readScalar("scope_trigger",scope_trigger);
  s->readScalar("scope_sent_l3",scope_sent_l3);
  s->readScalar("scope_has_event",scope_has_event);
  s->readScalar("scope_event_processed",scope_event_processed);
  s->readScalar("scope_image",scope_image);
  s->readScalar("scope_has_muon",scope_has_muon);
  scope_trigger_ev_num_hist.load(s->readStruct("scope_trigger_ev_num_hist"));
  scope_sent_l3_ev_num_hist.load(s->readStruct("scope_sent_l3_ev_num_hist"));
  scope_has_event_ev_num_hist.
    load(s->readStruct("scope_has_event_ev_num_hist"));

  cell = s->readCellVector("scope_tdc_diff");
  vsassert(cell);
  nscope = cell->dimensions();
  scope_tdc_diff.resize(nscope,VSLimitedHist<double>(0,0,0));
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!cell->isEmpty(iscope))
      {
	VSOctaveH5ReaderStruct* cs = cell->readStruct(iscope);
	scope_tdc_diff[iscope].load(cs);
	delete cs;
      }
  delete cell;
  scope_l3_ticks_veto_vdaq_hist.
    load(s->readStruct("scope_l3_ticks_veto_vdaq_hist"));

  // L2 PREFIX ----------------------------------------------------------------
  l2_rate_hist.load(s->readStruct("l2_rate_hist"));
  l2_high_res_rate_hist.load(s->readStruct("l2_high_res_rate_hist"));

  // GPS PREFIX ---------------------------------------------------------------
  if(s->isValid("gps_scope_event_dt"))
    gps_scope_event_dt.load(s->readStruct("gps_scope_event_dt"));
  else
    gps_scope_event_dt.load(s->readStruct("gps_scope_l3_dt"));

   if(s->isValid("gps_scope_diff_mean_ev_hist"))
    {
      s->readVector("gps_scope_diff_mean_ev_hist",gps_scope_diff_mean_ev_hist);
      s->readVector("gps_scope_diff_dev_ev_hist",gps_scope_diff_dev_ev_hist);
      s->readVector("gps_scope_diff_rms_ev_hist",gps_scope_diff_rms_ev_hist);
    }

  // CHANNEL PREFIX -----------------------------------------------------------
  s->readVector("channel_is_raw_max1",channel_is_raw_max1);
  s->readVector("channel_is_sig_max1",channel_is_sig_max1);
  s->readVector("channel_is_raw_max1_no_l2",channel_is_raw_max1_no_l2);
  s->readVector("channel_is_sig_max1_no_l2",channel_is_sig_max1_no_l2);
  s->readVector("channel_is_raw_top3",channel_is_raw_top3);
  s->readVector("channel_is_sig_top3",channel_is_sig_top3);
  s->readVector("channel_is_raw_top3_no_l2",channel_is_raw_top3_no_l2);
  s->readVector("channel_is_sig_top3_no_l2",channel_is_sig_top3_no_l2);
  s->readMatrix("channel_mean_higain_sample",
		m_nsample, m_nchan, channel_mean_hitrace_mtx);
  s->readMatrix("channel_mean_logain_sample",
		m_nsample, m_nchan, channel_mean_lotrace_mtx);
  s->readMatrix("channel_mean_higain_sample_no_l2",
		m_nsample, m_nchan, channel_mean_hitrace_mtx_no_l2);
  s->readMatrix("channel_mean_logain_sample_no_l2",
		m_nsample, m_nchan, channel_mean_lotrace_mtx_no_l2);
  s->readMatrix("channel_pedestal_covariance",
		m_nsample, m_nchan, channel_pedestal_covariance_mtx);
  s->readVector("channel_cfd_threshold",channel_cfd_threshold);

  // CAMERA PREFIX ------------------------------------------------------------
  camera_ntrigger.load(s->readStruct("camera_ntrigger"));
  camera_ntrigger_largest_region.
    load(s->readStruct("camera_ntrigger_largest_region"));
  camera_nimage_minus_ntrigger_largest_region.
    load(s->readStruct("camera_nimage_minus_ntrigger_largest_region"));
  camera_ntrigger_no_l2.load(s->readStruct("camera_ntrigger_no_l2"));
  camera_ntrigger_largest_region_no_l2.
    load(s->readStruct("camera_ntrigger_largest_region_no_l2"));
  camera_nimage_minus_ntrigger_largest_region_no_l2.
    load(s->readStruct("camera_nimage_minus_ntrigger_largest_region_no_l2"));
  camera_log10_charge.load(s->readStruct("camera_log10_charge"));
  camera_log10_charge_no_l2.load(s->readStruct("camera_log10_charge_no_l2"));
  camera_log10_charge_logain.
    load(s->readStruct("camera_log10_charge_logain"));
  camera_log10_max1.load(s->readStruct("camera_log10_max1"));
  camera_log10_max1_no_l2.load(s->readStruct("camera_log10_max1_no_l2"));
  camera_log10_sij_raw.load(s->readStruct("camera_log10_sij_raw"));
  camera_log10_sij_image.load(s->readStruct("camera_log10_sij_image"));
  camera_log10_sij_logain.load(s->readStruct("camera_log10_sij_logain"));
  camera_nimage.load(s->readStruct("camera_nimage"));
  camera_nimage_no_l2.load(s->readStruct("camera_nimage_no_l2"));
  camera_xc.load(s->readStruct("camera_xc"));
  camera_yc.load(s->readStruct("camera_yc"));
  camera_length.load(s->readStruct("camera_length"));
  camera_width.load(s->readStruct("camera_width"));
  camera_psi.load(s->readStruct("camera_psi"));
  intrinsic_length.load(s->readStruct("intrinsic_length"));
  intrinsic_width.load(s->readStruct("intrinsic_width"));

  // SCOPE PREFIX -------------------------------------------------------------
  VSOctaveH5ReaderStruct* smp = s->readStruct("scope_mean_position");
  vsassert(smp);
  smp->readVector("t",scope_pos_t);
  smp->readVector("az",scope_pos_az);
  smp->readVector("zn",scope_pos_zn);  
  smp->readVector("ra",scope_pos_ra);
  smp->readVector("dec",scope_pos_dec);  
  smp->readVector("ra_dev",scope_pos_ra_dev);
  smp->readVector("dec_dev",scope_pos_dec_dev);  
  smp->readVector("l",scope_pos_l);
  smp->readVector("b",scope_pos_b);  
  smp->readVector("l_dev",scope_pos_l_dev);
  smp->readVector("b_dev",scope_pos_b_dev);  
  smp->readVector("moon_separation",scope_moon_separation);
  delete smp;

  // MEDIAN PREFIX ------------------------------------------------------------
  if(s->isValid("median_l1_rate"))
    readXY(s,"median_l1_rate",median_l1_rate_time,median_l1_rate);
  if(s->isValid("median_current"))
    readXY(s,"median_current",median_current_time,median_current);
  if(s->isValid("median_pedvar"))
    readXY(s,"median_pedvar",median_pedvar_time,median_pedvar);  
}

void VSScopeDiagnosticsBase::
save(VSOctaveH5WriterStruct* s, bool slow_diagnostics) const
{
  VSOctaveH5WriterCellVector* cell;

  // SCOPE PREFIX -------------------------------------------------------------
  s->writeScalar("scope_trigger",scope_trigger);
  s->writeScalar("scope_sent_l3",scope_sent_l3);
  s->writeScalar("scope_has_event",scope_has_event);
  s->writeScalar("scope_event_processed",scope_event_processed);
  s->writeScalar("scope_image",scope_image);
  s->writeScalar("scope_has_muon",scope_has_muon);
  scope_trigger_ev_num_hist.save(s->writeStruct("scope_trigger_ev_num_hist"));
  scope_sent_l3_ev_num_hist.save(s->writeStruct("scope_sent_l3_ev_num_hist"));
  scope_has_event_ev_num_hist.
    save(s->writeStruct("scope_has_event_ev_num_hist"));
  cell = s->writeCellVector("scope_tdc_diff",scope_tdc_diff.size());
  vsassert(cell);
  for(unsigned iscope=0;iscope<scope_tdc_diff.size();iscope++)
    if(!scope_tdc_diff[iscope].empty())
      {
	VSOctaveH5WriterStruct* cs = cell->writeStruct(iscope);
	scope_tdc_diff[iscope].save(cs);
	delete cs;
      }
  delete cell;
  scope_l3_ticks_veto_vdaq_hist.
    save(s->writeStruct("scope_l3_ticks_veto_vdaq_hist"));

  // L2 PREFIX ----------------------------------------------------------------
  l2_rate_hist.save(s->writeStruct("l2_rate_hist"));
  l2_high_res_rate_hist.save(s->writeStruct("l2_high_res_rate_hist"));

  // GPS PREFIX ---------------------------------------------------------------
  gps_scope_event_dt.save(s->writeStruct("gps_scope_event_dt"));
  s->writeVector("gps_scope_diff_mean_ev_hist",gps_scope_diff_mean_ev_hist);
  s->writeVector("gps_scope_diff_dev_ev_hist",gps_scope_diff_dev_ev_hist);
  s->writeVector("gps_scope_diff_rms_ev_hist",gps_scope_diff_rms_ev_hist);

  // CHANNEL PREFIX -----------------------------------------------------------
  s->writeVector("channel_is_raw_max1",channel_is_raw_max1);
  s->writeVector("channel_is_sig_max1",channel_is_sig_max1);
  s->writeVector("channel_is_raw_max1_no_l2",channel_is_raw_max1_no_l2);
  s->writeVector("channel_is_sig_max1_no_l2",channel_is_sig_max1_no_l2);
  s->writeVector("channel_is_raw_top3",channel_is_raw_top3);
  s->writeVector("channel_is_sig_top3",channel_is_sig_top3);
  s->writeVector("channel_is_raw_top3_no_l2",channel_is_raw_top3_no_l2);
  s->writeVector("channel_is_sig_top3_no_l2",channel_is_sig_top3_no_l2);

  if(slow_diagnostics)
    {
      s->writeMatrix("channel_mean_higain_sample",
		     m_nsample, m_nchan, channel_mean_hitrace_mtx);
      s->writeMatrix("channel_mean_logain_sample",
		     m_nsample, m_nchan, channel_mean_lotrace_mtx);
      s->writeMatrix("channel_mean_higain_sample_no_l2",
		     m_nsample, m_nchan, channel_mean_hitrace_mtx_no_l2);
      s->writeMatrix("channel_mean_logain_sample_no_l2",
		     m_nsample, m_nchan, channel_mean_lotrace_mtx_no_l2);
      s->writeMatrix("channel_pedestal_covariance",
		     m_nsample, m_nchan, channel_pedestal_covariance_mtx);
      s->writeVector("channel_cfd_threshold",channel_cfd_threshold);
    }

  // CAMERA PREFIX ------------------------------------------------------------
  camera_ntrigger.save(s->writeStruct("camera_ntrigger"));
  camera_ntrigger_largest_region.
    save(s->writeStruct("camera_ntrigger_largest_region"));
  camera_nimage_minus_ntrigger_largest_region.
    save(s->writeStruct("camera_nimage_minus_ntrigger_largest_region"));
  camera_ntrigger_no_l2.save(s->writeStruct("camera_ntrigger_no_l2"));
  camera_ntrigger_largest_region_no_l2.
    save(s->writeStruct("camera_ntrigger_largest_region_no_l2"));
  camera_nimage_minus_ntrigger_largest_region_no_l2.
    save(s->writeStruct("camera_nimage_minus_ntrigger_largest_region_no_l2"));
  camera_log10_charge.save(s->writeStruct("camera_log10_charge"));
  camera_log10_charge_no_l2.save(s->writeStruct("camera_log10_charge_no_l2"));
  camera_log10_charge_logain.
    save(s->writeStruct("camera_log10_charge_logain"));
  camera_log10_max1.save(s->writeStruct("camera_log10_max1"));
  camera_log10_max1_no_l2.save(s->writeStruct("camera_log10_max1_no_l2"));
  camera_log10_sij_raw.save(s->writeStruct("camera_log10_sij_raw"));
  camera_log10_sij_image.save(s->writeStruct("camera_log10_sij_image"));
  camera_log10_sij_logain.save(s->writeStruct("camera_log10_sij_logain"));
  camera_nimage.save(s->writeStruct("camera_nimage"));
  camera_nimage_no_l2.save(s->writeStruct("camera_nimage_no_l2"));
  camera_xc.save(s->writeStruct("camera_xc"));
  camera_yc.save(s->writeStruct("camera_yc"));
  camera_length.save(s->writeStruct("camera_length"));
  camera_width.save(s->writeStruct("camera_width"));
  camera_psi.save(s->writeStruct("camera_psi"));
  intrinsic_length.save(s->writeStruct("intrinsic_length"));
  intrinsic_width.save(s->writeStruct("intrinsic_width"));

  // SCOPE PREFIX -------------------------------------------------------------
  VSOctaveH5WriterStruct* smp = s->writeStruct("scope_mean_position");
  smp->writeVector("t",scope_pos_t);
  smp->writeVector("az",scope_pos_az);
  smp->writeVector("zn",scope_pos_zn);  
  smp->writeVector("ra",scope_pos_ra);
  smp->writeVector("dec",scope_pos_dec);  
  smp->writeVector("ra_dev",scope_pos_ra_dev);
  smp->writeVector("dec_dev",scope_pos_dec_dev);  
  smp->writeVector("l",scope_pos_l);
  smp->writeVector("b",scope_pos_b);  
  smp->writeVector("l_dev",scope_pos_l_dev);
  smp->writeVector("b_dev",scope_pos_b_dev);  
  smp->writeVector("moon_separation",scope_moon_separation);
  delete smp;

  // MEDIAN PREFIX ------------------------------------------------------------
  if(!median_l1_rate_time.empty())
    writeXY(s,"median_l1_rate",median_l1_rate_time,median_l1_rate);
  if(!median_current_time.empty())
    writeXY(s,"median_current",median_current_time,median_current);
  if(!median_pedvar_time.empty())
    writeXY(s,"median_pedvar",median_pedvar_time,median_pedvar);  
}

// ============================================================================
//
// SCOPE DIAGNOSTICS DATA
//
// ============================================================================

VSScopeDiagnosticsData::
VSScopeDiagnosticsData(unsigned nchan, unsigned nsample, unsigned nscope):
  VSPartialScopeDiagnosticsBase(nchan), 
  VSScopeDiagnosticsBase(nchan, nsample, nscope)
{
  // nothing to see here
}

VSScopeDiagnosticsData::~VSScopeDiagnosticsData()
{
  // nothing to see here
}

void VSScopeDiagnosticsData::clear()
{
  VSPartialScopeDiagnosticsBase::clear();
  VSScopeDiagnosticsBase::clear();
}

void VSScopeDiagnosticsData::load(VSOctaveH5ReaderStruct* s)
{
  VSPartialScopeDiagnosticsBase::load(s);
  VSScopeDiagnosticsBase::load(s);
}

void VSScopeDiagnosticsData::
save(VSOctaveH5WriterStruct* s, bool slow_diagnostics) const
{
  VSPartialScopeDiagnosticsBase::save(s, slow_diagnostics);
  VSScopeDiagnosticsBase::save(s, slow_diagnostics);
}

// ============================================================================
//
// ARRAY DIAGNOSTICS DATA
//
// ============================================================================

VSArrayDiagnosticsBase::VSArrayDiagnosticsBase(unsigned nscope):
  packets_found(0), 
  events_found(0),
  events_processed(), 
  events_reconstructed(), 
  events_written(),
  events_ped(),
  events_l2(),
  events_failed_software_trigger(),
  gps_ticks_elapsed(), 
  gps_elaptime_sec(),
  gps_livetime_sec(),
  gps_l3_event_dt(100 /* nano-seconds */, -10000, 10000),
  gps_l3_diff_mean_ev_hist(),
  gps_l3_diff_dev_ev_hist(),
  gps_l3_diff_rms_ev_hist(),
  has_l3(), 
  has_l3_ev_num_hist(1000),
  l3_dt(0.0005,0.0,10.0), 
  l3_log10_dt(0.1,-10.0,10.0), 
  l3_rate_hist(HIST_PERIOD,-5,120), 
  l3_events_missing(),
  l3_ticks_elapsed(), 
  l3_ticks_veto_vdaq(), 
  l3_ticks_veto_lev3(), 
  l3_ticks_veto_both(),
  l3_ticks_elapsed_hist(HIST_PERIOD,-5,120), 
  l3_ticks_veto_vdaq_hist(HIST_PERIOD,-5,120),
  l3_ticks_veto_lev3_hist(HIST_PERIOD,-5,120),
  l3_ticks_veto_both_hist(HIST_PERIOD,-5,120),
  l3_elaptime_sec(),
  l3_livetime_sec(),
  scope_nimage(1), 
  scope_nquality(1), 
  scope_ntrigger(1), 
  scope_nhas_event(1), 
  scope_nsent_l3(1),
  moon_el(),
  moon_az(),
  moon_ra(),
  moon_dec(),
  moon_ra_string(),
  moon_dec_string(),
  moon_phase(),
  moon_angle(),
  moon_dphase_dt(),
  fir(),
  m_nmask(1<<nscope)
{ 

}

VSArrayDiagnosticsBase::~VSArrayDiagnosticsBase()
{

}

void VSArrayDiagnosticsBase::clear()
{
  // PACKETS AND EVENTS -------------------------------------------------------
  packets_found                                    = 0;
  events_found                                     = 0;
  events_processed                                 = 0;
  events_reconstructed                             = 0;
  events_written                                   = 0;
  events_ped                                       = 0;
  events_l2                                        = 0;
  events_failed_software_trigger                   = 0;

  // GPS PREFIX ---------------------------------------------------------------
  gps_ticks_elapsed                                = 0;
  gps_elaptime_sec                                 = 0;
  gps_livetime_sec                                 = 0;
  gps_l3_event_dt                                  .clear();
  gps_l3_diff_mean_ev_hist                         .clear();
  gps_l3_diff_dev_ev_hist                          .clear();
  gps_l3_diff_rms_ev_hist                          .clear();

  // L3 PREFIX ----------------------------------------------------------------
  has_l3                                           = 0;
  has_l3_ev_num_hist                               .clear();
  l3_dt                                            .clear();
  l3_log10_dt                                      .clear();
  l3_rate_hist                                     .clear();
  l3_events_missing                                .clear();
  l3_ticks_elapsed                                 = 0;
  l3_ticks_veto_vdaq                               = 0;
  l3_ticks_veto_lev3                               = 0;
  l3_ticks_veto_both                               = 0;
  l3_ticks_elapsed_hist                            .clear();
  l3_ticks_veto_vdaq_hist                          .clear();
  l3_ticks_veto_lev3_hist                          .clear();
  l3_ticks_veto_both_hist                          .clear();
  l3_elaptime_sec                                  = 0;
  l3_livetime_sec                                  = 0;

  // SCOPE PREFIX -------------------------------------------------------------
  scope_nimage                                     .clear();
  scope_nquality                                   .clear();
  scope_ntrigger                                   .clear();
  scope_nhas_event                                 .clear();
  scope_nsent_l3                                   .clear();

  // MOON PREFIX --------------------------------------------------------------
  moon_el                                          = 0;
  moon_az                                          = 0;
  moon_ra                                          = 0;
  moon_dec                                         = 0;
  moon_ra_string                                   .clear();
  moon_dec_string                                  .clear();
  moon_phase                                       = 0;
  moon_angle                                       = 0;
  moon_dphase_dt                                   = 0;

  // FIR PREFIX ---------------------------------------------------------------
  fir                                              .clear();
}

void VSArrayDiagnosticsBase::load(VSOctaveH5ReaderStruct* s)
{
  clear();

  // PACKETS AND EVENTS -------------------------------------------------------
  s->readScalar("packets_found",packets_found);
  s->readScalar("events_found",events_found);
  s->readScalar("events_processed",events_processed);
  s->readScalar("events_reconstructed",events_reconstructed);
  s->readScalar("events_written",events_written);
  s->readScalar("events_ped",events_ped);
  s->readScalar("events_l2",events_l2);
  s->readScalar("events_failed_software_trigger",
		events_failed_software_trigger);

  unsigned nmask = 0;
  unsigned mmask = 0;
  vsassert(nmask==mmask);
  m_nmask = nmask;

  // GPS PREFIX ---------------------------------------------------------------
  s->readScalar("gps_ticks_elapsed",gps_ticks_elapsed);
  s->readScalar("gps_elaptime_sec",gps_elaptime_sec);
  s->readScalar("gps_livetime_sec",gps_livetime_sec);
  if(s->isValid("gps_l3_event_dt"))
    gps_l3_event_dt.load(s->readStruct("gps_l3_event_dt"));
  if(s->isValid("gps_l3_diff_mean_ev_hist"))
    {
      s->readVector("gps_l3_diff_mean_ev_hist",gps_l3_diff_mean_ev_hist);
      s->readVector("gps_l3_diff_dev_ev_hist",gps_l3_diff_dev_ev_hist);
      s->readVector("gps_l3_diff_rms_ev_hist",gps_l3_diff_rms_ev_hist);
    }

  // L3 PREFIX ----------------------------------------------------------------
  s->readScalar("has_l3",has_l3);
  has_l3_ev_num_hist.load(s->readStruct("has_l3_ev_num_hist"));
  l3_dt.load(s->readStruct("l3_dt"));
  l3_log10_dt.load(s->readStruct("l3_log10_dt"));
  l3_rate_hist.load(s->readStruct("l3_rate_hist"));
  s->readVector("l3_events_missing",l3_events_missing);
  s->readScalar("l3_ticks_elapsed",l3_ticks_elapsed);
  s->readScalar("l3_ticks_veto_vdaq", l3_ticks_veto_vdaq);
  s->readScalar("l3_ticks_veto_lev3", l3_ticks_veto_lev3);
  s->readScalar("l3_ticks_veto_both",l3_ticks_veto_both);
  l3_ticks_elapsed_hist.load(s->readStruct("l3_ticks_elapsed_hist"));
  l3_ticks_veto_vdaq_hist.load(s->readStruct("l3_ticks_veto_vdaq_hist"));
  l3_ticks_veto_lev3_hist.load(s->readStruct("l3_ticks_veto_lev3_hist"));
  l3_ticks_veto_both_hist.load(s->readStruct("l3_ticks_veto_both_hist"));
  s->readScalar("l3_elaptime_sec",l3_elaptime_sec);
  s->readScalar("l3_livetime_sec",l3_livetime_sec);

  // SCOPE PREFIX -------------------------------------------------------------
  scope_nimage.load(s->readStruct("scope_nimage"));
  scope_nquality.load(s->readStruct("scope_nquality"));
  scope_ntrigger.load(s->readStruct("scope_ntrigger"));
  scope_nsent_l3.load(s->readStruct("scope_nsent_l3"));
  scope_nhas_event.load(s->readStruct("scope_nhas_event"));

  // MOON PREFIX --------------------------------------------------------------
  s->readScalar("moon_el",moon_el);
  s->readScalar("moon_az",moon_az);
  s->readScalar("moon_ra",moon_ra);
  s->readScalar("moon_dec",moon_dec);
  s->readString("moon_ra_string",moon_ra_string);
  s->readString("moon_dec_string",moon_dec_string);
  s->readScalar("moon_phase",moon_phase);
  s->readScalar("moon_angle",moon_angle);
  s->readScalar("moon_dphase_dt",moon_dphase_dt);

  // FIR PREFIX ---------------------------------------------------------------
  if(s->isValid("fir"))
    {
      VSOctaveH5ReaderCellVector* cv = s->readCellVector("fir");
      unsigned nfir = cv->dimensions();
      fir.resize(nfir);
      for(unsigned ifir=0;ifir<nfir;ifir++)
	cv->readCompositeVector(ifir,fir[ifir]);
      delete cv;
    }
}

void VSArrayDiagnosticsBase::
save(VSOctaveH5WriterStruct* s, bool slow_diagnostics) const
{
  // PACKETS AND EVENTS -------------------------------------------------------
  s->writeScalar("packets_found",packets_found);
  s->writeScalar("events_found",events_found);
  s->writeScalar("events_processed",events_processed);
  s->writeScalar("events_reconstructed",events_reconstructed);
  s->writeScalar("events_written",events_written);
  s->writeScalar("events_ped",events_ped);
  s->writeScalar("events_l2",events_l2);
  s->writeScalar("events_failed_software_trigger",
		 events_failed_software_trigger);

  // GPS PREFIX ---------------------------------------------------------------
  s->writeScalar("gps_ticks_elapsed",gps_ticks_elapsed);
  s->writeScalar("gps_elaptime_sec",gps_elaptime_sec);
  s->writeScalar("gps_livetime_sec",gps_livetime_sec);
  gps_l3_event_dt.save(s->writeStruct("gps_l3_event_dt"));
  s->writeVector("gps_l3_diff_mean_ev_hist",gps_l3_diff_mean_ev_hist);
  s->writeVector("gps_l3_diff_dev_ev_hist",gps_l3_diff_dev_ev_hist);
  s->writeVector("gps_l3_diff_rms_ev_hist",gps_l3_diff_rms_ev_hist);
  
  // L3 PREFIX ----------------------------------------------------------------
  s->writeScalar("has_l3",has_l3);
  has_l3_ev_num_hist.save(s->writeStruct("has_l3_ev_num_hist"));
  l3_dt.save(s->writeStruct("l3_dt"));
  l3_log10_dt.save(s->writeStruct("l3_log10_dt"));
  l3_rate_hist.save(s->writeStruct("l3_rate_hist"));
  s->writeVector("l3_events_missing",l3_events_missing);
  s->writeScalar("l3_ticks_elapsed",l3_ticks_elapsed);
  s->writeScalar("l3_ticks_veto_vdaq", l3_ticks_veto_vdaq);
  s->writeScalar("l3_ticks_veto_lev3", l3_ticks_veto_lev3);
  s->writeScalar("l3_ticks_veto_both",l3_ticks_veto_both);
  l3_ticks_elapsed_hist.save(s->writeStruct("l3_ticks_elapsed_hist"));
  l3_ticks_veto_vdaq_hist.save(s->writeStruct("l3_ticks_veto_vdaq_hist"));
  l3_ticks_veto_lev3_hist.save(s->writeStruct("l3_ticks_veto_lev3_hist"));
  l3_ticks_veto_both_hist.save(s->writeStruct("l3_ticks_veto_both_hist"));
  s->writeScalar("l3_elaptime_sec",l3_elaptime_sec);
  s->writeScalar("l3_livetime_sec",l3_livetime_sec);

  // SCOPE PREFIX -------------------------------------------------------------
  scope_nimage.save(s->writeStruct("scope_nimage"));
  scope_nquality.save(s->writeStruct("scope_nquality"));
  scope_ntrigger.save(s->writeStruct("scope_ntrigger"));
  scope_nsent_l3.save(s->writeStruct("scope_nsent_l3"));
  scope_nhas_event.save(s->writeStruct("scope_nhas_event"));

  // MOON PREFIX --------------------------------------------------------------
  s->writeScalar("moon_el",moon_el);
  s->writeScalar("moon_az",moon_az);
  s->writeScalar("moon_ra",moon_ra);
  s->writeScalar("moon_dec",moon_dec);
  s->writeString("moon_ra_string",moon_ra_string);
  s->writeString("moon_dec_string",moon_dec_string);
  s->writeScalar("moon_phase",moon_phase);
  s->writeScalar("moon_angle",moon_angle);
  s->writeScalar("moon_dphase_dt",moon_dphase_dt);

  // FIR PREFIX ---------------------------------------------------------------
  if(!fir.empty())
    {
      unsigned nfir = fir.size();
      VSOctaveH5WriterCellVector* cv = s->writeCellVector("fir",nfir);
      for(unsigned ifir=0;ifir<nfir;ifir++)
	cv->writeCompositeVector(ifir,fir[ifir]);
      delete cv;
    }
}

VSArrayDiagnosticsData::VSArrayDiagnosticsData(unsigned nscope):
  VSArrayDiagnosticsBase(nscope), scope()
{
  // nothing to see here
}

VSArrayDiagnosticsData::
VSArrayDiagnosticsData(const std::vector<bool>& config_mask,
		       const std::vector<unsigned>& nchan, 
		       const std::vector<unsigned>& nsample):
  VSArrayDiagnosticsBase(nchan.size()), scope()
{ 
  unsigned nscope = nchan.size();
  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned ns = 0;
      if(iscope < nsample.size())ns=nsample[iscope];
      if((config_mask[iscope])||(nchan[iscope]))
	scope[iscope]=
	  new VSScopeDiagnosticsData(nchan[iscope], ns, nscope);
    }
}

VSArrayDiagnosticsData::~VSArrayDiagnosticsData()
{
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)delete scope[iscope];
}

void VSArrayDiagnosticsData::clear()
{
  VSArrayDiagnosticsBase::clear();
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    delete scope[iscope];
  scope.clear();
}

void VSArrayDiagnosticsData::load(VSOctaveH5ReaderStruct* s)
{
  clear();
  VSArrayDiagnosticsBase::load(s);
  VSOctaveH5ReaderCellVector* c = s->readCellVector("t");
  vsassert(c);
  unsigned nscope = c->dimensions();
  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!c->isEmpty(iscope))
      {
	VSOctaveH5ReaderStruct* ss = c->readStruct(iscope);
	vsassert(ss);
	scope[iscope] = new VSScopeDiagnosticsData;
	scope[iscope]->load(ss);
	delete ss;
      }
  delete c;
}

void VSArrayDiagnosticsData::
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
