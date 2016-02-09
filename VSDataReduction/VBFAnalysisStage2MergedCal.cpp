//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFAnalysisStage2MergedCal.cpp

  Stage 2 analysis (merged calibration data)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 3.12 $
  \date       05/10/2006

  $Id: VBFAnalysisStage2MergedCal.cpp,v 3.12 2009/11/05 22:53:31 matthew Exp $

*/

#include <VSAnalysisStage2.hpp>

using namespace VERITAS;
using namespace SEphem;

// ----------------------------------------------------------------------------
//  __  __                     _  ___      _
// |  \/  |___ _ _ __ _ ___ __| |/ __|__ _| |
// | |\/| / -_) '_/ _` / -_) _` | (__/ _` | |
// |_|  |_\___|_| \__, \___\__,_|\___\__,_|_|
//                |___/
// ----------------------------------------------------------------------------

VBFAnalysisStage2::ArrayMergedCal::
ArrayMergedCal(const VSRunInfoData& run_info,
	       const VSTimeDepSupData* sup, const VSTimeDepPedData* ped,
	       const std::vector<unsigned> hi_sample_zero,
	       const std::vector<unsigned> lo_sample_zero,
	       CleaningScale cs, unsigned window_size,
	       const VSChannelMap* channel_map,
	       const VSAReconstruction::ArrayInfo& array_info,
	       const VSLaserData* laser, bool permissive_laser, 
	       bool unity_gain, const VSHiLoData* hilo,
	       const std::vector<double>& hi_lo_gain_ratio,
	       unsigned ped_event_count_min,
	       const std::vector<bool>& sup_scope,
	       const std::vector<double>& gain,
	       const bool auto_suppress_l2_chan,
	       const ChanList& sup_chan, const ChanList& no_pmt_chan,
	       const VSTimeDepSupData* pad_sup, 
	       const VSTimeDepPedData* pad_ped,
	       const VSMiscellaneousDBData* db_data,
	       const std::vector<LinearCoefficients>& dev_scaling):
  scope(), m_sup(sup), m_ped(ped), m_cs(cs), m_window_size(window_size),
  m_pad_sup(pad_sup), m_pad_ped(pad_ped), 
  m_sup_islice(), m_ped_islice(),
  m_ped_event_count_min(ped_event_count_min)
{ 
  static const double nochan = nan("");

  unsigned nscope=0;
  for(unsigned iscope=0;iscope<run_info.nchan.size();iscope++)
    if(run_info.nchan[iscope])
      nscope = iscope+1;

  const unsigned nslice = m_sup->nslice();

  // Merge calibration info ---------------------------------------------------
  
  vsassert(array_info.scopes.size() >= nscope);

  double nscaled_dev = 0;
  double scaled_dev_sum = 0;
  double scaled_dev_sumsq = 0;

  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const bool sup_iscope = 
	(iscope<sup_scope.size()) ? sup_scope[iscope] : false;
      const unsigned nchan = sup_iscope ? 0 : run_info.nchan[iscope];

      scope[iscope].nchan              = nchan;
      scope[iscope].channel            .resize(nchan,nchan);
      scope[iscope].suppress           = sup_iscope;
      scope[iscope].gain               = 1.0;
      scope[iscope].hi_sample_zero     = 0;
      scope[iscope].lo_sample_zero     = 0;
      scope[iscope].position_E_m       = array_info.scopes[iscope].ri.x();
      scope[iscope].position_N_m       = array_info.scopes[iscope].ri.y();
      scope[iscope].position_U_m       = array_info.scopes[iscope].ri.z();
      scope[iscope].median_median_dev  = 0;
      scope[iscope].scaled_dev         = 0;

      if(nchan == 0)continue;

      vsassert(iscope < hi_sample_zero.size());
      vsassert(iscope < lo_sample_zero.size());
      vsassert(iscope < m_ped->nscope());
      vsassert(iscope < m_sup->nscope());
      vsassert((!laser)||(permissive_laser)||(iscope < laser->scope.size()));
      vsassert((!m_pad_sup)||(iscope < m_pad_sup->nscope()));
      vsassert((!m_pad_ped)||(iscope < m_pad_ped->nscope()));

      scope[iscope].hi_sample_zero     = hi_sample_zero[iscope];
      scope[iscope].lo_sample_zero     = lo_sample_zero[iscope];
      if(iscope<gain.size())
	scope[iscope].gain             = gain[iscope];

      vsassert(m_sup->nchan(iscope) >= nchan);
      vsassert(m_ped->nchan(iscope) >= nchan);
      vsassert((!laser)||(permissive_laser)
	       ||(laser->scope[iscope].chan.size() >= nchan));
      vsassert((!hilo)||(iscope>=hilo->scope.size())
	       ||(hilo->scope[iscope].chan.empty())
	       ||(hilo->scope[iscope].chan.size() >= nchan));
      vsassert((!m_pad_sup)||(m_pad_sup->nchan(iscope) >= nchan));
      vsassert((!m_pad_ped)||(m_pad_ped->nchan(iscope) >= nchan));

      std::vector<double> chan_median_dev(nchan);

      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  const unsigned icrate = channel_map->crateForChannel(iscope,ichan);
	  const unsigned il2 = channel_map->l2ChannelForCrate(iscope, icrate);
	  ChanMergedCal& chan(scope[iscope].channel[ichan]);
	  chan.ped_hi                       = m_ped->hiPed(iscope,ichan);
	  chan.ped_lo                       = m_ped->loPed(iscope,ichan);
	  chan.median_dev                   = 0;
	  chan.hi_lo_gain_ratio             = 6.0;
	  chan.gain                         = 1.0;
	  chan.has_pmt                      = true;
	  chan.suppress_all_events          = false;
	  chan.suppress_lo_gain             = false;
	  chan.suppress_nslice              = 0;
	  chan.suppress_reason              = 0;
	  chan.crate                        = icrate;
	  chan.l2chan                       = il2;
	  chan.l2time                       = 0.0;
	  chan.chantime                     = 0.0;
	  chan.cratetime                    = 0.0;
	  chan.median_current               = 0;
	  chan.median_l1_rate               = 0;
	  chan.pmt_multiplier_gain          = 1.0;
	  chan.pmt_efficiency               = 1.0;
	  chan.pmt_coord_E_deg              = nochan;
	  chan.pmt_coord_U_deg              = nochan;

	  if(ichan < array_info.scopes[iscope].pij.size())
	    {
	      chan.pmt_coord_E_deg = 
		array_info.scopes[iscope].pij[ichan].x() * 180.0/M_PI;
	      chan.pmt_coord_U_deg = 
		array_info.scopes[iscope].pij[ichan].y() * 180.0/M_PI;
	    }

	  chan.median_dev                   =
	    m_ped->medianDev(iscope,ichan,m_window_size,m_ped_event_count_min);

	  scope[iscope].median_dev_hist.accumulate(chan.median_dev);
	  chan_median_dev[ichan]            = chan.median_dev;

	  if(iscope < hi_lo_gain_ratio.size())
	    chan.hi_lo_gain_ratio = hi_lo_gain_ratio[iscope];

	  if((hilo)&&(iscope<hilo->scope.size())
	     &&(!hilo->scope[iscope].chan.empty()))
	    {
	      const VSHiLoData::ChanData& C(hilo->scope[iscope].chan[ichan]);
	      if(C.has_lo_gain_ped)
		chan.ped_lo                 = C.lo_gain_ped;
	      else
		chan.suppress_lo_gain       = true;
	      if(C.has_hi_lo_gain_ratio)
		chan.hi_lo_gain_ratio       = C.hi_lo_gain_ratio;
	    }

	  if(std::isnan(chan.ped_hi) || chan.ped_hi<0)
	    {
	      chan.has_pmt                  = false;
	      chan.suppress_all_events      = true;
	      chan.suppress                 = true;
	      chan.suppress_reason         |= ChanMergedCal::SR_ILLEGAL_PEDVAL;
	    }

	  if(std::isnan(chan.ped_lo) || chan.ped_lo<0)
	    {
	      chan.suppress_lo_gain         = true;
	      chan.suppress_reason         |= ChanMergedCal::SR_LO_GAIN_PEDVAL;
	    }

	  if(m_sup->isAlwaysSuppressed(iscope,ichan))
	    {
	      chan.suppress_all_events      = true;
	      chan.suppress_reason         |= ChanMergedCal::SR_PED;
	    }

	  if((m_pad_sup)&&(m_pad_sup->isAlwaysSuppressed(iscope,ichan)))
	    {
	      chan.suppress_all_events      = true;
	      chan.suppress_reason         |= ChanMergedCal::SR_PAD_PED;
	    }

	  if((laser)&&(iscope<laser->scope.size())&&
	     (ichan<laser->scope[iscope].chan.size()))
	    {
	      const VSLaserData::ChanData& 
		C(laser->scope[iscope].chan[ichan]);
	      chan.crate                    = C.crate;
	      chan.l2chan                   = C.l2chan;
	      chan.chantime                 = C.chantime;
	      chan.l2time                   = C.l2time;
	      chan.cratetime                = C.cratetime;
	      chan.pmt_multiplier_gain      = C.absgain;
	      chan.pmt_efficiency           = C.eff;

	      if(!unity_gain)
		chan.gain                   = 1.0/C.gain;

	      if(C.uncalibrated)
		{
		  chan.suppress_all_events  = true;
		  chan.suppress_reason     |= ChanMergedCal::SR_LASER;
		}
	    }

	  if((db_data)
	     &&(iscope<db_data->scope.size())
	     &&(db_data->scope[iscope].has_scope))
	    {
	      if((!db_data->scope[iscope].hv_status.empty())
		 &&(ichan < db_data->scope[iscope].hv_status.front().
		    chan.size()))
		{
		  unsigned nmeas = db_data->scope[iscope].hv_status.size();
		  std::vector<double> current(nmeas);
		  for(unsigned imeas=0;imeas<nmeas;imeas++)
		    current[imeas]          = 
		      db_data->scope[iscope].hv_status[imeas].
		      chan[ichan].current;
		  chan.median_current       = median(current);
		}

	      if((!db_data->scope[iscope].l1_rate.empty())
		 &&(ichan < db_data->scope[iscope].l1_rate.front().
		    rate.size()))
		{
		  unsigned nmeas = db_data->scope[iscope].l1_rate.size();
		  std::vector<double> l1_rate(nmeas);
		  for(unsigned imeas=0;imeas<nmeas;imeas++)
		    l1_rate[imeas]          =
		      db_data->scope[iscope].l1_rate[imeas].rate[ichan];
		  chan.median_l1_rate = median(l1_rate);
		}
	    }


	  chan.clean_multiplier        = 1.0;
	  chan.suppress                = chan.suppress_all_events;
	  chan.dopad                   = false;
	  chan.dev                     = 0.0;
	  chan.pad                     = 0.0;
	}

      if(auto_suppress_l2_chan)
	{
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    {
	      unsigned l2chan = scope[iscope].channel[ichan].l2chan;
	      if(l2chan<nchan)
		{
		  ChanMergedCal& chan(scope[iscope].channel[l2chan]);
		  chan.has_pmt              = false;
		  chan.suppress_all_events  = true;
		  chan.suppress             = true;
		  chan.suppress_reason     |= ChanMergedCal::SR_NO_PMT;
		}
	    }
	}

      scope[iscope].median_median_dev  = median(chan_median_dev);

      if(iscope<dev_scaling.size() && dev_scaling[iscope].first>0)
	{
	  double s = 
	    scope[iscope].median_median_dev * dev_scaling[iscope].first
	    + dev_scaling[iscope].second;

	  scope[iscope].scaled_dev     = s;

	  scaled_dev_sum += s;
	  scaled_dev_sumsq += s*s;
	  nscaled_dev++;
	}
    }

  mean_scaled_dev = scaled_dev_sum/double(nscaled_dev);
  rms_scaled_dev = sqrt(scaled_dev_sumsq/double(nscaled_dev) 
			- mean_scaled_dev*mean_scaled_dev);

  // Demand suppressed channels -----------------------------------------------

  for(unsigned isup=0;isup<sup_chan.size();isup++)
    {
      unsigned iscope = sup_chan[isup].first;
      unsigned ichan = sup_chan[isup].second;
      if((iscope<scope.size())&&(ichan<scope[iscope].nchan))
	{
	  ChanMergedCal& chan(scope[iscope].channel[ichan]);
	  chan.suppress_all_events     = true;
      	  chan.suppress                = true;
	  chan.suppress_reason        |= ChanMergedCal::SR_DEMAND;
	}
    }

  // Demand no PMT channels ---------------------------------------------------

  for(unsigned inopmt=0;inopmt<no_pmt_chan.size();inopmt++)
    {
      unsigned iscope = no_pmt_chan[inopmt].first;
      unsigned ichan = no_pmt_chan[inopmt].second;
      if((iscope<scope.size())&&(ichan<scope[iscope].nchan))
	{
	  ChanMergedCal& chan(scope[iscope].channel[ichan]);
	  chan.has_pmt                 = false;
	  chan.suppress_all_events     = true;
      	  chan.suppress                = true;
	  chan.suppress_reason        |= ChanMergedCal::SR_NO_PMT;
	}
    }

  // Copy to cache variables prior to calling merge functions -----------------

  setFixedCacheVariables();

  // Merge suppression info ---------------------------------------------------

  mergeSupSlice();
  if((m_cs==CS_PEDRMS)||(m_pad_ped))mergePedSlice();	

  // Count number of suppressed slices ----------------------------------------

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      ScopeMergedCal& S(scope[iscope]);
      VSChannelSuppressedCount& SC(S.suppressed_count);
      const unsigned nchan = S.nchan;
      if(nchan==0)continue;

      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  ChanMergedCal& C(S.channel[ichan]);

	  if(C.suppress_all_events)
	    C.suppress_nslice          = nslice;

	  unsigned reason              = C.suppress_reason;

	  unsigned nsupslice_ped = 0;
	  unsigned nsupslice_pad = 0;
	  unsigned nsupslice = 0;

	  if(((reason & ChanMergedCal::SR_PED) == 0)
	     &&((reason & ChanMergedCal::SR_PAD_PED) == 0))
	    for(unsigned islice=0;islice<nslice;islice++)
	      if(m_sup->isSuppressed(islice, iscope, ichan))
		{
		  nsupslice_ped++;
		  nsupslice++;
		  if((m_pad_sup)
		     &&(m_pad_sup->isSuppressed(islice, iscope, ichan)))
		    nsupslice_pad++;
		}
	      else if((m_pad_sup)
		      &&(m_pad_sup->isSuppressed(islice, iscope, ichan)))
		{
		  nsupslice_pad++;
		  nsupslice++;
		}

	  if(nsupslice_ped*2 >= nslice)
	    reason                    |= ChanMergedCal::SR_PED,
	      C.suppress_reason       |= ChanMergedCal::SR_PED;
	  if(nsupslice_pad*2 >= nslice)
	    reason                    |= ChanMergedCal::SR_PAD_PED,
	      C.suppress_reason       |= ChanMergedCal::SR_PAD_PED;
	  if(!C.suppress_all_events)
	    C.suppress_nslice          = nsupslice;

	  switch(reason&(ChanMergedCal::SR_PAD_PED
			 |ChanMergedCal::SR_PED|ChanMergedCal::SR_LASER))
	    {
	    case ChanMergedCal::SR_PED:
	      SC.ped++;
	      break;
	    case ChanMergedCal::SR_PAD_PED:
	      SC.pad++;
	      break;
	    case ChanMergedCal::SR_LASER:
	      SC.laser++;
	      break;
	    case (ChanMergedCal::SR_PED
		  |ChanMergedCal::SR_PAD_PED):
	      SC.ped_pad++;
	      break;
	    case (ChanMergedCal::SR_PED
		  |ChanMergedCal::SR_LASER):
	      SC.ped_laser++;
	      break;
	    case (ChanMergedCal::SR_PAD_PED
		  |ChanMergedCal::SR_LASER):
	      SC.pad_laser++;
	      break;
	    case (ChanMergedCal::SR_PAD_PED
		  |ChanMergedCal::SR_PED
		  |ChanMergedCal::SR_LASER):
	      SC.all++;
	      break;
	    default:
	      break;
	    }

	  if(reason)SC.any++;
	  else SC.none++;
	}
    }
}

void VBFAnalysisStage2::ArrayMergedCal::print(const std::ostream& stream)
{

}

void VBFAnalysisStage2::ArrayMergedCal::save(VSOctaveH5WriterStruct* s) const
{
  s->writeComposite("array_info",*this);
  const unsigned nscope = scope.size();  
  VSOctaveH5WriterCellVector* wc_scope = 
    s->writeCellVector("scope_info", nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope].nchan)
      {
	VSOctaveH5WriterStruct* ws = wc_scope->writeStruct(iscope);
	vsassert(ws);
	scope[iscope].save(ws);
	delete ws;
      }
  delete wc_scope;

  VSOctaveH5WriterCellVector* wc_chan = 
    s->writeCellVector("channel_info", nscope);
  vsassert(wc_chan);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope].nchan)
      {
	VSOctaveH5WriterStruct* ws = wc_chan->writeStruct(iscope);	
	vsassert(ws);
	ws->writeComposite("suppressed_count",scope[iscope].suppressed_count);
	ws->writeCompositeVectorHere(scope[iscope].channel);
	delete ws;
      }
  delete wc_chan;
}

void VBFAnalysisStage2::ArrayMergedCal::mergeSupSlice()
{
  unsigned iword;
  unsigned mask;
  m_sup->getFastSliceFindInfo(m_sup_islice, iword, mask);
    
  const unsigned nscope=scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const unsigned nchan = scope[iscope].nchan;
      for(unsigned ichan=0;ichan<nchan;ichan++)
	if(!scope[iscope].ch_suppress_all_events[ichan])
	  {
	    bool s = m_sup->isSuppressedFast(iword, mask, iscope, ichan);
	    scope[iscope].channel[ichan].suppress = s;
	    scope[iscope].ch_suppress[ichan] = s;
	  }
    }

  if(m_pad_sup)
    {
      unsigned nslice = m_pad_sup->nslice();
      if(nslice)
	{
	  unsigned islice = m_sup_islice;
	  if(islice >= nslice)islice=nslice-1;
	  unsigned pad_iword;
	  unsigned pad_mask;
	  m_pad_sup->getFastSliceFindInfo(islice, pad_iword, pad_mask);
	  
	  for(unsigned iscope=0;iscope<nscope;iscope++)
	    {
	      const unsigned nchan = scope[iscope].nchan;
	      for(unsigned ichan=0;ichan<nchan;ichan++)
		if((!scope[iscope].ch_suppress_all_events[ichan])
		   &&(m_pad_sup->isSuppressedFast(pad_iword, pad_mask, 
						  iscope, ichan)))
		  {
		    scope[iscope].channel[ichan].suppress = true;
		    scope[iscope].ch_suppress[ichan] = true;
		  }
	    }
	}
    }
}

void VBFAnalysisStage2::ArrayMergedCal::mergePedSlice()
{
  const unsigned nscope=scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const unsigned nchan = scope[iscope].nchan;
      for(unsigned ichan=0;ichan<nchan;ichan++)
	scope[iscope].channel[ichan].dev =
	  m_ped->dev(m_ped_islice, iscope, ichan, m_window_size);
    }

  if(m_pad_ped)
    {
      const unsigned nslice = m_pad_ped->nslice();
      if(nslice)
	{
	  unsigned islice = m_ped_islice;
	  if(islice >= nslice)islice=nslice-1;

	  for(unsigned iscope=0;iscope<nscope;iscope++)
	    {
	      const unsigned nchan = scope[iscope].nchan;
	      for(unsigned ichan=0;ichan<nchan;ichan++)
		{
		  const double dev = 
		    scope[iscope].channel[ichan].dev;
		  const double pad_dev = 
		    m_pad_ped->dev(islice, iscope, ichan, m_window_size);
		  if(pad_dev > dev)
		    {
		      const double pad = sqrt(pad_dev*pad_dev - dev*dev);
		      scope[iscope].channel[ichan].dopad = true;
		      scope[iscope].channel[ichan].dev   = pad_dev;
		      scope[iscope].channel[ichan].pad   = pad;
		      scope[iscope].ch_dopad[ichan]      = true;
		      scope[iscope].ch_pad[ichan]        = pad;
		    }
		  else
		    {
		      scope[iscope].channel[ichan].dopad = false;
		      scope[iscope].channel[ichan].pad   = 0.0;
		      scope[iscope].ch_dopad[ichan]      = false;
		      scope[iscope].ch_pad[ichan]        = 0.0;
		    }
		}
	    }
	}
    }

  if(m_cs==CS_PEDRMS)
    for(unsigned iscope=0;iscope<nscope;iscope++)
      {
	const unsigned nchan = scope[iscope].nchan;
	for(unsigned ichan=0;ichan<nchan;ichan++)
	  {
	    double cm = 1.0/scope[iscope].channel[ichan].dev;
	    scope[iscope].channel[ichan].clean_multiplier = cm;
	    scope[iscope].ch_clean_multiplier[ichan] = cm;
	  }	      
      }
}

void VBFAnalysisStage2::ArrayMergedCal::setFixedCacheVariables()
{
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      VBFAnalysisStage2::ScopeMergedCal& scal(scope[iscope]);

      unsigned nchan = scal.channel.size();
      if(nchan==0)continue;

      scal.ch_has_pmt.resize(nchan);
      scal.ch_ped_hi.resize(nchan);
      scal.ch_ped_lo.resize(nchan);
      scal.ch_hi_lo_gain_ratio.resize(nchan);
      scal.ch_l2chan.resize(nchan);
      scal.ch_l2time.resize(nchan);
      scal.ch_chantime.resize(nchan);
      scal.ch_suppress.resize(nchan);
      scal.ch_suppress_all_events.resize(nchan);
      scal.ch_suppress_lo_gain.resize(nchan);
      scal.ch_dopad.resize(nchan);
      scal.ch_pad.resize(nchan);
      scal.ch_clean_multiplier.resize(nchan);
      scal.ch_gain.resize(nchan);

      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  VBFAnalysisStage2::ChanMergedCal& ccal(scal.channel[ichan]);
	  scal.ch_has_pmt[ichan]             = ccal.has_pmt;
	  scal.ch_ped_hi[ichan]              = ccal.ped_hi;
	  scal.ch_ped_lo[ichan]              = ccal.ped_lo;
	  scal.ch_hi_lo_gain_ratio[ichan]    = ccal.hi_lo_gain_ratio;
	  scal.ch_l2chan[ichan]              = ccal.l2chan;
	  scal.ch_l2time[ichan]              = ccal.l2time;
	  scal.ch_chantime[ichan]            = ccal.chantime;
	  scal.ch_suppress[ichan]            = ccal.suppress;
	  scal.ch_suppress_all_events[ichan] = ccal.suppress_all_events;
	  scal.ch_suppress_lo_gain[ichan]    = ccal.suppress_lo_gain;
	  scal.ch_dopad[ichan]               = ccal.dopad;
	  scal.ch_pad[ichan]                 = ccal.pad;
	  scal.ch_clean_multiplier[ichan]    = ccal.clean_multiplier;
	  scal.ch_gain[ichan]                = ccal.gain;
	}
    }
}

