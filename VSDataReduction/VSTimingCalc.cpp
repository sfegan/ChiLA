//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTimingCalc.cpp

  Calculate timing and integrated charge information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/19/2006

  $Id: VSTimingCalc.cpp,v 3.5 2008/02/23 22:41:34 sfegan Exp $

*/

#include <vsassert>

#include "fast_alloc.hpp"
#include "VBFSimplePeds.hpp" // only needed for definition of uint32_t from VBF
#include "VSTimingCalc.hpp"

using namespace VERITAS;

// ============================================================================
// VSTimingCalc
// ============================================================================

VSTimingCalc::~VSTimingCalc()
{
  // nothing to see here
}

// ============================================================================
// VSIntegralTimingCalc
// ============================================================================

VSIntegralTimingCalc::~VSIntegralTimingCalc()
{
  // nothing to see here
}

bool VSIntegralTimingCalc::
calc(double& signal, double& time,
     bool lo_gain, double lg_multiplier, unsigned nsample,
     const uint32_t* samples, const uint32_t* integrated, double ped)
{
  // --------------------------------------------------------------------------
  // Integrate signal
  // --------------------------------------------------------------------------

  if((nsample>m_sample_N)&&(m_sample_N>0))nsample=m_sample_N;
  std::vector<double> sample_signal(nsample-m_sample_0);
  double my_signal = 0;
  for(unsigned isample = m_sample_0; isample<nsample; isample++)
    {
      sample_signal[isample-m_sample_0] = my_signal;
      my_signal += double(samples[isample]) - ped;
    }
  signal = my_signal;
  if(lo_gain)signal*=lg_multiplier;

  // --------------------------------------------------------------------------
  // Interpolate point at which integral reaches threshold
  // --------------------------------------------------------------------------
  double signal_threshold = my_signal*m_threshold_frac;
  unsigned isample = nsample-m_sample_0;
  while((isample>1)
	&&(sample_signal[isample-1]>signal_threshold))
    isample--;
  double y1 = sample_signal[isample-1];
  double y2 = (isample==nsample)?signal:sample_signal[isample];
  time = (signal_threshold-y1)/(y2-y1)+isample-1+m_sample_0;

  return true;
}

// ============================================================================
// VSPeakTimingCalc
// ============================================================================

VSPeakTimingCalc::~VSPeakTimingCalc()
{
  // nothing to see here
}

bool VSPeakTimingCalc::
calc(double& signal, double& time,
     bool lo_gain, double lg_multiplier, unsigned nsample,
     const uint32_t* samples, const uint32_t* integrated, double ped)
{
  // --------------------------------------------------------------------------
  // Integrate signal and find peak
  // --------------------------------------------------------------------------

  if((nsample>m_sample_N)&&(m_sample_N>0))nsample=m_sample_N;
  std::vector<double> my_samples(nsample);
  double my_signal = 0;
  double peak = 0;
  unsigned ipeak = 0;
  for(unsigned isample = m_sample_0; isample<nsample; isample++)
    {
      unsigned sigsample = isample-m_sample_0;
      double s = double(samples[isample]) - ped;
      my_samples[sigsample] = s;
      my_signal += s;
      if((isample==0)||(s>peak))ipeak=sigsample, peak=s;
    }
  signal = my_signal;
  if(lo_gain)signal*=lg_multiplier;

  // --------------------------------------------------------------------------
  // Interpolate point at which signal reaches threshold
  // --------------------------------------------------------------------------

  double signal_threshold = peak*m_threshold_frac;
  unsigned isample = ipeak;
  while((isample>0)&&(my_samples[isample-1]>=signal_threshold))
    isample--;
  double y1 = (isample==0)?0:my_samples[isample-1];
  double y2 = my_samples[isample];
  time = (signal_threshold-y1)/(y2-y1)+isample+m_sample_0-1;

  // --------------------------------------------------------------------------
  // DEBUGGING
  // --------------------------------------------------------------------------

#ifdef DEBUG
  if(signal > 100)
    {
      for(unsigned isamp=0;isamp<samples.size();isamp++)
	{
	  if(isamp!=0)std::cerr << ' ';
	  if(isamp==m_sample_0)std::cerr << '[';
	  else if((isamp==m_sample_N)&&(m_sample_N!=0))std::cerr << ']';
	  else std::cerr << ' ';
	  std::cerr << ' ';
	  std::cerr << std::setw(3) << unsigned(samples[isamp]);
	}
      std::cerr << ' ';
      if((samples.size()==m_sample_N)||(m_sample_N==0))std::cerr << ']';
      else std::cerr << ' ';
      std::cerr << ' ' << signal << std::endl;
    }
#endif

  return true;
}

// ============================================================================
// VSConstantFractionCalc
// ============================================================================

VSConstantFractionCalc::~VSConstantFractionCalc()
{
  // nothing to see here
}

bool VSConstantFractionCalc::
calc(double& signal, double& time,
     bool lo_gain, double lg_multiplier, unsigned nsample,
     const uint32_t* samples, const uint32_t* integrated, double ped)
{
  unsigned window_width;
  double window_start;
  double total_signal;
  const uint32_t* peak_ptr;
  const unsigned sample_zero=lo_gain?m_lo_sample_zero:m_hi_sample_zero;
  return calc(signal, time, window_width, window_start, total_signal, peak_ptr,
	      lo_gain, lg_multiplier, nsample, samples, integrated, ped,
	      sample_zero);
}

bool VSConstantFractionCalc::
calc(double& signal, double& time, 
     unsigned& window_width, double& window_start,
     double& total_signal, const uint32_t*& peak_ptr,
     const bool lo_gain, const double lg_multiplier, unsigned nsample, 
     const uint32_t*const samples, const uint32_t*const integrated, 
     const double ped, const unsigned sample_zero)
{
  const uint32_t* int_nsample_ptr = integrated+nsample;
  const uint32_t* int_isample_ptr = integrated+sample_zero;

  const uint32_t integrated_0 = sample_zero ? *(int_isample_ptr-1) : 0;
  const uint32_t total_digital = *(int_nsample_ptr-1) - integrated_0;

  // --------------------------------------------------------------------------
  // Find peak
  // --------------------------------------------------------------------------

#ifdef VECTORIZE_FRIENDLY
  unsigned peak = samples[sample_zero];
  unsigned ipeak = sample_zero;
  for(unsigned isample=sample_zero+1;isample<nsample;isample++)
    if(samples[isample]>peak)peak=samples[isample],ipeak=isample;
  peak_ptr = samples+ipeak;
#else
  const uint32_t* nsample_ptr = samples+nsample;
  const uint32_t* isample_ptr = samples+sample_zero;

  unsigned peak = *isample_ptr;
  peak_ptr = isample_ptr;

  while(++isample_ptr<nsample_ptr)
    {
      const unsigned sample = *isample_ptr;
      if(sample>peak)peak_ptr=isample_ptr, peak=sample;
    }

  unsigned ipeak = peak_ptr - samples;
#endif

  // --------------------------------------------------------------------------
  // Calculate window width
  // --------------------------------------------------------------------------

  window_start = m_window_start;
  window_width = m_window_width;

  unsigned max_sample_width = nsample-sample_zero;
  total_signal = double(total_digital)-ped*max_sample_width;
  if(lo_gain)total_signal*=lg_multiplier;
  
  if((m_apply_increase)&&(total_signal > m_threshold_charge))
    {
      double increase_factor = log10(total_signal/m_threshold_charge);
      window_start += m_window_start_increase*increase_factor;
      window_width += 
	lround(m_window_width_increase*increase_factor);
      //unsigned(round(m_window_width_increase*increase_factor));
      if(window_width > max_sample_width)
	window_width=max_sample_width;
    }

  // --------------------------------------------------------------------------
  // Interpolate point at which signal reaches threshold
  // --------------------------------------------------------------------------

  double threshold = (double(peak)-ped)*m_threshold_frac;
  double t2 = threshold+ped;
  unsigned u_threshold = (t2>0)?unsigned(ceil(t2)):0;
  
  unsigned isample = ipeak;
  while((isample>sample_zero)&&(samples[isample-1]>=u_threshold))isample--;

  double y1 = (isample==sample_zero)?0:double(samples[isample-1])-ped;
  double y2 = double(samples[isample])-ped;
  time = (threshold-y1)/(y2-y1)+isample-1;

#if defined (DEBUG)
  unsigned rememberme = isample;
#endif 

  // --------------------------------------------------------------------------
  // Finally, calculate charge
  // --------------------------------------------------------------------------

  isample = 
    (time>window_start)?lround(time-window_start):0;
  //(time>window_start)?unsigned(round(time-window_start)):0;
  if(isample<sample_zero)
    {
      isample=sample_zero;
      if(isample+window_width > nsample)
	{
	  time = 0;
	  signal = 0;
	  return false;
	}
    }
  else if(isample+window_width > nsample)
    {
      isample=nsample-window_width;
      if(isample<sample_zero)
	{
	  time = 0;
	  signal = 0;
	  return false;
	}
    }

  if(isample==sample_zero)
    signal = double(integrated[isample+window_width-1]
		    - integrated_0);
  else 
    signal = double(integrated[isample+window_width-1] 
		    - integrated[isample-1]);
     
  signal -= ped*window_width;
  if(lo_gain)signal*=lg_multiplier;

  // --------------------------------------------------------------------------
  // DEBUGGING
  // --------------------------------------------------------------------------

#if defined(DEBUG)
  if(signal>100)
    {
      std::cerr << threshold << ' ' << t2 << ' ' << u_threshold << ' '
		<< y1 << ' ' << y2 << ' ' << (threshold-y1)/(y2-y1) << ' '
		<< ' ' << ipeak << ' ' << time
		<< std::endl;

      for(unsigned isamp=0;isamp<samples.size();isamp++)
	{
	  if(isamp!=0)std::cerr << ' ';
	  if(isamp==isample)std::cerr << '[';
	  else if(isamp==isample+window_width)std::cerr << ']';
	  else if(isamp==unsigned(floor(time))+1)std::cerr << '/';
	  else if(isamp==sample_zero)std::cerr << 'X';
	  else std::cerr << ' ';
	  std::cerr << ' ';
	  std::cerr << std::setw(3) << unsigned(samples[isamp]);
	}
      std::cerr << ' ';
      if(samples.size()==isample+window_width)std::cerr << ']';
      else std::cerr << ' ';
      std::cerr << ' ' << signal << std::endl;
    }
#endif

  return true;
}

#if 0
// ============================================================================
// VSMaxChargeCalc
// ============================================================================

VSMaxChargeCalc::~VSMaxChargeCalc()
{
  // nothing to see here
}

bool VSMaxChargeCalc::
calc(double& signal, double& time,
     bool lo_gain, double lg_multiplier, unsigned nsample,
     const uint32_t* samples, const uint32_t* integrated, double ped)
{
  unsigned window_width;
  double window_start;
  double total_signal;
  const uint32_t* peak_ptr;
  const unsigned sample_zero=lo_gain?m_lo_sample_zero:m_hi_sample_zero;
  return calc(signal, time, window_width, window_start, total_signal, peak_ptr,
	      lo_gain, lg_multiplier, nsample, samples, integrated, 
	      ped, sample_zero);
}

bool VSMaxChargeCalc::
calc(double& signal, double& time, 
     unsigned& window_width, double& window_start,
     double& total_signal, const uint32_t*& peak_ptr,
     const bool lo_gain, double lg_multiplier, 
     unsigned nsample, const uint32_t*const samples,
     const uint32_t*const integrated, const double ped,
     const unsigned sample_zero)
{
  const uint32_t* int_nsample_ptr = integrated+nsample;
  const uint32_t* int_isample_ptr = integrated+sample_zero;

  const uint32_t integrated_0 = sample_zero ? *(int_isample_ptr-1) : 0;
  const uint32_t total_digital = *(int_nsample_ptr-1) - integrated_0;

  // --------------------------------------------------------------------------
  // Find peak
  // --------------------------------------------------------------------------

#ifdef VECTORIZE_FRIENDLY
  unsigned peak = samples[sample_zero];
  unsigned ipeak = sample_zero;
  for(unsigned isample=sample_zero+1;isample<nsample;isample++)
    if(samples[isample]>peak)peak=samples[isample],ipeak=isample;
  peak_ptr = samples+ipeak;
#else
  const uint32_t* nsample_ptr = samples+nsample;
  const uint32_t* isample_ptr = samples+sample_zero;

  unsigned peak = *isample_ptr;
  peak_ptr = isample_ptr;

  while(++isample_ptr<nsample_ptr)
    {
      const unsigned sample = *isample_ptr;
      if(sample>peak)peak_ptr=isample_ptr, peak=sample;
    }

  unsigned ipeak = peak_ptr - samples;
#endif

  // --------------------------------------------------------------------------
  // Calculate window width
  // --------------------------------------------------------------------------

  window_start = m_window_start;
  window_width = m_window_width;

  unsigned max_sample_width = nsample-sample_zero;
  total_signal = double(total_digital)-ped*max_sample_width;
  if(lo_gain)total_signal*=lg_multiplier;
  
  if((m_apply_increase)&&(total_signal > m_threshold_charge))
    {
      double increase_factor = log10(total_signal/m_threshold_charge);
      window_start += m_window_start_increase*increase_factor;
      window_width += 
	lround(m_window_width_increase*increase_factor);
      //unsigned(round(m_window_width_increase*increase_factor));
      if(window_width > max_sample_width)
	window_width=max_sample_width;
    }

  // --------------------------------------------------------------------------
  // Interpolate point at which signal reaches threshold
  // --------------------------------------------------------------------------

  double threshold = (double(peak)-ped)*m_threshold_frac;
  double t2 = threshold+ped;
  unsigned u_threshold = (t2>0)?unsigned(ceil(t2)):0;
  
  unsigned isample = ipeak;
  while((isample>sample_zero)&&(samples[isample-1]>=u_threshold))isample--;

  double y1 = (isample==sample_zero)?0:double(samples[isample-1])-ped;
  double y2 = double(samples[isample])-ped;
  time = (threshold-y1)/(y2-y1)+isample-1;

#if defined (DEBUG)
  unsigned rememberme = isample;
#endif 

  // --------------------------------------------------------------------------
  // Finally, calculate charge
  // --------------------------------------------------------------------------

  isample = 
    (time>window_start)?lround(time-window_start):0;
  //(time>window_start)?unsigned(round(time-window_start)):0;
  if(isample<sample_zero)
    {
      isample=sample_zero;
      if(isample+window_width > nsample)
	{
	  time = 0;
	  signal = 0;
	  return false;
	}
    }
  else if(isample+window_width > nsample)
    {
      isample=nsample-window_width;
      if(isample<sample_zero)
	{
	  time = 0;
	  signal = 0;
	  return false;
	}
    }

  if(isample==sample_zero)
    signal = double(integrated[isample+window_width-1]
		    - integrated_0);
  else 
    signal = double(integrated[isample+window_width-1] 
		    - integrated[isample-1]);
     
  signal -= ped*window_width;
  if(lo_gain)signal*=lg_multiplier;

  // --------------------------------------------------------------------------
  // DEBUGGING
  // --------------------------------------------------------------------------

#if defined(DEBUG)
  if(signal>100)
    {
      std::cerr << threshold << ' ' << t2 << ' ' << u_threshold << ' '
		<< y1 << ' ' << y2 << ' ' << (threshold-y1)/(y2-y1) << ' '
		<< ' ' << ipeak << ' ' << time
		<< std::endl;

      for(unsigned isamp=0;isamp<samples.size();isamp++)
	{
	  if(isamp!=0)std::cerr << ' ';
	  if(isamp==isample)std::cerr << '[';
	  else if(isamp==isample+window_width)std::cerr << ']';
	  else if(isamp==unsigned(floor(time))+1)std::cerr << '/';
	  else if(isamp==sample_zero)std::cerr << 'X';
	  else std::cerr << ' ';
	  std::cerr << ' ';
	  std::cerr << std::setw(3) << unsigned(samples[isamp]);
	}
      std::cerr << ' ';
      if(samples.size()==isample+window_width)std::cerr << ']';
      else std::cerr << ' ';
      std::cerr << ' ' << signal << std::endl;
    }
#endif

  return true;
}

#endif
