//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTimingCalc.hpp

  Calculate timing and integrated charge information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/19/2006

  $Id: VSTimingCalc.hpp,v 3.2 2008/02/23 22:41:34 sfegan Exp $

*/

#ifndef VSTOFFSETCALC_HPP
#define VSTOFFSETCALC_HPP

#include<VSSimpleVBF.hpp>

#include<vector>

namespace VERITAS
{

  // ==========================================================================
  // VSTimingCalc
  // ==========================================================================

  class VSTimingCalc
  {
  public:
    virtual ~VSTimingCalc();
    virtual bool calc(double& signal, double& time,
		      bool lo_gain, double lg_multiplier, unsigned nsample, 
		      const uint32_t* samples, const uint32_t* integrated,
		      double ped) = 0;
  };
  
  // ==========================================================================
  // VSIntegralTimingCalc
  // ==========================================================================

  class VSIntegralTimingCalc: public VSTimingCalc
  {
  public:
    VSIntegralTimingCalc(double threshold_frac,
			  unsigned sample_0=0, unsigned sample_N=0):
      m_threshold_frac(threshold_frac),
      m_sample_0(sample_0), m_sample_N(sample_N)
    { /* nothing to see here */ }
    virtual ~VSIntegralTimingCalc();
    virtual bool calc(double& signal, double& time,
		      bool lo_gain, double lg_multiplier, unsigned nsample, 
		      const uint32_t* samples, const uint32_t* integrated, 
		      double ped);
    void setWindow(unsigned s0, unsigned sN) { m_sample_0=s0; m_sample_N=sN; }
  private:
    double                  m_threshold_frac;
    unsigned                m_sample_0;
    unsigned                m_sample_N;
  };

  // ==========================================================================
  // VSPeakTimingCalc
  // ==========================================================================
  
  class VSPeakTimingCalc: public VSTimingCalc
  {
  public:
    VSPeakTimingCalc(double threshold_frac,
		     unsigned sample_0=0, unsigned sample_N=0):
      m_threshold_frac(threshold_frac),
      m_sample_0(sample_0), m_sample_N(sample_N)
    { /* nothing to see here */ }
    virtual ~VSPeakTimingCalc();
    virtual bool calc(double& signal, double& time,
		      bool lo_gain, double lg_multiplier, unsigned nsample,
		      const uint32_t* samples, const uint32_t* integrated, 
		      double ped);
    void setWindow(unsigned s0, unsigned sN) { m_sample_0=s0; m_sample_N=sN; }
  private:
    double                  m_threshold_frac;
    unsigned                m_sample_0;
    unsigned                m_sample_N;
  };

  // ==========================================================================
  // VSConstantFractionCalc
  // ==========================================================================

  class VSConstantFractionCalc: public VSTimingCalc
  {
  public:
    VSConstantFractionCalc(double threshold_frac,
			   double window_start, unsigned window_width,
			   unsigned hi_sample_zero, unsigned lo_sample_zero,
			   double threshold_charge = 0,
			   double window_start_increase = 0,
			   double window_width_increase = 0):
      m_threshold_frac(threshold_frac),
      m_window_start(window_start), m_window_width(window_width),
      m_hi_sample_zero(hi_sample_zero), m_lo_sample_zero(lo_sample_zero),
      m_apply_increase(threshold_charge > 0),
      m_threshold_charge(threshold_charge), 
      m_window_start_increase(window_start_increase),
      m_window_width_increase(window_width_increase)
    { /* nothing to see here */ }
    virtual ~VSConstantFractionCalc();
    virtual bool calc(double& signal, double& time,
		      bool lo_gain, double lg_multiplier, unsigned nsample, 
		      const uint32_t* samples, const uint32_t* integrated, 
		      double ped);

    // This function presents an "extended" API that can return
    // additional information and allows for different window_zero to
    // be used (for example for different telescopes)

    bool calc(double& signal, double& time, 
	      unsigned& window_width, double& window_start,
	      double& total_signal, const uint32_t*& peak_ptr,
	      const bool lo_gain, double lg_multiplier, unsigned nsample, 
	      const uint32_t*const samples, const uint32_t*const integrated,
	      const double ped, const unsigned sample_zero);

  private:
    double                  m_threshold_frac;
    double                  m_window_start;
    unsigned                m_window_width;
    unsigned                m_hi_sample_zero;
    unsigned                m_lo_sample_zero;
    bool                    m_apply_increase;
    double                  m_threshold_charge;
    double                  m_window_start_increase;
    double                  m_window_width_increase;
  };
#if 0  
  // ==========================================================================
  // VSMaxchargeCalc
  // ==========================================================================

  class VSMaxChargeCalc: public VSTimingCalc
  {
  public:
    VSMaxChargeCalc(unsigned window_width,
		    unsigned hi_sample_zero, unsigned lo_sample_zero,
		    double threshold_charge = 0,
		    double window_start_increase = 0,
		    double window_width_increase = 0):
      m_window_width(window_width),
      m_hi_sample_zero(hi_sample_zero), m_lo_sample_zero(lo_sample_zero),
      m_apply_increase(threshold_charge > 0),
      m_threshold_charge(threshold_charge), 
      m_window_start_increase(window_start_increase),
      m_window_width_increase(window_width_increase)
    { /* nothing to see here */ }
    virtual ~VSMaxChargeCalc();
    virtual bool calc(double& signal, double& time,
		      bool lo_gain, double lg_multiplier, unsigned nsample,
		      const uint32_t* samples, const uint32_t* integrated, 
		      double ped);

    // This function presents an "extended" API that can return
    // additional information and allows for different window_zero to
    // be used (for example for different telescopes)

    bool calc(double& signal, double& time, 
	      unsigned& window_width, double& window_start,
	      double& total_signal, const uint32_t*& peak_ptr,
	      const bool lo_gain, const double lg_multiplier, 
	      unsigned nsample, const uint32_t*const samples,
	      const uint32_t*const integrated, const double ped,
	      const unsigned sample_zero);

  private:
    unsigned                m_window_width;
    unsigned                m_hi_sample_zero;
    unsigned                m_lo_sample_zero;
    bool                    m_apply_increase;
    double                  m_threshold_charge;
    double                  m_window_start_increase;
    double                  m_window_width_increase;
  };
#endif
}

#endif // VSTOFFSETCALC_HPP
