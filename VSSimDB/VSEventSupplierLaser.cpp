//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierLaser.cpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#include <cassert>
#include <vector>
#include <algorithm>

#include "VSEventSupplierLaser.hpp"

using namespace VERITAS;

VSEventSupplierLaser::
VSEventSupplierLaser(RandomNumbers* rng, unsigned nevents,
		     double mean_pe_per_pixel, 
		     double mean_pe_per_pixel_dev,
		     double mean_pulse_fwhm_ns,
		     unsigned nscope, unsigned nchan)
  : VSEventSupplier(), fRNG(rng), fNEvent(nevents),
    fMeanPEPerPixel(mean_pe_per_pixel), 
    fMeanPEPerPixelDev(mean_pe_per_pixel_dev),
    fMeanPulseFWHM(mean_pulse_fwhm_ns),
    fNScope(nscope), fNChannel(nchan), fIEvent()
{
  // nothing to see here
}

VSEventSupplierLaser::~VSEventSupplierLaser()
{
  // nothing to see here
}

std::vector<VSEventSupplier::SimParam> VSEventSupplierLaser::getSimParam()
{
  std::vector<SimParam> tp_vec;
  return tp_vec;
}

bool VSEventSupplierLaser::getNextEvent(Event& e)
{
  static const double fwhm = sqrt(log(2.0)*2.0);

  if((fNEvent >0)&&(fIEvent >= fNEvent))return false;
  
  e.fEventID                 = fIEvent+1;
  e.fTargetZenithRad         = 0;
  e.fTargetAzimuthRad        = 0;
  e.fPrimarySpeciesCORSIKA   = 0;
  e.fPrimaryEnergyTeV        = 0;
  e.fPrimaryZenithRad        = 0;
  e.fPrimaryAzimuthRad       = 0;
  e.fPrimaryCoreEastM        = 0;
  e.fPrimaryCoreNorthM       = 0;
  e.fPrimaryCoreUpASLM       = 0;
  e.fSamplingRadiusM         = 0;
  e.fTableIndex              = 0;
  e.fEventComplete           = 0;

  fIEvent++;

  e.fScopes.resize(fNScope);
  for(unsigned iscope=0; iscope<fNScope; iscope++)
    {
      e.fScopes[iscope].fScopeID         = iscope;
      e.fScopes[iscope].fScopeZenithRad  = 0;
      e.fScopes[iscope].fScopeAzimuthRad = 0;
      e.fScopes[iscope].fPixels.clear();
      e.fScopes[iscope].fPixels.reserve(fNChannel);

      // ----------------------------------------------------------------------
      // Randomize the amplitude of each laser flash
      // ----------------------------------------------------------------------
      double mean_pe = fMeanPEPerPixel + fMeanPEPerPixelDev*fRNG->Normal();
      while(mean_pe < 0)
	mean_pe = fMeanPEPerPixel + fMeanPEPerPixelDev*fRNG->Normal();

      unsigned islot=0;
      for(unsigned ichan=0; ichan<fNChannel; ichan++)
	{
	  unsigned npe = unsigned(fRNG->Poisson(mean_pe));
	  if(npe)
	    {
	      std::vector<double> pe(npe);
	      for(unsigned ipe=0;ipe<npe;ipe++)
		pe[ipe]=fRNG->Uniform()*fMeanPulseFWHM*fwhm;
	      std::sort(pe.begin(), pe.end());
	      e.fScopes[iscope].fPixels.resize(islot+1);
	      e.fScopes[iscope].fPixels[islot].fPixelID = ichan;
	      e.fScopes[iscope].fPixels[islot].fPEs.resize(npe);
	      for(unsigned ipe=0;ipe<npe;ipe++)
		e.fScopes[iscope].fPixels[islot].fPEs[ipe].fTimeNS = pe[ipe];
	      islot++;
	    }
	}
    }
  
  return true;
}

uint32_t VSEventSupplierLaser::getNumEvents()
{
  return fNEvent;
}
