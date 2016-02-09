//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierLaser.hpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#ifndef VSEVENTSUPPLIERLASER_HPP
#define VSEVENTSUPPLIERLASER_HPP

#include<RandomNumbers.hpp>

#include<VSEventSupplier.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierLaser: public VSEventSupplier
  {
  public:    
    VSEventSupplierLaser(RandomNumbers* rng, unsigned nevents,
			 double mean_pe_per_pixel = 100.0,
			 double mean_pe_per_pixel_dev = 10.0,
			 double mean_pulse_fwhm_ns = 5.0,
			 unsigned nscope = 4, unsigned nchan = 500);
    virtual ~VSEventSupplierLaser();
    virtual std::vector<SimParam> getSimParam();
    virtual bool getNextEvent(Event& e);
    virtual uint32_t getNumEvents();
  private:
    RandomNumbers*                 fRNG;
    unsigned                       fNEvent;
    double                         fMeanPEPerPixel;
    double                         fMeanPEPerPixelDev;
    double                         fMeanPulseFWHM;
    unsigned                       fNScope;
    unsigned                       fNChannel;
    unsigned                       fIEvent;
  };

}

#endif // defined VSEVENTSUPPLIERLASER_HPP
