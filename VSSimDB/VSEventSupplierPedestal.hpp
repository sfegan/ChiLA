//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierPedestal.hpp

  Generate a pedestal event.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       08/25/2006
  \note
*/

#ifndef VSEVENTSUPPLIERPEDESTAL_HPP
#define VSEVENTSUPPLIERPEDESTAL_HPP

#include<RandomNumbers.hpp>

#include<VSEventSupplier.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierPedestal: public VSEventSupplier
  {
  public:    
    VSEventSupplierPedestal(unsigned nevents,
			    unsigned nscope = 4, unsigned nchan = 500);
    virtual ~VSEventSupplierPedestal();
    virtual std::vector<SimParam> getSimParam();
    virtual bool getNextEvent(Event& e);
    virtual uint32_t getNumEvents();
  private:
    unsigned                       fNEvent;
    unsigned                       fNScope;
    unsigned                       fNChannel;
    unsigned                       fIEvent;
  };

}

#endif // defined VSEVENTSUPPLIERPEDESTAL_HPP
