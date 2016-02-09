//-*-mode:c++; mode:font-lock;-*-

/*! \file VSChannelReMapping.hpp

  Base class for channel remapping class

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       09/12/2006
  \note
*/

#ifndef VSCHANNELREMAPPING_HPP
#define VSCHANNELREMAPPING_HPP

#include<vector>
#include<stdint.h>

//! VERITAS namespace
namespace VERITAS 
{

  //! Provide mapping of channels used in analysis to those used in
  //simulations. Analysis channels care denoted "A-channels", while
  //those in simulations are denoted "S-channels". The mapping
  //basically allows some simulation channels to be skipped, so that
  //the simulations can be run with a larger field of view. Some of
  //the channels in this large camera can then be skipped to get back
  //the physical camera.

  class VSChannelReMapping
  {
  public:
    virtual ~VSChannelReMapping();
    virtual unsigned numScopes() = 0;
    virtual unsigned numAnlChannels(unsigned iscope) = 0;
    virtual unsigned numSimChannels(unsigned iscope) = 0;
    virtual unsigned mapAnlToSim(unsigned iscope, unsigned ianalysis) = 0;
    virtual unsigned mapSimToAnl(unsigned iscope, unsigned isim) = 0;
    virtual bool isPresentInAnl(unsigned iscope, unsigned isim) = 0;
  };

};

#endif // defined VSCHANNELREMAPPING_HPP
