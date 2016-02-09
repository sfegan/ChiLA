//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLTLibraryCalc.hpp

  Base class for lookup table library calculators.

  \author     Matthew Wood                \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/09/2010

  $Id: VSLTLibraryCalc.hpp,v 3.1 2010/06/22 00:00:15 matthew Exp $

*/

#ifndef VSLTLIBRARYCALC_HPP
#define VSLTLIBRARYCALC_HPP

#include<VSTime.hpp>

namespace VERITAS
{

  class VSLTLibraryCalc
  {
  public:
    VSLTLibraryCalc(): m_has_lt(false) { }
    virtual ~VSLTLibraryCalc() { }

    virtual void clear() = 0;
    
    virtual bool load(const std::string& sp_lookup_file,
		      const std::vector<unsigned>& nchan,
		      double zn_deg, double az_deg, double ped_rms) = 0;

    bool load(const std::string& file,
	      const VSTime& date,
	      const std::vector<unsigned>& nchan, double zn_deg, 
	      double az_deg, double ped_rms);

    

  protected:
    bool                           m_has_lt;

  };

}

#endif // VSLTLIBRARYCALC_HPP
