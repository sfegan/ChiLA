//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimplePedData.hpp
  Simple pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VSSimplePedData.hpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#ifndef VSSIMPLEPEDDATA_HPP
#define VSSIMPLEPEDDATA_HPP

#include<string>
#include<vector>

namespace VERITAS
{

  // ==========================================================================
  // SIMPLE PEDESTAL DATA
  // ==========================================================================

  class VSSimplePedData
  {
  public:

    struct ChanData
    {
      ChanData(): ped(), dev(), suppressed() { }
      double ped;
      std::vector<double> dev;
      bool suppressed;
    };

    typedef std::vector<ChanData> ScopeData;
  
    std::vector<ScopeData> m_data;

    VSSimplePedData(): m_data() { }
    bool load(const std::string& filename);
    bool save(const std::string& filename = "");
    void suppress(const double lo, const double hi, const unsigned window=0);
  };

};

#endif // VSSIMPLEPEDDATA_HPP
