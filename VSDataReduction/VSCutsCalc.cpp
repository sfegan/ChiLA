//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutsCalc.cpp

  Base class for applying cuts.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/09/2007

  $Id: VSCutsCalc.cpp,v 3.3 2010/10/20 03:11:02 matthew Exp $

*/

#include<VSCutsCalc.hpp>

using namespace VERITAS;

VSCutsCalc::VSCutsCalc():
  m_lo_date(VSTime::perpetual_past()),
  m_hi_date(VSTime::perpetual_future())
{

}

VSCutsCalc::~VSCutsCalc()
{

}

void VSCutsCalc::setDateRange(const std::string& lo_date,
			      const std::string& hi_date)
{
  if(lo_date != "*" && lo_date != "-")
    {
      uint64_t lo_date_uint64 = 0;
      VSDatumConverter< uint64_t >::fromString(lo_date_uint64,lo_date.c_str());
      lo_date_uint64*=1000000;
      m_lo_date.setFromDBTimeStamp(lo_date_uint64);
    }      
  
  if(hi_date != "*" && hi_date != "-")
    {
      uint64_t hi_date_uint64 = 0;
      VSDatumConverter< uint64_t >::fromString(hi_date_uint64,hi_date.c_str());
      hi_date_uint64*=1000000;
      m_hi_date.setFromDBTimeStamp(hi_date_uint64);
    }   
}

void VSCutsCalc::setDateRange(const VSTime& lo_date,
			      const VSTime& hi_date)
{
  m_lo_date = lo_date;
  m_hi_date = hi_date;
}
