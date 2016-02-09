//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutsCalc.hpp

  Base class for calculating whether an event passed a set of cuts.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/08/2007

  $Id: VSCutsCalc.hpp,v 3.7 2010/10/20 03:11:02 matthew Exp $

*/

#ifndef VSCUTSCALC_HPP
#define VSCUTSCALC_HPP

#include<VSCutsData.hpp>
#include<VSEventData.hpp>

namespace VERITAS
{
  class VSCutsCalc
  {
  public:
    VSCutsCalc();
    virtual ~VSCutsCalc();

    virtual void getCutResults(VSArrayCutsDatum& cut_results, 
			       const VSEventArrayDatum& event) = 0;    
    virtual void getScopeParamSet(std::set<std::string>& scope_set) = 0;
    virtual void getArrayParamSet(std::set<std::string>& array_set) = 0;
    virtual void print() = 0;
    virtual bool load(VSOctaveH5ReaderStruct* reader)
    {
      reader->readComposite("lo_date",m_lo_date);
      reader->readComposite("hi_date",m_hi_date);
      return true;
    }

    virtual bool save(VSOctaveH5WriterStruct* writer)
    {
      writer->writeComposite("lo_date",m_lo_date);
      writer->writeComposite("hi_date",m_hi_date);
      return true;
    }

    const VSTime& loDate() const { return m_lo_date; }
    const VSTime& hiDate() const { return m_hi_date; }

    void setDateRange(const std::string& lo_date,
		      const std::string& hi_date);
    void setDateRange(const VSTime& lo_date, const VSTime& hi_date);


    void setLoDate(const VSTime& lo_date) { m_lo_date = lo_date; }
    void setHiDate(const VSTime& hi_date) { m_hi_date = hi_date; }

  private:

    VSTime                                  m_lo_date;
    VSTime                                  m_hi_date;
  };
}

#endif // VSCUTSCALC_HPP
