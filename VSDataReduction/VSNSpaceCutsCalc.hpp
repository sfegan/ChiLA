//-*-mode:c++; mode:font-lock;-*-

/*! \file VSNSpaceCutsCalc.hpp

  Calculate whether event is in the array and scope filter volume

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/09/2007

  $Id: VSNSpaceCutsCalc.hpp,v 3.10 2008/11/24 02:02:02 matthew Exp $

*/

#ifndef VSNSPACECUTSCALC_HPP
#define VSNSPACECUTSCALC_HPP

#include<VSNSpaceFilterDatum.hpp>
#include<VSCutsData.hpp>
#include<VSCutsCalc.hpp>
#include<VSEventData.hpp>
#include<VSNSpace.hpp>
#include<VSOptions.hpp>

namespace VERITAS
{
  template<> class VSNSpaceFilterDatum<VSEventArrayDatum>
  {
  public:
    VSNSpaceFilterDatum(const VSNSpace::Volume& vol);
    ~VSNSpaceFilterDatum();
    bool isDatumInVolume(const VSEventArrayDatum& x) const;
    bool isDatumInVolume(const VSEventArrayDatum& x, unsigned index) const;
    const VSNSpace::Volume& getVolume() const;
  private:
    VSNSpace::Volume                                  m_volume;
    std::vector<VSH5DatumElement<VSEventArrayDatum>*> m_array_datums;
    std::vector<VSH5DatumElement<VSEventScopeDatum>*> m_scope_datums;
    std::vector<bool>                                 m_is_scope_datum;
    unsigned                                          m_ndatum;
  };

  class VSNSpaceCutsCalc : public VSCutsCalc
  {
  public:
    class Options
    {
    public:
      Options(): 
	nspace_cuts_file(), nscope(2)
      {}
      
      std::string           nspace_cuts_file;  
      unsigned              nscope;
    };

    VSNSpaceCutsCalc(const Options& opt = s_default_options);
    VSNSpaceCutsCalc(const VSNSpace::Volume* array_vol,
		     const VSNSpace::Volume* scope_vol_default,
		     const std::vector< const VSNSpace::Volume* >& scope_vol);

    ~VSNSpaceCutsCalc();

    void getCutResults(VSArrayCutsDatum& cut_results, 
		       const VSEventArrayDatum& event);
	
    virtual void print() {}

    virtual void clear();
    virtual bool load(VSOctaveH5ReaderStruct* reader);
    virtual bool save(VSOctaveH5WriterStruct* writer);

    virtual void getScopeParamSet(std::set< std::string >& scope_set)
    {

    }

    virtual void getArrayParamSet(std::set< std::string >& array_set)
    {

    }

    static void configure(VSOptions& options,
			  const std::string& profile = "", 
			  const std::string& opt_prefix = "");

  private:
    typedef VSNSpaceFilterDatum<VSEventArrayDatum> ArrayFilter;
    typedef VSNSpaceFilterDatum<VSEventScopeDatum> ScopeFilter;

    Options                   m_options;

    ArrayFilter*              m_array_filter;
    ArrayFilter*              m_scope_filter_default;
    std::vector<ScopeFilter*> m_scope_filter;

    static Options            s_default_options;
  };
}

#endif // defined VSNSPACECUTSCALC_HPP
