/*! \file VSNSpaceEventData.hpp

  Class for filling event data into nspace points.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/08/2007

*/

#include<VSNSpace.hpp>
#include<VSNSpaceDataSource.hpp>
#include<VSOctaveH5Reader.hpp>
#include<VSDatumElementExtractor.hpp>
#include<VSScaledParameterCalc.hpp>

#ifndef VSNSPACEEVENTDATA_HPP
#define VSNSPACEEVENTDATA_HPP

namespace VERITAS 
{
  class VSNSpaceEventData: public VSNSpaceDataSource
  {
  public:
    VSNSpaceEventData(const std::string& filename,
		      const VSNSpace::Space& space,
		      VSScaledParameterCalc* sp_calc);

    virtual ~VSNSpaceEventData();
    virtual bool getData(std::vector< VSNSpace::Point >& points, 
			 double& theta0, double& theta1);
  private:
    VSOctaveH5Reader*                            m_reader;
    VSEventDataReader*                           m_event_reader;
    VSScaledParameterCalc*                       m_sp_calc;
    VSEventArrayDatum                            m_event_data;
    unsigned                                     m_nevents;
    unsigned                                     m_ievent;
    std::vector<VSH5DatumElement<VSEventArrayDatum>*> m_array_datums;
    std::vector<VSH5DatumElement<VSEventScopeDatum>*> m_scope_datums;
    std::vector< unsigned >                      m_scope_index;
    unsigned                                     m_ndatum;
    std::vector< bool >                          m_is_scope_datum;
    bool                                         m_has_scope_datum;
  };
}

#endif // VSNSPACEEVENTDATA_HPP
