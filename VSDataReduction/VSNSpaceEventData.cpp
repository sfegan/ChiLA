/*! \file VSNSpaceEventData.cpp

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

#include "VSNSpaceDataSource.hpp"
#include "VSNSpaceOctaveH5IO.hpp"
#include "VSNSpaceEventData.hpp"

using namespace VERITAS;

VSNSpaceEventData::
VSNSpaceEventData(const std::string& filename, const VSNSpace::Space& space,
		  VSScaledParameterCalc* sp_calc):
  VSNSpaceDataSource(), 
  m_reader(),
  m_sp_calc(sp_calc),
  m_nevents(),
  m_ievent(),
  m_ndatum(),
  m_is_scope_datum(),
  m_has_scope_datum()
{
  m_reader = new VSOctaveH5Reader(filename);
  vsassert(m_reader);
  VSOctaveH5ReaderStruct* event_struct = m_reader->readStruct("events");
  vsassert(event_struct);
  m_event_reader = new VSEventDataReader(event_struct);
  vsassert(m_event_reader);

  m_nevents = m_event_reader->rows();

  const unsigned ndim = space.ndim;
  for(unsigned idim = 0; idim < ndim; idim++)
    {
      if(VSH5DatumElement<VSEventArrayDatum>::
	 hasElement(space.axes[idim].name))
	{
	  VSH5DatumElement<VSEventArrayDatum>* datum = 
	    VSH5DatumElement<VSEventArrayDatum>::
	    createDatumElement(space.axes[idim].name);

	  m_array_datums.push_back(datum);
	  m_is_scope_datum.push_back(false);
	}
      else if(VSH5DatumElement<VSEventScopeDatum>::
	      hasElement(space.axes[idim].name))
	{
	  VSH5DatumElement<VSEventScopeDatum>* datum = 
	    VSH5DatumElement<VSEventScopeDatum>::
	    createDatumElement(space.axes[idim].name);

	  m_scope_datums.push_back(datum);
	  std::string index = VSH5DatumElementParser::
	    getIndex(space.axes[idim].name);

	  m_is_scope_datum.push_back(true);
	  m_has_scope_datum = true;
	}
      else
	{
	  std::cerr << "Error parsing datum name: " 
		    << space.axes[idim].name << std::endl;
	  exit(EXIT_FAILURE);
	}
     
      m_ndatum++;
    }
}

VSNSpaceEventData::~VSNSpaceEventData()
{
  delete m_event_reader;
  delete m_reader;
}

bool VSNSpaceEventData::getData(std::vector<VSNSpace::Point>& points, 
				double& theta0, double& theta1)
{
  if(m_ievent >= m_nevents)
    return false;

  m_event_reader->element(m_event_data, m_ievent);
  points.clear();

  if(m_sp_calc)
    m_sp_calc->calcSP(m_event_data);

  theta0 = m_event_data.theta0;
  theta1 = m_event_data.theta1;

  unsigned np;
  if(m_has_scope_datum) np = m_event_data.scope.size();
  else np = 1;

  for(unsigned ipoint = 0; ipoint < np; ipoint++)
    {
      if(!m_event_data.scope[ipoint])
	continue;

      VSNSpace::Point p(m_ndatum);

      unsigned iscope = 0;
      unsigned iarray = 0;

      bool is_valid = true;

      for(unsigned idim = 0; idim < m_ndatum; idim++)
 	{
 	  if(m_is_scope_datum[idim])
	    {
	      is_valid =
		m_scope_datums[iscope]->getValue(*m_event_data.scope[ipoint],
						 p.x[idim]);
	      iscope++;
	    }
	  else
	    {
	      is_valid = 
		m_array_datums[iarray]->getValue(m_event_data,p.x[idim]);

	      iarray++;
	    }

	  if(!is_valid) break;
 	}

      if(is_valid) points.push_back(p);
    }

  m_ievent++;
  return true;
}
