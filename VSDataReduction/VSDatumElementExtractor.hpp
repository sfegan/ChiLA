/*! \file VSDatumElementExtractor.hpp

  Extract the the value of an element in a datum.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       07/08/2007

*/

#ifndef VSDATUMELEMENTEXTRACTOR_HPP
#define VSDATUMELEMENTEXTRACTOR_HPP

#include<iostream>
#include<cmath>
#include<VSOctaveIO.hpp>
#include<VSDataConverter.hpp>
#include<VSH5DatumElement.hpp>
#include<VSEventData.hpp>

namespace VERITAS
{
  class VSEventDatumElement
  {
  public:
    VSEventDatumElement() {}
    virtual ~VSEventDatumElement() {}

    virtual bool getValue(const VSEventArrayDatum& datum, double& value) = 0;
  };

  class VSEventArrayDatumElement : public VSEventDatumElement
  {
  public:
    VSEventArrayDatumElement(const std::string& element_name) 
    {
      m_datum_element = 
	VSH5DatumElement<VSEventArrayDatum>::createDatumElement(element_name);
    }
    
    ~VSEventArrayDatumElement() 
    {
      delete m_datum_element;
    }

    virtual bool getValue(const VSEventArrayDatum& datum, double& value)
    {
      return m_datum_element->getValue(datum,value);
    }
    

  private:    
    VSH5DatumElement<VSEventArrayDatum>* m_datum_element;
  };

  class VSEventScopeDatumElement : public VSEventDatumElement
  {
  public:
    VSEventScopeDatumElement(const std::string& element_name):
      m_scope_id()
    {
      std::string scope_index = 
	VSH5DatumElementParser::getIndex(element_name);
      vsassert(!scope_index.empty());
      VSDatumConverter< unsigned >::fromString(m_scope_id,scope_index.c_str());
      m_datum_element = 
	VSH5DatumElement<VSEventScopeDatum>::createDatumElement(element_name);
    }

    ~VSEventScopeDatumElement() 
    {
      delete m_datum_element;
    }

    virtual bool getValue(const VSEventArrayDatum& datum, double& value)
    {
      return m_datum_element->getValue(*(datum.scope[m_scope_id]),value);
    }
    
  private:    
    unsigned m_scope_id;
    VSH5DatumElement<VSEventScopeDatum>* m_datum_element;
  };
}

#endif // VSDATUMELEMENTEXTRACTOR_HPP
