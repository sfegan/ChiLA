/*! \file VSH5DatumElement.hpp

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

#ifndef VSH5DATUMELEMENT_HPP
#define VSH5DATUMELEMENT_HPP

#include<iostream>
#include<cmath>
#include<VSOctaveIO.hpp>
#include<VSDataConverter.hpp>
#include<VSEventData.hpp>

namespace VERITAS
{
  class VSH5DatumElementParser
  {
  public:
    static std::string getTransform(const std::string& name)
    {
      std::string transform_name = "";
      std::string::size_type pos = name.find_first_of(":");

      if(pos != std::string::npos)
	transform_name = name.substr(0,pos);

      return transform_name;
    }

    static std::string getName(const std::string& name)
    {
      std::string tmp = "";
      std::string::size_type pos = name.find_first_of(":");

      if(pos != std::string::npos)
	tmp = name.substr(pos+1,name.size()-pos-1);
      else
	tmp = name;

      return tmp;
    }

    static std::string getIndex(const std::string& name)
    {
      std::string index = "";
      
      std::string::size_type pos1 = name.find_first_of("{");
      std::string::size_type pos2 = name.find_first_of("}");
      
      if(pos1 != std::string::npos && pos2 != std::string::npos)
	index = name.substr(pos1+1,pos2-pos1-1);
      
      return index;
    }
    
    static std::string getElement(const std::string& name)
    {
      std::string element_name = "";
      std::string::size_type pos = name.find_last_of(".:");
      
      if(pos != std::string::npos)
	element_name = name.substr(pos+1,name.size()-pos-1);
      else
	element_name = name;

      return element_name;
    }

    static std::string getOctaveName(const std::string& name)
    {
      std::string octave_name = "";

      if(getTransform(name) != "")
	octave_name += getTransform(name) + ":";
	
      octave_name += getElement(name);

      //       if(getIndex(name) != "")
      // 	octave_name += "_" + getIndex(name);
      
      return octave_name;
    }
  };

  template<typename T>
  class VSH5DatumElement
  {
  public:
    VSH5DatumElement() {}
    virtual ~VSH5DatumElement() {}
    virtual VSH5DatumElement* clone() = 0;
    virtual bool getValue(const T& datum, double& value) = 0;
    virtual std::string getName() = 0;
    static VSH5DatumElement* 
    createDatumElement(const std::string& datum_element_name);

    static bool hasElement(const std::string& datum_element_name)
    {	  
      std::string name = VSH5DatumElementParser::getName(datum_element_name);

      VSOctaveH5CompositeDefinition c;
      VSOctaveH5CompositeCompose<T>::compose(c);
      std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
      std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator im;
      for(im = m.begin(); im!=m.end(); im++) 
	if(im->name == VSH5DatumElementParser::getElement(name))
	  return true;

      return false;
    }
  };
  
  template<typename T, typename Type>
  class VSH5DatumElementBasicType : public VSH5DatumElement<T>
  {      
  public:
    typedef bool (*TransformPtr) (double&);

    VSH5DatumElementBasicType(size_t offset, 
			      const std::string& element_name,
			      const std::string& transform_name);
    virtual ~VSH5DatumElementBasicType() {}
    virtual VSH5DatumElementBasicType* clone() 
    {
      return new VSH5DatumElementBasicType(m_offset,
					   m_element_name,
					   m_transform_name);
    }

    virtual bool getValue(const T& datum, double& value) 
    { 
      value = (double)*(reinterpret_cast<Type*>((char*)&datum + m_offset));
      return m_transform_ptr(value);
    }
    
    virtual std::string getName()
    {
      return m_element_name;
    }

    static bool transformNone(double& x) 
    { 
      if(std::isfinite(x))
	return true;
      else
	return false;
    }

    static bool transformSq(double& x) 
    { 
      if(std::isfinite(x))
	{
	  x = x*x;
	  return true;
	}
      else
	return false;
    }

    static bool transformLog(double& x) 
    { 
      if(x > 0 && std::isfinite(x))
	{
	  x = log10(x); 
	  return true;
	}
      else
	return false;
    } 

    static bool transformNegLog(double& x) 
    { 
      if(x > 0 && std::isfinite(x))
	{
	  x = -log10(x); 
	  return true;
	}
      else
	return false;
    } 

  private:
    size_t       m_offset;
    std::string  m_element_name;
    std::string  m_transform_name;
    TransformPtr m_transform_ptr;
  };

  template<typename T>
  VSH5DatumElement<T>*  VSH5DatumElement<T>::
  createDatumElement(const std::string& datum_element_name)
  {
    std::string transform_name = 
      VSH5DatumElementParser::getTransform(datum_element_name);
    std::string element_name = 
      VSH5DatumElementParser::getElement(datum_element_name);

    VSH5DatumElement* datum_element = NULL;
    VSOctaveH5CompositeDefinition c;
    VSOctaveH5CompositeCompose<T>::compose(c);
    const std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
    std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator im;
    for(im = m.begin(); im!=m.end(); im++) 
      {
	if(im->name == element_name)
	  break;
      }

    if(im == m.end())
      {
	std::cerr << "Unrecognized datum_element name: " 
		  << element_name << std::endl;
	exit(EXIT_FAILURE);
      }
    else if(im->type == VSOctaveH5Type<double>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,double>(im->offset,
						element_name,
						transform_name);
    else if(im->type == VSOctaveH5Type<bool>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,bool>(im->offset,
					      element_name,
					      transform_name);
    else if(im->type == VSOctaveH5Type<int8_t>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,int8_t>(im->offset,
						element_name,
						transform_name);
    else if(im->type == VSOctaveH5Type<int16_t>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,int16_t>(im->offset,
						 element_name,
						 transform_name);
    else if(im->type == VSOctaveH5Type<int32_t>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,int32_t>(im->offset,
						 element_name,
						 transform_name);
    else if(im->type == VSOctaveH5Type<uint8_t>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,uint8_t>(im->offset,
						 element_name,
						 transform_name);
    else if(im->type == VSOctaveH5Type<uint16_t>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,uint16_t>(im->offset,
						  element_name,
						  transform_name);
    else if(im->type == VSOctaveH5Type<uint32_t>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,uint32_t>(im->offset,
						  element_name,
						  transform_name);
    else if(im->type == VSOctaveH5Type<float>::int_type())
      datum_element = 
	new VSH5DatumElementBasicType<T,float>(im->offset,
					       element_name,
					       transform_name);
    else 
      {
	std::cerr << "Unrecognized datum_element type: " 
		  << element_name << std::endl;
	exit(EXIT_FAILURE);
      }

    return datum_element;
  }

  template<typename T,typename Type>
  VSH5DatumElementBasicType<T,Type>::
  VSH5DatumElementBasicType(size_t offset, 
			    const std::string& element_name,
			    const std::string& transform_name):
    VSH5DatumElement<T>(), m_offset(offset), 
    m_element_name(element_name), m_transform_name(transform_name) 
  {
    if(transform_name == "")
      m_transform_ptr = VSH5DatumElementBasicType<T,Type>::transformNone;
    else if(transform_name == "log")
      m_transform_ptr = VSH5DatumElementBasicType<T,Type>::transformLog;
    else if(transform_name == "neglog")
      m_transform_ptr = VSH5DatumElementBasicType<T,Type>::transformNegLog;
    else if(transform_name == "sq")
      m_transform_ptr = VSH5DatumElementBasicType<T,Type>::transformSq;
    else
      {
	std::cerr << "Unrecognized transform name: " 
		  << transform_name << std::endl;
	exit(EXIT_FAILURE);
      }
  }

}

#endif // VSH5DATUMELEMENT_HPP
