/*! \file VSNSpaceFilterDatum.hpp

  Apply NSpace multi-dimensional volume to a datum

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/08/2007

*/


#ifndef VSNSPACEFILTERDATUM_HPP
#define VSNSPACEFILTERDATUM_HPP

#include<VSNSpace.hpp>
#include<VSOctaveIO.hpp>
#include<VSH5DatumElement.hpp>

namespace VERITAS 
{
  template<typename T> class VSNSpaceFilterDatum
  {
  public:
    VSNSpaceFilterDatum(const VSNSpace::Volume& vol);
    ~VSNSpaceFilterDatum();
    bool isDatumInVolume(const T& x) const;
    bool isDatumInVolume(const T& x, unsigned index) const;
    const VSNSpace::Volume& getVolume() const;
  private:
    VSNSpace::Volume                  m_volume;
    std::vector<VSH5DatumElement<T>*> m_datums;
  };

  template<typename T>
  VSNSpaceFilterDatum<T>::VSNSpaceFilterDatum(const VSNSpace::Volume& vol):
    m_volume(vol), m_datums()
  {
    for(std::vector<VSNSpace::Axis>::const_iterator itr = 
	  m_volume.space().axes.begin();
	itr != m_volume.space().axes.end(); ++itr) 
      m_datums.push_back(VSH5DatumElement<T>::createDatumElement(itr->name));
  }

  template<typename T>
  VSNSpaceFilterDatum<T>::~VSNSpaceFilterDatum()
  {
    for(typename std::vector<VSH5DatumElement<T>*>::iterator itr = 
	  m_datums.begin(); itr != m_datums.end(); ++itr)
      delete (*itr);
  }

  template<typename T>
  bool VSNSpaceFilterDatum<T>::isDatumInVolume(const T& x) const
  {
    VERITAS::VSNSpace::Point point(m_datums.size());    
    for(unsigned i = 0; i < m_datums.size(); i++)
      if(!m_datums[i]->getValue(x,point.x[i])) return false;

    return m_volume.isPointInVolume(point);
  }

  template<typename T>
  bool VSNSpaceFilterDatum<T>::isDatumInVolume(const T& x, unsigned index) 
    const
  {
    VERITAS::VSNSpace::Point point(m_datums.size());    
    for(unsigned i = 0; i < m_datums.size(); i++)
      point.x[i] = m_datums[i]->getValue(x);
      
    return m_volume.isPointInVolume(point);
  }

  template<typename T>
  const VSNSpace::Volume& VSNSpaceFilterDatum<T>::getVolume() const
  {
    return m_volume;
  }
}

#endif // VSNSPACEFILTERDATUM_HPP
