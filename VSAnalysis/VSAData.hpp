//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAData.hpp

  Data classes.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2008
*/


#ifndef VSADATA_HPP
#define VSADATA_HPP

#include <vector>
#include <VSAMath.hpp>

namespace VERITAS
{
  namespace VSAMath
  {
    template< typename T = double >
    struct DataPoint
    {
      DataPoint(): x(), y(), sigma() { /* nothing to see here */ }
      DataPoint(T _x, double _y, double _s=1.0): 
	x(_x), y(_y), sigma(_s) { /* nothing to see here */ }

      T      x;
      double y;
      double sigma;
    };

    template< typename T = double >
    class Data
    {
    public:

      typedef typename std::vector< DataPoint<T> >::iterator iterator;
      typedef typename std::vector< DataPoint<T> >::const_iterator 
      const_iterator;

      Data(): m_data(), m_binSize() 
      { }

      // Accessors ------------------------------------------------------------
      unsigned size() const { return m_data.size(); }
      double binSize(unsigned idim) const 
      { 
	vsassert(idim < m_binSize.size());
	return m_binSize[idim]; 
      }

      const DataPoint<T>& operator[](unsigned i) const { return m_data[i]; }
      const std::vector< DataPoint<T> >& data() const { return m_data; }

      iterator begin() { return m_data.begin(); }
      iterator end() { return m_data.end(); }
      const_iterator begin() const { return m_data.begin(); }
      const_iterator end() const { return m_data.end(); }

      // Setters --------------------------------------------------------------
      void setBinSize(unsigned idim, double binsize)
      {
	m_binSize.resize(idim+1);
	m_binSize[idim] = binsize;
      }

      DataPoint<T>& operator[](unsigned i) { return m_data[i]; }
      std::vector< DataPoint<T> >& data() { return m_data; }

      void insert(const Data& data)
      {
	m_data.insert(m_data.end(),data.data().begin(),data.data().end());
      }

      void insert(const DataPoint<T>& data_point)
      {
	m_data.push_back(data_point);
      }

      void clear() { m_data.clear(); }

      iterator erase(iterator itr) { return m_data.erase(itr); }

    private:
      std::vector< DataPoint<T> > m_data;
      std::vector< double >       m_binSize;
    };
  }
}

#endif // VSADATA_HPP
