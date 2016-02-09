//-*-mode:c++; mode:font-lock;-*-

/*! \file VSACoord.hpp
  Coordinate classes.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSACOORD_HPP
#define VSACOORD_HPP

#include <typeinfo>

#include <VSAAlgebra.hpp>
#include <VSAData.hpp>

namespace VERITAS
{
  namespace VSACoord
  {
    // ------------------------------------------------------------------------
    // CoordND
    // ------------------------------------------------------------------------
    class CoordND
    {
    public:
      CoordND(): m_x() { }
      CoordND(unsigned n): m_x(n) { }
      CoordND(const std::vector<double>& x): m_x(x) { }
      CoordND(unsigned n, const double* x): m_x(n)
      {
	std::copy(x, x+n, m_x.begin());
      }

      virtual ~CoordND() { }

      // Accessors ------------------------------------------------------------
      unsigned size() const { return m_x.size(); }

      virtual const double& operator[] (unsigned idim) const 
      { return m_x[idim]; }

      // Setters --------------------------------------------------------------
      virtual void set(const std::vector<double>& x)
      {
	m_x = x;
      }

      virtual double& operator[] (unsigned idim) 
      { return m_x[idim]; }

    protected:
      std::vector<double> m_x;
    };    

    inline std::ostream& operator<<(std::ostream& o, const CoordND& c)
    {
      for(unsigned i = 0; i < c.size(); i++) o << std::setw(15) <<  c[i];
      return o;
    }

    // ------------------------------------------------------------------------
    // Coord1D
    // ------------------------------------------------------------------------
    class Coord1D
    {
    public:
      Coord1D(): m_x() { }
      Coord1D(double x): m_x(x) 
      { 

      }

      virtual ~Coord1D() { }

      // Accessors ------------------------------------------------------------
      unsigned size() const { return 1; }

      virtual double x() const { return m_x; }
      virtual double jacobian() const { return 1; }
      virtual double jacobian(unsigned idim) const { return 1; }
      virtual const double& operator[] (unsigned idim) const 
      { return m_x; }

      // Setters --------------------------------------------------------------
      virtual void set(const std::vector<double>& x)
      {
	m_x = x[0];
      }

      // Static Methods -------------------------------------------------------
      static unsigned ndim() { return 1; }

    protected:
      double m_x;
    };    

    inline std::ostream& operator<<(std::ostream& o, const Coord1D& c)
    {
      o << std::setw(15) << c.x();      
      return o;
    }

    // ------------------------------------------------------------------------
    // Coord2D
    // ------------------------------------------------------------------------
    class Coord2D
    {
    public:
      Coord2D(): m_x(2,0), m_polar(2,0) { }
      Coord2D(double x, double y): m_x(2,0), m_polar(2,0) 
      { 
	setCartesian(x,y);
      }

      Coord2D(const VSAAlgebra::Vec2D& xy): m_x(2,0), m_polar(2,0)
      {
	setCartesian(xy.x(),xy.y());
      }

      virtual ~Coord2D() { }

      // Accessors ------------------------------------------------------------
      unsigned size() const { return 2; }

      virtual double x() const { return m_x[0]; }
      virtual double y() const { return m_x[1]; }
      virtual double r() const { return m_polar[0]; }
      virtual double phi() const { return m_polar[1]; }

      const double& polar(unsigned idim) const { return m_polar[idim]; }
      const std::vector<double>& polar() const { return m_polar; }

      virtual double jacobian() const { return 1; }
      virtual double jacobian(unsigned idim) const { return 1; }
      virtual const double& operator[] (unsigned idim) const 
      { return m_x[idim]; }

      // Setters --------------------------------------------------------------
      void setPolar(double r, double phi)
      {
	m_polar[0] = r;
	m_polar[1] = phi;
	m_x[0] = r*cos(phi);
	m_x[1] = r*sin(phi);
      }

      void setCartesian(double x, double y)
      {
	m_x[0] = x;
	m_x[1] = y;
	m_polar[0] = sqrt(std::pow(m_x[0],2) + std::pow(m_x[1],2));
	m_polar[1] = atan2(y,x);
      }

      double& polar(unsigned idim) { return m_polar[idim]; }


      virtual void set(const std::vector<double>& x)
      {
	setCartesian(x[0],x[1]);
      }

      virtual void set(unsigned idim, double x)
      {
	m_x[idim] = x;
	setCartesian(m_x[0],m_x[1]);
      }
            
      // Static Methods -------------------------------------------------------
      static unsigned ndim() { return 2; }

      static Coord2D makeCartesian(double x, double y)
      {
	return Coord2D(x,y);
      }

      static Coord2D makePolar(double r, double phi)
      {
	Coord2D c;
	c.setPolar(r,phi);
	return c;
      }

    protected:
      std::vector< double > m_x;
      std::vector< double > m_polar;
    };    

    inline std::ostream& operator<<(std::ostream& o, const Coord2D& c)
    {
      o << std::setw(15) << c.x()
	<< std::setw(15) << c.y();
      return o;
    }

    template< typename T=Coord2D >
    class Polar : public T
    {
    public:
      Polar(const VSAAlgebra::Vec2D& xy, const T& rphi):
	T(rphi), m_xy(xy), m_r(), m_phi() 
      { 
	double x = m_xy.x()+T::x();
	double y = m_xy.y()+T::y();
	m_phi = atan2(y,x);
	m_r = sqrt(std::pow(x,2) + std::pow(y,2));	
      }

      Polar(const VSAAlgebra::Vec2D& xy, double r, double phi):
	T(), m_xy(xy), m_r(), m_phi() 
      { 
	T::setPolar(r,phi);
	double x = m_xy.x()+T::x();
	double y = m_xy.y()+T::y();
	m_phi = atan2(y,x);
	m_r = sqrt(std::pow(x,2) + std::pow(y,2));	
      }

      Polar(double r, double phi):
	T(), m_xy(), m_r(r), m_phi(phi)
      { 
	T::setPolar(r,phi);
      }

      virtual double x() const
      { 
	return m_r*cos(m_phi);
      }

      virtual double y() const
      { 
	return m_r*sin(m_phi);
      }

      virtual double r() const
      { 
	return m_r;
      }

      virtual double phi() const
      { 
	return m_phi;
      }

      virtual double jacobian() const { return T::r(); }
      virtual double jacobian(unsigned idim) const 
      { 
	if(idim == 0) return T::r();
	else return 1.;
      }

      virtual const double& operator[] (unsigned idim) const 
      { return T::polar(idim); }

      virtual void set(const std::vector<double>& _x)
      {
	T::setPolar(_x[0],_x[1]);
	double x = m_xy.x()+T::x();
	double y = m_xy.y()+T::y();
	m_phi = atan2(y,x);
	m_r = sqrt(std::pow(x,2) + std::pow(y,2));
      }

      virtual void set(unsigned idim, double x)
      {
	std::vector< double > xp = T::polar();
	xp[idim] = x;
	set(xp);
      }

    private:
      VSAAlgebra::Vec2D m_xy;
      double            m_r;
      double            m_phi;
    };

    template< typename T >
    inline std::ostream& operator<<(std::ostream& o, const Polar<T>& c)
    {
      o << std::setw(20) << c.x()
	<< std::setw(20) << c.y()
	<< std::setw(20) << c.r()
	<< std::setw(20) << c.phi();
	
      return o;
    }
  }
}

#endif // VSACOORD_HPP
