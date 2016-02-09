#ifndef VSMODELCOORD_HPP
#define VSMODELCOORD_HPP

#include <VSACoord.hpp>

namespace VERITAS
{
  class VSModelCoord : public VSACoord::Coord2D
  {
  public:
    VSModelCoord(): VSACoord::Coord2D() { }
    VSModelCoord(double x, double y, unsigned iobs): 
      VSACoord::Coord2D(x,y), m_iobs(iobs) { }
    VSModelCoord(const VSAAlgebra::Vec2D& xy, unsigned iobs): 
      VSACoord::Coord2D(xy), m_iobs(iobs) { }
    ~VSModelCoord() { }

    void setObs(unsigned iobs) { m_iobs = iobs; }

    unsigned iobs() const { return m_iobs; }

  private:

    unsigned m_iobs;

   };

  inline std::ostream& operator<<(std::ostream& o, const VSModelCoord& c)
  {
    o << std::setw(5) << c.iobs() 
      << std::setw(15) <<  c.x()
      << std::setw(15) <<  c.y();
    return o;
  }

}

#endif // VSMODELCOORD_HPP
