#ifndef VSDATAMODEL_HPP
#define VSDATAMODEL_HPP

#include <VSACoord.hpp>
#include <VSOctaveIO.hpp>
#include <VSAFunction.hpp>

#include <VSSourceModel.hpp>
#include <VSBkgndModel.hpp>
#include <VSModelCoord.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSOffRegion
  // ==========================================================================
  class VSOffRegion
  {
  public:
    static void getRegions(const VSAAlgebra::Vec2D& src_xy,
			   const VSAAlgebra::Vec2D& obs_xy,
			   double theta, unsigned max_nregion,
			   std::vector< VSAAlgebra::Vec2D >& off_coords);
  };

  // ==========================================================================
  // VSExclusionRegion
  // ==========================================================================
  class VSExclusionRegion
  {
  public:

    class Region
    {
    public:

      Region():
	m_name(), m_type(), m_xy(), m_ra_J2000_rad(), m_dec_J2000_rad(),
	m_radius_deg(), m_vmag()
      { }

      Region(const VSAAlgebra::Vec2D& xy, 
	     double radius_deg, 
	     const std::string type = "star",
	     double vmag = 0, const std::string& name = ""):
	m_name(name), m_type(type), m_xy(xy),
	m_ra_J2000_rad(),
	m_dec_J2000_rad(),
	m_radius_deg(radius_deg), m_vmag(vmag)
      { }

      Region(const VSAAlgebra::Vec2D& xy, 
	     const SEphem::SphericalCoords& radec_j2000,
	     double radius_deg, 
	     const std::string type = "star",
	     double vmag = 0, const std::string& name = ""):
	m_name(name), m_type(type), m_xy(xy),
	m_ra_J2000_rad(radec_j2000.longitude().radPM()), 
	m_dec_J2000_rad(radec_j2000.latitude().radPM()),
	m_radius_deg(radius_deg), m_vmag(vmag)
      { }

      bool isExcluded(const VSAAlgebra::Vec2D& xy) const
      {
	if(m_xy.d(xy) < m_radius_deg) return true;
	else return false;
      }
      
      double radius() const { return m_radius_deg; }
      const VSAAlgebra::Vec2D& xy() const { return m_xy; }

      void save(VSOctaveH5WriterStruct* writer) const
      {
	writer->writeCompositeHere(*this);
	m_xy.save(writer,"xy");
      }

      bool load(VSOctaveH5ReaderStruct* reader) 
      {
	reader->readCompositeHere(*this);
	m_xy.load(reader,"xy");
	return true;
      }

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDNAMEDMEMBER(c,Region,m_name,"name");
	H5_ADDNAMEDMEMBER(c,Region,m_type,"type");
	H5_ADDNAMEDMEMBER(c,Region,m_ra_J2000_rad,"ra_J2000_rad");
	H5_ADDNAMEDMEMBER(c,Region,m_dec_J2000_rad,"dec_J2000_rad");
	H5_ADDNAMEDMEMBER(c,Region,m_radius_deg,"radius_deg");
	H5_ADDNAMEDMEMBER(c,Region,m_vmag,"vmag");
      }

    private:
      std::string       m_name;
      std::string       m_type;
      VSAAlgebra::Vec2D m_xy;
      double            m_ra_J2000_rad;
      double            m_dec_J2000_rad;
      double            m_radius_deg;
      double            m_vmag;
    };

    VSExclusionRegion(): m_regions() { }

    typedef std::vector< Region >::iterator iterator;

//     void add(const std::pair<double,double>& xy, double radius_deg)
//     {
//       m_regions.
// 	push_back(Region(VSAAlgebra::Vec2D(xy.first,xy.second),radius_deg));
//     }

//     void add(const VSAAlgebra::Vec2D& coord, double radius_deg)
//     {
//       m_regions.push_back(Region(coord,radius_deg));
//     }

    void addStar(double x, double y, 
		 const SEphem::SphericalCoords& radec_j2000,
		 double radius_deg, double vmag)
    {
      Region r(VSAAlgebra::Vec2D(x,y), radec_j2000, radius_deg,"star",vmag);
      m_regions.push_back(r);
    }

    void addSource(double x, double y, 
		   const SEphem::SphericalCoords& radec_j2000,
		   double radius_deg)
    {
      Region r(VSAAlgebra::Vec2D(x,y), radec_j2000, radius_deg,"source",0.);
      m_regions.push_back(r);
    }

    void clear()
    {
      m_regions.clear();
    }

    bool isExcluded(const VSAAlgebra::Vec2D& xy)
    {
      for(std::vector< Region >::iterator itr = m_regions.begin(); itr != 
	    m_regions.end(); ++itr)
	if(itr->isExcluded(xy)) return true;

      return false;
    }

    iterator begin() { return m_regions.begin(); }
    iterator end() { return m_regions.end(); }

    std::vector< Region >& regions() { return m_regions; }

    void save(VSOctaveH5WriterStruct* writer) const;
    bool load(VSOctaveH5ReaderStruct* reader);

  private:

    std::vector< Region >   m_regions;
  };

  // ==========================================================================
  // VSDataModel
  // ==========================================================================
  class VSDataModel : public VSAFunction::CompositeSumFn<VSModelCoord>
  {
  public:
    VSDataModel(VSSourceModel* source_fn, 
		VSAcceptanceModel* bkgnd_fn);
    ~VSDataModel();

    // Accessors --------------------------------------------------------------
    using VSAFunction::CompositeSumFn<VSModelCoord>::val;
    double val(const VSAAlgebra::Vec2D& xy) const 
    {
      return m_source_fn->val(xy) + m_bkgnd_fn->val(xy);
    }

    VSSourceModel* getSourceModel() 
    { 
      return m_source_fn; 
    }

    VSAcceptanceModel* getBkgndModel() 
    { 
      return m_bkgnd_fn; 
    }

    unsigned nBkgndParm() const { return m_bkgnd_fn->nparm(); }
    unsigned nSourceParm() const { return m_source_fn->nparm(); }

    VSAAlgebra::VecND bkgndParam() const { return m_bkgnd_fn->param(); }
    VSAAlgebra::VecND sourceParam() const { return m_source_fn->param(); }

    
    double integrate(const VSAAlgebra::Vec2D& xy, double R1, double R2) const
    {
      return m_source_fn->integrate(xy,R1,R2) + 
	m_bkgnd_fn->integrate(xy,R1,R2);
    }

    // Setters ----------------------------------------------------------------
    void setBkgndFn(VSAcceptanceModel* bkgnd_fn)
    {
      m_bkgnd_fn = bkgnd_fn->clone();
      setFn1(m_bkgnd_fn);
    }

    void setSourceFn(VSSourceModel* source_fn)
    {
      m_source_fn = source_fn->clone();
      setFn2(m_source_fn);
    }

    void fixBkgndParam(bool fixed = true) {  m_fn1->fixParam(fixed); }
    void fixSourceParam(bool fixed = true) {  m_fn2->fixParam(fixed); }

    // Virtual Constructor --------------------------------------------------
    virtual VSDataModel* clone() const
    {
      return new VSDataModel(*this);      
    }

    // Copy-Constructor -----------------------------------------------------
    VSDataModel(const VSDataModel& o);

  private:
    VSSourceModel*                                 m_source_fn;
    VSAcceptanceModel*                             m_bkgnd_fn;
  };
}

#endif // VSDATAMODEL_HPP
