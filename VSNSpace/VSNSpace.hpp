/*! \file VSNSpace.hpp

  NSpace multi-dimensional cutting

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/16/2006

  Based on Utah version, rewritten for simplicity and to agree with
  conventions of VS package

  N dimensional space class header file.  This class represents an
  NDimensional space.  This class will also query data in the Whipple
  mysql database.

  \author Jeter Hall \n
          Physics Department \n
          University of Utah \n
          E-Mail: jeter@physics.utah.edu

  \author V.V. Vassiliev \n
          Physics Department \n
          University of Utah \n
          E-Mail: vvv@physics.utah.edu

  \author T. Nagai \n
          Physics Department \n
          University of Utah \n
          E-mail: tn68@physics.utah.edu

  OLD_date   May 27, 2003

  OLD_version 0.0

*/

#include<ostream>
#include<string>
#include<set>
#include<vector>
#include<algorithm>
#include<cmath>
#include<cassert>
#include<stdexcept>

#ifndef VSNSPACE_HPP
#define VSNSPACE_HPP

#include <VSSimple2DHist.hpp>

namespace VERITAS 
{

  class VSNSpace
  {
  public:
    typedef double Coord;
    typedef double Weight;
    typedef unsigned Index;
    
    class Point
    {
    public:
      Point(): ndim(), x(), x_base(new std::vector<Coord>), maxdim()
      { setx(); }
      Point(unsigned dim): 
	ndim(dim), x(), x_base(new std::vector<Coord>(ndim)), maxdim()
      { setx(); }
      Point(unsigned dim, Coord* _x): ndim(dim), x(_x), x_base(), maxdim(dim)
      { /* nothing to see here */ }
      Point(unsigned dim, Coord* _x, unsigned maxdim):
	ndim(dim), x(_x), x_base(), maxdim(maxdim) { assert(dim<=maxdim); }
      Point(const Point& o): 
	ndim(o.ndim), x(), x_base(new std::vector<Coord>(o.begin(), o.end())),
	maxdim()
      { setx(); }
      Point& operator=(const Point& o)
      { if(x_base) { *x_base = *o.x_base; ndim=o.ndim; setx(); }
	else { assert(o.ndim<=maxdim); ndim=o.ndim; 
	  std::copy(o.begin(), o.end(), x); } return *this; }

      ~Point() { delete x_base; }

      Coord* begin() { return x; }
      Coord* end() { return x+ndim; }
      const Coord* begin() const { return x; }
      const Coord* end() const { return x+ndim; }
      inline void resize(unsigned dim)
      { if(dim!=ndim) { if(x_base) { ndim=dim; x_base->resize(dim); setx(); }
	  else { assert(dim<=maxdim); ndim=dim; } } }

      inline bool removeDimension(unsigned dim);
      inline bool addDimension(unsigned dim, Coord _x);

      inline bool marginalize(unsigned dim) { return removeDimension(dim); }
      inline bool slice(unsigned dim) { return removeDimension(dim); }
      inline bool project(unsigned dim);
      inline bool project(const std::set<unsigned>& dims);

      inline Point& operator += (const Point& o);
      inline Point& operator -= (const Point& o);
      inline Point& operator *= (const Point& o);
      inline Point& operator /= (const Point& o);
      inline Point& operator *= (const Coord& w);

      unsigned                ndim;
      Coord*                  x;
    private:
      void setx() { x=&x_base->front(); }

      std::vector<Coord>*     x_base;
      unsigned                maxdim;
    };
    
    class Cell
    {
    public:
      Cell(): ndim(), i(), i_base(new std::vector<Index>), maxdim()
      { seti(); }
      Cell(unsigned dim):
	ndim(dim), i(), i_base(new std::vector<Index>(ndim,0)), maxdim()
      { seti(); }
      Cell(unsigned dim, Index* _i): ndim(dim), i(_i), i_base(), maxdim(dim)
      { /* nothing to see here */ }
      Cell(unsigned dim, Index* _i, unsigned maxdim): 
	ndim(dim), i(_i), i_base(), maxdim(maxdim) { assert(dim<=maxdim); }
      Cell(const Cell& o): 
	ndim(o.ndim), i(), i_base(new std::vector<Index>(o.begin(), o.end())),
	maxdim()
      { seti(); }
      Cell& operator=(const Cell& o)
      { if(i_base) { *i_base = *o.i_base; ndim=o.ndim; seti(); }
	else { assert(o.ndim<=maxdim); ndim=o.ndim; 
	  std::copy(o.begin(), o.end(), i); } return *this; }

      ~Cell() { delete i_base; }

      Index* begin() { return i; }
      Index* end() { return i+ndim; }
      const Index* begin() const { return i; }
      const Index* end() const { return i+ndim; }
      inline void resize(unsigned dim)
      { if(dim!=ndim) { if(i_base) { ndim=dim; i_base->resize(dim); seti(); }
	  else { assert(dim<=maxdim); ndim=dim; } } }

      inline bool isNeighbor(const Cell& c);

      inline bool removeDimension(unsigned dim);
      inline bool addDimension(unsigned dim, Index _i);

      inline bool marginalize(unsigned dim) { return removeDimension(dim); }
      inline bool slice(unsigned dim) { return removeDimension(dim); }
      inline bool project(unsigned dim);
      inline bool project(const std::set<unsigned>& dims);

      inline bool operator == (const Cell& c) const;

      inline Cell& operator += (const Cell& c);

      unsigned                ndim;
      Index*                  i;
    private:
      void seti() { i=&i_base->front(); }

      std::vector<Index>*     i_base;
      unsigned                maxdim;
    };
    
    class Axis
    {
    public:
      Axis(): lo_bound(), hi_bound(), nbin(), name(), bin_size(), bin_factor()
      { /* nothing to see here */ }
      Axis(Coord lo, Coord hi, Index n, const std::string _name = ""):
	lo_bound(lo), hi_bound(hi), nbin(n), name(_name), 
	bin_size((hi-lo)/Coord(n)), bin_factor(1.0/bin_size)
      { vsassert(lo<=hi); vsassert(n>=1); }
      Axis(Coord lo, Coord hi, Coord delta, Index unused, 
	   const std::string _name):
	lo_bound(lo), hi_bound(lo+ceil((hi-lo)/delta)*delta),
	nbin(unsigned(ceil((hi-lo)/delta))), name(_name), 
	bin_size(delta), bin_factor(1.0/delta) 
      { vsassert(lo<=hi); vsassert(delta>0); }
      Axis(Coord lo, Coord hi, Coord delta, const std::string _name = ""):
	lo_bound(lo), hi_bound(lo+ceil((hi-lo)/delta)*delta),
	nbin(unsigned(ceil((hi-lo)/delta))), name(_name), 
	bin_size(delta), bin_factor(1.0/delta) 
      { vsassert(lo<=hi); vsassert(delta>0); }
      Coord lo_bound;
      Coord hi_bound;
      Index nbin;
      std::string name;
      Coord bin_size;
      Coord bin_factor;
      inline bool isIndexCompatible(const Index& i) const;
      inline bool isCoordCompatible(const Coord x) const;
      inline bool isAxisEqual(const Axis& a) const;
      inline bool isRangeCompatible(const Coord lo, const Coord hi) const;
      inline Index indexUnchecked(const Coord x) const;
      inline Coord midCoordUnchecked(const Index i) const;
      inline Coord maxCoordUnchecked(const Index i) const;
      inline Coord minCoordUnchecked(const Index i) const;
    };
      
    class Space
    {
    public:
      Space(): ndim(), axes() { }
      Space(unsigned dim): ndim(dim), axes(dim) { }
      unsigned ndim;
      std::vector<Axis> axes;
      inline void resize(unsigned dim)
      { if(dim!=ndim) { ndim=dim; axes.resize(dim); } }
      inline bool marginalize(unsigned dim);
      inline bool slice(unsigned dim) { return marginalize(dim); }
      inline bool project(unsigned dim);
      inline bool project(const std::set<unsigned>& dims);
      Cell cell() const { return Cell(ndim); }
      Point point() const { return Point(ndim); }
      const Axis& axis(unsigned idim) const { return axes[idim]; }
      inline Point loPoint() const;
      inline Point hiPoint() const;
      inline bool isPointCompatible(const Point& p) const;
      inline bool isCellCompatible(const Cell& c) const;
      inline bool isSpaceEqual(const Space& s) const;
      inline bool isSpaceCompatible(const Space& s) const;
      inline Index size() const;
      inline Index indexOfCellUnchecked(const Cell& c) const;
      inline Index indexOfPointUnchecked(const Point& p) const;
      inline void cellOfIndexUnchecked(const Index i, Cell& c) const;
      inline void cellOfPointUnchecked(const Point& p, Cell& c) const;
      inline void cellAndResidualOfPointUnchecked(const Point& p, Cell& c, Point& res) const;
      inline void midPointOfIndexUnchecked(const Index i, Point& p) const;
      inline void midPointOfCellUnchecked(const Cell& c, Point& p) const;
      inline void maxPointOfIndexUnchecked(const Index i, Point& p) const;
      inline void maxPointOfCellUnchecked(const Cell& c, Point& p) const;
      inline void minPointOfIndexUnchecked(const Index i, Point& p) const;
      inline void minPointOfCellUnchecked(const Cell& c, Point& p) const;
    };
      
    class Volume
    {
    public:
      Volume():
	m_space(), m_volume() { }
      Volume(const Space& space): 
	m_space(space), m_volume(space.size()) { }

      bool load(const Space& space, const std::vector<bool>& volume);
      bool loadSparse(const Space& space, std::vector<Index>& index);

      // Simple accessors
      const Space& space() const { return m_space; }
      const std::vector<bool>& volume() const { return m_volume; }
      void sparse(std::vector<Index>& index) const;
      unsigned volumeSize() const;

      // Test whether index/cell/point is in volume
      inline bool isIndexInVolumeUnchecked(const Index i) const;
      inline bool isCellInVolume(const Cell& c) const;
      inline bool isPointInVolume(const Point& p) const;

      inline bool operator[](const Index i) const 
      { return isIndexInVolumeUnchecked(i); }
      inline bool operator[](const Cell& c) const
      { return isCellInVolume(c); }
      inline bool operator[](const Point& p) const
      { return isPointInVolume(p); }

      inline bool isIndexOnEdge(const Index i) const;
      inline bool isIndexAdjacent(const Index i) const;
      inline bool isIndexIsolated(const Index i) const;

      // Setters
      inline void setIndexUnchecked(const Index i, bool set = true);
      inline bool setCell(const Cell& c, bool set = true);
      inline bool setPoint(const Point& p, bool set = true);

      // Set operators
      void setInvert();
      bool setIntersect(const Volume& o);
      bool setUnion(const Volume& o);
      bool setLess(const Volume& o);
      bool clearInsideRange(unsigned dim, const Index lo, const Index hi);

      // Get total number of cells in volume
      unsigned countCellsInVolume() const;

      //Increase volume by adding cells on the boundary
      void expandFromEdge();

    private:
      friend class VSNSpace;
      Space m_space;
      std::vector<bool> m_volume;
    };

    class Ordering
    {
    public:
      Ordering():
	m_space(),m_index(), m_counts_on(),m_counts_off(),m_excess(),
	m_Q(),m_sigma() {}
      Ordering(const Space& space):
	m_space(space),m_index(), m_counts_on(),m_counts_off(),m_excess(),
	m_Q(),m_sigma() {}

      Space m_space;
      std::vector<Index>  m_index;
      std::vector<double> m_counts_on;
      std::vector<double> m_counts_off;
      std::vector<double> m_excess;
      std::vector<double> m_Q;
      std::vector<double> m_sigma;
    };

    // ------------------------------------------------------------------------
    // CONSTRUCTOR
    // ------------------------------------------------------------------------

    //! Construct empty histogram
    VSNSpace(): m_comment(), m_space(), m_weight() 
    { /* nothing to see here */ }

    //! Construct histogram over given space, with default occupancy weights
    VSNSpace(const Space& space, const std::string& comment = "",
	     const Weight zero_weight = 0);

    //! Construct histogram from a Volume, weight is either zero or one
    VSNSpace(const Volume& volume, const std::string& comment = "",
	     const Weight zero_weight = 0, const Weight one_weight = 1);

    //! Construct an nspace of dimension ndim
    VSNSpace(unsigned ndim, Coord lo, Coord hi, Index nbin);

    //! Construct a 1D nspace
    VSNSpace(Coord lo, Coord hi, Index nbin);

    //! Construct a 2D nspace
    VSNSpace(Coord lox, Coord hix, Index nbinx,
	     Coord loy, Coord hiy, Index nbiny);

    //! Construct nspace from 1D histogram
    VSNSpace(const VSLimitedHist<double,double>& h);

    //! Construct nspace from 2D histogram
    VSNSpace(const VSSimple2DHist<double,double>& h);

    bool load(const Space& space, const std::vector<Weight>& weight,
	      const std::string& comment = "");

    void hist(const std::set<unsigned>& dims, 
	      const std::vector< double >& coords,
	      VSLimitedHist<double,double>& h);

    void hist(const std::set<unsigned>& dims, 
	      const std::vector< unsigned >& indices,
	      VSLimitedHist<double,double>& h);

    // ------------------------------------------------------------------------
    // PRIMARY FUNCTIONALITY -- ACCESSORS
    // ------------------------------------------------------------------------

    //! Get total number of elements in space
    unsigned size() const { return m_weight.size(); }

    //! Get weight in histogram at point
    Weight getWeightUnchecked(Index i) const { return m_weight[i]; }
    inline bool getWeight(Index i, Weight& weight) const;
    inline bool getWeight(const Cell& c, Weight& weight) const;
    inline bool getWeight(const Point& p, Weight& weight) const;

    inline const Weight& operator[](Index i) const;
    inline const Weight& operator[](const Cell& c) const;
    inline const Weight& operator[](const Point& p) const;

    //! Interpolate (linearly) the weight at a point in the space
    //  using the values of the point
    bool interpolateWeight(const Point& p, Weight& weight) const;

    //! Comment
    const std::string& comment() const { return m_comment; }

    //! Access all areas
    const std::vector<Weight>& weight() const { return m_weight; }

    // ------------------------------------------------------------------------
    // PRIMARY FUNCTIONALITY -- SETTERS
    // ------------------------------------------------------------------------

    //! Accumulate weight into histogram
    inline bool accumulate(const Point& p, const Weight weight = 1);

    //! Set weight in histogram at point
    void setWeightUnchecked(Index i, const Weight& w) { m_weight[i]=w; }
    inline bool setWeight(Index i, const Weight& weight);
    inline bool setWeight(const Cell& c, const Weight& weight);
    inline bool setWeight(const Point& p, const Weight& weight);

    inline Weight& operator[](Index i);
    inline Weight& operator[](const Cell& c);
    inline Weight& operator[](const Point& p);

    //! Comment
    void setComment(const std::string& comment) { m_comment=comment; }

    // ------------------------------------------------------------------------
    // ACCESSORS
    // ------------------------------------------------------------------------

    //! Access to space
    const Space& space() const { return m_space; }

    //! Access to space
    const unsigned ndim() const { return m_space.ndim; }

    const Axis& axis(unsigned idim) const 
    { return m_space.axes[idim]; }

    //! Get total weight in histogram
    Weight totalWeight() const;

    //! Get maximum weight in histogram
    Weight maxWeight() const;

    //! Get maximum weight in histogram, also return index of maximum
    Weight maxWeight(Index& index) const;

    //! Get minimum weight in histogram
    Weight minWeight() const;

    //! Get minimum weight in histogram, also return index of minimum
    Weight minWeight(Index& index) const;

    //! Get total weight in histogram region defined by volume
    bool totalWeight(const Volume& volume, Weight& total) const;
    
    //! Get maximum weight in histogram region defined by volume
    bool maxWeight(const Volume& volume, Weight& max) const;

    //! Get minimum weight in histogram region defined by volume
    bool minWeight(const Volume& volume, Weight& min) const;

    //! Get maximum weight in histogram, also return index of maximum
    bool maxWeight(const Volume& volume, Weight& max, Index& index) const;

    //! Get minimum weight in histogram, also return index of minimum
    bool minWeight(const Volume& volume, Weight& max, Index& index) const;

    //! Calculate ln likelihood over space
    Weight logLikelihood() const;

    //! Calculate ln likelihood over volume
    bool logLikelihood(const Volume& volume, Weight& loglike) const;

    //! General integration over the space. The Functional type should
    //! have a member: Weight operator() (const VSNSpace& p) whic returns
    //! some value at each point in the space. This is multiplied by the
    //! weight and accumulated over all cells in the space
    template<typename Functional, typename T>
    void integrate(Functional& functional, T& integral) const;

    //! General integration over the space in the volume
    template<typename Functional, typename T>
    bool integrate(Functional& functional, const Volume& volume, 
		   T& integral) const;

    //! General integration over the space (with const functor)
    template<typename Functional, typename T>
    void integrate(const Functional& functional, T& integral) const;
    
    //! General integration over the space in the volume (with const functor)
    template<typename Functional, typename T>
    bool integrate(const Functional& functional, const Volume& volume, 
		   T& integral) const;

    //! Make a volume to include all cells
    void volumeAll(Volume& volume) const;

    //! Make a volume to include all cells within subspace
    bool volumeSubSpace(const Space& space, Volume& volume) const;    

    //! Make a volume to include those cells above some threshold
    void volumeAbove(Weight threshold, Volume& volume) const;

    //! Make a volume to include those cells below some threshold
    void volumeBelow(Weight threshold, Volume& volume) const;

    //! Make a volume to include those cells above/equal to some threshold
    void volumeAboveEqual(Weight threshold, Volume& volume) const;

    //! Make a volume to include those cells below/equal to some threshold
    void volumeBelowEqual(Weight threshold, Volume& volume) const;

    //! Make a volume to include those cells above some threshold
    VSNSpace::Volume operator>(Weight threshold) const 
    { VSNSpace::Volume vol; volumeAbove(threshold,vol); return vol; }

    //! Make a volume to include those cells below some threshold
    VSNSpace::Volume operator<(Weight threshold) const
    { VSNSpace::Volume vol; volumeBelow(threshold,vol); return vol; }
      
    //! Make a volume to include those cells above some threshold
    VSNSpace::Volume operator>=(Weight threshold) const 
    { VSNSpace::Volume vol; volumeAboveEqual(threshold,vol); return vol; }

    //! Make a volume to include those cells below some threshold
    VSNSpace::Volume operator<=(Weight threshold) const
    { VSNSpace::Volume vol; volumeBelowEqual(threshold,vol); return vol; }

    // ------------------------------------------------------------------------
    // ORDERING
    // ------------------------------------------------------------------------

    static bool
    forwardOrderingSimple(std::vector<Index>& ordering,
			  const VSNSpace& excess, const VSNSpace& variance,
			  Weight min_variance = 0);

    static bool
    cumulativeOrderingSimple(std::vector<Index>& ordering,
			     const VSNSpace& excess, const VSNSpace& variance,
			     Weight min_variance = 0);
    
    // ------------------------------------------------------------------------
    // SETTERS AND OPERATORS
    // ------------------------------------------------------------------------

    //! Super-access to comment
    std::string& nonConstComment() { return m_comment; }

    //! Super-access all areas - be careful!
    std::vector<Weight>& nonConstWeight() { return m_weight; }

    //! Reset the space
    void clear(const Weight zero_weight = 0);

    //! Zero elements in space outside of volume
    bool clearOutsideVolume(const Volume& volume,
			    const Weight zero_weight = 0);

    //! Zero elements below threshold
    void clearKeepOnlyAbove(Weight threshold, const Weight zero_weight = 0);

    //! Zero elements below threshold
    void clearKeepOnlyBelow(Weight threshold, const Weight zero_weight = 0);

    //! Modify space so that weight_j = sqrt(weight_j)
    VSNSpace& sqrt();

    //! Modify space so that weight_j = ln(weight_j)
    VSNSpace& log();

    //! Normalize to total weight in space
    VSNSpace& normalize();

    //! Normalize to total weight in volume
    bool normalize(const Volume& volume);

    //! Transform to cumulative space. New weight (W) is given in terms
    //  of old weight (w) as: W(I,J,K) = Sum(i<=I,j<=J,j<=k) w(i,j,k)
    bool partiallyIntegrate();

    //! Marginalize over given dimension
    bool marginalize(unsigned dim);

    //! Marginalize cells within volume over given dimension
    bool marginalize(unsigned dim, const Volume& volume);

    //! Marginalize over given dimension weighting cells using supplied
    //  1-D space "dim_weight"
    bool marginalize(unsigned dim, const VSNSpace& dim_weight);

    bool slice(const std::vector<std::pair<unsigned,Index> >& dim, 
	       VSNSpace& o);

    bool slice(const std::vector<std::pair<unsigned,Coord> >& dim, 
	       VSNSpace& o);

    //! Extract slice in one dimension
    bool slice(unsigned dim, Index index);

    //! Extract slice in one dimension
    bool slice(unsigned dim, Coord coord);

    //! Project histogram onto given dimension
    bool project(unsigned dim);

    //! Project histogram volume onto given dimension
    bool project(unsigned dim, const Volume& volume);

    //! Project histogram onto set of dimensions
    bool project(const std::set<unsigned>& dims);

    //! Project histogram volume onto set of dimensions
    bool project(const std::set<unsigned>& dims, const Volume& volume);

    //! Scale all histogram elements by w
    VSNSpace& operator *= (Weight w);
    VSNSpace operator * (Weight w) const;
    
    //! Element-wise addition of two spaces (note: returns a bool)
    bool operator += (const VSNSpace& s);

    //! Element-wise subtraction of two spaces (note: returns a bool)
    bool operator -= (const VSNSpace& s);

    //! Element-wise multiplication of two spaces (note: returns a bool)
    bool operator *= (const VSNSpace& s);

    //! Element-wise division of two spaces (note: returns a bool)
    bool operator /= (const VSNSpace& s);

#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

  private:
    static double gammaln(double); 

    std::string         m_comment;
    Space               m_space;
    std::vector<Weight> m_weight;
  };

  class VSNSpaceIO
  {
  public:
    virtual ~VSNSpaceIO();
    virtual bool readSpace(const std::string& filename,
			   VSNSpace::Space& space) = 0;
    virtual bool writeSpace(const std::string& filename,
			    const VSNSpace::Space& space) = 0;
    virtual bool readVolume(const std::string& filename, 
			    VSNSpace::Volume& vol) = 0;
    virtual bool writeVolume(const std::string& filename,
			     const VSNSpace::Volume& vol,
			     bool sparse = false) = 0;
    virtual bool readOrdering(const std::string& filename,
			      VSNSpace::Ordering& ordering) = 0;
    virtual bool writeOrdering(const std::string& filename,
			       const VSNSpace::Ordering& ordering) = 0;
    virtual bool readHistogram(const std::string& filename,
			       VSNSpace& vol) = 0;
    virtual bool writeHistogram(const std::string& filename,
				const VSNSpace& vol) = 0;
  };

  std::ostream& operator<<(std::ostream& stream, const VSNSpace::Cell& c);
  std::ostream& operator<<(std::ostream& stream, const VSNSpace::Axis& a);
  std::ostream& operator<<(std::ostream& stream, const VSNSpace::Space& a);
  std::ostream& operator<<(std::ostream& stream, const VSNSpace& n);

  // ==========================================================================
  //
  // INLINED MEMBER FUNCTIONS FOR --- VSNSpace::Point
  //
  // ==========================================================================
  
  inline bool VSNSpace::Point::
  addDimension(unsigned dim, Coord _x)
  {
    if(dim>ndim)return false;
    resize(ndim+1);
    unsigned idim=ndim;
    while(--idim>dim)x[idim]=x[idim-1];
    x[dim] = _x;
    return true;
  }

  inline bool VSNSpace::Point::removeDimension(unsigned dim)
  {
    if(dim>=ndim)return false;
    while(++dim<ndim)x[dim-1]=x[dim];
    resize(ndim-1);
    return true;
  }

  inline bool VSNSpace::Point::project(unsigned dim)
  {
    if(dim>=ndim)return false;
    Coord _x = x[dim];
    resize(1);
    x[0]=_x;
    return true;
  }
  
  inline bool VSNSpace::Point::project(const std::set<unsigned>& dims)
  {
    if(dims.empty())return false;
    if(*(--dims.end()) >= ndim)return false;
    unsigned new_ndim = 0;
    for(std::set<unsigned>::const_iterator idim = dims.begin(); 
	idim!=dims.end(); idim++)x[new_ndim++] = x[*idim];
    resize(new_ndim);
    return true;
  }
  
  inline VSNSpace::Point& 
  VSNSpace::Point::operator += (const VSNSpace::Point& o)
  {
    for(unsigned idim=0;idim<ndim;idim++)x[idim]+=o.x[idim];
    return *this;
  }

  inline VSNSpace::Point& 
  VSNSpace::Point::operator -= (const VSNSpace::Point& o)
  {
    for(unsigned idim=0;idim<ndim;idim++)x[idim]-=o.x[idim];
    return *this;
  }

  inline VSNSpace::Point& 
  VSNSpace::Point::operator *= (const VSNSpace::Point& o)
  {
    for(unsigned idim=0;idim<ndim;idim++)x[idim]*=o.x[idim];
    return *this;
  }

  inline VSNSpace::Point&
  VSNSpace::Point::operator /= (const VSNSpace::Point& o)
  {
    for(unsigned idim=0;idim<ndim;idim++)x[idim]/=o.x[idim];
    return *this;
  }

  inline VSNSpace::Point& 
  VSNSpace::Point::operator *= (const VSNSpace::Weight& w)
  {
    for(unsigned idim=0;idim<ndim;idim++)x[idim]*=w;
    return *this;
  }

  // ==========================================================================
  //
  // INLINED MEMBER FUNCTIONS FOR --- VSNSpace::Cell
  //
  // ==========================================================================

  inline bool VSNSpace::Cell::
  addDimension(unsigned dim, Index _i)
  {
    if(dim>ndim)return false;
    resize(ndim+1);
    unsigned idim=ndim;
    while(--idim>dim)i[idim]=i[idim-1];
    i[dim] = _i;
    return true;
  }

  inline bool VSNSpace::Cell::removeDimension(unsigned dim)
  {
    if(dim>=ndim)return false;
    while(++dim<ndim)i[dim-1]=i[dim];
    resize(ndim-1);
    return true;
  }

  inline bool VSNSpace::Cell::project(unsigned dim)
  {
    if(dim>=ndim)return false;
    Index _i = i[dim];
    resize(1);
    i[0]=_i;
    return true;
  }

  inline bool VSNSpace::Cell::project(const std::set<unsigned>& dims)
  {
    if(dims.empty())return false;
    if(*(--dims.end()) >= ndim)return false;
    unsigned new_ndim = 0;
    for(std::set<unsigned>::const_iterator idim = dims.begin();
	idim!=dims.end(); idim++)i[new_ndim++] = i[*idim];
    resize(new_ndim);
    return true;
  }

  inline bool VSNSpace::Cell::operator == (const Cell& c) const
  {
    for(unsigned idim = 0; idim < ndim; idim++)
      if(c.i[idim] != i[idim]) return false;
    
    return true;
  }
  
  inline VSNSpace::Cell& VSNSpace::Cell::operator += (const Cell& c)
  {
    for(unsigned idim = 0; idim < ndim; idim++) i[idim] += c.i[idim];
    return *this;
  }

  // ==========================================================================
  //
  // INLINED MEMBER FUNCTIONS FOR --- VSNSpace::Axis
  //
  // ==========================================================================

  inline bool VSNSpace::Axis::isIndexCompatible(const Index& i) const
  {
    return i<nbin;
  }

  inline bool VSNSpace::Axis::isCoordCompatible(const Coord x) const
  {
    return (x>=lo_bound)&&(x<hi_bound);
  }

  inline bool VSNSpace::Axis::isAxisEqual(const Axis& a) const
  {
    return
      (a.lo_bound == lo_bound)&&(a.hi_bound == hi_bound)&&(a.nbin == nbin);
  }

  inline bool VSNSpace::Axis::isRangeCompatible(const Coord lo, const Coord hi)
    const
  {
    return (lo>=lo_bound)&&(hi<=hi_bound);
  }

  inline VSNSpace::Index VSNSpace::Axis::indexUnchecked(const Coord x) const
  {
    return Index(floor((x-lo_bound)*bin_factor));
  }

  inline VSNSpace::Coord VSNSpace::Axis::midCoordUnchecked(const Index i) const
  {
    return (Coord(i)+0.5)*bin_size+lo_bound;
  }

  inline VSNSpace::Coord VSNSpace::Axis::maxCoordUnchecked(const Index i) const
  {
    return (Coord(i)+1.0)*bin_size+lo_bound;
  }

  inline VSNSpace::Coord VSNSpace::Axis::minCoordUnchecked(const Index i) const
  {
    return (Coord(i))*bin_size+lo_bound;
  }
  
  // ==========================================================================
  //
  // INLINED MEMBER FUNCTIONS FOR --- VSNSpace::Space
  //
  // ==========================================================================

  inline bool VSNSpace::Space::marginalize(unsigned dim)
  {
    if(dim>=ndim)return false;
    while(++dim<ndim)axes[dim-1]=axes[dim];
    axes.resize(--ndim);
    return true;
  }

  inline bool VSNSpace::Space::project(unsigned dim)
  {
    if(dim>=ndim)return false;
    Axis _a=axes[dim];
    ndim=1;
    axes.resize(1);
    axes[0]=_a;
    return true;
  }

  inline bool VSNSpace::Space::project(const std::set<unsigned>& dims)
  {
    if(dims.empty())return false;
    if(*(--dims.end()) >= ndim)return false;
    ndim = 0;
    for(std::set<unsigned>::const_iterator idim = dims.begin(); 
	idim!=dims.end(); idim++)axes[ndim++] = axes[*idim];
    axes.resize(ndim);
    return true;
  }

  inline VSNSpace::Point VSNSpace::Space::loPoint() const
  {
    Point p(ndim);
    const Axis* iaxis      = &(*axes.begin());
    const Axis*const naxis = &(*axes.end());
    Coord* ix = p.begin();
    while(iaxis!=naxis)*(ix++) = (iaxis++)->lo_bound;
    return p;
  }

  inline VSNSpace::Point VSNSpace::Space::hiPoint() const
  {
    Point p(ndim);
    const Axis* iaxis      = &(*axes.begin());
    const Axis*const naxis = &(*axes.end());
    Coord* ix = p.begin();
    while(iaxis!=naxis)*(ix++) = (iaxis++)->hi_bound;
    return p;  
  }

  inline bool VSNSpace::Space::
  isPointCompatible(const Point& p) const
  {
    if(p.ndim != ndim)return false;
    const Coord* ix = p.begin();
    const Coord* nx = p.end();
    const Axis* iaxis = &(*axes.begin());
    while(ix!=nx)
      {
	if(!iaxis->isCoordCompatible(*ix))return false;
	ix++;
	iaxis++;
      }
    return true;
  }
  
  inline bool VSNSpace::Space::
  isCellCompatible(const Cell& c) const
  {
    if(c.ndim != ndim)return false;
    const Index* iindex = c.begin();
    const Index* nindex = c.end();
    const Axis* iaxis = &(*axes.begin());
    while(iindex!=nindex)
      {
	if(!iaxis->isIndexCompatible(*iindex))return false;
	iindex++;
	iaxis++;
      }
    return true;
  }

  inline bool VSNSpace::Space::
  isSpaceEqual(const Space& s) const
  {
    if(s.ndim != ndim)return false;
    const Axis* iaxis_s = &(*s.axes.begin());
    const Axis* naxis_s = &(*s.axes.end());
    const Axis* iaxis_this = &(*axes.begin());
    while(iaxis_s!=naxis_s)
      {
	if(!iaxis_this->isAxisEqual(*iaxis_s))return false;
	iaxis_s++;
	iaxis_this++;
      }
    return true;
  }

  inline bool VSNSpace::Space::
  isSpaceCompatible(const Space& s) const
  {
    if(s.ndim != ndim)return false;
    const Axis* iaxis_s = &(*s.axes.begin());
    const Axis* naxis_s = &(*s.axes.end());
    const Axis* iaxis_this = &(*axes.begin());
    while(iaxis_s!=naxis_s)
      {
	if(!iaxis_this->
	   isRangeCompatible(iaxis_s->lo_bound, iaxis_s->hi_bound))
	  return false;
	iaxis_s++;
	iaxis_this++;
      }
    return true;
  }

  inline VSNSpace::Index VSNSpace::Space::size() const
  {
    if(ndim==0) return 0;
    const Axis* iaxis = &(*axes.begin());
    const Axis* naxis = &(*axes.end());
    Index nbin = iaxis->nbin;
    while(++iaxis!=naxis)nbin *= iaxis->nbin;
    return nbin;
  }
      
  inline VSNSpace::Index 
  VSNSpace::Space::indexOfCellUnchecked(const Cell& c) const
  {
    const Index* iindex = c.begin();
    const Index* nindex = c.end();
    const Axis* iaxis = &(*axes.begin());
    Index index = *iindex;
    while(++iindex!=nindex)
      {
	iaxis++; 
	index = index*iaxis->nbin + *iindex; 
      }
    return index;
  }

  inline VSNSpace::Index
  VSNSpace::Space::indexOfPointUnchecked(const Point& p) const
  {
    const Coord* ix = p.begin();
    const Coord* nx = p.end();
    const Axis* iaxis = &(*axes.begin());
    Index index = iaxis->indexUnchecked(*ix);
    while(++ix!=nx)
      {
	iaxis++; 
	index = index*iaxis->nbin + iaxis->indexUnchecked(*ix); 
      }
    return index;
  }
  
  inline void 
  VSNSpace::Space::cellOfIndexUnchecked(const Index i, Cell& c) const
  {
    c.resize(ndim);
    const Axis* iaxis = &(*axes.begin());
    const Axis* naxis = &(*axes.end());
    Index* nindex = c.end();
    Index iremain=i;
    while(--naxis!=iaxis)
      {
	--nindex;
	*nindex = iremain%naxis->nbin;
	iremain /= naxis->nbin;
      }
    --nindex;
    *nindex = iremain;
  }

  inline void 
  VSNSpace::Space::cellOfPointUnchecked(const Point& p, Cell& c) const
  {
    c.resize(ndim);
    const Coord* ix = p.begin();
    const Coord* nx = p.end();
    Index* iindex = c.begin();
    const Axis* iaxis = &(*axes.begin());
    while(ix!=nx)
      {
	*iindex = iaxis->indexUnchecked(*ix);
	ix++;
	iindex++;
	iaxis++;
      }
  }

  inline void VSNSpace::Space::
  cellAndResidualOfPointUnchecked(const Point& p, Cell& c, Point& res) const
  {
    c.resize(ndim);
    res.resize(ndim);
    const Coord* ix = p.begin();
    const Coord* nx = p.end();
    Index* iindex = c.begin();
    Coord* ires = res.begin();
    const Axis* iaxis = &(*axes.begin());
    while(ix!=nx)
      {
	*iindex = iaxis->indexUnchecked(*ix);
	*ires = *ix - iaxis->minCoordUnchecked(*iindex);
	ix++;
	iindex++;
	ires++;
	iaxis++;
      }
  }

  inline void 
  VSNSpace::Space::midPointOfIndexUnchecked(const Index i, Point& p) const
  {
    p.resize(ndim);
    const Axis* iaxis = &(*axes.begin());
    const Axis* naxis = &(*axes.end());
    Coord* nx = p.end();
    Index iremain=i;
    while(--naxis!=iaxis)
      {
	--nx;
	*nx = naxis->midCoordUnchecked(iremain%naxis->nbin);
	iremain /= naxis->nbin;
      }
    --nx;
    *nx = naxis->midCoordUnchecked(iremain);
  }

  inline void 
  VSNSpace::Space::midPointOfCellUnchecked(const Cell& c, Point& p) const
  {
    p.resize(ndim);
    const Index* iindex = c.begin();
    const Index* nindex = c.end();
    Coord* ix = p.begin();
    const Axis* iaxis = &(*axes.begin());
    while(iindex!=nindex)
      {
	*ix = iaxis->midCoordUnchecked(*iindex);
	iindex++;
	ix++;
	iaxis++;
      }
  }

  inline void 
  VSNSpace::Space::maxPointOfIndexUnchecked(const Index i, Point& p) const
  {
    p.resize(ndim);
    const Axis* iaxis = &(*axes.begin());
    const Axis* naxis = &(*axes.end());
    Coord* nx = p.end();
    Index iremain=i;
    while(--naxis!=iaxis)
      {
	--nx;
	*nx = naxis->maxCoordUnchecked(iremain%naxis->nbin);
	iremain /= naxis->nbin;
      }
    --nx;
    *nx = naxis->maxCoordUnchecked(iremain);
  }

  inline void 
  VSNSpace::Space::maxPointOfCellUnchecked(const Cell& c, Point& p) const
  {
    p.resize(ndim);
    const Index* iindex = c.begin();
    const Index* nindex = c.end();
    Coord* ix = p.begin();
    const Axis* iaxis = &(*axes.begin());
    while(iindex!=nindex)
      {
	*ix = iaxis->maxCoordUnchecked(*iindex);
	iindex++;
	ix++;
	iaxis++;
      }
  }

  inline void 
  VSNSpace::Space::minPointOfIndexUnchecked(const Index i, Point& p) const
  {
    p.resize(ndim);
    const Axis* iaxis = &(*axes.begin());
    const Axis* naxis = &(*axes.end());
    Coord* nx = p.end();
    Index iremain=i;
    while(--naxis!=iaxis)
      {
	--nx;
	*nx = naxis->minCoordUnchecked(iremain%naxis->nbin);
	iremain /= naxis->nbin;
      }
    --nx;
    *nx = naxis->minCoordUnchecked(iremain);
  }

  inline void 
  VSNSpace::Space::minPointOfCellUnchecked(const Cell& c, Point& p) const
  {
    p.resize(ndim);
    const Index* iindex = c.begin();
    const Index* nindex = c.end();
    Coord* ix = p.begin();
    const Axis* iaxis = &(*axes.begin());
    while(iindex!=nindex)
      {
	*ix = iaxis->minCoordUnchecked(*iindex);
	iindex++;
	ix++;
	iaxis++;
      }
  }

  // ==========================================================================
  //
  // INLINED MEMBER FUNCTIONS FOR --- VSNSpace::Volume
  //
  // ==========================================================================

  inline bool VSNSpace::Volume::isIndexInVolumeUnchecked(const Index i) const
  {    
    return m_volume[i];
  }

  inline bool VSNSpace::Volume::isPointInVolume(const Point& p) const
  {    
    if(!m_space.isPointCompatible(p))return false;
    Index i = m_space.indexOfPointUnchecked(p);
    return m_volume[i];
  }

  inline bool VSNSpace::Volume::isCellInVolume(const Cell& c) const
  {    
    if(!m_space.isCellCompatible(c))return false;
    Index i = m_space.indexOfCellUnchecked(c);
    return m_volume[i];
  }

  inline bool VSNSpace::Volume::isIndexIsolated(const Index i) const
  {    
    VSNSpace::Cell c;
    m_space.cellOfIndexUnchecked(i,c);
    for(unsigned idim = 0; idim < m_space.ndim; idim++)
      {
	VSNSpace::Cell c2 = c;
	c2.i[idim]=c.i[idim]+1;
	if(isCellInVolume(c2)) return false;
	c2.i[idim]=c.i[idim]-1;
	if(isCellInVolume(c2)) return false;	    
      }

    return true;
  }

  inline bool VSNSpace::Volume::isIndexOnEdge(const Index i) const
  {    
    VSNSpace::Cell c;
    m_space.cellOfIndexUnchecked(i,c);

    if(!isCellInVolume(c)) return false;

    for(unsigned idim = 0; idim < m_space.ndim; idim++)
      {
	if(c.i[idim] == 0 || c.i[idim] == m_space.axes[idim].nbin-1) 
	  return true;

	VSNSpace::Cell c2 = c;
	c2.i[idim]=c.i[idim]+1;
	if(!isCellInVolume(c2)) return true;
	c2.i[idim]=c.i[idim]-1;
	if(!isCellInVolume(c2)) return true;
      }

    return false;
  }

  inline bool VSNSpace::Volume::isIndexAdjacent(const Index i) const
  {    
    VSNSpace::Cell c;
    m_space.cellOfIndexUnchecked(i,c);

    if(isCellInVolume(c)) return false;

    unsigned icmax = lround(std::pow(3.0,(int)m_space.ndim));
    VSNSpace::Cell cc = m_space.cell();

    // std::cout << std::endl;
    // std::cout << "----------------------------------" << std::endl;
    // std::cout << c << std::endl;
    // std::cout << "icmax " << icmax << std::endl;
    unsigned icell = 0;
    while(icell++ < icmax)
      {
	bool z = 0;
	VSNSpace::Cell ctmp = cc;
	ctmp += c;
	for(unsigned idim = 0; idim < m_space.ndim; idim++)
	  {
	    if(ctmp.i[idim] == 0) z = true;
	    else if(cc.i[idim] == m_space.axes[idim].nbin+1) z = true;
	    else ctmp.i[idim]--;
	  }

	for(unsigned idim = 0; idim < m_space.ndim; idim++)
	  {
	    cc.i[idim]++;
	    cc.i[idim] = cc.i[idim] % 3;
	    if(cc.i[idim]) break;
	  }

	if(z || ctmp == c) continue;
	if(isCellInVolume(ctmp)) return true;
      }

    return false;
  }

  inline void VSNSpace::Volume::setIndexUnchecked(const Index i, bool set)
  {
    m_volume[i]=set;
  }

  inline bool VSNSpace::Volume::setPoint(const Point& p, bool set)
  {
    if(!m_space.isPointCompatible(p))return false;
    Index i = m_space.indexOfPointUnchecked(p);
    m_volume[i] = set;
    return true;
  }

  inline bool VSNSpace::Volume::setCell(const Cell& c, bool set)
  {
    if(!m_space.isCellCompatible(c))return false;
    Index i = m_space.indexOfCellUnchecked(c);
    m_volume[i] = set;
    return true;
  }

  // ==========================================================================
  //
  // INLINED MEMBER FUNCTIONS FOR --- VSNSpace
  //
  // ==========================================================================

  inline bool VSNSpace::accumulate(const Point& p, const Weight weight)
  {
    if(!m_space.isPointCompatible(p))return false;
    m_weight[m_space.indexOfPointUnchecked(p)] += weight;
    return true;
  }

  inline bool VSNSpace::getWeight(const Index i, Weight& weight) const
  {
    if(i > m_weight.size())return false;
    weight = m_weight[i];
    return true;
  }

  inline bool VSNSpace::getWeight(const Cell& c, Weight& weight) const
  {
    if(!m_space.isCellCompatible(c))return false;
    weight = m_weight[m_space.indexOfCellUnchecked(c)];
    return true;
  }

  inline bool VSNSpace::getWeight(const Point& p, Weight& weight) const
  {
    if(!m_space.isPointCompatible(p))return false;
    weight = m_weight[m_space.indexOfPointUnchecked(p)];
    return true;
  }

  inline const VSNSpace::Weight& VSNSpace::operator[](Index i) const
  {
    if(i > m_weight.size())throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[i];
  }

  inline const VSNSpace::Weight& VSNSpace::operator[](const Cell& c) const
  {
    if(!m_space.isCellCompatible(c))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfCellUnchecked(c)];
  }

  inline const VSNSpace::Weight& VSNSpace::operator[](const Point& p) const
  {
    if(!m_space.isPointCompatible(p))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfPointUnchecked(p)];
  }

  bool VSNSpace::setWeight(Index i, const Weight& weight)
  {
    if(i > m_weight.size())return false;
    m_weight[i] = weight;
    return true;
  }

  bool VSNSpace::setWeight(const Cell& c, const Weight& weight)
  {
    if(!m_space.isCellCompatible(c))return false;
    m_weight[m_space.indexOfCellUnchecked(c)] = weight;
    return true;
  }

  bool VSNSpace::setWeight(const Point& p, const Weight& weight)
  {
    if(!m_space.isPointCompatible(p))return false;
    m_weight[m_space.indexOfPointUnchecked(p)] = weight;
    return true;
  }

  VSNSpace::Weight& VSNSpace::operator[](Index i)
  {
    if(i > m_weight.size())throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[i];
  }

  VSNSpace::Weight& VSNSpace::operator[](const Cell& c)
  {
    if(!m_space.isCellCompatible(c))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfCellUnchecked(c)];
  }

  VSNSpace::Weight& VSNSpace::operator[](const Point& p)
  {
    if(!m_space.isPointCompatible(p))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfPointUnchecked(p)];
  }

  template<typename Functional, typename T> void
  VSNSpace::integrate(Functional& functional, T& integral) const
  {
    const Weight* iweight = &(*m_weight.begin());
    const Weight* nweight = &(*m_weight.end());
    unsigned iindex = 0;    
    while(iweight!=nweight)
      {
	Point p;
	m_space.midPointOfIndexUnchecked(iindex,p);
	T x;
	functional(p, x);
	x *= (*iweight);
	if(iindex)integral += x;
	else integral = x;
	iweight++;
	iindex++;
      }
    return;
  }
  
  template<typename Functional, typename T> bool
  VSNSpace::integrate(Functional& functional, const Volume& volume,
		      T& integral) const
  {
    integral = T();
    if(!m_space.isSpaceEqual(volume.m_space))return false;
    const Weight* iweight = &(*m_weight.begin());
    const Weight* nweight = &(*m_weight.end());
    unsigned iindex = 0;    
    bool got_one = false;
    while(iweight!=nweight)
      {
	if(volume.m_volume[iindex])
	  {
	    Point p;
	    m_space.midPointOfIndexUnchecked(iindex,p);
	    T x;
	    functional(p, x);
	    x *= (*iweight);
	    if(got_one)integral += x;
	    else integral = x, got_one=true;
	  }
	iweight++;
	iindex++;
      }
    return true;
  }

  template<typename Functional, typename T> void
  VSNSpace::integrate(const Functional& functional, T& integral) const
  {
    integral = T();
    const Weight* iweight = &(*m_weight.begin());
    const Weight* nweight = &(*m_weight.end());
    unsigned iindex = 0;    
    while(iweight!=nweight)
      {
	Point p;
	m_space.midPointOfIndexUnchecked(iindex,p);
	T x;
	functional(p, x);
	x *= (*iweight);
	if(iindex)integral += x;
	else integral = x;
	iweight++;
	iindex++;
      }
    return;
  }
  
  template<typename Functional, typename T> bool
  VSNSpace::integrate(const Functional& functional, const Volume& volume,
		      T& integral) const
  {
    integral = T();
    if(!m_space.isSpaceEqual(volume.m_space))return false;
    const Weight* iweight = &(*m_weight.begin());
    const Weight* nweight = &(*m_weight.end());
    unsigned iindex = 0;    
    bool got_one = false;
    while(iweight!=nweight)
      {
	if(volume.m_volume[iindex])
	  {
	    Point p;
	    m_space.midPointOfIndexUnchecked(iindex,p);
	    T x;
	    functional(p, x);
	    x *= (*iweight);
	    if(got_one)integral += x;
	    else integral = x, got_one=true;
	  }
	iweight++;
	iindex++;
      }
    return true;
  }
  
}

#endif // defined VSNSPACE_HPP
