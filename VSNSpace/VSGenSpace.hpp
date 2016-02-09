/*! \file VSGenSpace.hpp

  GenSpace generalizeed multi-dimensional space

  Multi dimensional space whose cells are arbitrary types. Has less
  functionality than normal NSpace but can be useful in some
  circumstances.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/26/2007

*/

#include<algorithm>
#include<vector>

#include<VSSimpleStat.hpp>
#include<VSNSpace.hpp>

#ifndef VSGENSPACE_HPP
#define VSGENSPACE_HPP

namespace VERITAS 
{

  template<typename GenWeight> class VSGenSpace
  {
  public:    
    typedef VSNSpace::Coord Coord;
    typedef VSNSpace::Index Index;
    typedef VSNSpace::Weight Weight;
    typedef VSNSpace::Point Point;
    typedef VSNSpace::Cell Cell;
    typedef VSNSpace::Axis Axis;
    typedef VSNSpace::Space Space;

    // ------------------------------------------------------------------------
    // CONSTRUCTOR
    // ------------------------------------------------------------------------

    //! Construct empty histogram
    VSGenSpace(): m_comment(), m_space(), m_weight() 
    { /* nothing to see here */ }

    //! Construct histogram over given space, with default occupancy weights
    VSGenSpace(const Space& space, const std::string& comment = "",
	       const GenWeight& zero_weight = GenWeight());

    // ------------------------------------------------------------------------
    // ACCESSORS
    // ------------------------------------------------------------------------

    //! Comment
    const std::string& comment() const { return m_comment; }

    //! Access to space
    const Space& space() const { return m_space; }

    //! Access to space
    const unsigned ndim() const { return m_space.ndim; }

    //! Get total number of elements in space
    unsigned size() const { return m_weight.size(); }

    //! Get weight in histogram at point
    Weight getWeightUnchecked(Index i) const { return m_weight[i]; }
    inline bool getWeight(Index i, Weight& weight) const;
    inline bool getWeight(const Cell& c, Weight& weight) const;
    inline bool getWeight(const Point& p, Weight& weight) const;

    //! Get gen weight in histogram at point
    GenWeight getGenWeightUnchecked(Index i) const { return m_weight[i]; }
    inline bool getGenWeight(Index i, GenWeight& genweight) const;
    inline bool getGenWeight(const Cell& c, GenWeight& genweight) const;
    inline bool getGenWeight(const Point& p, GenWeight& genweight) const;

    inline const GenWeight& operator[](Index i) const;
    inline const GenWeight& operator[](const Cell& c) const;
    inline const GenWeight& operator[](const Point& p) const;

    //! Access all areas
    const std::vector<GenWeight>& genweight() const { return m_weight; }

    // ------------------------------------------------------------------------
    // SETTERS
    // ------------------------------------------------------------------------

    //! Accumulate weight into histogram
    inline bool accumulate(const Point& p, const Weight weight = 1);

    //! Add an entire NSpace into thos
    bool operator += (const VSNSpace& s);

    //! Comment
    void setComment(const std::string& comment) { m_comment=comment; }

    //! Super-access to comment
    std::string& nonConstComment() { return m_comment; }

    //! Reset the space
    void clear(const GenWeight& zero_weight = GenWeight());

    //! Set data elements
    inline GenWeight& operator[](Index i);
    inline GenWeight& operator[](const Cell& c);
    inline GenWeight& operator[](const Point& p);

    //! Super-access all areas - be careful!
    std::vector<GenWeight>& nonConstGenWeight() { return m_weight; }

    // ------------------------------------------------------------------------
    // GENERATE NSPACE
    // ------------------------------------------------------------------------

    inline VSNSpace* nspace() const;    
    inline void nspace(VSNSpace& ns) const;    

  private:
    std::string            m_comment;
    Space                  m_space;
    std::vector<GenWeight> m_weight;
  };

  // ==========================================================================
  //
  // Weight classes
  //
  // ==========================================================================

  class OneSidedIntervalWeight
  {
  public:
    OneSidedIntervalWeight(double frac = 0.5): 
      m_frac(frac), m_weights(), m_sum() 
    { }

    OneSidedIntervalWeight& operator+=(VSNSpace::Weight x)
    {
      m_weights.push_back(std::make_pair(x,1.0)); 
      m_sum++;
      return *this;
    }

    void accumulate(VSNSpace::Weight x, VSNSpace::Weight w)
    {
      m_weights.push_back(std::make_pair(x,w)); 
      m_sum += w;
    }

    operator VSNSpace::Weight() const
    {
      std::vector<std::pair<VSNSpace::Weight,VSNSpace::Weight> > 
	weights(m_weights);
      std::sort(weights.begin(), weights.end());
      return interpolateArrayFraction(weights, m_sum, m_frac);
    }

    static VSNSpace::Weight 
    interpolateArrayFraction(const std::vector<VSNSpace::Weight>& w, 
			     double frac)
    {
      unsigned nw = w.size();
      if(nw==0)return 0;
      else if((nw==1)||(frac<0.0))return w.front();
      else if(frac>=1.0)return w.back();
      const double dindex = double(nw-1)*frac;
      const double x = dindex - floor(dindex);
      const unsigned iindex = unsigned(floor(dindex));
      return w[iindex]*(1.0-x) + w[iindex+1]*x;
    }

    static VSNSpace::Weight 
    interpolateArrayFraction(const std::vector<std::pair<VSNSpace::Weight,
			     VSNSpace::Weight> >& w, double sum,
			     double frac);

  private:
    double                        m_frac;
    std::vector< std::pair<VSNSpace::Weight,VSNSpace::Weight> > m_weights;
    double m_sum;
  };

  class TwoSidedIntervalWeight
  {
  public:
    TwoSidedIntervalWeight(double frac = 0.68269): 
      m_frac(frac), m_weights(), m_sum() 
    { }

    TwoSidedIntervalWeight& operator+=(VSNSpace::Weight x)
    { 
      m_weights.push_back(std::make_pair(x,1.0)); 
      m_sum++;
      return *this;
    }

    void accumulate(VSNSpace::Weight x, VSNSpace::Weight w)
    {
      m_weights.push_back(std::make_pair(x,w)); 
      m_sum += w;
    }

    operator VSNSpace::Weight() const
    {
      std::vector<std::pair<VSNSpace::Weight,VSNSpace::Weight> > 
	weights(m_weights);
      std::sort(weights.begin(), weights.end());
      return 0.5*(OneSidedIntervalWeight::
 		  interpolateArrayFraction(weights, m_sum,0.5*(1.0+m_frac)) -
 		  OneSidedIntervalWeight::
 		  interpolateArrayFraction(weights, m_sum,0.5*(1.0-m_frac)));
    }

  private:
    double                        m_frac;
    std::vector< std::pair<VSNSpace::Weight,VSNSpace::Weight> > m_weights;
    double m_sum;
  };

  class MeanWeight
  {
  public:
    MeanWeight(): m_stat() { }

    MeanWeight& operator+=(VSNSpace::Weight w)
    { 
      m_stat.accumulate(w); 
      return *this;
    }

    void accumulate(VSNSpace::Weight x, VSNSpace::Weight w)
    {
      m_stat.accumulate(x,w); 
    }

    operator VSNSpace::Weight() const { return m_stat.mean(); }

  private:
    VSSimpleStat1<VSNSpace::Weight,VSNSpace::Weight> m_stat;
  };

  class StandardDevWeight
  {
  public:
    StandardDevWeight(): m_stat() { }

    StandardDevWeight& operator+=(VSNSpace::Weight w)
    {
      m_stat.accumulate(w); 
      return *this;
    }
    
    void accumulate(VSNSpace::Weight x, VSNSpace::Weight w)
    {
      m_stat.accumulate(x,w); 
    }

    operator VSNSpace::Weight() const { return m_stat.dev(); }

  private:
    VSSimpleStat2<VSNSpace::Weight,VSNSpace::Weight> m_stat;
  };

  // ==========================================================================
  //
  // CONSTRUCTOR
  //
  // ==========================================================================

  template<typename GenWeight> VSGenSpace<GenWeight>::
  VSGenSpace(const Space& space, const std::string& comment,
	     const GenWeight& zero_weight):
    m_comment(comment), m_space(space), m_weight()
  {
    clear(zero_weight);
  }

  // ==========================================================================
  //
  // ACCESSORS
  //
  // ==========================================================================

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  getWeight(const Index i, Weight& weight) const
  {
    if(i > m_weight.size())return false;
    weight = m_weight[i];
    return true;
  }

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  getWeight(const Cell& c, Weight& weight) const
  {
    if(!m_space.isCellCompatible(c))return false;
    weight = m_weight[m_space.indexOfCellUnchecked(c)];
    return true;
  }

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  getWeight(const Point& p, Weight& weight) const
  {
    if(!m_space.isPointCompatible(p))return false;
    weight = m_weight[m_space.indexOfPointUnchecked(p)];
    return true;
  }

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  getGenWeight(const Index i, GenWeight& genweight) const
  {
    if(i > m_weight.size())return false;
    genweight = m_weight[i];
    return true;
  }

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  getGenWeight(const Cell& c, GenWeight& genweight) const
  {
    if(!m_space.isCellCompatible(c))return false;
    genweight = m_weight[m_space.indexOfCellUnchecked(c)];
    return true;
  }

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  getGenWeight(const Point& p, GenWeight& genweight) const
  {
    if(!m_space.isPointCompatible(p))return false;
    genweight = m_weight[m_space.indexOfPointUnchecked(p)];
    return true;
  }

  template<typename GenWeight> inline const GenWeight& 
  VSGenSpace<GenWeight>::operator[](Index i) const
  {
    if(i > m_weight.size())throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[i];
  }

  template<typename GenWeight> inline const GenWeight& 
  VSGenSpace<GenWeight>::operator[](const Cell& c) const
  {
    if(!m_space.isCellCompatible(c))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfCellUnchecked(c)];
  }

  template<typename GenWeight> inline const GenWeight& 
  VSGenSpace<GenWeight>::operator[](const Point& p) const
  {
    if(!m_space.isPointCompatible(p))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfPointUnchecked(p)];
  }

  // ==========================================================================
  //
  // SETTERS
  //
  // ==========================================================================

  template<typename GenWeight> inline bool VSGenSpace<GenWeight>::
  accumulate(const Point& p, const Weight weight)
  {
    if(!m_space.isPointCompatible(p))return false;
    m_weight[m_space.indexOfPointUnchecked(p)] += weight;
    return true;
  }

  template<typename GenWeight> bool VSGenSpace<GenWeight>::
  operator += (const VSNSpace& s)
  {
    if(!m_space.isSpaceEqual(s.space()))return false;
    GenWeight* iweight = &(*m_weight.begin());
    GenWeight*const nweight = &(*m_weight.end());
    const Weight* oweight = &(*s.weight().begin());
    while(iweight!=nweight)*(iweight++) += (*oweight++);
    return true;
  }

  template<typename GenWeight> void VSGenSpace<GenWeight>::
  clear(const GenWeight& zero_weight)
  {
    m_weight.clear();
    m_weight.resize(m_space.size(),zero_weight);    
  }

  template<typename GenWeight> inline GenWeight& 
  VSGenSpace<GenWeight>::operator[](Index i)
  {
    if(i > m_weight.size())throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[i];
  }

  template<typename GenWeight> inline GenWeight& 
  VSGenSpace<GenWeight>::operator[](const Cell& c)
  {
    if(!m_space.isCellCompatible(c))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfCellUnchecked(c)];
  }

  template<typename GenWeight> inline GenWeight& 
  VSGenSpace<GenWeight>::operator[](const Point& p)
  {
    if(!m_space.isPointCompatible(p))
      throw std::out_of_range(__PRETTY_FUNCTION__);
    return m_weight[m_space.indexOfPointUnchecked(p)];
  }

  // ==========================================================================
  //
  // GENERATE NSPACE
  //
  // ==========================================================================

  template<typename GenWeight> inline VSNSpace* VSGenSpace<GenWeight>::
  nspace() const
  {
    VSNSpace* ns = new VSNSpace(m_space, m_comment);
    unsigned nel = size();
    assert(ns->size() == nel);
    for(unsigned iel = 0; iel<nel; iel++)(*ns)[iel] = m_weight[iel];
    return ns;
  }

  template<typename GenWeight> inline void VSGenSpace<GenWeight>::
  nspace(VSNSpace& ns) const
  {
    ns = VSNSpace(m_space, m_comment);
    unsigned nel = size();
    assert(ns.size() == nel);
    for(unsigned iel = 0; iel<nel; iel++)ns[iel] = m_weight[iel];
  }

}

#endif // defined VSNSPACE_HPP
