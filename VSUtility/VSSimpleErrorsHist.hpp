//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleErrorsHist.hpp

  Simple histogram class with errors.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005
*/

#include <string>
#include <vector>
#include <cmath>
#include <vsassert>

#ifndef VSSIMPLEERRORSHIST_HPP
#define VSSIMPLEERRORSHIST_HPP

#include<VSSimpleHist.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSLimitedErrorsHist
  // ==========================================================================
  template<typename T, typename count_type=double, 
	   typename BIN=VSBinCalcLinear<T> >
  class VSLimitedErrorsHist
  {
  public:
    class iterator
    {
    public:
      iterator(const iterator& i): fHist(i.fHist), fBin(i.fBin) { }
      count_type count() const { return fHist->count(fBin); }
      T val() const { return fHist->val(fBin); }
      T center() const { return fHist->center(fBin); }
      count_type err() const { return fHist->err(fBin); }
      count_type var() const { return fHist->var(fBin); }
      int bin() const { return fBin; }
      const iterator* operator->() { return this; }
      iterator& operator++() { fBin++; return *this; }
      iterator operator++(int) { iterator i=*this; fBin++; return i; }
      iterator& operator--() { fBin--; return *this; }
      iterator operator--(int) { iterator i=*this; fBin--; return i; }
      iterator& operator=(const iterator& i) 
      { fHist=i.fHist; fBin=i.fBin; return *this; }
      iterator& operator+=(int b) { fBin+=b; return *this; }
      iterator& operator-=(int b) { fBin-=b; return *this; }
      iterator operator+(int b) const { iterator i(*this); i+=b; return i; }
      iterator operator-(int b) const { iterator i(*this); i-=b; return i; }
      int operator-(const iterator& i) { return fBin-i.fBin; }
      bool operator!=(const iterator& o) const { return fBin != o.fBin; }
      bool operator==(const iterator& o) const { return fBin == o.fBin; }
      bool operator<(const iterator& o) const { return fBin < o.fBin; }
      bool operator<=(const iterator& o) const { return fBin <= o.fBin; }
      bool operator>(const iterator& o) const { return fBin > o.fBin; }
      bool operator>=(const iterator& o) const { return fBin >= o.fBin; }
    private:
      friend class VSLimitedErrorsHist;
      iterator(const VSLimitedErrorsHist* hist, int bin): 
	fHist(hist), fBin(bin) { }
      const VSLimitedErrorsHist* fHist;
      int                 fBin;
    };

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------
    VSLimitedErrorsHist(const std::string& name = ""):
      fName(name), fHistContents("contents"), fHistVariance("variance")
    { }

    VSLimitedErrorsHist(const T& binsize, 
			const T& lo_limit, const T& hi_limit,
			const std::string& name = ""):
	fName(name),
	fHistContents(binsize,lo_limit,hi_limit,"contents"), 
	fHistVariance(binsize,lo_limit,hi_limit,"variance")
    { }

    // ------------------------------------------------------------------------
    // Simple accessors
    // ------------------------------------------------------------------------
    const T& binSize() const { return fHistContents.binSize(); }
    const T& loLimit() const { return fHistContents.loLimit(); }
    const T& hiLimit() const { return fHistContents.hiLimit(); }
    const unsigned nBins() const { return fHistContents.nBins(); }
    const count_type& underflow() const { return fHistContents.underflow(); }
    const count_type& overflow() const { return fHistContents.overflow(); } 
    const std::string& name() const{ return fName; }
    inline count_type count(int bin) const { return fHistContents.count(bin); }
    inline count_type countForVal(const T& val) const
    { return fHistContents.countForVal(val); }
    inline count_type err(int bin) const 
    { return sqrt(fHistVariance.count(bin)); }
    inline count_type var(int bin) const 
    { return fHistVariance.count(bin); }
    inline count_type varForVal(const T& val) const 
    { return fHistVariance.countForVal(val); }
    inline T val(int bin) const { return fHistContents.val(bin); }
    inline T center(int bin) const { return fHistContents.center(bin); }
    inline int bin(const T& val) const { return fHistContents.bin(val); }

    T mean() const { return fHistContents.mean(); }
    count_type sum() const { return fHistContents.sum(); }

    // ------------------------------------------------------------------------
    // Iterator
    // ------------------------------------------------------------------------
    inline int firstBin() const { return fHistContents.firstBin(); }
    inline int onePastLastBin() const 
    { return fHistContents.onePastLastBin(); }

    inline unsigned size() const { return fHistContents.size(); }
    inline iterator begin() const;
    inline iterator end() const;

    inline bool empty() const { return fHistContents.empty(); }
        
    // ------------------------------------------------------------------------
    // Enter a value into the histogram
    // ------------------------------------------------------------------------
    void accumulate(const T& x)
    {
      fHistContents.accumulate(x);
      fHistVariance.accumulate(x);
    }

    void accumulate(const T& x, count_type count)
    {
      fHistContents.accumulate(x,count);
      fHistVariance.accumulate(x,count);
    }

    void accumulate(const T& x, count_type count, count_type var)
    {
      fHistContents.accumulate(x,count);
      fHistVariance.accumulate(x,var);
    }

    void setBin(int bin, count_type count, count_type var)
    {
      fHistContents.setBin(bin,count);
      fHistVariance.setBin(bin,var);
    }

    void fill(count_type count = 0, count_type var = 0)
    {
      fHistContents.fill(count);
      fHistVariance.fill(var);
    }

    void set(const std::vector<double>& count)
    {
      vsassert(count.size() == nBins());
      clear();
      fill();
      unsigned n = 0;
      for(int ibin = firstBin(); ibin < onePastLastBin(); ibin++)
	setBin(ibin,count[n++],0);
    }

    void set(const std::vector<double>& count, const std::vector<double>& var)
    {
      vsassert(count.size() == nBins());
      clear();
      fill();
      unsigned n = 0;
      for(int ibin = firstBin(); ibin < onePastLastBin(); ibin++)
	{
	  setBin(ibin,count[n],var[n]);
	  n++;
	}
    }

    // ------------------------------------------------------------------------
    // Operators
    // ------------------------------------------------------------------------
    inline VSLimitedErrorsHist& operator+= (const VSLimitedErrorsHist& o);
    inline VSLimitedErrorsHist& operator-= (const VSLimitedErrorsHist& o);
    inline VSLimitedErrorsHist& operator*= (const VSLimitedErrorsHist& o);
    inline VSLimitedErrorsHist& operator/= (const VSLimitedErrorsHist& o);

    inline VSLimitedErrorsHist operator+ (const VSLimitedErrorsHist& o) const;
    inline VSLimitedErrorsHist operator- (const VSLimitedErrorsHist& o) const;
    inline VSLimitedErrorsHist operator* (const VSLimitedErrorsHist& o) const;
    inline VSLimitedErrorsHist operator/ (const VSLimitedErrorsHist& o) const;

    inline VSLimitedErrorsHist& operator+= (count_type x);
    inline VSLimitedErrorsHist& operator-= (count_type x);
    inline VSLimitedErrorsHist& operator*= (count_type x);

    inline VSLimitedErrorsHist operator+ (count_type x) const;
    inline VSLimitedErrorsHist operator- (count_type x) const;
    inline VSLimitedErrorsHist operator* (count_type x) const;

    // ------------------------------------------------------------------------
    // Utility functions 
    // ------------------------------------------------------------------------
    VSLimitedErrorsHist getNormalizedHist() const;
    VSLimitedErrorsHist getCumulativeHist() const;
    const VSLimitedHist<T,count_type,BIN>& contentsHist() const
    { return fHistContents; }
    const VSLimitedHist<T,count_type,BIN>& varianceHist() const
    { return fHistVariance; }

    VSLimitedHist<T,count_type,BIN>& contentsHist() { return fHistContents; }
    VSLimitedHist<T,count_type,BIN>& varianceHist() { return fHistVariance; }

    void merge(const VSLimitedErrorsHist& o)
    {
      fHistContents.merge(o.contentsHist());
      fHistVariance.merge(o.varianceHist());
    }

    // ------------------------------------------------------------------------
    // Load/Save/Reset Methods
    // ------------------------------------------------------------------------
    void clear() { fHistContents.clear(); fHistVariance.clear(); }
#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

  private:
    std::string                     fName;
    VSLimitedHist<T,count_type,BIN> fHistContents;
    VSLimitedHist<T,count_type,BIN> fHistVariance;
  };

  template<typename T, typename count_type, typename BIN> 
  inline typename VSLimitedErrorsHist<T,count_type,BIN>::iterator 
  VSLimitedErrorsHist<T,count_type,BIN>::begin() const
  {
    return iterator(this,firstBin());
  }

  template<typename T, typename count_type, typename BIN> 
  inline typename VSLimitedErrorsHist<T,count_type,BIN>::iterator 
  VSLimitedErrorsHist<T,count_type,BIN>::end() const
  {
    return iterator(this,onePastLastBin());
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator+= 
  (const VSLimitedErrorsHist& o)
  {
    fHistContents += o.contentsHist();
    fHistVariance += o.varianceHist();
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator-= 
  (const VSLimitedErrorsHist& o)
  {
    fHistContents -= o.contentsHist();
    fHistVariance += o.varianceHist();
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator*= 
  (const VSLimitedErrorsHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      {
	count_type count = countForVal(ibin.center());
	count_type var = varForVal(ibin.center());

	fHistContents.setBin(bin(ibin.center()),count*ibin.count());
	fHistVariance.
	  setBin(bin(ibin.center()),
		 std::pow(count,2)*ibin.var()+std::pow(ibin.count(),2)*var);
      }

    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator/= 
  (const VSLimitedErrorsHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      {
	count_type count = countForVal(ibin.center());
	count_type var = varForVal(ibin.center());

	if(ibin.count() == 0 && count == 0)
	  {
	    fHistContents.setBin(bin(ibin.center()),1);
	    fHistVariance.setBin(bin(ibin.center()),0);
	  }
	else
	  {
	    fHistContents.setBin(bin(ibin.center()),count/ibin.count());
	    fHistVariance.
	      setBin(bin(ibin.center()),
		     std::pow(count/ibin.count(),2)*
		     (std::pow(count,2)/var+
		      std::pow(ibin.count(),2)/ibin.var()));
	  }
      }

    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator+ 
  (const VSLimitedErrorsHist& o) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) += o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator- 
  (const VSLimitedErrorsHist& o) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) -= o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator* 
  (const VSLimitedErrorsHist& o) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) *= o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator/ 
  (const VSLimitedErrorsHist& o) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) /= o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator+= (count_type x)
  {
    fHistContents+=x;
    fHistVariance+=x*x;
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator-= (count_type x)
  {
    fHistContents-=x;
    fHistVariance+=x*x;
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN>& 
  VSLimitedErrorsHist<T,count_type,BIN>::operator*= (count_type x)
  {
    fHistContents*=x;
    fHistVariance*=x*x;
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator+ (count_type x) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) += x);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator- (count_type x) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) -= x);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::operator* (count_type x) const
  {
    return (VSLimitedErrorsHist<T,count_type,BIN>(*this) *= x);
  }

  template<typename T, typename count_type, typename BIN> 
  VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::getNormalizedHist() const
  {
    VSLimitedErrorsHist<T,count_type,BIN> h(*this);
    count_type sum = VSLimitedErrorsHist<T,count_type,BIN>::sum();
    if(sum) h*=1./sum;
    return h;
  }

  template<typename T, typename count_type, typename BIN> 
  VSLimitedErrorsHist<T,count_type,BIN> 
  VSLimitedErrorsHist<T,count_type,BIN>::getCumulativeHist() const
  {
    VSLimitedErrorsHist<T,count_type,BIN> h(*this);
    count_type sum = 0;
    count_type var = 0;
    for(iterator ibin = h.begin(); ibin != h.end(); ibin++) 
      {
	sum += ibin->count();
	var += ibin->var();
	h.setBin(ibin->bin(),sum,var);
      }
    return h;
  }

#ifndef NOHDF5
  template<typename T, typename count_type, typename BIN> 
  bool VSLimitedErrorsHist<T,count_type,BIN>::
  load(VSOctaveH5ReaderStruct* reader)
  {
    if(!reader) return false;
    reader->readString("name",fName);
    if(!fHistContents.load(reader->readStruct("contents"))) return false;
    if(!fHistVariance.load(reader->readStruct("variance"))) return false;
    return true;
  }

  template<typename T, typename count_type, typename BIN> 
  bool VSLimitedErrorsHist<T,count_type,BIN>::
  save(VSOctaveH5WriterStruct* writer) const
  {
    writer->writeString("name",fName);
    fHistContents.save(writer->writeStruct("contents"));
    fHistVariance.save(writer->writeStruct("variance"));
    return true;
  }
#endif

}

#endif // VSSIMPLEERRORSHIST_HPP
