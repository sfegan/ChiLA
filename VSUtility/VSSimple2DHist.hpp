//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleHist.hpp

  Simple histogram class.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/30/2007
*/

#include <vector>

#ifndef VSSIMPLE2DHIST_HPP
#define VSSIMPLE2DHIST_HPP

#include<VSSimpleHist.hpp>

namespace VERITAS
{

  template<typename T, typename count_type=unsigned, 
	   typename BIN=VSBinCalcLinear<T> >
  class VSSimple2DHist
  {
  public:

    class iterator
    {
    public:
      iterator(const iterator& i): fHist(i.fHist), fBin(i.fBin) { }
      count_type count() const { return fHist->countForIndex(fBin); }
      bool val(T&x, T&y) const { return fHist->val(fBin,x,y); }
      bool center(T&x, T&y) const { return fHist->center(fBin,x,y); }
      T x() const { return fHist->binToCenterX(fBin); }
      T y() const { return fHist->binToCenterY(fBin); }
      int bin() const { return fBin; }
      const iterator* operator->() { return this; }
      iterator& operator++() { fBin++; return *this; }
      iterator operator++(int) { iterator i=*this; fBin++; return i; }
      iterator& operator--() { fBin--; return *this; }
      iterator operator--(int) { iterator i=*this; fBin--; return i; }
      iterator& operator=(const iterator& i) { fHist=i.fHist; fBin+=i.fBin; return *this; }
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
      friend class VSSimple2DHist;
      iterator(const VSSimple2DHist* hist, int bin): fHist(hist), fBin(bin) { }
      const VSSimple2DHist* fHist;
      int                   fBin;
    };

    VSSimple2DHist(const std::string& name = ""):
      fXBinCalc(0), fXLoLimit(), fXHiLimit(), fXZero(), fXNBins(), 
      fYBinCalc(0), fYLoLimit(), fYHiLimit(), fYZero(), fYNBins(),
      fNBins(), fData(), fOverflow(), fName(name)
    {
      // nothing to see here
    }

    VSSimple2DHist(const T& xybinsize, const T& xylo,  const T& xyhi,
		   const std::string& name = ""):
      fXBinCalc(xybinsize,xylo), fXLoLimit(xylo), fXHiLimit(xyhi), 
      fXZero(), fXNBins(), 
      fYBinCalc(xybinsize,xylo), fYLoLimit(xylo), fYHiLimit(xyhi), 
      fYZero(), fYNBins(),
      fNBins(), fData(), fOverflow(), fName(name)
    {
      fXZero = fXBinCalc.valToBin(xylo+xybinsize/2,0);
      fYZero = fYBinCalc.valToBin(xylo+xybinsize/2,0);
      fXNBins = fYNBins = lround((xyhi-xylo)/xybinsize);
      fNBins = fXNBins*fYNBins;
      fData = new count_type[fNBins];
      clear();      
    }

    VSSimple2DHist(const T& xbinsize, const T& xlo,  const T& xhi,
		   const T& ybinsize, const T& ylo,  const T& yhi,
		   const std::string& name = ""): 
      fXBinCalc(xbinsize,xlo),  fXLoLimit(xlo), fXHiLimit(xhi), 
      fXZero(), fXNBins(), 
      fYBinCalc(ybinsize,ylo), fYLoLimit(ylo), fYHiLimit(yhi), 
      fYZero(), fYNBins(),
      fNBins(), fData(), fOverflow(), fName(name)
    {
      fXZero = fXBinCalc.valToBin(xlo+xbinsize/2,0);
      fYZero = fYBinCalc.valToBin(ylo+ybinsize/2,0);
      fXNBins = lround((xhi-xlo)/xbinsize);
      fYNBins = lround((yhi-ylo)/ybinsize);
      fNBins = fXNBins*fYNBins;
      fData = new count_type[fNBins];
      clear();      
    }

    ~VSSimple2DHist()
    {
      delete[] fData;
    }
    
    // ------------------------------------------------------------------------
    // Reset
    // ------------------------------------------------------------------------

    void clear()
    {
      for(unsigned iel=0;iel<fNBins;iel++)fData[iel]=0;
      fOverflow=0;
    }

#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

    // ------------------------------------------------------------------------
    // Simple accessors
    // ------------------------------------------------------------------------

    const T& xBinSize() const { return fXBinCalc.binWidth(); }
    const T& xLoLimit() const { return fXLoLimit; }
    const T& xHiLimit() const { return fXHiLimit; }
    const T& yBinSize() const { return fYBinCalc.binWidth(); }
    const T& yLoLimit() const { return fYLoLimit; }
    const T& yHiLimit() const { return fYHiLimit; }
    const unsigned nXBins() const { return fXNBins; }
    const unsigned nYBins() const { return fYNBins; }
    const unsigned nBins() const { return fNBins; }
    const std::string& name() const { return fName; }

    // ------------------------------------------------------------------------
    // Enter a value into the histogram
    // ------------------------------------------------------------------------
    
    void accumulate(const T& x, const T& y)
    {
      int ibin;
      if(!bin(x,y,ibin))fOverflow++;
      else fData[ibin]++;
    }

    void accumulate(const T& x, const T& y, count_type count)
    {
      int ibin;
      if(!bin(x,y,ibin))fOverflow+=count;
      else fData[ibin]+=count;
    }

    void accumulateBin(int ibin)
    {
      if(ibin >= (int)fNBins)fOverflow++;
      else fData[ibin]++;
    }

    void accumulateBin(int ibin, count_type count)
    {
      if(ibin >= (int)fNBins)fOverflow+=count;
      else fData[ibin]+=count;
    }

    void setBin(int ibin, count_type count)
    {
      vsassert(ibin < (int)fNBins);
      fData[ibin] = count;
    }

    void setBin(int ixbin, int iybin, count_type count)
    {
      int ibin = xyBinToIndex(ixbin,iybin);
      vsassert(ibin < (int)fNBins);
      fData[ibin] = count;
    }

    void set(const std::vector<double>& x)
    {
      vsassert(x.size() == fNBins);
      for(unsigned ibin = 0; ibin < fNBins; ibin++) setBin(ibin,x[ibin]);
    }

    // ------------------------------------------------------------------------
    // Operators
    // ------------------------------------------------------------------------
    inline VSSimple2DHist& operator+= (const VSSimple2DHist& o);
    inline VSSimple2DHist& operator-= (const VSSimple2DHist& o);
    inline VSSimple2DHist& operator*= (const VSSimple2DHist& o);

    inline VSSimple2DHist operator+ (const VSSimple2DHist& o) const;
    inline VSSimple2DHist operator- (const VSSimple2DHist& o) const;

    inline VSSimple2DHist& operator+= (count_type x);
    inline VSSimple2DHist& operator-= (count_type x);
    inline VSSimple2DHist& operator*= (count_type x);

    inline VSSimple2DHist operator+ (count_type x) const;
    inline VSSimple2DHist operator- (count_type x) const;
    inline VSSimple2DHist operator* (count_type x) const;

    // ------------------------------------------------------------------------
    // Assignment operator
    // ------------------------------------------------------------------------

    VSSimple2DHist& operator= (const VSSimple2DHist& o)
    {
      if(this == &o)
	return *this;

      delete[] fData;
      fData = new count_type[o.fNBins];
      
      fName     = o.fName;
      fXBinCalc = o.fXBinCalc;
      fXLoLimit = o.fXLoLimit;
      fXHiLimit = o.fXHiLimit;
      fXZero    = o.fXZero;
      fXNBins   = o.fXNBins;
      fYBinCalc = o.fYBinCalc;
      fYLoLimit = o.fYLoLimit;
      fYHiLimit = o.fYHiLimit;
      fYZero    = o.fYZero;
      fYNBins   = o.fYNBins;
      fNBins    = o.fNBins;         

      clear();
  
      fOverflow = o.fOverflow;
      unsigned nbin = o.size();
      for(unsigned ibin=0;ibin<nbin;ibin++)
	fData[ibin] = o.fData[ibin];

      return *this;
    }

    // ------------------------------------------------------------------------
    // Copy Constructor
    // ------------------------------------------------------------------------
    VSSimple2DHist(const VSSimple2DHist& o)
    {
      fData = new count_type[o.fNBins];
      
      fName     = o.fName;
      fXBinCalc = o.fXBinCalc;
      fXLoLimit = o.fXLoLimit;
      fXHiLimit = o.fXHiLimit;
      fXZero    = o.fXZero;
      fXNBins   = o.fXNBins;
      fYBinCalc = o.fYBinCalc;
      fYLoLimit = o.fYLoLimit;
      fYHiLimit = o.fYHiLimit;
      fYZero    = o.fYZero;
      fYNBins   = o.fYNBins;
      fNBins    = o.fNBins;         

      clear();
  
      fOverflow = o.fOverflow;
      unsigned nbin = o.size();
      for(unsigned ibin=0;ibin<nbin;ibin++)
	fData[ibin] = o.fData[ibin];
    }

    // ------------------------------------------------------------------------
    // Counts for bins within the histogram or zero otherwise
    // ------------------------------------------------------------------------

    count_type countForIndex(int bin) const
    {
      if((bin<0)||(bin>=(int)fNBins))return 0;
      return fData[bin];
    }

    count_type count(int ixbin, int iybin) const
    {
      if((ixbin<0)||(iybin<0)||(ixbin>=(int)fXNBins)||(iybin>=(int)fYNBins))
	return 0;
      return fData[xyBinToIndex(ixbin,iybin)];
    }

    count_type countForVal(const T& x, const T& y) const
    {
      int ibin;
      if(!bin(x,y,ibin))return 0;
      return fData[ibin];
    }

    count_type overflow() const { return fOverflow; }

    // ------------------------------------------------------------------------
    // Simple conversions which do not require that the actual bin exists
    // ------------------------------------------------------------------------

    T xBinToVal(int ixbin) const { return fXBinCalc.binToVal(ixbin,fXZero); }
    T yBinToVal(int iybin) const { return fYBinCalc.binToVal(iybin,fYZero); }

    T xBinToCenter(int ixbin) const 
    { return fXBinCalc.binToCenter(ixbin,fXZero); }
    T yBinToCenter(int iybin) const 
    { return fYBinCalc.binToCenter(iybin,fYZero); }

    T binToCenterX(int ibin) const 
    { return fXBinCalc.binToCenter(indexToXBin(ibin),fXZero); }
    T binToCenterY(int ibin) const 
    { return fYBinCalc.binToCenter(indexToYBin(ibin),fYZero); }

    int xValToBin(const T& x) const { return fXBinCalc.valToBin(x,fXZero); }
    int yValToBin(const T& y) const { return fYBinCalc.valToBin(y,fYZero); }

    int indexToXBin(int ibin) const 
    { return (ibin<0)?(fXNBins-abs(ibin)%fXNBins):(ibin%fXNBins); }
    int indexToYBin(int ibin) const 
    { return (ibin<0)?(-(abs(ibin+1)/fXNBins+1)):(ibin/fXNBins); }

    int xyBinToIndex(int ixbin, int iybin) const 
    { return iybin*fXNBins+ixbin; }

    // ------------------------------------------------------------------------
    // Conversions that require the bin to be valid
    // ------------------------------------------------------------------------

    bool val(int ibin, T& x, T& y) const
    {
      if((ibin<0)||(ibin>=(int)fNBins))return false;
      x = xBinToVal(indexToXBin(ibin));
      y = yBinToVal(indexToYBin(ibin));
      return true;
    }

    bool center(int ibin, T& x, T& y) const
    {
      if((ibin<0)||(ibin>=(int)fNBins))return false;
      x = xBinToCenter(indexToXBin(ibin));
      y = yBinToCenter(indexToYBin(ibin));
      return true;
    }

    bool bin(const T& x, const T& y, int& bin) const
    {
      int ixbin = xValToBin(x); 
      if((ixbin<0)||(ixbin>=(int)fXNBins))return false;
      int iybin = yValToBin(y); 
      if((iybin<0)||(iybin>=(int)fYNBins))return false;
      bin = xyBinToIndex(ixbin,iybin);
      return true;
    }

    // ------------------------------------------------------------------------
    // Utility functions
    // ------------------------------------------------------------------------
    
    void merge(const VSSimple2DHist& o)
    {
      if(nBins() == 0) *this = o;
      else if(xLoLimit() > o.xLoLimit() || xHiLimit() < o.xHiLimit() ||
	      yLoLimit() > o.yLoLimit() || yHiLimit() < o.yHiLimit())
	{
	  double xlo = std::min(xLoLimit(),o.xLoLimit());
	  double xhi = std::max(xHiLimit(),o.xHiLimit());
	  double ylo = std::min(yLoLimit(),o.yLoLimit());
	  double yhi = std::max(yHiLimit(),o.yHiLimit());
	  VSSimple2DHist h(xBinSize(),xlo,xhi,yBinSize(),ylo,yhi);
	  h += *this;
	  *this = h;
	  *this += o;
	}
      else *this += o;
    }

    void normalizeX()
    {
      for(unsigned iybin = 0; iybin < fYNBins; iybin++)
	{
	  count_type sum = 0;

	  for(unsigned ixbin = 0; ixbin < fXNBins; ixbin++)
	    sum += count(ixbin,iybin);

	  for(unsigned ixbin = 0; ixbin < fXNBins; ixbin++)
	    setBin(ixbin,iybin,count(ixbin,iybin)/sum);
	}
    }

    void normalizeY()
    {
      for(unsigned ixbin = 0; ixbin < fXNBins; ixbin++)
	{
	  count_type sum = 0;

	  for(unsigned iybin = 0; iybin < fYNBins; iybin++)
	    sum += count(ixbin,iybin);

	  for(unsigned iybin = 0; iybin < fYNBins; iybin++)
	    setBin(ixbin,iybin,count(ixbin,iybin)/sum);
	}
    }

    // X-Projections to 1D Hisograms:
    VSLimitedHist<T, count_type> ProjectX(unsigned ybin_lo, unsigned ybin_hi)
    {
      vsassert(ybin_lo <= ybin_hi && ybin_lo >= 0 && ybin_hi <= fYNBins);

      //fXBinCalc.binWidth();
      VSLimitedHist<T, count_type> 
	projX(fXBinCalc.binWidth(), fXLoLimit, fXHiLimit);

      for(unsigned ixbin = 0; ixbin < fXNBins; ixbin++)
	{
	  count_type xcount = 0;	  
	  for(unsigned iybin = ybin_lo; iybin <= ybin_hi; iybin++)
	    xcount += count(ixbin,iybin);

	  projX.setBin(ixbin, xcount);
	}

      return projX;
    }

    VSLimitedHist<T, count_type> ProjectX(T yval, unsigned num_bins)
    {
      unsigned ybin_mean = fYBinCalc.valToBin(yval,fYZero);
      unsigned ylo = ybin_mean - (num_bins + 1)/2 + 1;
      unsigned yhi = ybin_mean + (num_bins - num_bins/2);
      
      return ProjectX(ylo,yhi);
      
    }

    // Y-Projections to 1D Hisograms:
    VSLimitedHist<T, count_type> ProjectY(unsigned xbin_lo, unsigned xbin_hi)
    {
      vsassert(xbin_lo <= xbin_hi && xbin_lo >= 0 && xbin_hi <= fXNBins);

      //fYBinCalc.binWidth();
      VSLimitedHist<T, count_type> 
	projY(fYBinCalc.binWidth(), fYLoLimit, fYHiLimit);

      for(unsigned iybin = 0; iybin < fYNBins; iybin++)
	{
	  count_type ycount = 0;	  
	  for(unsigned ixbin = xbin_lo; ixbin <= xbin_hi; ixbin++)
	    ycount += count(ixbin,iybin);

	  projY.setBin(iybin, ycount);
	}

      return projY;
    }

    VSLimitedHist<T, count_type> ProjectY(T xval, unsigned num_bins)
    {
      unsigned xbin_mean = fXBinCalc.valToBin(xval,fXZero);
      unsigned xlo = xbin_mean - (num_bins + 1)/2 + 1;
      unsigned xhi = xbin_mean + (num_bins - num_bins/2);
      
      return ProjectY(xlo,xhi);
      
    }


    // ------------------------------------------------------------------------
    // Iterator
    // ------------------------------------------------------------------------

    int firstIndex() const { return 0; }
    int onePastLastIndex() const { return fNBins; }
    
    inline iterator begin() const { return iterator(this, 0); }
    inline iterator end() const { return iterator(this, fNBins); }

    unsigned size() const { return fNBins; }
    inline bool empty() const { return fNBins==0; }
    
  private:
    BIN         fXBinCalc;
    T           fXLoLimit;
    T           fXHiLimit;
    int         fXZero;
    unsigned    fXNBins;
    BIN         fYBinCalc;
    T           fYLoLimit;
    T           fYHiLimit;
    int         fYZero;
    unsigned    fYNBins;
    unsigned    fNBins;
    count_type* fData;
    count_type  fOverflow;
    std::string fName;
  };

  // ==========================================================================
  // Operators
  // ==========================================================================
  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN>& 
  VSSimple2DHist<T,count_type,BIN>::operator+= (const VSSimple2DHist& o)
  {
    unsigned nbin = o.size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      {
	T x;
	T y;
	vsassert(o.val(ibin,x,y));
	accumulate(x+xBinSize()/2,y+yBinSize()/2,o.fData[ibin]);
      }
    fOverflow+=o.fOverflow;
    return *this;
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN>& 
  VSSimple2DHist<T,count_type,BIN>::operator-= (const VSSimple2DHist& o)
  {
    unsigned nbin = o.size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      {
	T x;
	T y;
	vsassert(o.val(ibin,x,y));
	accumulate(x+xBinSize()/2,y+yBinSize()/2,-o.fData[ibin]);
      }
    fOverflow+=o.fOverflow;
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN>& 
  VSSimple2DHist<T,count_type,BIN>::operator*= (const VSSimple2DHist& o)
  {
    unsigned nbin = o.size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      {
	T x = o.binToCenterX(ibin);
	T y = o.binToCenterY(ibin);

	int ibin2;
	if(!bin(x,y,ibin2)) continue;

	count_type c = countForIndex(ibin2)*o.fData[ibin];
	setBin(ibin2,c);
      }
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN> 
  VSSimple2DHist<T,count_type,BIN>::operator+ (const VSSimple2DHist& o) const
  {
    return (VSSimple2DHist<T,count_type,BIN>(*this)+=o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN> 
  VSSimple2DHist<T,count_type,BIN>::operator- (const VSSimple2DHist& o) const
  {
    return (VSSimple2DHist<T,count_type,BIN>(*this)-=o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN>& 
  VSSimple2DHist<T,count_type,BIN>::operator+= (count_type x)
  {
    unsigned nbin = size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      setBin(ibin,fData[ibin]+x);
    
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN>& 
  VSSimple2DHist<T,count_type,BIN>::operator-= (count_type x)
  {
    unsigned nbin = size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      setBin(ibin,fData[ibin]-x);
    
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN>& 
  VSSimple2DHist<T,count_type,BIN>::operator*= (count_type x)
  {
    unsigned nbin = size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      setBin(ibin,fData[ibin]*x);
    
    return *this;
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN> 
  VSSimple2DHist<T,count_type,BIN>::operator+ (count_type x) const
  {
    return (VSSimple2DHist<T,count_type,BIN>(*this) += x);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN> 
  VSSimple2DHist<T,count_type,BIN>::operator- (count_type x) const
  {
    return (VSSimple2DHist<T,count_type,BIN>(*this) -= x);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSSimple2DHist<T,count_type,BIN> 
  VSSimple2DHist<T,count_type,BIN>::operator* (count_type x) const
  {
    return (VSSimple2DHist<T,count_type,BIN>(*this) *= x);
  }

#ifndef NOHDF5
  template<typename T, typename count_type, typename BIN> 
  bool VSSimple2DHist<T,count_type,BIN>::load(VSOctaveH5ReaderStruct* reader)
  {
    if(!reader) return false;
    if(!fXBinCalc.load(reader->readStruct("xbincalc"))) return false;
    if(!fYBinCalc.load(reader->readStruct("ybincalc"))) return false;

    std::vector<T> x;
    std::vector<T> y;
    reader->readString("name",fName);
    reader->readScalar("xlo",fXLoLimit);
    reader->readScalar("xhi",fXHiLimit);
    reader->readScalar("ylo",fYLoLimit);
    reader->readScalar("yhi",fYHiLimit);
    reader->readVector("x",x);
    reader->readVector("y",y);

    fXNBins = x.size();
    fYNBins = y.size();
    fNBins = fXNBins*fYNBins;
    
    if(fNBins > 0)
      {
	fXZero = fXBinCalc.valToBin(x[0],0);
	fYZero = fYBinCalc.valToBin(y[0],0);
	delete[] fData;
	reader->readMatrix("z",fYNBins,fXNBins,fData);
	reader->readScalar("overflow",fOverflow);
      }
    return true;
  }

  template<typename T, typename count_type, typename BIN> 
  bool VSSimple2DHist<T,count_type,BIN>::
  save(VSOctaveH5WriterStruct* writer) const
  {
    writer->writeString("name",fName);
    fXBinCalc.save(writer->writeStruct("xbincalc"));
    fYBinCalc.save(writer->writeStruct("ybincalc"));
    writer->writeScalar("xlo",fXLoLimit);
    writer->writeScalar("xhi",fXHiLimit);
    writer->writeScalar("ylo",fYLoLimit);
    writer->writeScalar("yhi",fYHiLimit);
    
    std::vector<T> x;
    std::vector<T> y;
    count_type* z;

    x.resize(fXNBins);
    for(unsigned ixbin=0;ixbin<fXNBins;ixbin++)
      x[ixbin] = fXBinCalc.binToVal(ixbin,fXZero);
    y.resize(fYNBins);
    for(unsigned iybin=0;iybin<fYNBins;iybin++)
      y[iybin] = fYBinCalc.binToVal(iybin,fYZero);
    z=new count_type[fNBins];
    for(unsigned ibin=0;ibin<fNBins;ibin++)
      z[ibin] = fData[ibin];
    
    writer->writeVector("x",x);
    writer->writeVector("y",y);
    writer->writeMatrix("z",y.size(),x.size(),z);
    writer->writeScalar("overflow",fOverflow);
    delete[] z;
    return true;
  }
#endif

}

template<typename T, typename count_type, typename BIN> 
inline typename VERITAS::VSSimple2DHist<T,count_type,BIN>::iterator
operator+(int b, 
	  const typename VERITAS::VSSimple2DHist<T,count_type,BIN>::iterator& i)
{ 
  typename VERITAS::VSSimple2DHist<T,count_type,BIN>::iterator j(i); 
  j+=b; 
  return j; 
}

template<typename T, typename count_type, typename BIN> 
inline typename VERITAS::VSSimple2DHist<T,count_type,BIN>::iterator
operator-(int b, 
	  const typename VERITAS::VSSimple2DHist<T,count_type,BIN>::iterator& i)
{ 
  typename VERITAS::VSSimple2DHist<T,count_type,BIN>::iterator j(i); 
  j-=b; 
  return j; 
}

#endif // defined VSSIMPLE2DHIST_HPP
