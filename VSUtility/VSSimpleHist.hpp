//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleHist.hpp

  Simple histogram class.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005
*/

#include <string>
#include <vector>
#include <cmath>
#include <vsassert>

#ifndef NOHDF5
#include <VSOctaveIO.hpp>
#endif

#ifndef VSSIMPLEHIST_HPP
#define VSSIMPLEHIST_HPP

namespace VERITAS
{
  // ==========================================================================
  // VSBinCalcLinear
  // ==========================================================================
  template<typename T> class VSBinCalcLinear
  {
  public:
    typedef struct { } Options;
    static inline Options defOptions() { return Options(); }
    
    VSBinCalcLinear(const T& bin_width = T(), const T& bin_start = T(), 
		    const Options& o = defOptions()):
      m_bin_width(bin_width), m_bin_start(bin_start)
    { }

    inline int valToBin(T val, int zero=0) const
    {
      val -= m_bin_start;
      if(val>=0)return int(val/m_bin_width)-zero;
      else 
	{
	  int offset = int(val/m_bin_width);
	  if((val-m_bin_width*offset)<0)offset-=1;
	  return offset-zero;
	}
    }

    inline T binToVal(int bin, int zero=0) const
    {
      return T(bin+zero)*m_bin_width+m_bin_start;
    }

    inline T binToCenter(int bin, int zero=0) const
    {
      return T(bin+zero)*m_bin_width+m_bin_start + m_bin_width/2;
    }

    inline std::string scheme() const { return "linear"; }

    const T& binWidth() const { return m_bin_width; }
    const T& binStart() const { return m_bin_start; }

#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

  private:
    T m_bin_width;
    T m_bin_start;
  };

#ifndef NOHDF5
  template<typename T> inline bool VSBinCalcLinear<T>::
  load(VSOctaveH5ReaderStruct* reader)
  {
    if(!reader) return false;
    std::string bin_scheme;   
    if((!reader->readString("bin_scheme",bin_scheme))
       ||(bin_scheme == "linear")
       ||(bin_scheme == "left_aligned"))
      {
	if((!reader->readScalar("bin_width",m_bin_width))
	   &&(!reader->readScalar("dx",m_bin_width)))return false;
	if(!reader->readScalar("bin_start",m_bin_start))m_bin_start = T();
      }
    else if(bin_scheme == "center_aligned")
      {
	if(!reader->readScalar("dx",m_bin_width))return false;
	m_bin_start = -m_bin_width/2;
      }
    else return false;
    return true;
  }

  template<typename T> inline bool VSBinCalcLinear<T>::
  save(VSOctaveH5WriterStruct* writer) const
  {
    //writer->writeString("bin_scheme","linear");
    writer->writeScalar("bin_width",m_bin_width);
    if(m_bin_start != T())writer->writeScalar("bin_start",m_bin_start);
    return true;
  }
#endif
  
  // ==========================================================================
  // VSSimpleHist
  // ==========================================================================
  template<typename T, typename count_type=unsigned, 
	   typename BIN=VSBinCalcLinear<T> >
  class VSSimpleHist
  {
  public:
    class iterator
    {
    public:
      iterator(const iterator& i): fHist(i.fHist), fBin(i.fBin) { }
      count_type count() const { return fHist->count(fBin); }
      T val() const { return fHist->val(fBin); }
      T center() const { return fHist->center(fBin); }
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
      friend class VSSimpleHist;
      iterator(const VSSimpleHist* hist, int bin): fHist(hist), fBin(bin) { }
      const VSSimpleHist* fHist;
      int                 fBin;
    };
    
    typedef std::vector<unsigned>::size_type size_type;
    
    // ------------------------------------------------------------------------
    // Constructors 
    // ------------------------------------------------------------------------
    VSSimpleHist(const std::string& name = "",
		 const typename BIN::Options& opt = BIN::defOptions()):
      fBinCalc(0,0,opt), fNegStore(), fPosStore(), fZeroSet(false), fZero(),
      fName(name) { }

    VSSimpleHist(const T& binsize, const std::string& name = "",
		 const typename BIN::Options& opt = BIN::defOptions()):
      fBinCalc(binsize,0,opt), 
      fNegStore(), fPosStore(), fZeroSet(false), fZero(),
      fName(name) { }

    VSSimpleHist(const T& binsize, const T& binstart, 
		 const std::string& name = "",
		 const typename BIN::Options& opt = BIN::defOptions()):
      fBinCalc(binsize,binstart,opt), 
      fNegStore(), fPosStore(), fZeroSet(false), fZero(),
      fName(name) { }

    // ------------------------------------------------------------------------
    // Simple accessors
    // ------------------------------------------------------------------------
    const T& binSize() const { return fBinCalc.binWidth(); }
    const std::string& name() const { return fName; }
    inline count_type count(int bin) const;
    inline count_type countForVal(const T& val) const;

    inline T val(int bin) const;
    inline T center(int bin) const;
    inline int bin(const T& val) const;

    T mean() const;
    count_type sum() const;

    // ------------------------------------------------------------------------
    // Iterator
    // ------------------------------------------------------------------------
    inline int firstBin() const;
    inline int onePastLastBin() const;
    
    inline size_type size() const { return fPosStore.size()+fNegStore.size(); }
    inline iterator begin() const;
    inline iterator end() const;

    inline bool empty() const { return !fZeroSet; }

    // ------------------------------------------------------------------------
    // Enter a value into the histogram
    // ------------------------------------------------------------------------
    inline void accumulate(const T& x);
    inline void accumulate(const T& x, count_type count);
    inline void setBin(int bin, count_type count);

    // ------------------------------------------------------------------------
    // Operators
    // ------------------------------------------------------------------------
    inline VSSimpleHist& operator+= (const VSSimpleHist& o);
    inline VSSimpleHist& operator-= (const VSSimpleHist& o);
    inline VSSimpleHist& operator*= (const VSSimpleHist& o);
    inline VSSimpleHist& operator/= (const VSSimpleHist& o);
    
    inline VSSimpleHist operator+ (const VSSimpleHist& o) const;
    inline VSSimpleHist operator- (const VSSimpleHist& o) const;
    inline VSSimpleHist operator* (const VSSimpleHist& o) const;
    inline VSSimpleHist operator/ (const VSSimpleHist& o) const;

    // ------------------------------------------------------------------------
    // Load/Save/Reset Methods
    // ------------------------------------------------------------------------
    void clear() { fNegStore.clear(); fPosStore.clear(); fZeroSet=false; }
#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

  private:
    BIN                     fBinCalc;
    std::vector<count_type> fNegStore;
    std::vector<count_type> fPosStore;
    bool                    fZeroSet;
    int                     fZero;
    std::string             fName;
  };

#ifndef NOHDF5
  template<typename T, typename count_type, typename BIN> 
  bool VSSimpleHist<T,count_type,BIN>::load(VSOctaveH5ReaderStruct* reader)
  {
    if(!reader) return false;
    if(!fBinCalc.load(reader)) return false;

    std::vector<T> x;
    std::vector<count_type> y;
    if(!reader->readString("name",fName))fName.clear();
    reader->readVector("x",x);
    reader->readVector("y",y);
    const unsigned nbin = x.size();
    for(unsigned ibin=0;ibin<nbin;ibin++)
      accumulate(x[ibin]+fBinCalc.binWidth()/2,y[ibin]);
    return true;
  }

  template<typename T, typename count_type, typename BIN> 
  bool VSSimpleHist<T,count_type,BIN>::
  save(VSOctaveH5WriterStruct* writer) const
  {
    fBinCalc.save(writer);

    std::vector<T> x;
    std::vector<count_type> y;
    
    const unsigned nbin = size();
    x.resize(nbin);
    y.resize(nbin);
    unsigned iibin = 0;
    for(iterator ibin=begin(); ibin!=end();ibin++)
      {
 	x[iibin] = ibin->val();
 	y[iibin] = ibin->count();
 	iibin++;
      }

    if(!fName.empty())writer->writeString("name",fName);
    writer->writeVector("x",x);
    writer->writeVector("y",y);
    return true;
  }
#endif

  template<typename T, typename count_type, typename BIN> 
  inline void VSSimpleHist<T,count_type,BIN>::
  accumulate(const T& x)
  {  
    int bin = 0;
    if(!fZeroSet)
      {
	fZero = fBinCalc.valToBin(x,0);
	fZeroSet = true;
	bin = 0;
      }
    else bin = fBinCalc.valToBin(x,fZero);
    
    if(bin>=0)
      {
	unsigned ubin = unsigned(bin);
	if(fPosStore.size()<=ubin)fPosStore.resize(ubin+1);
	fPosStore[ubin]++;
    }
    else
      {
	unsigned ubin=unsigned(-bin)-1;
	if(fNegStore.size()<=ubin)fNegStore.resize(ubin+1);
	fNegStore[ubin]++;      
      }
  }

  template<typename T, typename count_type, typename BIN> 
  inline void VSSimpleHist<T,count_type,BIN>::
  accumulate(const T& x, count_type count)
  {  
    int bin = 0;
    if(!fZeroSet)
      {
	fZero = fBinCalc.valToBin(x,0);
	fZeroSet = true;
	bin = 0;
      }
    else bin = fBinCalc.valToBin(x,fZero);
    
    if(bin>=0)
      {
	unsigned ubin = unsigned(bin);
	if(fPosStore.size()<=ubin)fPosStore.resize(ubin+1);
	fPosStore[ubin]+=count;
    }
    else
      {
	unsigned ubin=unsigned(-bin)-1;
	if(fNegStore.size()<=ubin)fNegStore.resize(ubin+1);
	fNegStore[ubin]+=count;
      }
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline void VSSimpleHist<T,count_type,BIN>::
  setBin(int bin, count_type count)
  {  
    if(bin>=0)
      {
	unsigned ubin = unsigned(bin);
	if(fPosStore.size()<=ubin)fPosStore.resize(ubin+1);
	fPosStore[ubin] = count;
      }
    else
      {
	unsigned ubin=unsigned(-bin)-1;
	if(fNegStore.size()<=ubin)fNegStore.resize(ubin+1);
	fNegStore[ubin] = count;
      }
  }

  // ==========================================================================
  // Operators
  // ==========================================================================
  template<typename T, typename count_type, typename BIN> 
  inline VSSimpleHist<T,count_type,BIN>& VSSimpleHist<T,count_type,BIN>::
  operator+= (const VSSimpleHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      accumulate(ibin.center(), ibin.count());
    return *this;
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline VSSimpleHist<T,count_type,BIN>& VSSimpleHist<T,count_type,BIN>::
  operator-= (const VSSimpleHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      accumulate(ibin.center(), -ibin.count());
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline count_type VSSimpleHist<T,count_type,BIN>::
  count(int bin) const
  {
    if(bin>=0)
      {
	unsigned ubin = unsigned(bin);
	if(ubin<fPosStore.size())return fPosStore[ubin];
	else return 0;
      }
    else
      {
	unsigned ubin = unsigned((-bin)-1);
	if(ubin<fNegStore.size())return fNegStore[ubin];
	else return 0;
      }
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline count_type VSSimpleHist<T,count_type,BIN>::
  countForVal(const T& val) const
  {
    return count(bin(val));
  }

  template<typename T, typename count_type, typename BIN> 
  inline T VSSimpleHist<T,count_type,BIN>::
  val(int bin) const
  {
    return fBinCalc.binToVal(bin,fZero);
  }

  template<typename T, typename count_type, typename BIN> 
  inline T VSSimpleHist<T,count_type,BIN>::
  center(int bin) const
  {
    return fBinCalc.binToCenter(bin,fZero);
  }

  template<typename T, typename count_type, typename BIN> 
  inline int VSSimpleHist<T,count_type,BIN>::
  bin(const T& val) const
  {
    return fBinCalc.valToBin(val,fZero);
  }

  template<typename T, typename count_type, typename BIN> 
  inline int VSSimpleHist<T,count_type,BIN>::
  firstBin() const
  {
    if(fNegStore.empty())return 0;
    else return -int(fNegStore.size());
  }

  template<typename T, typename count_type, typename BIN> 
  inline int VSSimpleHist<T,count_type,BIN>::
  onePastLastBin() const
  {
    return int(fPosStore.size());
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline typename VSSimpleHist<T,count_type,BIN>::iterator 
  VSSimpleHist<T,count_type,BIN>::
  begin() const
  {
    if(fNegStore.empty())return iterator(this,0);
    else return iterator(this,-int(fNegStore.size()));
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline typename VSSimpleHist<T,count_type,BIN>::iterator 
  VSSimpleHist<T,count_type,BIN>::
  end() const
  {
    return iterator(this,fPosStore.size());
  }

  template<typename T, typename count_type, typename BIN> 
  inline T VSSimpleHist<T,count_type,BIN>::
  mean() const
  {
    T sum(0);
    unsigned count(0);
    for(iterator ibin=begin(); ibin!=end(); ibin++)
      sum += ibin->val()*ibin->count(), count += ibin->count();
    if(count)return sum/count;
    else return sum;
  }

  template<typename T, typename count_type, typename BIN> 
  inline count_type VSSimpleHist<T,count_type,BIN>::sum() const
  {
    count_type count(0);
    for(iterator ibin=begin(); ibin!=end(); ibin++)
      count += ibin->count();
    return count;
  }

  // ==========================================================================
  // VSLimitedHist
  // ==========================================================================
  template<typename T, typename count_type, typename BIN> 
  class VSLimitedHist;

  template<typename T, typename count_type, typename BIN> 
  VSLimitedHist<T,count_type,BIN> operator+ 
  (count_type x, const VSLimitedHist<T,count_type,BIN> &o);
  template<typename T, typename count_type, typename BIN> 
  VSLimitedHist<T,count_type,BIN> operator- 
  (count_type x, const VSLimitedHist<T,count_type,BIN> &o);
  template<typename T, typename count_type, typename BIN> 
  VSLimitedHist<T,count_type,BIN> operator* 
  (count_type x, const VSLimitedHist<T,count_type,BIN> &o);
  template<typename T, typename count_type, typename BIN> 
  VSLimitedHist<T,count_type,BIN> operator/ 
  (count_type x, const VSLimitedHist<T,count_type,BIN> &o);

  template<typename T, typename count_type=unsigned, 
	   typename BIN=VSBinCalcLinear<T> >
  class VSLimitedHist: protected VSSimpleHist<T,count_type,BIN>
  {
  private:
    typedef VSSimpleHist<T,count_type,BIN> BASE;

  public:
    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------
    VSLimitedHist(const std::string& name = ""):
      VSSimpleHist<T,count_type,BIN>(name),
      fLoLimit(), fHiLimit(), fNBins(), fUnderflow(), fOverflow() { }

    VSLimitedHist(const T& binsize, const T& lo_limit, const T& hi_limit,
		  bool auto_fill = false, const std::string& name = ""):
      VSSimpleHist<T,count_type,BIN>(binsize,lo_limit,name), 
      fLoLimit(lo_limit), fHiLimit(), fNBins(),
      fUnderflow(), fOverflow() 
    { 
      fNBins = lround((hi_limit-lo_limit)/binsize);
      fHiLimit = fLoLimit + binsize*(T)(fNBins);
      if(auto_fill) fill(0);
    }

    VSLimitedHist(const T& binsize, const T& lo_limit, const T& hi_limit,
		  const std::vector<double>& v,
		  bool auto_fill = false, const std::string& name = ""):
      VSSimpleHist<T,count_type,BIN>(binsize,lo_limit,name), 
      fLoLimit(lo_limit), fHiLimit(), fNBins(),
      fUnderflow(), fOverflow() 
    { 
      fNBins = lround((hi_limit-lo_limit)/binsize);
      fHiLimit = fLoLimit + binsize*(T)(fNBins);
      vsassert(v.size() == fNBins);
      fill(0);
      for(unsigned ibin = 0; ibin < fNBins; ibin++) setBin(ibin,v[ibin]);
    }

    typedef typename BASE::iterator iterator;
    typedef typename BASE::size_type size_type;

    // ------------------------------------------------------------------------
    // Simple accessors
    // ------------------------------------------------------------------------
    const T& binSize() const { return BASE::binSize(); }
    const T& loLimit() const { return fLoLimit; }
    const T& hiLimit() const { return fHiLimit; }
    unsigned nBins() const { return fNBins; }
    const count_type& underflow() const { return fUnderflow; }
    const count_type& overflow() const { return fOverflow; }
    const std::string& name() const{ return BASE::name(); }
    inline count_type count(int bin) const { return BASE::count(bin); }
    inline count_type countForVal(const T& val) const
    { return BASE::countForVal(val); }

    inline T val(int bin) const { return BASE::val(bin); }
    inline T center(int bin) const { return BASE::center(bin); }
    inline int bin(const T& val) const { return BASE::bin(val); }

    T mean() const { return BASE::mean(); }
    count_type sum() const { return BASE::sum(); }

    // ------------------------------------------------------------------------
    // Iterator
    // ------------------------------------------------------------------------
    inline int firstBin() const { return BASE::firstBin(); }
    inline int onePastLastBin() const { return BASE::onePastLastBin(); }
    
    inline size_type size() const { return BASE::size(); }
    inline iterator begin() const { return BASE::begin(); }
    inline iterator end() const { return BASE::end(); }

    inline bool empty() const { return BASE::empty(); }

    // ------------------------------------------------------------------------
    // Enter a value into the histogram
    // ------------------------------------------------------------------------
    void accumulate(const T& x)
    {
      if(!std::isfinite(x)) return;
      if(x>fHiLimit)fOverflow++;
      else if(x<fLoLimit)fUnderflow++;
      else BASE::accumulate(x);
    }

    void accumulate(const T& x, count_type count)
    {
      if(!std::isfinite(x) || !std::isfinite(count)) return;
      if(x>fHiLimit)fOverflow+=count;
      else if(x<fLoLimit)fUnderflow+=count;
      else BASE::accumulate(x,count);
    }

    void setBin(int bin, count_type count) { BASE::setBin(bin,count); }

    void fill(count_type count = 0)
    {
      T x = loLimit() + binSize()/(T)2;

      while(x < hiLimit())
	{
	  accumulate(x,count);
	  x += binSize();
	}
    }

    // void set(const std::vector<double>& x)
    // {
    //   fill();
    //   unsigned n = 0;
    //   for(unsigned ibin = firstBin(); ibin < onePastLastBin(); ibin++)
    // 	setBin(bin,x[n++]);
    // }

    // ------------------------------------------------------------------------
    // Operators
    // ------------------------------------------------------------------------
    inline VSLimitedHist& operator+= (const VSLimitedHist& o);
    inline VSLimitedHist& operator-= (const VSLimitedHist& o);
    inline VSLimitedHist& operator*= (const VSLimitedHist& o);
    inline VSLimitedHist& operator/= (const VSLimitedHist& o);
    
    inline VSLimitedHist operator+ (const VSLimitedHist& o) const;
    inline VSLimitedHist operator- (const VSLimitedHist& o) const;
    inline VSLimitedHist operator* (const VSLimitedHist& o) const;
    inline VSLimitedHist operator/ (const VSLimitedHist& o) const;
    
    inline VSLimitedHist& operator+= (count_type x);
    inline VSLimitedHist& operator-= (count_type x);
    inline VSLimitedHist& operator*= (count_type x);
    inline VSLimitedHist& operator/= (count_type x);

    inline VSLimitedHist operator+ (count_type x) const;
    inline VSLimitedHist operator- (count_type x) const;
    inline VSLimitedHist operator* (count_type x) const;
    inline VSLimitedHist operator/ (count_type x) const;

    friend VSLimitedHist<T,count_type,BIN> 
    operator+ (count_type x, const VSLimitedHist<T,count_type,BIN> &i)
    {
      return i+x;
    }

    friend VSLimitedHist<T,count_type,BIN> 
    operator- (count_type x, const VSLimitedHist<T,count_type,BIN> &i)
    {
      VSLimitedHist<T,count_type,BIN> o(i);
      for(typename VSLimitedHist<T,count_type,BIN>::iterator ibin = o.begin(); 
	  ibin!=o.end(); ibin++)
	o.setBin(ibin.bin(),x-ibin.count());
      return o;
    }

    // ------------------------------------------------------------------------
    // Utility functions 
    // ------------------------------------------------------------------------
    VSLimitedHist getNormalizedHist();
    VSLimitedHist getCumulativeHist();

    void merge(const VSLimitedHist& o)
    {
      if(this == &o) return;
      else if(nBins() == 0) *this = o;
      else if(loLimit() > o.loLimit() || hiLimit() < o.hiLimit())
	{
	  T lo = std::min(loLimit(),o.loLimit());
	  T hi = std::max(hiLimit(),o.hiLimit());
	  VSLimitedHist h(o.binSize(),lo,hi);
	  h += *this;
	  h += o;
	  *this = h;
	}
      else *this += o;
    }

    // ------------------------------------------------------------------------
    // Load/Save/Reset Methods
    // ------------------------------------------------------------------------
    void clear() 
    { 
      BASE::clear(); 
      fUnderflow=count_type(); fOverflow=count_type(); 
    }
#ifndef NOHDF5
    bool load(VSOctaveH5ReaderStruct* reader);
    bool save(VSOctaveH5WriterStruct* writer) const;    
#endif

  private:
    T                       fLoLimit;
    T                       fHiLimit;
    unsigned                fNBins;
    count_type              fUnderflow;
    count_type              fOverflow;
  };

  // ==========================================================================
  // Operators
  // ==========================================================================
  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& VSLimitedHist<T,count_type,BIN>::
  operator+= (const VSLimitedHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      accumulate(ibin.center(), ibin.count());
    fUnderflow += o.fUnderflow; // what else can we do with these?
    fOverflow += o.fOverflow;
    return *this;
  }
  
  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& VSLimitedHist<T,count_type,BIN>::
  operator-= (const VSLimitedHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      accumulate(ibin.center(), -ibin.count());
    fUnderflow -= o.fUnderflow; // what else can we do with these?
    fOverflow -= o.fOverflow;
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& VSLimitedHist<T,count_type,BIN>::
  operator*= (const VSLimitedHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      setBin(bin(ibin.center()), countForVal(ibin.center())*ibin.count());
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& VSLimitedHist<T,count_type,BIN>::
  operator/= (const VSLimitedHist& o)
  {
    for(iterator ibin = o.begin(); ibin!=o.end(); ibin++)
      {
	vsassert(ibin.count()>0 || countForVal(ibin.center()) == 0);
	if(ibin.count() == 0 && countForVal(ibin.center()) == 0)
	  setBin(bin(ibin.center()), 1);
	else
	  setBin(bin(ibin.center()), 
		 countForVal(ibin.center())/ibin.count());
      }
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN> 
  VSLimitedHist<T,count_type,BIN>::operator+ 
  (const VSLimitedHist& o) const
  {
    return (VSLimitedHist<T,count_type,BIN>(*this) += o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN> 
  VSLimitedHist<T,count_type,BIN>::operator- 
  (const VSLimitedHist& o) const
  {
    return (VSLimitedHist<T,count_type,BIN>(*this) -= o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN> 
  VSLimitedHist<T,count_type,BIN>::operator* 
  (const VSLimitedHist& o) const
  {
    return (VSLimitedHist<T,count_type,BIN>(*this) *= o);
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN> 
  VSLimitedHist<T,count_type,BIN>::operator/ 
  (const VSLimitedHist& o) const
  {
    return (VSLimitedHist<T,count_type,BIN>(*this) /= o);
  }

    template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& 
  VSLimitedHist<T,count_type,BIN>::operator+= (count_type x)
  {
    for(iterator ibin = begin(); ibin!=end(); ibin++)
      setBin(ibin.bin(),ibin.count()+x);
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& 
  VSLimitedHist<T,count_type,BIN>::operator-= (count_type x)
  {
    for(iterator ibin = begin(); ibin!=end(); ibin++)
      setBin(ibin.bin(),ibin.count()-x);
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& 
  VSLimitedHist<T,count_type,BIN>::operator*= (count_type x)
  {
    for(iterator ibin = begin(); ibin!=end(); ibin++)
      setBin(ibin.bin(),ibin.count()*x);
    return *this;
  }

  template<typename T, typename count_type, typename BIN> 
  inline VSLimitedHist<T,count_type,BIN>& 
  VSLimitedHist<T,count_type,BIN>::operator/= (count_type x)
  {
    for(iterator ibin = begin(); ibin!=end(); ibin++)
      setBin(ibin.bin(),ibin.count()/x);
    return *this;
  }

#ifndef NOHDF5
  template<typename T, typename count_type, typename BIN> 
  bool VSLimitedHist<T,count_type,BIN>::load(VSOctaveH5ReaderStruct* reader)
  {
    if(!reader) return false;
    if(!VSSimpleHist<T,count_type,BIN>::load(reader))
      return false;
    
    if(reader->isValid("lo_limit"))
      {
	reader->readScalar("lo_limit",fLoLimit);
	reader->readScalar("hi_limit",fHiLimit);
	if(!reader->readScalar("underflow",fUnderflow))fUnderflow=0;
	if(!reader->readScalar("overflow",fOverflow))fOverflow=0;
      }
    else
      {
	fLoLimit = val(firstBin());
	fHiLimit = val(onePastLastBin());
	fUnderflow = count_type(0);
	fOverflow = count_type(0);
      }

    fNBins = lround((fHiLimit-fLoLimit)/binSize());
    return true;
  }

  template<typename T, typename count_type, typename BIN> 
  bool VSLimitedHist<T,count_type,BIN>::
  save(VSOctaveH5WriterStruct* writer) const
  {
    if(!VSSimpleHist<T,count_type,BIN>::save(writer))
      return false;

    writer->writeScalar("lo_limit",fLoLimit);
    writer->writeScalar("hi_limit",fHiLimit);
    if(fUnderflow)writer->writeScalar("underflow",fUnderflow);
    if(fOverflow)writer->writeScalar("overflow",fOverflow);
    return true;
  }
#endif

}

template<typename T, typename count_type, typename BIN> 
inline typename VERITAS::VSSimpleHist<T,count_type,BIN>::iterator
operator+(int b, 
	  const typename VERITAS::VSSimpleHist<T,count_type,BIN>::iterator& i)
{ 
  typename VERITAS::VSSimpleHist<T,count_type,BIN>::iterator j(i); 
  j+=b; 
  return j; 
}

template<typename T, typename count_type, typename BIN> 
inline typename VERITAS::VSSimpleHist<T,count_type,BIN>::iterator
operator-(int b, 
	  const typename VERITAS::VSSimpleHist<T,count_type,BIN>::iterator& i)
{ 
  typename VERITAS::VSSimpleHist<T,count_type,BIN>::iterator j(i); 
  j-=b; 
  return j; 
}

#endif // VSSIMPLEHIST_HPP
