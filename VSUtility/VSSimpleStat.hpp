//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleStat.hpp
  Simple statistics class

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005
*/

#include<vector>
#include<algorithm>

#include<cmath>

#ifndef VSSIMPLESTAT_HPP
#define VSSIMPLESTAT_HPP

namespace VERITAS
{

  // ==========================================================================
  // VSSimpleStat1
  // ==========================================================================

  // Calculator of first moment: mean

  template<typename T, typename count_type = unsigned> class VSSimpleStat1
  {
  public:
    VSSimpleStat1(): fCount(0), fSum(0) { /* nothing to see here */ }

    void accumulate(const T& x) { fCount++; fSum+=x; }
    void accumulate(const T& x, count_type n) { fCount+=n, fSum+=x*T(n); }
    void clear() { fCount=0; fSum=0; }
    unsigned count() const { return fCount; }
    T sum() const { return fSum; }
    T mean() const { return fSum/fCount; }
    VSSimpleStat1& operator+= (const VSSimpleStat1& o);
  private:
    count_type fCount;
    T fSum;
  };

  template<typename T, typename count_type> 
  VERITAS::VSSimpleStat1<T,count_type>& VSSimpleStat1<T,count_type>::
  operator +=(const VSSimpleStat1& o)
  {
    fCount += o.fCount;
    fSum   += o.fSum;
    return *this;
  }

  // ==========================================================================
  // VSSimpleStat2
  // ==========================================================================

  // Calculator of first two moments: mean, variance

  template<typename T, typename count_type = unsigned> class VSSimpleStat2
  {
  public:
    VSSimpleStat2(): fCount(0), fSum(0), fSumSq(0) { /* NTSH */ }
    void accumulate(const T& x) { fCount++; fSum+=x; fSumSq+=x*x; }
    void accumulate(const T& x, count_type n) 
    { fCount+=n; const T xn(x*T(n)); fSum+=xn; fSumSq+=x*xn; }
    void clear() { fCount=0; fSum=0; fSumSq=0; }
    count_type count() const { return fCount; }
    T sum() const { return fSum; }
    T sumsq() const { return fSumSq; }
    T mean() const { return fSum/fCount; }
    T mean_var() const { return var()/fCount; }
    T var() const  { T m=mean(); T v=fSumSq/fCount-m*m; return (v<T())?T():v; }
    T var_var() const  { return 2./(fCount-1)*std::pow(var(),2); }
    T dev() const { return sqrt(var()); }
    T dev_var() const { return 0.5*var()/(fCount-1); }
    T chi2(const T& model) const;
    VSSimpleStat2& operator+= (const VSSimpleStat2& o);
  private:
    count_type fCount;
    T fSum;
    T fSumSq;
  };

  template<typename T, typename count_type> 
  VSSimpleStat2<T,count_type>& VERITAS::VSSimpleStat2<T,count_type>::
  operator +=(const VSSimpleStat2& o)
  {
    fCount += o.fCount;
    fSum   += o.fSum;
    fSumSq += o.fSumSq;
    return *this;
  }

  template<typename T, typename count_type> 
  T VERITAS::VSSimpleStat2<T,count_type>::
  chi2(const T& model) const
  {
    // 1/N sum ( x_i - model )^2 
    // = 1/N ( sum x_i^2   -   2*model*sum x_i   +    model^2 )
    // = sumsq/N   -   2*model*mean   +   model^2

    T m = mean(); 
    T c = fSumSq/fCount - 2*model*m + model*model; 
    return (c<T())?T():c;
  }

  // ==========================================================================
  // VSSimpleStat3
  // ==========================================================================
  
  template<typename T> class VSSimpleStat3
  {
  public:
    VSSimpleStat3(): fCount(0), fSum(0), fSumSq(0), fSumCu(0) { /* NTSH */ }
    void accumulate(const T& x) 
    { fCount++; fSum+=x; const T x2=x*x; fSumSq+=x2; fSumCu+=x2*x; }
    void accumulate(const T& x, unsigned n) 
    { fCount+=n; const T xn(x*T(n)); fSum+=xn; const T xxn(x*xn);
      fSumSq+=xxn; fSumCu+=x*xxn; }
    void clear() { fCount=0; fSum=0; fSumSq=0; fSumCu=0; }
    unsigned count() const { return fCount; }
    T sum() const { return fSum; }
    T sumsq() const { return fSumSq; }
    T mean() const { return fSum/fCount; }
    T var() const { T m=mean(); T v=fSumSq/fCount-m*m; return (v<T())?T():v; }
    T mom3() const 
    { T m=mean(); return (fSumCu-T(3)*m*fSumSq)/fCount+T(2)*m*m*m; }
    T dev() const { return sqrt(var()); }
    T chi2(const T& model) const;
    VSSimpleStat3& operator+= (const VSSimpleStat3& o);
  private:
    unsigned fCount;
    T fSum;
    T fSumSq;
    T fSumCu;
  };

  template<typename T> VSSimpleStat3<T>& VERITAS::VSSimpleStat3<T>::
  operator +=(const VSSimpleStat3& o)
  {
    fCount += o.fCount;
    fSum   += o.fSum;
    fSumSq += o.fSumSq;
    fSumCu += o.fSumCu;
    return *this;
  }

  template<typename T> T VERITAS::VSSimpleStat3<T>::
  chi2(const T& model) const
  {
    // 1/N sum ( x_i - model )^2 
    // = 1/N ( sum x_i^2   -   2*model*sum x_i   +    model^2 )
    // = sumsq/N   -   2*model*mean   +   model^2

    T m = mean(); 
    T c = fSumSq/fCount - 2*model*m + model*model; 
    return (c<T())?T():c;
  }

  // ==========================================================================
  // VSSimpleStatN
  // ==========================================================================

  template<typename T> class VSSimpleStatN
  {
  public:
    inline VSSimpleStatN(unsigned num_moments=3);
    inline VSSimpleStatN(const VSSimpleStatN& o);
    inline VSSimpleStatN& operator= (const VSSimpleStatN& o);
    ~VSSimpleStatN() { delete[] fRawMoments; }
    inline void accumulate(const T& x);
    inline void accumulate(const T& x, unsigned n);
    unsigned count() const { return fCount; }
    T sum() const { return fRawMoments[0]; }
    T sumsq() const { return fRawMoments[1]; }
    T mean() const { return fRawMoments[0]/fCount; }
    T var() const { T m=mean(); return fRawMoments[1]/fCount-m*m; }
    T dev() const { return sqrt(var()); }
    T rawMoment(unsigned m) const { return fRawMoments[m]/fCount; }
    inline T moment(unsigned m) const;
    VSSimpleStatN& operator+= (const VSSimpleStatN& o);
  private:
    unsigned fNumMoments;
    unsigned fCount;
    T*       fRawMoments;
  };
  
  template<typename T> inline VERITAS::VSSimpleStatN<T>::
  VSSimpleStatN(unsigned num_moments): 
    fNumMoments(num_moments), fCount(0), fRawMoments(new T[num_moments])
  {
    for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]=T();
  }
  
  template<typename T> inline VERITAS::VSSimpleStatN<T>::
  VSSimpleStatN(const VSSimpleStatN& o):
    fNumMoments(o.fNumMoments), fCount(o.fCount), 
    fRawMoments(new T[o.fNumMoments])
  {
    for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]=o.fRawMoments[i];
  }
  
  template<typename T> inline VERITAS::VSSimpleStatN<T>& 
  VERITAS::VSSimpleStatN<T>::operator= (const VSSimpleStatN& o)
  {
    if(fNumMoments!=o.fNumMoments)
      {
	delete[] fRawMoments;
	fRawMoments = new T[o.fNumMoments];
      }
    fNumMoments = o.fNumMoments;
    fCount = o.fCount;
    for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]=o.fRawMoments[i];
    return *this;
  }
  
  template<typename T> inline void VERITAS::VSSimpleStatN<T>::
  accumulate(const T& x)
  { 
    T y(x);
    fCount++; 
    for(unsigned i=0; i<fNumMoments; i++)
      { 
	if(i>0)y*=x; 
	fRawMoments[i]+=y;
      }
  }

  template<typename T> inline void VERITAS::VSSimpleStatN<T>::
  accumulate(const T& x, unsigned n)
  { 
    T y(x);
    y*=T(n);
    fCount+=n; 
    for(unsigned i=0; i<fNumMoments; i++)
      { 
	if(i>0)y*=x; 
	fRawMoments[i]+=y;
      }
  }
  
  template<typename T> inline T VERITAS::VSSimpleStatN<T>::
  moment(unsigned m) const
  {
    if(m==0)return 0;
    else if(m==1)return var();
    else if(m==2)
      {
	T m1=rawMoment(0);
	T m2=rawMoment(1);
	T m3=rawMoment(2);
	return m3-3.0*m2*m1+2*m1*m1*m1;
      }
    else
      {
	T sum = T();
	
	sum += fRawMoments[m];
	
	T m1n(fRawMoments[0]);
	for(unsigned i=1; i<m; i++)
	  {
	    double iCn=1;
	    for(unsigned j=0;j<i;j++)iCn *= double(m+1-j)/double(j+1);
	    T sss(iCn*fRawMoments[m-i]*iCn);
	    if(i%2==0)sum += sss;
	    else sum -= sss;
	    m1n*=fRawMoments[0];
	  }
	m1n*=fRawMoments[0];
	
	if(m%2==0)sum += double(m)*m1n;
	else sum -= sum -= double(m)*m1n;

	return sum;
      }
  }

  template<typename T> VERITAS::VSSimpleStatN<T>& VERITAS::VSSimpleStatN<T>::
  operator +=(const VSSimpleStatN& o)
  {
    fCount += o.fCount;
    for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]+=o.fRawMoments[i];
    return *this;
  }

  // ==========================================================================
  // MEDIAN
  // ==========================================================================

  template<typename T> static
  inline T median(const std::vector<std::pair<bool, T> >& data)
  {
    std::vector<T> good_data;
    good_data.reserve(data.size());
    for(typename std::vector<std::pair<bool, T> >::const_iterator 
	  idatum = data.begin(); idatum != data.end(); idatum++)
      if(idatum->first)good_data.push_back(idatum->second);
    unsigned ngood_data = good_data.size();
    if(ngood_data == 0)return T();
    std::sort(good_data.begin(), good_data.end());
    if(ngood_data%2 == 1)return good_data[ngood_data/2];
    else return (good_data[ngood_data/2-1]+good_data[ngood_data/2])/2;
  }

  template<typename T> static
  inline T median(const std::vector<T>& data)
  {
    std::vector<T> good_data(data);
    unsigned ngood_data = good_data.size();
    if(ngood_data == 0)return T();
    std::sort(good_data.begin(), good_data.end());
    if(ngood_data%2 == 1)return good_data[ngood_data/2];
    else return (good_data[ngood_data/2-1]+good_data[ngood_data/2])/2;
  }

  // ==========================================================================
  // CONTAINMENT INTERVAL
  // ==========================================================================

  template<typename T> static
  inline bool containmentInterval(const std::vector<std::pair<bool, T> >& data,
				  double fraction, const double scale,
				  T& lo, T& hi)
  {
    std::vector<T> good_data;
    good_data.reserve(data.size());
    for(typename std::vector<std::pair<bool, T> >::const_iterator 
	  idatum = data.begin(); idatum != data.end(); idatum++)
      if(idatum->first)good_data.push_back(idatum->second);
    unsigned ngood_data = good_data.size();
    if(ngood_data == 0)
      {
	lo = T();
	hi = T();
	return false;
      }
    if(fraction<0.0)fraction=0.0;
    else if(fraction>1.0)fraction=1.0;
    std::sort(good_data.begin(), good_data.end());
    T med;
    if(ngood_data%2 == 1)med = good_data[ngood_data/2];
    else med = (good_data[ngood_data/2-1]+good_data[ngood_data/2])/2;
    unsigned ihi = unsigned(round(double(ngood_data-1)*(1+fraction)/2.0));
    unsigned ilo = unsigned(round(double(ngood_data-1)*(1-fraction)/2.0));
    hi = med + (good_data[ihi] - med)*scale;
    lo = med + (good_data[ilo] - med)*scale;
    return true;
  }
  
  template<typename T> static
  inline bool containmentIntervalLower(const std::vector<T>& data,
				       double fraction, T& lo, T& lo_err)
  {
    std::vector< double > d = data;

    const unsigned ndata = d.size();
    if(ndata == 0) 
      {
	lo = T();
	lo_err = T();
	return false;
      }

    std::sort(d.begin(), d.end());
    unsigned ilo;

    if(fraction<=0.0) ilo = ndata-1;
    else if(fraction>=1.0) ilo = 0;
    else ilo = unsigned(floor(ndata*(1-fraction)));
      
    unsigned err = (unsigned)sqrt(ndata*fraction*(1-fraction));
    
    unsigned ilo_err1 = std::max((int)0,(int)(ilo-err));
    unsigned ilo_err2 = std::min((int)(ndata-1),(int)(ilo+err));

    lo = d[ilo];
    lo_err = (d[ilo_err2] - d[ilo_err1])/2.;

    return true;
  }

  template<typename T> static
  inline bool containmentIntervalUpper(const std::vector<T>& data,
				       double fraction, T& hi, T& hi_err)
  {
    std::vector< double > d = data;

    const unsigned ndata = d.size();
    if(ndata == 0) 
      {
	hi = T();
	hi_err = T();
	return false;
      }

    std::sort(d.begin(), d.end());
    unsigned ihi;

    if(fraction<=0.0) ihi = 0;
    else if(fraction>=1.0) ihi = ndata-1;
    else ihi = unsigned(floor(ndata*fraction));

    unsigned err = (unsigned)sqrt(ndata*fraction*(1-fraction));
    
    unsigned ihi_err1 = std::max((int)0,(int)(ihi-err));
    unsigned ihi_err2 = std::min((int)(ndata-1),(int)(ihi+err));

    hi = d[ihi];  
    hi_err = (d[ihi_err2] - d[ihi_err1])/2.;

    return true;
  }

}

#endif // VSSIMPLESTAT_HPP
