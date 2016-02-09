//-*-mode:c++; mode:font-lock;-*-

/*! \file VASimpleStat.hpp
  Simple statistice class

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005
*/

#ifndef VASIMPLESTAT_HPP
#define VASIMPLESTAT_HPP

#if 1

namespace VERITAS
{

  template<typename T> class VASimpleStat
  {
  public:
    VASimpleStat(): fCount(0), fSum(0), fSumSq(0) { /* nothing */ }
    void accumulate(const T& x) 
    {
      fCount++; 
      fSum+=x; 
      fSumSq+=x*x; 
#ifdef VASIMPLESTAT_MOM3
      fSumCu+=x*x*x;
#endif
    }
    void clear() { fCount=0; fSum=0; fSumSq=0; }
    unsigned count() const { return fCount; }
    T sum() const { return fSum; }
    T sumsq() const { return fSumSq; }
    T mean() const { return fSum/fCount; }
    T var() const { T m=mean(); return fSumSq/fCount-m*m; }
#ifdef VASIMPLESTAT_MOM3
    T mom3() const { T m=mean(); return (fSumCu-T(3)*m*fSumSq)/fCount+T(2)*m*m*m; }
#endif
    T dev() const { return sqrt(var()); }
    VASimpleStat& operator+= (const VASimpleStat& o);
  private:
    unsigned fCount;
    T fSum;
    T fSumSq;
#ifdef VASIMPLESTAT_MOM3
    T fSumCu;
#endif
  };

}  
 
template<typename T> VERITAS::VASimpleStat<T>& VERITAS::VASimpleStat<T>::
operator +=(const VASimpleStat& o)
{
  fCount += o.fCount;
  fSum   += o.fSum;
  fSumSq += o.fSumSq;
  return *this;
}

#else

namespace VERITAS
{

  template<typename T> class VASimpleStat
  {
  public:
    //inline VASimpleStat(): fNumMoments(), fCount(), fRawMoments() { }
    inline VASimpleStat(unsigned num_moments=3);
    inline VASimpleStat(const VASimpleStat& o);
    inline VASimpleStat& operator= (const VASimpleStat& o);
    ~VASimpleStat() { delete[] fRawMoments; }
    inline void accumulate(const T& x);
    unsigned count() const { return fCount; }
    T sum() const { return fRawMoments[0]; }
    T sumsq() const { return fRawMoments[1]; }
    T mean() const { return fRawMoments[0]/fCount; }
    T var() const { T m=mean(); return fRawMoments[1]/fCount-m*m; }
    T dev() const { return sqrt(var()); }
    T rawMoment(unsigned m) const { return fRawMoments[m]/fCount; }
    inline T moment(unsigned m) const;
    VASimpleStat& operator+= (const VASimpleStat& o);
  private:
    unsigned fNumMoments;
    unsigned fCount;
    T*       fRawMoments;
  };
  
}

template<typename T> inline VERITAS::VASimpleStat<T>::
VASimpleStat(unsigned num_moments): 
  fNumMoments(num_moments), fCount(0), fRawMoments(new T[num_moments])
{
  for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]=T();
}

template<typename T> inline VERITAS::VASimpleStat<T>::
VASimpleStat(const VASimpleStat& o):
  fNumMoments(o.fNumMoments), fCount(o.fCount), 
  fRawMoments(new T[o.fNumMoments])
{
  for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]=o.fRawMoments[i];
}

template<typename T> inline VERITAS::VASimpleStat<T>& 
VERITAS::VASimpleStat<T>::operator= (const VASimpleStat& o)
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

template<typename T> inline void VERITAS::VASimpleStat<T>::
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

template<typename T> inline T VERITAS::VASimpleStat<T>::
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

template<typename T> VERITAS::VASimpleStat<T>& VERITAS::VASimpleStat<T>::
operator +=(const VASimpleStat& o)
{
  fCount += o.fCount;
  for(unsigned i=0; i<fNumMoments; i++)fRawMoments[i]+=o.fRawMoments[i];
  return *this;
}

#endif

#endif // VASIMPLESTAT_HPP
