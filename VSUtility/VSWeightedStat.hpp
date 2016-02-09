//-*-mode:c++; mode:font-lock;-*-

/*! \file VSWeightedStat.hpp
  Weighted statistics class

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       09/16/2005
*/

#ifndef VSWEIGHTEDSTAT_HPP
#define VSWEIGHTEDSTAT_HPP

namespace VERITAS
{

  template<typename T> class VSWeightedStat
  {
  public:
    VSWeightedStat(): 
      fNorm(0), fSum(0), fSumSq(0) 
#ifdef VSWEIGHTEDSTAT_MOM3
      , fSumCu(0)
#endif
    {
      /* nothing to see here */
    }
    
    void accumulate(const T& x, const T& w=1.0) 
    {
      fNorm+=w;
      fSum+=w*x; 
      fSumSq+=w*x*x; 
#ifdef VSWEIGHTEDSTAT_MOM3
      fSumCu+=w*x*x*x;
#endif
    }

    void clear() 
    {
      fNorm=0;
      fSum=0;
      fSumSq=0;
#ifdef VSWEIGHTEDSTAT_MOM3
      fSumCu=0;
#endif
    }

    unsigned norm() const { return fNorm; }
    T sum() const { return fSum; }
    T sumsq() const { return fSumSq; }
    T mean() const { return fSum/fNorm; }
    T var() const { T m=mean(); return fSumSq/fNorm-m*m; }
#ifdef VSWEIGHTEDSTAT_MOM3
    T mom3() const 
    { 
      T m=mean(); 
      return (fSumCu-T(3)*m*fSumSq)/fNorm+T(2)*m*m*m; 
    }
#endif
    T dev() const { return sqrt(var()); }
    VSWeightedStat& operator+= (const VSWeightedStat& o);
  private:
    T fNorm;
    T fSum;
    T fSumSq;
#ifdef VSWEIGHTEDSTAT_MOM3
    T fSumCu;
#endif
  };

}  
 
template<typename T> VERITAS::VSWeightedStat<T>& VERITAS::VSWeightedStat<T>::
operator +=(const VSWeightedStat& o)
{
  fNorm += o.fNorm;
  fSum   += o.fSum;
  fSumSq += o.fSumSq;
#ifdef VSWEIGHTEDSTAT_MOM3
  fSumCu += o.fSumCu;
#endif
  return *this;
}

#endif // VSWEIGHTEDSTAT_HPP
