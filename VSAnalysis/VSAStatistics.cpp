//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAStatistics.cpp
  Various statistics functions

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#include <VSAMath.hpp>
#include <VSAStatistics.hpp>

using namespace VERITAS;

double VSAStatistics::heleneUL(double e, double err, double cl)
{
  double s = sqrt(2.)*err;
  double erf = VSAMath::erf(e/s);
  return e + s*VSAMath::erf_inverse(cl*(1+erf) - erf);
}

double VSAStatistics::limaSignificance(unsigned on, unsigned off, double alpha)
{
  if(on + off == 0) return 0;

  double sigma = 0;
  double a = 0;
  double b = 0;

  if(on != 0) 
    a = (double)on*log((1+alpha)/alpha*((double)on/(double)(on+off)));
  if(off != 0) 
    b = (double)off*log((1+alpha)*((double)off/(double)(on+off)));
    
  sigma = sqrt(2.)*sqrt( a + b );
  
  if((double)on - alpha*(double)off < 0) return -sigma;
  else return sigma;
}

double VSAStatistics::limaSignificance(std::vector< double >& on, 
				       std::vector< double >& off, 
				       std::vector< double >& alpha) 
{
  double x = 0;
  double y = 0;
  double stot = 0;
  double btot = 0;
  double etot = 0;

  for(unsigned i = 0; i < on.size(); i++) 
    {
      stot += on[i];
      btot += off[i];
      etot += on[i] - alpha[i]*off[i];

      x += alpha[i]/(1+alpha[i])*(on[i] + off[i]);
      y += 1./(1+alpha[i])*(on[i] + off[i]);
    }

  double sigma = 0;
  double a = 0;
  double b = 0;

  if(stot != 0) a = stot*log(stot/x);
  if(btot != 0) b = btot*log(btot/y);
    
  sigma = sqrt(2.)*sqrt( a + b );
  
  if(etot < 0) return -sigma;
  else return sigma;
}
