//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAAtmosphere.hpp
  Integrated atmospheric thickness profile

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/20/2006
*/

#ifndef VSAATMOSPHERE_HPP
#define VSAATMOSPHERE_HPP

#define VSA_ASSERT

#include<iostream>
#include<string>
#include<vector>

namespace VERITAS 
{

  class VSAAtmosphere
  {
  public:
    VSAAtmosphere(): fGood(false), fAtmo(), fAtmoLn()
    { /* nothing to see here */ }
    VSAAtmosphere(const std::string& filename);
    void interpolate(const double& height, double& rho, double& thickness, 
		     double& n_minus_one, double& delay) const;
    double rho(const double& height) const;
    double thickness(const double& height) const;
    double nMinusOne(const double& height) const;
    double delay(const double& height) const;
    
    bool good() const { return fGood; }

    static VSAAtmosphere usStandard();
    
  private:
    void integrate();

    struct AtmoSlice
    {
      AtmoSlice(double h=0, double r=0, double t=0, double n=0, double d=0):
	fHeight(h), fRho(r), fThickness(t), fNMinusOne(n), fDelay(d)
      { /* nothing to see here */ }
      double            fHeight;
      double            fRho;
      double            fThickness;
      double            fNMinusOne;
      double            fDelay;
      bool operator< (const AtmoSlice& o) const { return fHeight<o.fHeight; }
    };
    
    bool fGood;
    std::vector<AtmoSlice> fAtmo;
    std::vector<AtmoSlice> fAtmoLn;
  };
  
}


#endif // VSAATMOSPHERE_HPP
