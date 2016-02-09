//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRingFitter.cpp

  Fit a ring to an image

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/07/2007

  $Id: VSRingFitter.cpp,v 3.3 2007/12/03 00:49:05 sfegan Exp $

*/

#include<vector>

#include<VSRingFitter.hpp>
#include<VSAMinimizationWithDerivative.hpp>
#include<VSSimpleStat.hpp>

using namespace VERITAS;

class VSARingFitterMinimizationFunction
{
public:
  VSARingFitterMinimizationFunction(unsigned _ni,
				    const double* _xi,
				    const double* _yi,
				    const double* _si):
    ni(_ni), zxi(_xi), zyi(_yi), zsi(_si), nsi(_si+_ni), 
    sum_si(), sum_sixi_r(), sum_siyi_r(), sum_siri2_r(),
    sum_sixic(), sum_siyic(), sum_siric2(), sum_siric(), 
    sum_sixic_ric(), sum_siyic_ric(), r0(), count()
  { 
    const double* ixi=zxi;
    const double* iyi=zyi;
    const double* isi=zsi;
    while(isi<nsi)
      {
	const double si     = *isi++;
	sum_si             += si;
	const double xi     = *ixi++;
	const double sixi   = si * xi;
	sum_sixi_r         += sixi;
	const double yi     = *iyi++;
	const double siyi   = si * yi;
	sum_siyi_r         += siyi;
	sum_siri2_r        += sixi*xi + siyi*yi;
      }
  }

  void setXY(const double x, const double y)
  {
    const double sum_six = sum_si * x;
    const double sum_siy = sum_si * y;
    sum_siric     = 0;
    sum_siric2    = 
      sum_siri2_r - 2.0*(sum_sixi_r*x + sum_siyi_r*y) + sum_six*x + sum_siy*y;
    sum_sixic     = sum_sixi_r - sum_six;
    sum_siyic     = sum_siyi_r - sum_siy;
    sum_sixic_ric = 0;
    sum_siyic_ric = 0;

    const double* ixi=zxi;
    const double* iyi=zyi;
    const double* isi=zsi;
    while(isi<nsi)
      {
	const double si     = *(isi++);
	const double xic    = *ixi++ - x;
	const double yic    = *iyi++ - y;
	const double ric    = sqrt(xic*xic + yic*yic);
	const double si_ric = si/ric;
	sum_siric          += si*ric;
	sum_sixic_ric      += si_ric*xic;
	sum_siyic_ric      += si_ric*yic;

#if 0
	std::cout << xic << ' ' << yic << ' ' << ric2 << ' ' << ric << ' '
		  << sum_siric << ' ' << sum_siric2 << ' ' 
		  << sum_sixic << ' ' << sum_siyic << ' ' 
		  << sum_sixic_ric << ' ' << sum_siyic_ric << '\n';
#endif
	
      }
    r0 = sum_siric/sum_si;
    count++;
  }

  void eval(double& f) const
  {
    f = sum_siric2 - (2*sum_siric - sum_si*r0)*r0;
  }

  void deriv(double& dx, double& dy) const
  {
    dx = -2.0 * (sum_sixic - r0*sum_sixic_ric);
    dy = -2.0 * (sum_siyic - r0*sum_siyic_ric);
  }

  double signal() const { return sum_si; }
  double radius() const { return r0; }
  unsigned neval() const { return count; }
  
  ~VSARingFitterMinimizationFunction()
  { 
    // nothing to see here
  }

private:
  VSARingFitterMinimizationFunction(const VSARingFitterMinimizationFunction&);
  VSARingFitterMinimizationFunction& 
  operator=(const VSARingFitterMinimizationFunction&);

  const unsigned ni;
  const double*const zxi;
  const double*const zyi;
  const double*const zsi;
  const double*const nsi;

  double sum_si;
  double sum_sixi_r;
  double sum_siyi_r;
  double sum_siri2_r;

  double sum_sixic;
  double sum_siyic;
  double sum_siric2;
  double sum_siric;
  double sum_sixic_ric;
  double sum_siyic_ric;
  double r0;
  unsigned count;
};

VSRingFitter::
VSRingFitter(const VSAReconstruction::ScopeInfo& scope): 
  m_scope(scope), m_neval()
{
  // nothing to see here
}

VSRingFitter::
~VSRingFitter()
{
  std::cerr << "Ring fitter: " << m_neval << " function evaluations\n";
}

bool VSRingFitter::
fitRing(const VSAReconstruction::ScopeImage& image,
	const double x0, const double y0,
	double& chi2, double& signal, 
	double& xc, double& yc, double& r0,
	double& dx, double& dy) const
{
  unsigned ni = image.pixels.size();
  
  xc = x0;
  yc = y0;

  std::vector<double> xi(ni);
  std::vector<double> yi(ni);
  std::vector<double> si(ni);
  for(unsigned ii=0;ii<ni;ii++)
    {
      unsigned ipixel = image.pixels[ii].j;
      xi[ii] = m_scope.pij[ipixel].x();
      yi[ii] = m_scope.pij[ipixel].y();
      si[ii] = image.pixels[ii].nij;
    }

  VSARingFitterMinimizationFunction 
    function(ni,&(*xi.begin()),&(*yi.begin()),&(*si.begin()));
  signal = function.signal();
  try
    {
      VSAMath::minimizeWithDerivative2(xc, yc, chi2, 1e-4, function, 2e-3);
    }
  catch(const VSAMath::MinimizationFailed& x)
    {
      return false;
    }
    
  function.setXY(xc,yc);
  r0 = function.radius();

  dx = x0 - xc;
  dy = y0 - yc;

  m_neval += function.neval();

  return true;
}
