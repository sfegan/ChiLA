//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRingFitter.hpp

  Fit a ring to an image

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/07/2007

  $Id: VSRingFitter.hpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#ifndef VSARINGFITTER_HPP
#define VSARINGFITTER_HPP

#include<VSAReconstruction.hpp>

namespace VERITAS
{

  class VSRingFitter
  {
  public:
    VSRingFitter(const VSAReconstruction::ScopeInfo& scope);
    bool fitRing(const VSAReconstruction::ScopeImage& image,
		 const double x0, const double y0,
		 double& chi2, double& signal, 
		 double& xc, double& yc, double& r0,
		 double& dx, double& dy) const;
    ~VSRingFitter();
    unsigned getNEval() const { return m_neval; }
  private:
    const VSAReconstruction::ScopeInfo m_scope;
    mutable unsigned                   m_neval;
  };

} // namespace VERITAS

#endif // defined VSARINGFITTER_HPP
