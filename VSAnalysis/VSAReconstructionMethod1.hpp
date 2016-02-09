//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionMethod1.hpp
  Implementation of reconstruction 1, independent telescope images

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/11/2005
*/

#ifndef VSARECONSTRUCTIONMETHOD1_HPP
#define VSARECONSTRUCTIONMETHOD1_HPP

#define VSA_ASSERT

#include<vector>

#include<VSAAlgebra.hpp>
#include<VSAReconstructionCommon.hpp>

namespace VERITAS 
{
  namespace VSAReconstructionInternals
  {

    class Method1
    {
    public:
      Method1() { /* nothing to see here */ }
      virtual ~Method1();
      
      struct Results
      {
	Results():
	  fCoreLocation(), fArrivalDirection(),
	  fWeightedSignal(), fChi2R(), fChi2e(), fYabEigen()
	{ /* nothing to see here */ }

	VSAAlgebra::Vec3D                  fCoreLocation;
	VSAAlgebra::Vec3D                  fArrivalDirection;
	double                             fWeightedSignal;
	double                             fChi2R;
	double                             fChi2e;
	VSAAlgebra::Eigen3D                fYabEigen;
      };
      
      bool reconstruct(Results& results, const ArrayMoments& event);
      
    private:
      
    };
    
  } // namespace VSAReconstructionInternals
} // namespace VERITAS

#endif // VSARECONSTRUCTIONMETHOD1_HPP
