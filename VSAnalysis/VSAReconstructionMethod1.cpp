//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionMethod1.hpp
  Implementation of reconstruction 1, independent telescope images

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/11/2005
*/

#include<VSAReconstructionMethod1.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;
using namespace VERITAS::VSAReconstructionInternals;

Method1::~Method1()
{
  // nothing to see here
}

bool Method1::reconstruct(Results& results, const ArrayMoments& event)
{
  unsigned nscope=event.scopes.size();

  double N = 0;
  Symmetric3D Y;
  Vec3D FR;
  double SR = 0;
  Vec3D Fe;
  double Se = 0;

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const ScopeMoments& P(event.scopes[iscope]);

      // Skip this telescope completely
      if(!P.use_in_reconstruction)continue;

      const double WiNi = P.Ni*P.Wi;
      if(WiNi<=0)continue;

      const Vec3D& ri(P.ri);
      const Vec3D& Vi(P.Vi);
      const Vec3D& rhoi(P.HillasParameters.vec[1]);

      const double wi = WiNi;
      const double ri_rhoi = ri*rhoi;
      const double Vi_rhoi = Vi*rhoi;

      // Weighting
      N += wi;

      // Tensor
      Y += wi*Symmetric3D(rhoi);

      // Vector and scalar for core minimization
      FR += wi*ri_rhoi*rhoi;
      SR += wi*ri_rhoi*ri_rhoi;

      // Vector and scalar for array minimization
      Fe += wi*Vi_rhoi*rhoi;
      Se += wi*Vi_rhoi*Vi_rhoi;      
    }

  results.fWeightedSignal = N;

  if(N == 0)
    {
      results.fCoreLocation.set(0,0,0);
      results.fArrivalDirection.set(0,0,0);
      results.fChi2R = 0;
      results.fChi2e = 0;
      return false;
    }

  Eigen3D& eigen(results.fYabEigen);
  Y.eigen(eigen);

  // Transform into 2D space defined by the two eigen-vectors with the
  // largest eigen-values -- label the coordinates by "l,m,s"
  
  // Basis vectors
  Vec3D basis_l(eigen.vec[2]);
  Vec3D basis_m(eigen.vec[1]);
  Vec3D basis_s(eigen.vec[0]);

  // Check that the s basis vector points to positive elevations in
  // the earth frame
  if(event.U_earth_to_event.transBwd(basis_s).z()<0)
    {
      // rotate around "l" axis so that "s" axis points "up"
      basis_s = -basis_s;
      basis_m = -basis_m;
    }

  // Y is diagonalized in this basis
  const double Y_ll = eigen.val[2];
  const double Y_mm = eigen.val[1];
  //const double Y_ss = eigen.val[0];
  //const double Y_lm = 0;
  //const double Y_ls = 0;
  //const double Y_ms = 0;

  const double FR_l = FR*basis_l;
  const double FR_m = FR*basis_m;
  //const double FR_s = FR*basis_s;
  
  const double Fe_l = Fe*basis_l;
  const double Fe_m = Fe*basis_m;
  //const double Fe_s = Fe*basis_s;
  
#if 0
  std::cerr << "Yll: " << Y_ll << "  Ymm: " << Y_mm << "  Yss: " << Y_ss 
	    << std::endl
	    << "FRl: " << FR_l << "  FRm: " << FR_m << "  FRs: " << FR_s
	    << std::endl
	    << "Fel: " << Fe_l << "  Fem: " << Fe_m << "  Fes: " << Fe_s
	    << std::endl
	    << "SR: " << SR << std::endl
	    << "Se: " << Se << std::endl;
#endif

  // Solve equations for R and e in 2D sub-space using this basis
  // by multiplying by Y-inverse. Easy since Y-inverse is:
  // diag(1/Y_ll,1/Y_mm)
  
  const double R_l = FR_l/Y_ll;
  const double R_m = FR_m/Y_mm;
  //const double R_s = 0;
  
  const double e_l = Fe_l/Y_ll;
  const double e_m = Fe_m/Y_mm;
  const double e_s = -sqrt(1-e_l*e_l-e_m*e_m);
  
  results.fChi2R = SR-R_l*FR_l-R_m*FR_m;
  results.fChi2e = Se-e_l*Fe_l-e_m*Fe_m;
  
#warning what to do here ?? - in 2 telescope case chi^2 is likely exactly zero
  if(results.fChi2R < DBL_EPSILON)results.fChi2R=0;
  if(results.fChi2e < DBL_EPSILON)results.fChi2e=0;

  // Transform R and e back to global coordinates

  results.fCoreLocation     = R_l*basis_l + R_m*basis_m;
  results.fArrivalDirection = e_l*basis_l + e_m*basis_m + e_s*basis_s;

#if 0
  if(!basis_invert)
    {
      results.fCoreLocation     = -results.fCoreLocation;
      results.fArrivalDirection = -results.fArrivalDirection;
    }
#endif

  return true;
}
