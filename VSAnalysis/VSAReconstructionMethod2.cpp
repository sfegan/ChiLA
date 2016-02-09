//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionMethod2.hpp 

  Implementation of reconstruction 2, 2D simultaneously in focal
  planes of telescope

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/11/2005
*/

#include<VSAReconstructionMethod2.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;
using namespace VERITAS::VSAReconstructionInternals;

Method2::~Method2()
{
  // nothing to see here
}

double Method2::
calculateChi2ForHypothesis(VSAAlgebra::Vec3D& e, const VSAAlgebra::Vec3D& R,
			   VSAAlgebra::Symmetric3D& Y,
			   const ArrayMoments& event, bool minimize_e) const
{
  const unsigned nscope = event.scopes.size();

  Vec3D F;
  double S = 0;
  Y.clear();

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const ScopeMoments& P(event.scopes[iscope]);

      // Skip this telescope completely
      if(!P.use_in_reconstruction)continue; 

      const double NiWi = P.Ni * P.Wi;
      if(NiWi <= 0)continue;

      const Vec3D& ei(P.ei);
      const Vec3D deltari(P.ri-R);

      Symmetric3D Qi;

      if(deltari.norm2() < DBL_EPSILON*DBL_EPSILON)
	{
	  Qi = Symmetric3D::kronecker;
	  Qi -= Symmetric3D(ei);
	}
      else
	{
	  const double ei_deltari = ei*deltari;
	  Vec3D kappai = deltari - ei*ei_deltari;
	  kappai.normalize();
	  Vec3D rhoi = ei^kappai;

	  Qi = Symmetric3D(rhoi);
	}

      Y += NiWi*Qi;
      F += NiWi*Qi*P.Vi;
      S += NiWi*Qi.traceProduct(P.Ti);
    }

  double chi2 = 0;
  if(minimize_e)
    {
      Eigen3D eigen;
      Y.eigen(eigen);
      
      // Transform into 2D space defined by the two eigen-vectors with
      // the largest eigen-values -- label the coordinates by "l,m,s"
  
      // Basis vectors
      Vec3D& basis_l(eigen.vec[2]);
      Vec3D& basis_m(eigen.vec[1]);
      Vec3D& basis_s(eigen.vec[0]);

      if(basis_s.z()<0)
	{
	  // rotate around "l" axis so that "s" axis points "up"
	  basis_s = -basis_s;
	  basis_m = -basis_m;
	}
 
      // Y is diagonalized in this basis
      const double Y_ll = eigen.val[2];
      const double Y_mm = eigen.val[1];
      //const double Y_ss = eigen.val[0]; // Should be small
      //const double Y_lm = 0;
      //const double Y_ls = 0;
      //const double Y_ms = 0;
      
      const double F_l = F*basis_l;
      const double F_m = F*basis_m;
      //const double F_s = F*basis_s;
      
      // Solve equations for e in 2D sub-space using this basis by
      // multiplying by Y-inverse. Easy since Y-inverse is:
      // diag(1/Y_ll,1/Y_mm)
	  
      const double e_l = F_l/Y_ll;
      const double e_m = F_m/Y_mm;
      const double e_s = -sqrt(1-e_l*e_l-e_m*e_m);
      
      chi2 = S - F_l*e_l - F_m*e_m;
      
#if 0
      if(isnan(chi2))
	{
	  std::cerr << __PRETTY_FUNCTION__ << ": chi2 is NAN\n"
		    << S << '\n'
		    << e_l << ' ' << e_m << '\n'
		    << F_l << ' ' << F_m << '\n'
		    << Y_ll << ' ' << Y_mm << '\n';
	}
#endif

      // Transform e back to global coordinates
      
      e = e_l*basis_l + e_m*basis_m + e_s*basis_s;
    }
  else
    {
      chi2 = S - 2*F*e + Y.traceProduct(e);
    }

#if 0
  std::cout << "m2 " << e << ' ' << R << ' ' << chi2 << std::endl;
#endif

  return chi2;
}

double Method2::GridCalculator::operator () (double x, double y) const
{
  VSAAlgebra::Vec3D e;
  VSAAlgebra::Vec3D R;
  VSAAlgebra::Symmetric3D Y;
  return calculate(x,y,e,R,Y,true);
}

double Method2::GridCalculator::
calculate(double x, double y, VSAAlgebra::Vec3D& e, VSAAlgebra::Vec3D& R,
	  VSAAlgebra::Symmetric3D& Y, bool minimize) const
{
  R.set(x,y,0);
  return fReconstruction.calculateChi2ForHypothesis(e,R,Y,fEvent,minimize);
}

void Method2::gridSearch(MinimizationResults& results,
			 const ArrayMoments& event,
			 const std::vector<GridSearchSpecs>& grid)
{
  NOOPGridSearchElementOutputIterator noop_iterator;
  gridSearch(results,event,grid,noop_iterator);
}
