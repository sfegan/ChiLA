//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionMethod3.hpp 

  Implementation of reconstruction 3, 3D simultaneously in real
  space

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/09/2005
*/

#include<VSAReconstructionMethod3.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;
using namespace VERITAS::VSAReconstructionInternals;

Method3::~Method3()
{
  // nothing to see here
}

double Method3::
calculateChi2ForHypothesis(const VSAAlgebra::Vec3D& e, VSAAlgebra::Vec3D& R,
			   VSAAlgebra::Symmetric3D& Y,
			   const ArrayMoments& event, 
			   const std::vector<ScopeImage>& images,
			   const std::vector<ScopeInfo>& scopes,
			   bool minimize_R) const
{
  const unsigned nscope = event.scopes.size();

  Vec3D F;
  double S = 0;
  Y.clear();

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const ScopeMoments& scopeparam(event.scopes[iscope]);
      const ScopeImage& scopeimage(images[iscope]);
      const ScopeInfo& scopeinfo(scopes[iscope]);

      const Orthogonal3D& U_cam_to_event(scopeparam.U_cam_to_event);
      Vec3D e_cam(U_cam_to_event.transBwd(e));

      // Skip this telescope completely
      if(!scopeparam.use_in_reconstruction)continue; 

      const double Wi = scopeparam.Wi;
      const double Ni = scopeparam.Ni;
      const double NiWi = Ni*Wi;
      if(NiWi <= 0)continue;

      Symmetric3D QiNiWi;

      unsigned npixel=scopeimage.pixels.size();
      for(unsigned ipixel=0;ipixel<npixel;ipixel++)
	{
	  const unsigned pixelid = scopeimage.pixels[ipixel].j;

	  const double nij = scopeimage.pixels[ipixel].nij;
	  const Vec2D& pij(scopeinfo.pij[pixelid]);

	  // const Vec3D pij_3d = 
	  //   scopeparam.U_cam_to_event.transFwd(Vec3D(pij_3d.x(),pij_3d.y(),0));	  
	  // const Vec3D eij = pij_3d - sqrt(1-pij_3d.norm2())*ei;

	  // In the "cam" coordinate system
	  // pij is in x-y plane only and ei=(0,0,-1)
	  const Vec3D eij(pij.x(), pij.y(), sqrt(1-pij.norm2()));

	  const double eps_prod = eij*e_cam;
	  const double eps_divi = 1 - eps_prod*eps_prod;

	  Symmetric3D Qij(Symmetric3D::kronecker);
	  if(eps_divi < DBL_EPSILON)
	    {
	      Qij -= Symmetric3D(e_cam);
	    }
	  else
	    {
	      const Vec3D eps1 = (e_cam - eps_prod*eij)/eps_divi;
	      const Vec3D eps2 = (eij - eps_prod*e_cam)/eps_divi;
	      Qij -= Symmetric3D(eps1);
	      Qij -= Symmetric3D(eps2);
	      Qij -= Symmetric3D(eps1,eps2)*eps_prod;
	    }

	  QiNiWi += nij*Wi*Qij;
	}

      QiNiWi = U_cam_to_event.transFwd(QiNiWi);

      const Vec3D& ri(scopeparam.ri);
      const Vec3D Fi = QiNiWi*ri;

      Y += QiNiWi;
      F += Fi;
      S += ri*Fi;
    }

  double chi2 = 0;
  if(minimize_R)
    {
      Eigen3D eigen;
      Y.eigen(eigen);
      
      // Transform into 2D space defined by the two eigen-vectors with
      // the largest eigen-values -- label the coordinates by "l,m,s"
  
      // Basis vectors
      const Vec3D& basis_l(eigen.vec[2]);
      const Vec3D& basis_m(eigen.vec[1]);
      
      // Y is diagonalized in this basis
      const double Y_ll = eigen.val[2];
      const double Y_mm = eigen.val[1];
      
      const double F_l = F*basis_l;
      const double F_m = F*basis_m;
      
      // Solve equations for e in 2D sub-space using this basis by
      // multiplying by Y-inverse. Easy since Y-inverse is:
      // diag(1/Y_ll,1/Y_mm)
      
      const double R_l = F_l/Y_ll;
      const double R_m = F_m/Y_mm;
	  
      chi2 = S - F_l*R_l - F_m*R_m;
      
      // Transform R back to global coordinates
      
      R = R_l*basis_l + R_m*basis_m;
    }
  else
    {
      chi2 = S - 2*F*R + Y.traceProduct(R);
    }

  //std::cout << "m3 " << e << ' ' << R << ' ' << chi2 << std::endl;

  return chi2;
}

double Method3::GridCalculator::operator () (double x, double y) const
{
  VSAAlgebra::Vec3D e;
  VSAAlgebra::Vec3D R;
  VSAAlgebra::Symmetric3D Y;
  return calculate(x,y,e,R,Y,true);
}

double Method3::GridCalculator::
calculate(double x, double y, VSAAlgebra::Vec3D& e, VSAAlgebra::Vec3D& R,
	  VSAAlgebra::Symmetric3D& Y, bool minimize) const
{
  e.set(0,0,-1);
  e.rotate(Vec3D(0,-sqrt(x*x+y*y)/180*M_PI,0));
  e.rotate(Vec3D(0,0,atan2(y,x)));
  return 
    fReconstruction.calculateChi2ForHypothesis(e,R,Y,fEvent,fImages,fScopes,
					       minimize);
}

void Method3::gridSearch(MinimizationResults& results,
			 const ArrayMoments& event,
			 const std::vector<ScopeImage>& images,
			 const std::vector<ScopeInfo>& scopes,
			 const std::vector<GridSearchSpecs>& grid)
{
  NOOPGridSearchElementOutputIterator noop_iterator;
  gridSearch(results,event,images,scopes,grid,noop_iterator);
}
