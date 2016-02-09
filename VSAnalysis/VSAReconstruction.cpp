//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstruction.cpp
  Common interface to various reconstruction methods

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/12/2005
*/

//#define DEBUG

#include<VSAReconstruction.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;
using namespace VERITAS::VSAReconstructionInternals;

VSAReconstruction::~VSAReconstruction()
{
  // nothing to se ehere
}

bool VSAReconstruction::
reconstruct(Reconstruction& recon, 
	    Method method, ScopeWeighting weighting,
	    const std::vector<ScopeImage>& images) const
{
  bool good=false;

  calculateArrayMoments(recon.moments,weighting,images);

  if(recon.moments.nuse_in_reconstruction
     < fQualityCuts.nuse_in_reconstruction_min)return false;

  switch(method)
    {
    case M_1:
      {
	Method1::Results m1results;

	good = method1()->reconstruct(m1results, recon.moments);

	recon.e       = 
	  recon.moments.U_earth_to_event.transBwd(m1results.fArrivalDirection);
	recon.R       = 
	  recon.moments.U_earth_to_event.transBwd(m1results.fCoreLocation);
	recon.chi2e   = m1results.fChi2e;
	recon.chi2R   = m1results.fChi2R;
	recon.deltael = sqrt(m1results.fChi2e/m1results.fYabEigen.val[1]);
	recon.deltaew = sqrt(m1results.fChi2e/m1results.fYabEigen.val[2]);
	recon.deltaRl = sqrt(m1results.fChi2R/m1results.fYabEigen.val[1]);
	recon.deltaRw = sqrt(m1results.fChi2R/m1results.fYabEigen.val[2]);
      }
      break;

    case M_2:
    case M_2_GRIDSEARCH:
      {
	Method2::MinimizationResults m2results;
	std::vector<GridSearchSpecs> grid;
	grid.push_back(GridSearchSpecs(VSAR_M2_GRID*1.000,VSAR_M2_GRID*0.100));
	grid.push_back(GridSearchSpecs(VSAR_M2_GRID*0.100,VSAR_M2_GRID*0.010));
	grid.push_back(GridSearchSpecs(VSAR_M2_GRID*0.010,VSAR_M2_GRID*0.001));

	method2()->gridSearch(m2results, recon.moments, grid);

#ifdef DEBUG
	std::cerr << recon.moments.U_earth_to_event << std::endl;
	std::cerr << recon.moments.scopes[0].U_cam_to_event << std::endl;
	std::cerr << "e: " << m2results.e << std::endl;
	std::cerr << "R: " << m2results.R << std::endl;
#endif //DEBUG

	recon.e       = 
	  recon.moments.U_earth_to_event.transBwd(m2results.e);
	recon.R       = 
	  recon.moments.U_earth_to_event.transBwd(m2results.R);
	recon.chi2e   = m2results.chi2;
	recon.chi2R   = m2results.chi2;
	recon.deltael = m2results.chi2_e_curv_l;
	recon.deltaew = m2results.chi2_e_curv_w;
	recon.deltaRl = m2results.chi2_R_curv_l;
	recon.deltaRw = m2results.chi2_R_curv_w;

	good = !m2results.at_edge;
      }
      break;

    case M_3:
    case M_3_GRIDSEARCH:
      {
	Method3::MinimizationResults m3results;
	std::vector<GridSearchSpecs> grid;
	grid.push_back(GridSearchSpecs(VSAR_M3_GRID*1.000,VSAR_M3_GRID*0.100));
	grid.push_back(GridSearchSpecs(VSAR_M3_GRID*0.100,VSAR_M3_GRID*0.010));
	grid.push_back(GridSearchSpecs(VSAR_M3_GRID*0.010,VSAR_M3_GRID*0.001));

	method3()->gridSearch(m3results, recon.moments, images, 
			      fArray.scopes, grid);

	recon.e       = 
	  recon.moments.U_earth_to_event.transBwd(m3results.e);
	recon.R       = 
	  recon.moments.U_earth_to_event.transBwd(m3results.R);
	recon.chi2e   = m3results.chi2;
	recon.chi2R   = m3results.chi2;
	recon.deltael = m3results.chi2_e_curv_l;
	recon.deltaew = m3results.chi2_e_curv_w;
	recon.deltaRl = m3results.chi2_R_curv_l;
	recon.deltaRw = m3results.chi2_R_curv_w;

	good = !m3results.at_edge;
      }
      break;
    }

  return good;
}

/*! 
  Calculate Moments for each telescope in camera frame
*/
void VSAReconstruction::
calculateArrayMoments(ArrayMoments& eventmoments, ScopeWeighting weighting,
		      const std::vector<ScopeImage>& images,
		      const ArrayQualityCuts* qc) const
{
  unsigned nscopes = images.size();
#ifdef VSA_ASSERT
  assert(nscopes == fArray.scopes.size());
#endif

  if(qc==0)qc=&fQualityCuts;

  eventmoments.nuse_in_reconstruction = 0;
  eventmoments.N = 0;
  eventmoments.scopes.clear();
  eventmoments.scopes.resize(nscopes);

  Vec3D psi(0,0,0);
  Vec3D phi(0,0,0);

  // --------------------------------------------------------------------------
  // First run through the data to calculate relative weighting
  // between telescopes and to choose appropriate basis vectors
  // --------------------------------------------------------------------------

  unsigned QC_nimage_min = 0;
  double QC_Ni_min = 0;
  double QC_width_min = -1;
  double QC_dist2_max = std::numeric_limits<double>::infinity();

  if(qc->scopes.size()==1)
    QC_nimage_min = qc->scopes[0].nimage_min,
      QC_Ni_min = qc->scopes[0].N_min,
      QC_width_min = qc->scopes[0].width_min,
      QC_dist2_max = qc->scopes[0].dist2_max;

  if(QC_Ni_min<0)QC_Ni_min=0;

  for(unsigned iscope=0; iscope<nscopes; iscope++)
    {
      const ScopeImage& scopeimage(images[iscope]);
      ScopeMoments& scopemoments(eventmoments.scopes[iscope]);
      const ScopeInfo& scopeinfo(fArray.scopes[iscope]);

      double Ni = 0;

      unsigned npixel=scopeimage.pixels.size();
      for(unsigned ipixel=0;ipixel<npixel;ipixel++)
	{
	  const unsigned pixelid = scopeimage.pixels[ipixel].j;
#ifdef VSA_ASSERT
	  assert(pixelid < scopeinfo.pij.size());
#endif

	  const double nij = scopeimage.pixels[ipixel].nij;
	  Ni += nij;
	}

      if(qc->scopes.size() >= nscopes)
	{
	  QC_nimage_min = qc->scopes[iscope].nimage_min;
	  QC_Ni_min = qc->scopes[iscope].N_min;
	  if(QC_Ni_min<0)QC_Ni_min=0;
	}

      scopemoments.use_in_reconstruction = 
	scopeimage.use_in_reconstruction;
      scopemoments.Ni = Ni;

      if((npixel<QC_nimage_min)||(Ni <= QC_Ni_min))
	scopemoments.use_in_reconstruction = false;

      if(scopemoments.use_in_reconstruction)
	{
	  eventmoments.N += Ni;

	  Vec3D ei(sin(scopeimage.zenith)*sin(scopeimage.azimuth),
		   sin(scopeimage.zenith)*cos(scopeimage.azimuth),
		   cos(scopeimage.zenith));

	  psi += ei*Ni;
	  phi += scopeinfo.ri*Ni;
	}  
    }

  if(eventmoments.N <= 0)
    {
      for(unsigned iscope=0; iscope<nscopes; iscope++)
	eventmoments.scopes[iscope].use_in_reconstruction = false;
      return;
    }

#ifdef VSA_ASSERT
  assert(psi.norm2() >= DBL_EPSILON);
#endif

  if(1.0-fabs(phi*psi)/psi.norm()/phi.norm() < DBL_EPSILON)
    {
      Vec3D blah;
      psi.perp(blah,phi);
    }

  eventmoments.U_earth_to_event = Orthogonal3D(psi,phi,SS_31);

  // --------------------------------------------------------------------------
  // Second run through computes vector and tensor moments and
  // tranforms to the event frame
  // --------------------------------------------------------------------------
  
  for(unsigned iscope=0; iscope<nscopes; iscope++)
    {
      ScopeMoments& scopemoments(eventmoments.scopes[iscope]);
      const ScopeImage& scopeimage(images[iscope]);
      const ScopeInfo& scopeinfo(fArray.scopes[iscope]);

      Vec3D cam_rotation(0,0,-scopeimage.camera_rotation);
      // Subtle extra rotation of PI/2 to get from camera frame where
      // tubes are in x-y plane with z pointing to the telescope to
      // global frame. Then rotate again by PI/2-Zenith to get to ERF
      cam_rotation &= Vec3D(M_PI-scopeimage.zenith,0,0);
      cam_rotation &= Vec3D(0,0,-scopeimage.azimuth);

      Orthogonal3D U_cam_to_earth(cam_rotation);
      Orthogonal3D U = eventmoments.U_earth_to_event * U_cam_to_earth;
      
      scopemoments.U_cam_to_event = U;
      scopemoments.ei = U.transFwd(Vec3D(0,0,1));
      scopemoments.ri = eventmoments.U_earth_to_event.transFwd(scopeinfo.ri);
      
      const double Ni = scopemoments.Ni;

      if(Ni<=0)
	{
	  scopemoments.Wi = 0.0;	
	  scopemoments.Vi_cam.set(0,0);
	  scopemoments.Ti_cam.set(0,0,0);
	  scopemoments.Vi.set(0,0,0);
	  scopemoments.Ti.set(0,0,0,0,0,0);
	  scopemoments.HillasParameters.clear();
	  continue;
	}

      Vec2D Vi;
      Symmetric2D Ti;

      unsigned npixel=scopeimage.pixels.size();
      for(unsigned ipixel=0;ipixel<npixel;ipixel++)
	{
	  const unsigned pixelid = scopeimage.pixels[ipixel].j;

	  const double nij = scopeimage.pixels[ipixel].nij;
	  const Vec2D& pij(scopeinfo.pij[pixelid]);

	  Vi += nij*pij;
	  Ti += nij*Symmetric2D(pij);
	}
      
      Vi /= Ni;
      Ti /= Ni;

      scopemoments.Vi_cam = Vi;
      scopemoments.Ti_cam = Ti;
      
      scopemoments.Vi = U.transFwd(Vec3D(Vi.x(),Vi.y(),0));
      scopemoments.Ti = U.transFwd(Symmetric3D(Ti.xx(),Ti.yy(),0,Ti.xy(),0,0));

      const Symmetric3D m2(scopemoments.Ti - Symmetric3D(scopemoments.Vi));
      Eigen3D& eigeni(scopemoments.HillasParameters);
      const unsigned __attribute((unused)) nroot = m2.eigen(eigeni);

      if(qc->scopes.size() >= nscopes)
	{
	  QC_width_min = qc->scopes[iscope].width_min;
	  QC_dist2_max = qc->scopes[iscope].dist2_max;
	}

      const double li2 = eigeni.val[2];
      const double wi2 = eigeni.val[1];

      if((QC_width_min>=0)&&(wi2<=QC_width_min*QC_width_min))
	{
	  scopemoments.use_in_reconstruction = false;
	  scopemoments.Wi = 0.0;	
	  continue;
	}

      if((std::isnormal(QC_dist2_max))&&(QC_dist2_max>=0)&&
	 (scopemoments.Vi_cam.norm2()>QC_dist2_max))
	{
	  scopemoments.use_in_reconstruction = false;
	  scopemoments.Wi = 0.0;	
	  continue;
	}

      switch(weighting)
	{
	case SW_ONE:
	  scopemoments.Wi = 1.0;
	  break;

	case SW_SIZE:
	  scopemoments.Wi = Ni;
	  break;

	case SW_ELLIPTICITY_MEMO:
	  if(wi2<=0)
	    {
	      scopemoments.Wi = 0.0;
	      scopemoments.use_in_reconstruction = false;
	    }
	  else scopemoments.Wi = (li2/wi2-1);
	  break;

	case SW_SIZE_ELLIPTICITY_MEMO:
	  if(wi2<=0)
	    {
	      scopemoments.Wi = 0.0;
	      scopemoments.use_in_reconstruction = false;
	    }
	  else scopemoments.Wi = Ni*(li2/wi2-1);
	  break;

	case SW_ELLIPTICITY_TRAD:
	  if(wi2<=0)
	    {
	      scopemoments.Wi = 0.0;
	      scopemoments.use_in_reconstruction = false;
	    }
	  else scopemoments.Wi = (1-wi2/li2);
	  break;

	case SW_WIDTH:
	  if(wi2<=0)
	    {
	      scopemoments.Wi = 0.0;
	      scopemoments.use_in_reconstruction = false;
	    }
	  else scopemoments.Wi = 1.0/sqrt(wi2);
	  break;
	}

      if(scopemoments.use_in_reconstruction)
	eventmoments.nuse_in_reconstruction++;

#ifdef DEBUG
      std::cerr << U_cam_to_earth.transFwd(Vec3D(0,0,1))
		<< std::endl << std::endl;
      std::cerr << eventmoments.U_earth_to_event.transFwd(U_cam_to_earth.transFwd(Vec3D(0,0,1))) << std::endl << std::endl;
      
      std::cerr << scopemoments.ei << std::endl << std::endl;
      std::cerr << scopemoments.ri << std::endl << std::endl;

      std::cerr << Vi << std::endl << std::endl;
      std::cerr << Ti << std::endl << std::endl;
      
      std::cerr << scopemoments.Vi << std::endl << std::endl;
      std::cerr << scopemoments.Ti << std::endl << std::endl;

      Eigen2D eigen2;
      Ti.eigen(eigen2);
      std::cerr << eigen2.vec[0] << ' ' << eigen2.val[0] << std::endl;
      std::cerr << eigen2.vec[1] << ' ' << eigen2.val[1] << std::endl;
      std::cerr << std::endl;

      Eigen3D eigen3;
      scopemoments.Ti.eigen(eigen3);
      std::cerr << eigen3.vec[0] << ' ' << eigen3.val[0] << std::endl;
      std::cerr << eigen3.vec[1] << ' ' << eigen3.val[1] << std::endl;
      std::cerr << eigen3.vec[2] << ' ' << eigen3.val[2] << std::endl;
      std::cerr << std::endl;
#endif
    }
}

void VSAReconstruction::
calculateArrayParameters(ArrayParameters& param,
			 const std::vector<ScopeImage>& images,
			 const Reconstruction& reconstruction) const
{
  unsigned nscopes=images.size();
#ifdef VSA_ASSERT
  assert(nscopes == fArray.scopes.size());
  assert(reconstruction.moments.scopes.size() == fArray.scopes.size());
#endif

  param = ArrayParameters();
  param.scopes.resize(nscopes);

  const ArrayMoments& moments(reconstruction.moments);
  const Vec3D& e(reconstruction.e);
  const Vec3D& R(reconstruction.R);

  const Orthogonal3D& U_earth_to_event(moments.U_earth_to_event);

  double R_e = R*e;

  double N                 = 0;
  double d                 = 0;
  double d2                = 0;
  Vec3D Delta;
  Symmetric3D Delta2;
  double h                 = 0;
  double h2                = 0;
  double G                 = 0;
  double G2                = 0;
  double t0                = 0;
  double t02               = 0;

  double h0 = fArray.elevation;

  for(unsigned iscope=0; iscope<nscopes; iscope++)
    {
      static const double c = 2.9979246e+08 / 1e9;

      const ScopeImage& scopeimage(images[iscope]);
      const ScopeInfo& scopeinfo(fArray.scopes[iscope]);
      ScopeParameters& sp(param.scopes[iscope]);
      
      const Orthogonal3D& 
	U_cam_to_event(moments.scopes[iscope].U_cam_to_event);

      Orthogonal3D U_cam_to_earth = U_earth_to_event.inverse()*U_cam_to_event;

      Vec3D ei = U_cam_to_earth.transFwd(Vec3D(0,0,1));

      const Vec3D& ri(scopeinfo.ri);
      double D2i = scopeinfo.Di*scopeinfo.Di;

      const Vec3D Ri(R-ri);
      double Ri_e = Ri*e;

      double Ni            = 0;
      double di            = 0;
      double d2i           = 0;
      double li            = 0;
      double l2i           = 0;
      double thetai        = 0;
      double theta2i       = 0;
      Vec3D Deltai;
      Symmetric3D Delta2i;
      double hi            = 0;
      double h2i           = 0;
      double Gi            = 0;
      double G2i           = 0;

      double t0i           = 0;
      double t02i          = 0;

      unsigned npixel=scopeimage.pixels.size();
      for(unsigned ipixel=0;ipixel<npixel;ipixel++)
	{
	  const unsigned pixelid = scopeimage.pixels[ipixel].j;

	  const double nij = scopeimage.pixels[ipixel].nij;
	  const Vec2D& pij(scopeinfo.pij[pixelid]);

	  Vec3D eij = U_cam_to_earth.transFwd(Vec3D(pij.x(),pij.y(),
						    sqrt(1-pij.norm2())));
	  
	  double eij_e = eij*e;
	  double ri_eij = ri*eij;
	  double one_minus_eij_e_sq = 1-eij_e*eij_e;
	  if(one_minus_eij_e_sq < DBL_EPSILON)continue;

	  // Angle of ray to primary
	  double thetaij = asin(sqrt(one_minus_eij_e_sq));

	  // Distance to emission along primary trajectory
	  double dij = (e-eij_e*eij)*(ri-R)/one_minus_eij_e_sq + R_e;
	  Vec3D Dij = R+(dij-R_e)*e;

	  // Distance from telescope to emission point
	  double lij = -(eij-eij_e*e)*(ri-R)/one_minus_eij_e_sq + ri_eij;
	  Vec3D Lij = ri + (lij - ri_eij)*eij;

	  // Vector joining closest approach of photon/primary trajectory
	  Vec3D Deltaij = Dij-Lij;

	  // Depth in atmosphere and propagation delay
	  double hij = Dij.z() + h0;

	  double atm_rho;
	  double atm_thickness;
	  double atm_n_minus_one;
	  double atm_delay;
	  fAtmo.interpolate(hij/1000.0, atm_rho, atm_thickness,
			    atm_n_minus_one, atm_delay);

	  double Gij = atm_thickness/fabs(e.z());
	  double tauij = (atm_delay-fObsDelay)/fabs(e.z());

	  // Timing
	  double t0ij = scopeimage.pixels[ipixel].tij + lij/c - dij/c - tauij;

#ifdef DEBUG
	  std::cerr << "t0(" << iscope << ',' << pixelid << ") =  " 
		    << scopeimage.pixels[ipixel].tij << " + " 
		    << lij/c << " - " << dij/c << " - " << tauij
		    << " = " << t0ij << std::endl;
#endif

	  // Luminosity

	  Ni               += nij;
	  di               += nij*dij;
	  d2i              += nij*dij*dij;
	  li               += nij*lij;
	  l2i              += nij*lij*lij;
	  thetai           += nij*thetaij;
	  theta2i          += nij*thetaij*thetaij;
	  Deltai           += nij*Deltaij;
	  Delta2i          += nij*Symmetric3D(Deltaij);
	  hi               += nij*hij;
	  h2i              += nij*hij*hij;
	  Gi               += nij*Gij;
	  G2i              += nij*Gij*Gij;
	  t0i              += nij*t0ij;
	  t02i             += nij*t0ij*t0ij;
	}

      if(Ni > 0)
	{
	  if(scopeimage.use_in_reconstruction)
	    {
	      // Accumulate array parameters
	      N                += Ni;
	      d                += di;
	      d2               += d2i;
	      Delta            += Deltai;
	      Delta2           += Delta2i;
	      G                += Gi;
	      G2               += G2i;
	      h                += hi;
	      h2               += h2i;
	      t0               += t0i;
	      t02              += t02i;
	    }

	  // Transfer scope parameters to structure
	  sp.Ni            = Ni;
	  sp.Ri            = sqrt(Ri*Ri-Ri_e*Ri_e);
	  sp.d1i           = di/Ni;
	  sp.d2i           = sqrt(d2i/Ni - sp.d1i*sp.d1i);
	  sp.l1i           = li/Ni;
	  sp.l2i           = sqrt(l2i/Ni - sp.l1i*sp.l1i);
	  sp.theta1i       = thetai/Ni;
	  sp.theta2i       = sqrt(theta2i/Ni - sp.theta1i*sp.theta1i);
	  sp.delta1i       = Deltai/Ni;
	  sp.delta2i       = Delta2i/Ni - Symmetric3D(sp.delta1i);
	  sp.D1i           = R+(sp.d1i-R_e)*e;
	  sp.h1i           = hi/Ni;
	  sp.h2i           = sqrt(h2i/Ni - sp.h1i*sp.h1i);
	  sp.G1i           = Gi/Ni;
	  sp.G2i           = sqrt(G2i/Ni - sp.G1i*sp.G1i);

	  double thetac    = acos(1.0/(1.0+fAtmo.nMinusOne(sp.h1i/1000.0)));
	  double sinthetac = sin(thetac);
	  sp.thetac1       = thetac;

	  sp.t01i          = t0i/Ni;
	  sp.t02i          = sqrt(t02i/Ni - sp.t01i*sp.t01i);

	  sp.lambdadi      = 
	    fabs(8*Ni*sp.l1i*sp.l1i/sp.d2i/sinthetac/sinthetac/D2i);
	  sp.lambdaci      = fabs(8*Ni*sp.l1i/D2i);
	}
    }
}

std::auto_ptr<Method1> VSAReconstruction::sMethod1;
std::auto_ptr<Method2> VSAReconstruction::sMethod2;
std::auto_ptr<Method3> VSAReconstruction::sMethod3;

Method1* VSAReconstruction::method1()
{
  if(sMethod1.get() == 0)sMethod1.reset(new Method1);
  return sMethod1.get();
}

Method2* VSAReconstruction::method2()
{
  if(sMethod2.get() == 0)sMethod2.reset(new Method2);
  return sMethod2.get();
}

Method3* VSAReconstruction::method3()
{
  if(sMethod3.get() == 0)sMethod3.reset(new Method3);
  return sMethod3.get();
}
