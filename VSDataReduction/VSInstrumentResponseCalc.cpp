//-*-mode:c++; mode:font-lock;-*-

/*! \file VSInstrumentResponseCalc.cpp

  General purpose class for generating detector response functions.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSInstrumentResponseCalc.cpp,v 3.5 2010/10/20 17:51:38 matthew Exp $

*/

#include <VSFileUtility.hpp>
#include <VSInstrumentResponseCalc.hpp>
#include <VSPSFCalc.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSAQuadrature.hpp>

using namespace VERITAS;

// ============================================================================
// VSPSFFn
// ============================================================================
VSPSFFn::VSPSFFn():
  VSAFunction::ParamFn<double>(4)
{
  
}

double VSPSFFn::val(const double& x) const 
{ 
  return val(x,param());
}

double VSPSFFn::val(const double& x, const VSAAlgebra::VecND& a) const 
{ 
  double s1i = std::pow(a(1),-2);
  double s2i = std::pow(a(2),-2);

  return a(0)*std::pow(2*M_PI,-1)*
    (a(3)*s1i*exp(-x*s1i/2.) + (1-a(3))*s2i*exp(-x*s2i/2.));
}

void VSPSFFn::dyda(const double& x, VSAAlgebra::VecND& dyda) const 
{
  VSPSFFn::dyda(x,param(),dyda);
}

void VSPSFFn::dyda(const double& x, const VSAAlgebra::VecND& a, 
		   VSAAlgebra::VecND& dyda) const 
{
  dyda.resize(4);

  double s1i = std::pow(a(1),-2);
  double s2i = std::pow(a(2),-2);

  dyda(0) = std::pow(2*M_PI,-1)*
    (a(3)*s1i*exp(-x*s1i/2.) + (1-a(3))*s2i*exp(-x*s2i/2.));
  dyda(1) = a(0)*std::pow(2*M_PI,-1)*
    a(3)*s1i*exp(-x*s1i/2.)*(x-2*a(1)*a(1))/std::pow(a(1),3);
  dyda(2) = a(0)*std::pow(2*M_PI,-1)*
    (1-a(3))*s2i*exp(-x*s2i/2.)*(x-2*a(2)*a(2))/std::pow(a(2),3);
  dyda(3) = a(0)*std::pow(2*M_PI,-1)*
    (s1i*exp(-x*s1i/2.) - s2i*exp(-x*s2i/2.));
}

double VSPSFFn::integrate(double xlo, double xhi)
{
  return integrate(param(),xlo,xhi);
}

double VSPSFFn::integrate(const VSAAlgebra::VecND& a, 
			  double xlo, double xhi)
{
  double s1i = std::pow(a(1),-2);
  double s2i = std::pow(a(2),-2);

  return a(0)*(exp(-xlo*s2i/2.) - exp(-xhi*s2i/2.) + 
    a(3)*(exp(-xlo*s1i/2.) - exp(-xhi*s1i/2.) - 
	  exp(-xlo*s2i/2.) + exp(-xhi*s2i/2.)));
}

// ============================================================================
// VSEnergyResponseFn
// ============================================================================
VSEnergyResponseFn::VSEnergyResponseFn(double dx, double norm):
  VSAFunction::ParamFn<VSACoord::CoordND>(5), m_norm(norm), m_dx(dx)
{

}

double VSEnergyResponseFn::val(const VSACoord::CoordND& x) const 
{ 
  return val(x,param());
}

double VSEnergyResponseFn::val(const VSACoord::CoordND& x, 
			       const VSAAlgebra::VecND& a) const 
{ 
  return val(x,a,m_norm);
}

double VSEnergyResponseFn::val(const VSACoord::CoordND& x, 
			       const VSAAlgebra::VecND& a, 
			       double norm) const 
{ 
  double g1 = m_fn.val(x[0],norm*a(4),fabs(a(0)),a(1));
  double g2 = m_fn.val(x[0],norm*(1-a(4)),fabs(a(2)),a(3));

  if(!std::isfinite(g1)) g1 = 0;
  if(!std::isfinite(g2)) g2 = 0;

  double v = (g1+g2)*m_dx;

  return v + 1E-100;
}

void VSEnergyResponseFn::dyda(const VSACoord::CoordND& x, 
			      VSAAlgebra::VecND& dyda) const 
{
  VSEnergyResponseFn::dyda(x,param(),dyda);
}

void VSEnergyResponseFn::dyda(const VSACoord::CoordND& x, 
			      const VSAAlgebra::VecND& a, 
			      VSAAlgebra::VecND& dyda) const 
{
  dyda.resize(5);

  double g1 = m_fn.val(x[0],m_norm,fabs(a(0)),a(1));
  double g2 = m_fn.val(x[0],m_norm,fabs(a(2)),a(3));

  if(!std::isfinite(g1)) g1 = 0;
  if(!std::isfinite(g2)) g2 = 0;

  VSAAlgebra::VecND a1(3);
  a1(0) = m_norm*a(4);
  a1(1) = fabs(a(0));
  //  a1(2) = a(1) + x[1];
  a1(2) = a(1);

  VSAAlgebra::VecND a2(3);
  a2(0) = m_norm*(1-a(4));
  a2(1) = fabs(a(2));
  //  a2(2) = a(3) + x[1];
  a2(2) = a(3);

  VSAAlgebra::VecND d1(3);
  m_fn.dyda(x[0],a1,d1);

  VSAAlgebra::VecND d2(3);
  m_fn.dyda(x[0],a2,d2);

  dyda(0) = a(4)*d1(1);
  dyda(1) = a(4)*d1(2);
  dyda(2) = (1-a(4))*d2(1);
  dyda(3) = (1-a(4))*d2(2);
  dyda(4) = g1-g2;

  dyda *= m_dx;
}

// ============================================================================
// VSEnergyKernelEffareaFn
// ============================================================================
VSEnergyKernelEffareaFn::
VSEnergyKernelEffareaFn(const VSNSpace&  effarea,
			const VSNSpace&  krn_sigma1,
			const VSNSpace&  krn_bias1,
			const VSNSpace&  krn_sigma2,
			const VSNSpace&  krn_bias2,
			const VSNSpace&  krn_alpha,
			VSSpectrumFn* sp_fn,
			double offset,
			double theta_cut):
  m_effarea(effarea), 
  m_krn_sigma1(krn_sigma1), m_krn_bias1(krn_bias1), 
  m_krn_sigma2(krn_sigma2), m_krn_bias2(krn_bias2), m_krn_alpha(krn_alpha),
  m_dpdloge_fn(), m_sp_fn(sp_fn->clone()), 
  m_offset(offset), m_theta_cut(theta_cut)
{

}

double VSEnergyKernelEffareaFn::val(const VSACoord::Coord2D& x) const
{
  const double erec = x[0];
  const double emc = x[1];


  VSACoord::CoordND cnd(2);

  cnd[0] = erec;
  cnd[1] = emc;

  VSNSpace::Point p(1);
  p.x[0] = emc;

  VSAAlgebra::VecND param = m_dpdloge_fn.param();

  m_krn_sigma1.interpolateWeight(p,param[0]);
  m_krn_bias1.interpolateWeight(p,param[1]);
  m_krn_sigma2.interpolateWeight(p,param[2]);
  m_krn_bias2.interpolateWeight(p,param[3]);
  m_krn_alpha.interpolateWeight(p,param[4]);

  double dpdloge = m_dpdloge_fn.val(cnd,param);
  double dfde = m_sp_fn->val(emc);
  double effarea = 0;

  m_effarea.interpolateWeight(p,effarea);

  double w = effarea*dpdloge*dfde*std::pow(10,emc)*std::log(10);
  return w;
}

// ============================================================================
// VSInstrumentResponseCalc
// ============================================================================
VSInstrumentResponseCalc::Options 
VSInstrumentResponseCalc::s_default_options = 
  VSInstrumentResponseCalc::Options();

VSInstrumentResponseCalc::VSInstrumentResponseCalc(const Options& opt):
  VSLTLibraryCalc(),
  m_options(opt),
  m_effarea_nspace(),
  m_psf_sigma1_nspace(),
  m_psf_sigma2_nspace(),
  m_psf_alpha_nspace(),
  m_krn_sigma1_nspace(),
  m_krn_bias1_nspace(),
  m_krn_sigma2_nspace(),
  m_krn_bias2_nspace(),
  m_krn_alpha_nspace(),
  m_theta_max()
{

}

VSInstrumentResponseCalc::~VSInstrumentResponseCalc()
{

}

// VSInstrumentResponseCalc::
// VSInstrumentResponseCalc(const VSNSpace& effarea,
// 			     const VSNSpace& psf_sigma1,
// 			     const VSNSpace& psf_sigma2,
// 			     const VSNSpace& psf_alpha,
// 			     double theta_max):
//   m_effarea_nspace(effarea),
//   m_psf_sigma1_nspace(psf_sigma1),
//   m_psf_sigma2_nspace(psf_sigma2),
//   m_psf_alpha_nspace(psf_alpha),
//   m_theta_max(theta_max)
// {

// }

VSNSpace VSInstrumentResponseCalc::getPSF(const VSSpectrumFn* spectrum_fn,
					  double offset) const
{
  VSNSpace::Space space(1);

  space.axes[0] = VSNSpace::Axis(0, 0.15, 0.0005, "Theta-Squared"); 

  VSNSpace psf(space);

  if(offset > m_psf_sigma1_nspace.axis(1).hi_bound) return psf;

  VSNSpace::Cell c(2);

  double dloge = m_effarea_nspace.axis(0).bin_size;
  double sum_rate = 0;
  // Loop on Energy -------------------------------------------------------
  for(c.i[0] = 0; c.i[0] < m_effarea_nspace.axis(0).nbin; c.i[0]++)
    {
      double emc = m_effarea_nspace.axis(0).midCoordUnchecked(c.i[0]);

      VSNSpace::Point p(2);
      
      p.x[0] = emc;
      p.x[1] = offset;

      double effarea;

      m_effarea_nspace.interpolateWeight(p,effarea);

      double egy_tev = std::pow(10,p.x[0]);
      double dfde = spectrum_fn->diffFlux(p.x[0]);
      double rate = effarea*egy_tev*dfde*dloge*log(10.);
            
      double sigma1, sigma2, alpha;
      
      m_psf_sigma1_nspace.interpolateWeight(p,sigma1);
      m_psf_sigma2_nspace.interpolateWeight(p,sigma2);
      m_psf_alpha_nspace.interpolateWeight(p,alpha);

      sum_rate += rate;

      VSPSFFn psf_fn;
      psf_fn.setParam(0,1);
      psf_fn.setParam(1,sigma1);
      psf_fn.setParam(2,sigma2);
      psf_fn.setParam(3,alpha);
      
      VSNSpace::Cell c2(1);
      for(c2.i[0] = 0; c2.i[0] < psf.axis(0).nbin; c2.i[0]++)
	{
	  double thsq = psf.axis(0).midCoordUnchecked(c2.i[0]);
	  psf[c2] += psf_fn.val(thsq)*rate;
	}
    }

  psf *= (1./sum_rate)*M_PI;

  return psf;
}

VSNSpace VSInstrumentResponseCalc::
getEffectiveArea(double ebin, double elo, double ehi,
		 double offset, double theta) const
{
  VSNSpace::Space space(1);

  space.axes[0] = VSNSpace::Axis(elo, ehi, ebin, "MC Energy"); 

  VSNSpace effarea(space);

  if(offset > m_effarea_nspace.axis(1).hi_bound) return effarea;

  VSNSpace::Cell c(1);
  // Loop on Energy -----------------------------------------------------------
  for(c.i[0] = 0; c.i[0] < effarea.axis(0).nbin; c.i[0]++)
    {
      double emc = effarea.axis(0).midCoordUnchecked(c.i[0]);
      VSNSpace::Point p(2);



      p.x[0] = emc;
      p.x[1] = offset;
      
      double sigma1, sigma2, alpha;

      m_psf_sigma1_nspace.interpolateWeight(p,sigma1);
      m_psf_sigma2_nspace.interpolateWeight(p,sigma2);
      m_psf_alpha_nspace.interpolateWeight(p,alpha);

      VSPSFFn psf_fn;
      psf_fn.setParam(0,1);
      psf_fn.setParam(1,sigma1);
      psf_fn.setParam(2,sigma2);
      psf_fn.setParam(3,alpha);

      double v1 = psf_fn.integrate(0,0.3*0.3);
      double v2 = psf_fn.integrate(0,theta*theta);

      m_effarea_nspace.interpolateWeight(p,effarea[c]);      

      effarea[c] *= v2/v1;
    }

  return effarea;
}

VSNSpace VSInstrumentResponseCalc::
getEffectiveArea(double offset, double theta) const
{
  double dloge = m_effarea_nspace.axis(0).bin_size;

  return getEffectiveArea(dloge,
  			  m_effarea_nspace.axis(0).lo_bound,
  			  m_effarea_nspace.axis(0).hi_bound,
  			  offset, theta);

  // double dloge = m_effarea_nspace.axis(0).bin_size;

  // return getEffectiveArea(dloge,
  // 			  m_effarea_nspace.axis(0).lo_bound,
  // 			  m_effarea_nspace.axis(0).hi_bound,
  // 			  offset, theta);
}

double VSInstrumentResponseCalc::
getEffectiveArea(double offset, double theta_cut, double egy) const
{  
  if(offset > m_effarea_nspace.axis(1).hi_bound) return 0;

  VSNSpace::Point p(2);
  p.x[0] = egy;
  p.x[1] = offset;

  double sigma1, sigma2, alpha;

  m_psf_sigma1_nspace.interpolateWeight(p,sigma1);
  m_psf_sigma2_nspace.interpolateWeight(p,sigma2);
  m_psf_alpha_nspace.interpolateWeight(p,alpha);

  VSPSFFn psf_fn;
  psf_fn.setParam(0,1);
  psf_fn.setParam(1,sigma1);
  psf_fn.setParam(2,sigma2);
  psf_fn.setParam(3,alpha);
  
  double v1 = psf_fn.integrate(0,0.3*0.3);
  double v2 = psf_fn.integrate(0,theta_cut*theta_cut);

  double w = 0;
  m_effarea_nspace.interpolateWeight(p,w);      
  return w*v2/v1;
}

void VSInstrumentResponseCalc::
getEffectiveAreaAperture(const std::string& spectrum_fn,
			 double theta, VSNSpace& effarea) const
{
  VSSpectrumFn* sp_fn = VSSpectrumFn::create(spectrum_fn);
  sp_fn->setNormalization(1);


  effarea = VSNSpace(1,m_effarea_nspace.axis(1).lo_bound,
		     m_effarea_nspace.axis(1).hi_bound,
		     m_effarea_nspace.axis(1).nbin);

  VSNSpace::Cell c(2);
  // Loop on Offset Angle -----------------------------------------------------
  for(c.i[1] = 0; c.i[1] < m_effarea_nspace.axis(1).nbin; c.i[1]++)
    {
      double dloge = m_effarea_nspace.axis(0).bin_size;
      double effarea_sum = 0;
      // Loop on Energy -------------------------------------------------------
      for(c.i[0] = 0; c.i[0] < m_effarea_nspace.axis(0).nbin; c.i[0]++)
	{
	  VSNSpace::Point p;
	  m_effarea_nspace.space().midPointOfCellUnchecked(c,p);
	  double egy_tev = std::pow(10,p.x[0]);
	  double flux = sp_fn->diffFlux(p.x[0]);
	  double effarea = m_effarea_nspace[c]*egy_tev*flux*dloge*log(10.);
	  
	  if(effarea <= 0) continue;

	  double sigma1 = m_psf_sigma1_nspace[c];
	  double sigma2 = m_psf_sigma2_nspace[c];
	  double alpha  = m_psf_alpha_nspace[c];

	  VSPSFFn psf_fn;
	  psf_fn.setParam(0,1);
	  psf_fn.setParam(1,sigma1);
	  psf_fn.setParam(2,sigma2);
	  psf_fn.setParam(3,alpha);

	  double v1 = psf_fn.integrate(0,0.3*0.3);
	  double v2 = psf_fn.integrate(0,theta*theta);

	  effarea_sum += effarea*v2/v1;
	}

      effarea[c.i[1]] = effarea_sum;
    }

  delete sp_fn;
}

double VSInstrumentResponseCalc::
getEffectiveAreaAperture(const std::string& spectrum_fn,
			 double offset, double theta_cut) const
{
  VSSpectrumFn* sp_fn = VSSpectrumFn::create(spectrum_fn);
  sp_fn->setNormalization(1);

  const double emin = m_effarea_nspace.axis(0).lo_bound;
  const double emax = m_effarea_nspace.axis(0).hi_bound;
  const double ebin = m_effarea_nspace.axis(0).bin_size/2.;

  double effarea_sum = 0;

  for(double emc = emin; emc < emax; emc += ebin)
    {
      VSNSpace::Point p = m_effarea_nspace.space().point();

      p.x[0] = emc;
      p.x[1] = offset;

      double effarea = 0;
      m_effarea_nspace.interpolateWeight(p,effarea);

      double egy_tev = std::pow(10,p.x[0]);
      double flux = sp_fn->diffFlux(p.x[0]);
	  
      if(effarea <= 0) continue;

      effarea *= egy_tev*flux*ebin*log(10.);

      double sigma1 = 0;
      double sigma2 = 0;
      double alpha  = 0;
      
      m_psf_sigma1_nspace.interpolateWeight(p,sigma1);
      m_psf_sigma2_nspace.interpolateWeight(p,sigma2);
      m_psf_alpha_nspace.interpolateWeight(p,alpha);

      VSPSFFn psf_fn;
      psf_fn.setParam(0,1);
      psf_fn.setParam(1,sigma1);
      psf_fn.setParam(2,sigma2);
      psf_fn.setParam(3,alpha);
      
      double v1 = psf_fn.integrate(0,0.3*0.3);
      double v2 = psf_fn.integrate(0,std::pow(theta_cut,2));
      
      effarea_sum += effarea*v2/v1;
    }
  
  delete sp_fn;

  return effarea_sum;
}

void VSInstrumentResponseCalc::
getEffectiveAreaOffset(const std::string& spectrum_fn,
		       VSNSpace& effarea_nspace,
		       VSNSpace& effarea_domega_nspace) const
{
  VSSpectrumFn* sp_fn = VSSpectrumFn::create(spectrum_fn);
  sp_fn->setNormalization(1);

  getEffectiveAreaOffset(sp_fn,effarea_nspace,effarea_domega_nspace);

  delete sp_fn;  
}

void VSInstrumentResponseCalc::
getEffectiveAreaOffset(const VSSpectrumFn* sp_fn,
		       VSNSpace& effarea_nspace,
		       VSNSpace& effarea_domega_nspace) const
{



  effarea_nspace = VSNSpace(1,m_effarea_nspace.axis(1).lo_bound,
			    m_effarea_nspace.axis(1).hi_bound,
			    m_effarea_nspace.axis(1).nbin);

  effarea_domega_nspace = VSNSpace(m_effarea_nspace.axis(1).lo_bound,
				m_effarea_nspace.axis(1).hi_bound,
				m_effarea_nspace.axis(1).nbin,0,0.4,80);

  VSNSpace::Cell c(2);
  // Loop on Offset Angle -----------------------------------------------------
  for(c.i[1] = 0; c.i[1] < m_effarea_nspace.axis(1).nbin; c.i[1]++)
    {
      double dloge = m_effarea_nspace.axis(0).bin_size;
      double effarea_sum = 0;
      // Loop on Energy -------------------------------------------------------
      for(c.i[0] = 0; c.i[0] < m_effarea_nspace.axis(0).nbin; c.i[0]++)
	{
	  VSNSpace::Point p;
	  m_effarea_nspace.space().midPointOfCellUnchecked(c,p);
	  double egy_tev = std::pow(10,p.x[0]);
	  double flux = sp_fn->diffFlux(p.x[0]);
	  double effarea = m_effarea_nspace[c]*egy_tev*flux*dloge*log(10.);
	  
	  if(effarea <= 0) continue;

	  double sigma1 = m_psf_sigma1_nspace[c];
	  double sigma2 = m_psf_sigma2_nspace[c];
	  double alpha  = m_psf_alpha_nspace[c];

	  effarea_sum += effarea;

	  VSPSFFn psf_fn;
	  psf_fn.setParam(0,1);
	  psf_fn.setParam(1,sigma1);
	  psf_fn.setParam(2,sigma2);
	  psf_fn.setParam(3,alpha);

	  VSNSpace::Cell c2(2);
	  c2.i[0] = c.i[1];
	  for(c2.i[1] = 0; c2.i[1] < 
		effarea_domega_nspace.axis(1).nbin; c2.i[1]++)
	    {
	      double theta = 
		effarea_domega_nspace.axis(1).midCoordUnchecked(c2.i[1]);
	      double thsq = std::pow(theta,2);
	      effarea_domega_nspace[c2] += psf_fn.val(thsq)*effarea;
	    }
	}

      effarea_nspace[c.i[1]] = effarea_sum;
    }


}


VSNSpace VSInstrumentResponseCalc::
getKernel(double ebin, double elo, double ehi,
	  double offset, double theta_cut) const
{
  const double erec_lo = -1.34375;
  const double erec_hi = 2.15625;
  const double erec_bin = 0.0625;

  return getKernel(ebin,elo,ehi,erec_bin,erec_lo,erec_hi,offset,theta_cut);
}

VSNSpace VSInstrumentResponseCalc::
getKernel(double emc_bin, double emc_lo, double emc_hi,
	  double erec_bin, double erec_lo, double erec_hi,
	  double offset, double theta_cut) const
{
  VSNSpace::Space space(2);

  space.axes[0] = 
    VSNSpace::Axis(erec_lo, erec_hi, erec_bin, "Reconstructed Energy"); 
  space.axes[1] = VSNSpace::Axis(emc_lo, emc_hi, emc_bin, "MC Energy"); 

  VSNSpace krn(space);
	       
  if(offset > m_krn_sigma1_nspace.axis(1).hi_bound) return krn;

  VSNSpace::Cell c(2);
  
  VSEnergyResponseFn fn;
  VSSpectrumFn* sp_fn = new VSSpectrumFnPowerLaw(2.5);

  VSNSpace effarea = getEffectiveArea(offset,theta_cut);
  VSNSpace krn_bias1 = m_krn_bias1_nspace;
  VSNSpace krn_sigma1 = m_krn_sigma1_nspace;
  VSNSpace krn_bias2 = m_krn_bias2_nspace;
  VSNSpace krn_sigma2 = m_krn_sigma2_nspace;
  VSNSpace krn_alpha = m_krn_alpha_nspace;
  krn_bias1.slice(1,offset);
  krn_sigma1.slice(1,offset);
  krn_bias2.slice(1,offset);
  krn_sigma2.slice(1,offset);
  krn_alpha.slice(1,offset);

  // Loop on true energy ------------------------------------------------------
  for(c.i[1] = 0; c.i[1] < krn.axis(1).nbin; c.i[1]++)
    {
      double emc = krn.axis(1).midCoordUnchecked(c.i[1]); 
      double emc_lo = krn.axis(1).minCoordUnchecked(c.i[1]);
      double emc_hi = krn.axis(1).maxCoordUnchecked(c.i[1]);

      VSAAlgebra::VecND param = fn.param();
      VSNSpace::Point p(2);

      p.x[0] = emc;
      p.x[1] = offset;

      m_krn_sigma1_nspace.interpolateWeight(p,param(0));
      m_krn_bias1_nspace.interpolateWeight(p,param(1));
      m_krn_sigma2_nspace.interpolateWeight(p,param(2));
      m_krn_bias2_nspace.interpolateWeight(p,param(3));
      m_krn_alpha_nspace.interpolateWeight(p,param(4));

      // Loop on reconstructed energy -----------------------------------------
      for(c.i[0] = 0; c.i[0] < krn.axis(0).nbin; c.i[0]++)
	{	  
	  VSACoord::Coord2D clo;
	  VSACoord::Coord2D chi;
	  VSACoord::CoordND cmid(2);

	  clo.setCartesian(krn.axis(0).minCoordUnchecked(c.i[0]),
			   krn.axis(1).minCoordUnchecked(c.i[1]));

	  chi.setCartesian(krn.axis(0).maxCoordUnchecked(c.i[0]),
			   krn.axis(1).maxCoordUnchecked(c.i[1]));

	  cmid[0] = krn.axis(0).midCoordUnchecked(c.i[0]);
	  cmid[1] = krn.axis(1).midCoordUnchecked(c.i[1]);

	  double v4 = 0;

	  // ------------------------------------------------------------------
	  // Perform 2D integration over True and Reconstructed Energy
	  if(fn.val(cmid,param) > 1E-10)
	    {
	      VSEnergyKernelEffareaFn 
		kfn(effarea,
		    krn_sigma1,krn_bias1,krn_sigma2,krn_bias2,krn_alpha,
		    sp_fn,offset,theta_cut);
	      try
		{
		  v4 = VSAMath::Quadrature::integrate(kfn,clo,chi,1E-5);
		}
	      catch(const std::exception& e)
		{
		  std::cerr 
		    << std::string(__PRETTY_FUNCTION__) 
		    << ": Caught exception." << std::endl
		    << "emc_lo: " << emc_lo << std::endl
		    << "emc_hi: " << emc_hi << std::endl
		    << e.what() << std::endl;
		}
	    }
	  
// 	  std::cout << std::setw(15) << krn.axis(0).midCoordUnchecked(c.i[0])
// 		    << std::setw(15) << krn.axis(1).midCoordUnchecked(c.i[1])
// 		    << std::setw(15) << fn.val(cmid,param)
// 		    << std::setw(15) << v4
// 		    << std::endl;

	  // Normalize by integral flux
	  v4 /= sp_fn->integralFlux(emc_lo,emc_hi);

	  krn.setWeight(c,v4);
	}
    }

  delete sp_fn;

  return krn; 
}

VSNSpace VSInstrumentResponseCalc::
getKernel(double ebin, double elo, double ehi, double offset) const
{
  VSNSpace::Space space(2);

  space.axes[0] = VSNSpace::Axis(elo, ehi, ebin, "MC Energy"); 
  space.axes[1] = VSNSpace::Axis(elo, ehi, ebin, "Reconstructed Energy"); 

  VSNSpace krn(space);
	       
  VSNSpace::Cell c(2);
  
  VSEnergyResponseFn fn(ebin);

  for(c.i[0] = 0; c.i[0] < krn.space().axes[0].nbin; c.i[0]++)
    {
      double emc = krn.space().axes[0].midCoordUnchecked(c.i[0]);

      VSAAlgebra::VecND param = fn.param();

      VSNSpace::Point p(2);

      p.x[0] = emc;
      p.x[1] = offset;

      m_krn_sigma1_nspace.interpolateWeight(p,param(0));
      m_krn_bias1_nspace.interpolateWeight(p,param(1));
      m_krn_sigma2_nspace.interpolateWeight(p,param(2));
      m_krn_bias2_nspace.interpolateWeight(p,param(3));
      m_krn_alpha_nspace.interpolateWeight(p,param(4));

      double vsum = 0;
      for(c.i[1] = 0; c.i[1] < krn.space().axes[1].nbin; c.i[1]++)
	{
	  double erec = krn.space().axes[1].midCoordUnchecked(c.i[1]);
	  VSACoord::CoordND cnd(2);
	  cnd[0] = erec;
	  cnd[1] = emc;

// 	  double erec1 = krn.axis(0).minCoordUnchecked(c.i[1]);
// 	  double erec2 = krn.axis(0).maxCoordUnchecked(c.i[1]);
// 	  double erec_step = (erec2-erec1)/10.;
	  
// 	  VSAFunction::Spline erec_spline;

// 	  std::cout << erec << " " << emc << std::endl;

// 	  for(unsigned i = 0; i <= 10; i++)
// 	    {
// 	      double e = erec1+i*erec_step;
// 	      VSACoord::CoordND ec(2);
// 	      ec[0] = e;
// 	      ec[1] = emc;
// 	      erec_spline.setPoint(e,fn.val(ec,param)/ebin);
// 	    }

// 	  erec_spline.spline();


	  double v = fn.val(cnd,param);
	    //VSAMath::Quadrature::integrate(erec_spline,erec1,erec2);
	  //fn.val(cnd,param);
	  if(v < 1E-8) v = 0;

	  vsum += v;

// 	  std::cout << erec << " " << emc << " " << fn.val(cnd,param) 
// 		    << std::endl;

	  krn.setWeight(c,v);
	}

      for(c.i[1] = 0; c.i[1] < krn.space().axes[1].nbin; c.i[1]++)
	krn[c] /= vsum;

    }

  return krn;
}

double VSInstrumentResponseCalc::
getKernel(double offset, double emc, double erec) const
{
  VSEnergyResponseFn fn(1);
  VSAAlgebra::VecND param = fn.param();

  VSNSpace::Point p(2);
  
  p.x[0] = emc;
  p.x[1] = offset;
  
  m_krn_sigma1_nspace.interpolateWeight(p,param(0));
  m_krn_bias1_nspace.interpolateWeight(p,param(1));
  m_krn_sigma2_nspace.interpolateWeight(p,param(2));
  m_krn_bias2_nspace.interpolateWeight(p,param(3));
  m_krn_alpha_nspace.interpolateWeight(p,param(4));

  VSACoord::CoordND cnd(2);
  cnd[0] = erec;
  cnd[1] = emc;

  return fn.val(cnd,param);
}

// ============================================================================

bool VSInstrumentResponseCalc::load(const std::vector<unsigned>& nchan,
				 double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();

  if(m_options.irf_lookup_file.empty()) return false;
  else return load(m_options.irf_lookup_file,nchan,zn_deg,az_deg,ped_rms);
}

bool VSInstrumentResponseCalc::load(const VSTime& date,
				 const std::vector<unsigned>& nchan,
				 double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();

  return VSLTLibraryCalc::load(m_options.irf_lookup_file,date,
			       nchan,zn_deg,az_deg,ped_rms);
}

bool VSInstrumentResponseCalc::load(const std::string& lookup_file,
				    const std::vector<unsigned>& nchan,
				    double zn_deg, double az_deg, 
				    double ped_rms)
{
  clear();

  std::string name = lookup_file;
  VSFileUtility::expandFilename(name);
  
  if(!VSFileUtility::isFile(name))
    {
      std::cerr << std::string(__PRETTY_FUNCTION__) + ": " << std::endl
		<< "Failed to load lookup table: "
		<< name << std::endl;
      std::cerr << "File does not exist." << std::endl;
      return false;
    }
  else if(!VSOctaveH5ReaderStruct::isHDF5(name))
    {
      std::cerr << std::string(__PRETTY_FUNCTION__) + ": " << std::endl
		<< "Failed to load lookup table: "
		<< name << std::endl;
      std::cerr << "File is not HDF5." << std::endl;
      return false;
    }

  try
    {
      std::cout << "Loading Instrument Response Library: "
		<< name << std::endl;

      typedef VSScaledParameterLibraryWriter SPL;
      VSOctaveH5Reader* file = new VSOctaveH5Reader(name);
      VSScaledParameterLibraryReader* library = 
	new VSScaledParameterLibraryReader(file);

      load(library,zn_deg,az_deg,ped_rms);

      delete library;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr
	<< "Caught instance of VSOctaveH5Exception loading " << std::endl
	<< "instrument response lookup tables from:" << name << std::endl
	<< e.message() << std::endl;
      exit(EXIT_FAILURE);
    }

  m_has_lt = true;
  return true;
}

void VSInstrumentResponseCalc::load(VSScaledParameterLibraryReader* reader,
				    double zenith_deg,
				    double azimuth_deg,
				    double ped_dev)
{
  typedef VSScaledParameterLibraryWriter SPL;

  VSNSpace* effarea = 
    reader->readAndInterpolate(SPL::SPS_EFFECTIVE_AREA,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* psf_sigma1 = 
    reader->readAndInterpolate(SPL::SPS_PSF_SIGMA1,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* psf_sigma2 = 
    reader->readAndInterpolate(SPL::SPS_PSF_SIGMA2,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* psf_alpha = 
    reader->readAndInterpolate(SPL::SPS_PSF_ALPHA,
			       zenith_deg, azimuth_deg, 0, ped_dev);

  VSNSpace* krn_sigma1 = 
    reader->readAndInterpolate(SPL::SPS_KERNEL_SIGMA1,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* krn_bias1 = 
    reader->readAndInterpolate(SPL::SPS_KERNEL_BIAS1,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* krn_sigma2 = 
    reader->readAndInterpolate(SPL::SPS_KERNEL_SIGMA2,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* krn_bias2 = 
    reader->readAndInterpolate(SPL::SPS_KERNEL_BIAS2,
			       zenith_deg, azimuth_deg, 0, ped_dev);
  VSNSpace* krn_alpha = 
    reader->readAndInterpolate(SPL::SPS_KERNEL_ALPHA,
			       zenith_deg, azimuth_deg, 0, ped_dev);


  m_effarea_nspace = *effarea;
  m_psf_sigma1_nspace = *psf_sigma1;
  m_psf_sigma2_nspace = *psf_sigma2;
  m_psf_alpha_nspace = *psf_alpha;
  m_krn_sigma1_nspace = *krn_sigma1;
  m_krn_bias1_nspace  = *krn_bias1;
  m_krn_sigma2_nspace = *krn_sigma2;
  m_krn_bias2_nspace  = *krn_bias2;
  m_krn_alpha_nspace  = *krn_alpha;

  delete effarea;
  delete psf_sigma1;
  delete psf_sigma2;
  delete psf_alpha;
  delete krn_sigma1;
  delete krn_bias1;
  delete krn_sigma2;
  delete krn_bias2;
  delete krn_alpha;
}

bool VSInstrumentResponseCalc::load(VSOctaveH5ReaderStruct* reader)
{
  VSNSpaceOctaveH5IO io;

  bool status = true;

  status &= io.readHistogram(reader->readStruct("effarea"),m_effarea_nspace);
  status &= 
    io.readHistogram(reader->readStruct("psf_sigma1"),m_psf_sigma1_nspace);
  status &= 
    io.readHistogram(reader->readStruct("psf_sigma2"),m_psf_sigma2_nspace);
  status &= 
    io.readHistogram(reader->readStruct("psf_alpha"),m_psf_alpha_nspace);
  status &= 
    io.readHistogram(reader->readStruct("krn_sigma1"),m_krn_sigma1_nspace);
  status &= 
    io.readHistogram(reader->readStruct("krn_bias1"),m_krn_bias1_nspace);
  status &= 
    io.readHistogram(reader->readStruct("krn_sigma2"),m_krn_sigma2_nspace);
  status &= 
    io.readHistogram(reader->readStruct("krn_bias2"),m_krn_bias2_nspace);
  status &= 
    io.readHistogram(reader->readStruct("krn_alpha"),m_krn_alpha_nspace);

  return status;
}

void VSInstrumentResponseCalc::save(VSOctaveH5WriterStruct* writer) const
{
  VSNSpaceOctaveH5IO io;

  io.writeHistogram(writer->writeStruct("effarea"),m_effarea_nspace);
  io.writeHistogram(writer->writeStruct("psf_sigma1"),m_psf_sigma1_nspace);
  io.writeHistogram(writer->writeStruct("psf_sigma2"),m_psf_sigma2_nspace);
  io.writeHistogram(writer->writeStruct("psf_alpha"),m_psf_alpha_nspace);
  io.writeHistogram(writer->writeStruct("krn_sigma1"),m_krn_sigma1_nspace);
  io.writeHistogram(writer->writeStruct("krn_bias1"),m_krn_bias1_nspace);
  io.writeHistogram(writer->writeStruct("krn_sigma2"),m_krn_sigma2_nspace);
  io.writeHistogram(writer->writeStruct("krn_bias2"),m_krn_bias2_nspace);
  io.writeHistogram(writer->writeStruct("krn_alpha"),m_krn_alpha_nspace);
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSInstrumentResponseCalc::
configure(VSOptions& options, const std::string& profile, 
	  const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"irf_lookup"),
			s_default_options.irf_lookup_file,
			"Set the simulation library file for effective "
			"area, energy response matrix, and psf.");
}
