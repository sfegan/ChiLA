//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMuonAnalysis.cpp

  Try to fit a ring to an image, cut on the image parameters, do the muon
  analysis and return the muon paramaters

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/14/2007

  $Id: VSMuonAnalysis.cpp,v 3.3 2008/04/04 00:12:42 matthew Exp $

*/

#include<cmath>

#include<VSMuonAnalysis.hpp>
#include<VSAMath.hpp>
#include<VSABracketedMonotonicRootFinder.hpp>
#include<fast_alloc.hpp>

using namespace VERITAS;

VSMuonAnalysis::RawDataFetcher::~RawDataFetcher()
{
  // nothing to see here
}

static inline double E0(double xi)
{
  static const double C = 2.0/M_PI;
  if(fabs(xi-1.0)<DBL_EPSILON)
    return C;
  else if(xi<1)
    return C*VSAMath::elliptic2Complete(xi);
  else
    {
      const double xi_inv = 1.0/xi;
      return C*(xi*VSAMath::elliptic2Complete(xi_inv)
		- (xi-xi_inv)*VSAMath::elliptic1Complete(xi_inv));
    }
}

static inline double F0(double xi)
{
  static const double C = 2.0/M_PI;
  if(fabs(xi-1.0)<DBL_EPSILON)
    return std::numeric_limits<double>::infinity();
  else if(xi<1)
    return C*VSAMath::elliptic1Complete(xi);
  else
    {
      const double xi_inv = 1.0/xi;
      return C*VSAMath::elliptic1Complete(xi_inv)*xi_inv;
    }
}

static inline double E2(double xi)
{
  double xi2 = xi*xi;
  if(xi<DBL_EPSILON)
    return 1.0;
  else if(fabs(xi-1.0)<DBL_EPSILON)
    return 4.0/(3.0*M_PI);
  else 
    return ((4.0*xi2-2.0)*E0(xi)-2.0*(xi2-1.0)*F0(xi))/(3.0*xi2);
}

static inline double E4(double xi)
{
  double xi2 = xi*xi;
  double xi4 = xi2*xi2;
  if(xi<DBL_EPSILON)
    return 1.0;
  else if(fabs(xi-1.0)<DBL_EPSILON)
    return 16.0/(15.0*M_PI);
  else
    return ((8.0*xi4 - 3.0*xi2 - 2.0)*E0(xi) 
	    + 2.0*(-2.0*xi4 + xi2 + 1.0)*F0(xi))*(8.0/45.0)/xi4;
}

static inline
void E024(double xi, double& e0, double& e2, double& e4)
{
  vsassert(xi>=0);
  if(xi<DBL_EPSILON)
    {
      e0 = 1.0;
      e2 = 1.0;
      e4 = 1.0;
    }
  else if(fabs(xi-1.0)<DBL_EPSILON)
    {
      e0 = 2.0/M_PI;
      e2 = 4.0/(3.0*M_PI);
      e4 = 16.0/(15.0*M_PI);
    }
  else
    {
      static const double C = 2.0/M_PI;
      double f0;
      if(xi<1)
	{
	  f0 = C*VSAMath::elliptic1Complete(xi);
	  e0 = C*VSAMath::elliptic2Complete(xi);
	}
      else
	{
	  const double xi_inv = 1.0/xi;
	  f0 = C*VSAMath::elliptic1Complete(xi_inv)*xi_inv;
	  e0 = C*VSAMath::elliptic2Complete(xi_inv)*xi - (xi*xi-1)*f0;
	}
      const double xi2 = xi*xi;
      const double xi4 = xi2*xi2;
      e2 = ((4.0*xi2-2.0)*e0-2.0*(xi2-1.0)*f0)/(3.0*xi2);
      e4 = ((8.0*xi4 - 3.0*xi2 - 2.0)*e0 
	    + 2.0*(-2.0*xi4 + xi2 + 1.0)*f0)*(8.0/45.0)/xi4;
    }
}

class XiFinderFunction
{
public:
  XiFinderFunction(double d0_r0): y(2.0*d0_r0) {}
  double operator() (double xi) { return xi/E0(xi) - y; }
private:
  double y;
};

class XiFinderReciprocalFunction
{
public:
  XiFinderReciprocalFunction(double d0_r0): y(2.0*d0_r0) {}
  double operator() (double xi_inv) 
  { if(xi_inv == 0)return 2.0-y; else return xi_inv/E0(1.0/xi_inv) - y; }
private:
  double y;
};

static inline double
bracketToFindXi(const double dc_r0, const double tol = 1e-8)
{
  if(dc_r0 <= M_PI/4)
    {
      XiFinderFunction func(dc_r0);
      return VSAMath::findBracketedMonotonicRoot(func,0.0,1.0,tol);
    }
  else
    {
      XiFinderReciprocalFunction func(dc_r0);
      return 1.0/VSAMath::findBracketedMonotonicRoot(func,0.0,1.0,tol);
    }
}

#ifdef MUON_TEST_NON_UNIFORMITY

static inline double
Lambda(const double xi, const double beta, 
       const double e0, const double e2, const double e4)
{
  static const double TWO_THIRD = 2.0/3.0;
  return e0 + beta*TWO_THIRD*(e0 + xi*xi*e2);
}

static inline double
Phi(const double xi, const double beta, 
    const double e0, const double e2, const double e4)
{
  if(xi>1)return 0.5*(1+beta)/(xi*Lambda(xi,beta,e0,e2,e4));
  else return 0.5*(1+beta*xi*xi)*xi/Lambda(xi,beta,e0,e2,e4);
}

static inline double PhiXiOne(const double beta)
{
  return (9.0*M_PI/4.0)*(1.0+beta)/(9.0+10.0*beta);
}

static inline double 
Beta(const double Wtilde, const double xi,
     const double e0, const double e2, const double e4)
{
  const double xi2=xi*xi;
  return 1.5*(e2 - Wtilde*e0)/(Wtilde*(e0 + xi2*e2) - (e2 + 1.5*xi2*e4));
}

static inline double Q(const double phi, const double xi, const double beta)
{
  const double xx   = xi*xi;
  const double c    = cos(phi);
  const double s    = sin(phi);
  const double xxss = xx*s;
  return (xxss>1.0) ? 0.0 : (sqrt(1-xxss)*(1.0+2.0/3.0*beta*(1.0+2.0*xxss))
			     + xi*c*
			     + 2.0*beta*xx*xi*(1.0 - 2.0/3.0*c*c)*c);
}

class PhiFunction
{
public:
  PhiFunction(double d0_r0, double _beta): y(d0_r0), beta(_beta) {}
  double operator() (double xi) 
  { double e0,e2,e4; E024(xi,e0,e2,e4); return Phi(xi,beta,e0,e2,e4) - y; }
private:
  double y;
  double beta;
};

class PhiReciprocalFunction
{
public:
  PhiReciprocalFunction(double d0_r0, double _beta): y(d0_r0), beta(_beta) {}
  double operator() (double xi_inv) 
  { if(xi_inv == 0)return 1.0-y;
    else { double xi,e0,e2,e4; xi=1.0/xi_inv; E024(xi,e0,e2,e4);
      return Phi(xi,beta,e0,e2,e4) - y; } }
private:
  double y;
  double beta;
};

static inline double
bracketToFindXiBeta(const double beta, 
		    const double dc_r0, const double tol = 1e-8)
{
  if(dc_r0 <= PhiXiOne(beta))
    {
      PhiFunction func(dc_r0,beta);
      return VSAMath::findBracketedMonotonicRoot(func,0.0,1.0,tol);
    }
  else
    {
      PhiReciprocalFunction func(dc_r0,beta);
      return 1.0/VSAMath::findBracketedMonotonicRoot(func,0.0,1.0,tol);
    }
}

#endif

bool VSMuonAnalysis::analyze(VSMuonAnalysisDatum& datum,
			     unsigned ievent,
			     const VSAReconstruction::ScopeImage& image,
			     RawDataFetcher* raw_image_fetcher) const
{
  unsigned nimage = image.pixels.size();
  if(nimage <= m_nimage_cut)
    return false;

  // --------------------------------------------------------------------------
  // Fit pass one -- all channels included
  // --------------------------------------------------------------------------

  double sum_ni   = 0;
  double sum_nixi = 0;
  double sum_niyi = 0;
  for(unsigned iimage=0;iimage<nimage;iimage++)
    {
      const unsigned ipixel = image.pixels[iimage].j;
      const double xi = m_scope.pij[ipixel].x();
      const double yi = m_scope.pij[ipixel].y();
      const double ni = image.pixels[iimage].nij;
      sum_nixi += ni*xi;
      sum_niyi += ni*yi;
      sum_ni   += ni;
    }

  double fit_chi2;
  double fit_signal;
  double fit_x0;
  double fit_y0;
  double fit_r0;
  double fit_cx;
  double fit_cy;

  m_fitter->
    fitRing(image, sum_nixi/sum_ni, sum_niyi/sum_ni,
	    fit_chi2, fit_signal, fit_x0, fit_y0, fit_r0, 
	    fit_cx, fit_cy);

  double fit_dist0 = sqrt(fit_x0*fit_x0 + fit_y0*fit_y0);
  double fit_mom0  = sqrt(fit_cx*fit_cx + fit_cy*fit_cy);
  double fit_rms   = sqrt(fit_chi2/fit_signal);

  // --------------------------------------------------------------------------
  // Cuts
  // --------------------------------------------------------------------------

  if((fit_r0       < m_radius_min_cut)
     ||(fit_r0     > m_radius_max_cut)
     ||(fit_rms    > m_rms_max_cut)
     ||(fit_dist0  > m_ring_edge_dist_max_cut)
     ||(fit_mom0   > fit_r0*m_centroid_radius_ratio_max_cut))
    return false;

#if 1 // defined TWOPASS
  // --------------------------------------------------------------------------
  // Fit pass two -- exclude channels outside of the ring
  // --------------------------------------------------------------------------

  const double r_inner = fit_r0 - 0.3/180.0*M_PI;
  const double r_outer = fit_r0 + 0.3/180.0*M_PI;

  sum_nixi = 0;
  sum_niyi = 0;
  sum_ni   = 0;

  // Look through channels to see if any need to be elimiated
  unsigned iimage=0;
  while(iimage<nimage)
    {
      const unsigned ipixel = image.pixels[iimage].j;
      const double xi = m_scope.pij[ipixel].x() - fit_x0;
      const double yi = m_scope.pij[ipixel].y() - fit_y0;
      const double ri = sqrt(xi*xi+yi*yi);
      if(ri<r_inner || ri>r_outer)break;
      const double ni = image.pixels[iimage].nij;
      sum_ni   += ni;
      sum_nixi += ni*m_scope.pij[ipixel].x();
      sum_niyi += ni*m_scope.pij[ipixel].y();
      iimage++;
    }

  // If we get all the way to the end then we don't need to redo the fit
  if(iimage<nimage)
    {
      iimage++;

      VSAReconstruction::ScopeImage pass2_image;
      pass2_image.pixels.reserve(nimage);

      for(unsigned jimage=0;jimage<iimage;jimage++)
	pass2_image.pixels.push_back(image.pixels[jimage]);

      while(iimage<nimage)
	{
	  const unsigned ipixel = image.pixels[iimage].j;
	  const double xi = m_scope.pij[ipixel].x() - fit_x0;
	  const double yi = m_scope.pij[ipixel].y() - fit_y0;
	  const double ri = sqrt(xi*xi+yi*yi);
	  if(ri>=r_inner && ri<=r_outer)
	    {
	      const double ni = image.pixels[iimage].nij;
	      pass2_image.pixels.push_back(image.pixels[iimage]);
	      sum_ni   += ni;
	      sum_nixi += ni*m_scope.pij[ipixel].x();
	      sum_niyi += ni*m_scope.pij[ipixel].y();
	    }
	  iimage++;
	}

      nimage = pass2_image.pixels.size(); 
      if(nimage <= m_nimage_cut)
	return false;

      m_fitter->
	fitRing(pass2_image, fit_x0, fit_y0,
		fit_chi2, fit_signal, fit_x0, fit_y0, fit_r0, fit_cx, fit_cy);

      fit_cx   = sum_nixi/sum_ni - fit_x0;
      fit_cy   = sum_niyi/sum_ni - fit_y0;

      fit_mom0 = sqrt(fit_cx*fit_cx + fit_cy*fit_cy);
      fit_rms  = sqrt(fit_chi2/fit_signal);

      // ----------------------------------------------------------------------
      // Cut again
      // ----------------------------------------------------------------------

      if((fit_r0       < m_radius_min_cut)
	 ||(fit_r0     > m_radius_max_cut)
	 ||(fit_rms    > m_rms_max_cut)
	 ||(fit_dist0  > m_ring_edge_dist_max_cut)
	 ||(fit_mom0   > fit_r0*m_centroid_radius_ratio_max_cut))
	return false;
    }

#endif

  // --------------------------------------------------------------------------
  // Fill in the "common" parameters
  // --------------------------------------------------------------------------

  datum.ievent = ievent;
  datum.chi2     = fit_chi2;
  datum.x0       = fit_x0;
  datum.y0       = fit_y0;
  datum.r0       = fit_r0;

  // --------------------------------------------------------------------------
  // Fill in the "cleaned" parameters
  // --------------------------------------------------------------------------

  datum.c_nimage = nimage;
  datum.c_signal = fit_signal;
  datum.c_rms    = fit_rms;
  datum.c_cx     = fit_cx;
  datum.c_cy     = fit_cy;

  const double c_mom0_r0 = fit_mom0/fit_r0;

  datum.c_xi     = 1.0;
  datum.c_U0_r0  = 0;

  if(c_mom0_r0<1)
    {
      // Calculate XI and U0 --------------------------------------------------

      try
	{
	  datum.c_xi = bracketToFindXi(c_mom0_r0, 1e-8);
	}
      catch(const VSAMath::FunctionNotMonotonic& x)
	{
	  std::cout 
	    << ievent << ": non monotonic function exception d/theta (c)="
	    << c_mom0_r0 << std::endl;
	  return false;
	}
      datum.c_U0_r0 = datum.c_signal/E0(datum.c_xi);
    }

  // --------------------------------------------------------------------------
  // Fill in the "raw" parameters
  // --------------------------------------------------------------------------

  if(raw_image_fetcher)
    {
      const double r_inner = fit_r0 - m_raw_ring_width;
      const double r_outer = fit_r0 + m_raw_ring_width;

      nimage            = 0;
      sum_ni            = 0;
      sum_nixi          = 0;
      sum_niyi          = 0;
      double sum_niri2  = 0;
      double sum_nitot  = 0;

#ifdef MUON_TEST_NON_UNIFORMITY
      double sum_nixi2  = 0;
      double sum_niyi2  = 0;
      double sum_nixyi  = 0;
#endif

      unsigned npixel = m_scope.pij.size();
      double* raw_image = FASTCALLOC(double,npixel);
      double* raw_image_total = FASTCALLOC(double,npixel);
      unsigned* ring_ipixel = FASTCALLOC(unsigned,npixel);
      unsigned nring_ipixel = 0;

      // Use fetcher to get raw data ------------------------------------------

      if(raw_image_fetcher->fetchRawData(npixel, raw_image, raw_image_total))
	{
	  for(unsigned ipixel=0; ipixel<npixel; ipixel++)
	    {
	      // Loop through pixels in ring accumulating moments -------------
	      // and recording pixel IDs for later use

	      const double xi = m_scope.pij[ipixel].x() - fit_x0;
	      const double yi = m_scope.pij[ipixel].y() - fit_y0;
	      const double ri = sqrt(xi*xi+yi*yi);
	      if(ri>r_inner && ri<r_outer)
		{
		  const double ni = raw_image[ipixel];
		  if(!std::isnan(ni))
		    {
		      nimage++;
		      const double nixi = ni*xi;
		      sum_nixi  += nixi;
#ifdef MUON_TEST_NON_UNIFORMITY
		      sum_nixi2 += nixi*xi;
		      sum_nixyi += nixi*yi;
#endif
		      const double niyi = ni*yi;
		      sum_niyi  += niyi;
#ifdef MUON_TEST_NON_UNIFORMITY
		      sum_niyi2 += niyi*yi;
#endif
		      sum_niri2 += ni*(ri-fit_r0)*(ri-fit_r0);
		      sum_ni    += ni;
		      sum_nitot += raw_image_total[ipixel];
		      ring_ipixel[nring_ipixel++]=ipixel;
		    }
		  else
		    {
		      ring_ipixel[nring_ipixel++]=ipixel+npixel;
		    }
		}
	    }

	  datum.r_nimage                  = nimage;
	  datum.r_signal                  = sum_ni;
	  datum.r_rms                     = sqrt(sum_niri2/sum_ni);
	  datum.r_cx                      = sum_nixi/sum_ni;
	  datum.r_cy                      = sum_niyi/sum_ni;
	  datum.r_xi                      = -1.0;
	  datum.r_U0_r0                   = 0;
	  datum.r_U0_r0_corr              = 0;
#ifdef MUON_TEST_NON_UNIFORMITY
	  datum.nu_xi                     = -1.0;
	  datum.nu_U0_r0                  = 0;
	  datum.nu_U0_r0_corr             = 0;
#endif

	  const double r_mom0 = 
	    sqrt(datum.r_cx*datum.r_cx+datum.r_cy*datum.r_cy);
	  const double r_mom0_r0 = r_mom0/fit_r0;

	  if(r_mom0_r0<1)
	    {
	      // Calculate XI and U0 ------------------------------------------
	      
	      try
		{
		  datum.r_xi = bracketToFindXi(r_mom0_r0, 1e-8);
		}
	      catch(const VSAMath::FunctionNotMonotonic& x)
		{
		  std::cout 
		    << ievent 
		    << ": non monotonic function exception d/theta (r)="
		    << r_mom0_r0 << std::endl;
		  FASTFREE(ring_ipixel);
		  FASTFREE(raw_image_total);
		  FASTFREE(raw_image);
		  return false;
		}

	      double scale                = 1.0/E0(datum.r_xi);

	      datum.r_U0_r0               = datum.r_signal*scale;

	      // Attempt to correct for missing channels ----------------------

	      double present_weight        = 0;
	      double missing_weight        = 0;

	      for(unsigned imissing=0;imissing<nring_ipixel;imissing++)
		{
		  unsigned ipixel = ring_ipixel[imissing];
		  bool missing = false;
		  if(ipixel>=npixel)ipixel-=npixel, missing=true;

		  const double xi  = m_scope.pij[ipixel].x() - fit_x0;
		  const double yi  = m_scope.pij[ipixel].y() - fit_y0;
		  const double dri = sqrt(xi*xi+yi*yi) - fit_r0;
		  const double phi = atan2(yi,xi)-atan2(datum.r_cy,datum.r_cx);
		  const double cp  = cos(phi);
		  const double sp  = sin(phi);
		  const double xs2 = datum.r_xi*datum.r_xi*sp*sp;
		  const double phi_factor = 
		    xs2<1 ? sqrt(1-xs2)+datum.r_xi*cp : 0;
		  const double r_factor =
		    exp(-dri*dri/(2.0*datum.r_rms*datum.r_rms));
		  const double weight = phi_factor*r_factor;

		  if(missing)missing_weight += weight;
		  else present_weight += weight;
		}

	      datum.r_U0_r0_corr = 
		datum.r_U0_r0*(1.0+missing_weight/present_weight)
		* sum_nitot/sum_ni;

#ifdef MUON_TEST_NON_UNIFORMITY
	      // Estimate non`uniformty ---------------------------------------

	      datum.nu_xi                  = datum.r_xi;
	      datum.nu_U0_r0               = datum.r_U0_r0;
	      datum.nu_U0_r0_corr          = datum.r_U0_r0_corr;

	      if(fabs(m_nu_beta)>=std::numeric_limits<double>::epsilon())
		{
		  try
		    {
		      datum.nu_xi          =
			bracketToFindXiBeta(m_nu_beta, r_mom0_r0, 1e-8);
		    }
		  catch(const VSAMath::FunctionNotMonotonic& x)
		    {
		      std::cout 
			<< ievent 
			<< ": non monotonic function exception d/theta (nu)="
			<< r_mom0_r0 << std::endl;
		      FASTFREE(ring_ipixel);
		      FASTFREE(raw_image_total);
		      FASTFREE(raw_image);
		      datum.nu_xi          = -1;
		      datum.nu_U0_r0       = 0;
		      datum.nu_U0_r0_corr  = 0;
		      return true;
		    }

		  double e0;
		  double e2;
		  double e4;
		  E024(datum.nu_xi,e0,e2,e4);
		  scale                    = 
		    (1+m_nu_beta)/Lambda(datum.nu_xi,m_nu_beta,e0,e2,e4);

		  datum.nu_U0_r0           = datum.r_signal*scale;

		  // Attempt to correct for missing channels ------------------

		  present_weight           = 0;
		  missing_weight           = 0;

		  for(unsigned imissing=0;imissing<nring_ipixel;imissing++)
		    {
		      unsigned ipixel = ring_ipixel[imissing];
		      bool missing = false;
		      if(ipixel>=npixel)ipixel-=npixel, missing=true;

		      const double xi  = m_scope.pij[ipixel].x() - fit_x0;
		      const double yi  = m_scope.pij[ipixel].y() - fit_y0;
		      const double dri = sqrt(xi*xi+yi*yi) - fit_r0;
		      const double phi = 
			atan2(yi,xi)-atan2(datum.r_cy,datum.r_cx);
		      const double phi_factor = Q(phi,datum.nu_xi,m_nu_beta);
		      const double r_factor =
			exp(-dri*dri/(2.0*datum.r_rms*datum.r_rms));
		      const double weight = phi_factor*r_factor;

		      if(missing)missing_weight += weight;
		      else present_weight += weight;
		    }

		  datum.nu_U0_r0_corr = 
		    datum.nu_U0_r0*(1.0+missing_weight/present_weight)
		    * sum_nitot/sum_ni;
#endif
		}
	    }
	}
      else
	{
	  datum.r_nimage      = 0;
	  datum.r_signal      = 0;
	  datum.r_rms         = 0;
	  datum.r_cx          = 0;
	  datum.r_cy          = 0;
	  datum.r_xi          = -1.00;
	  datum.r_U0_r0       = 0;
	  datum.r_U0_r0_corr  = 0;
#ifdef MUON_TEST_NON_UNIFORMITY
	  datum.nu_xi         = -1.0;
	  datum.nu_U0_r0      = 0;
	  datum.nu_U0_r0_corr = 0;
#endif
	}

      FASTFREE(ring_ipixel);
      FASTFREE(raw_image_total);
      FASTFREE(raw_image);
    }
  else
    {
      datum.r_nimage          = 0;
      datum.r_signal          = 0;
      datum.r_rms             = 0;
      datum.r_cx              = 0;
      datum.r_cy              = 0;
      datum.r_xi              = -1.0;
      datum.r_U0_r0           = 0;
      datum.r_U0_r0_corr      = 0;
#ifdef MUON_TEST_NON_UNIFORMITY
      datum.nu_xi             = -1.0;
      datum.nu_U0_r0          = 0;
      datum.nu_U0_r0_corr     = 0;
#endif
    }

  return true;
}


#if 0
	      // Estimate non-uniformty -- old method -------------------------

	      const double Txx     = sum_nixi2/sum_ni - datum.r_cx*datum.r_cx;
	      const double Txy     = sum_nixyi/sum_ni - datum.r_cx*datum.r_cy;
	      const double Tyy     = sum_niyi2/sum_ni - datum.r_cy*datum.r_cy;

	      const double traceT  = Txx+Tyy;
	      const double detT    = Txx*Tyy-Txy*Txy;
	      
	      const double r02     = fit_r0*fit_r0;
	      const double r04     = r02*r02;
	      const double trT_r02 = traceT/r02;

	      const double Wdiff   = sqrt(trT_r02*trT_r02 - 4.0*detT/r04);
	      const double Wpos    = trT_r02 + Wdiff;

	      double last_xi       = -1;
	      double last_beta     = -1;
 	      double curr_xi       = datum.r_xi;
	      double curr_beta     = 0;

	      for(unsigned i=0;i<100;i++)
		{
		  last_xi = curr_xi;
		  last_beta = curr_beta;

		  double e0;
		  double e2;
		  double e4;
		  E024(curr_xi, e0, e2, e4);

		  curr_beta = Beta(Wpos, curr_xi, e0, e2, e4);
		  if(curr_beta<-0.5)curr_beta=-0.5;

		  try
		    {
		      curr_xi = bracketToFindXiBeta(curr_beta, r_mom0_r0);
		    }
		  catch(const VSAMath::FunctionNotMonotonic& x)
		    {
		      std::cout 
			<< "Caught non-monotonic function exception: beta="
			<< curr_beta << " r/theta=" << r_mom0_r0 << '\n';
		      break;
		    }

		  if((fabs(curr_xi - last_xi)<1e-5)
		     &&(fabs(curr_beta - last_beta)<1e-5))
		    {
		      datum.nu_xi = curr_xi;
		      datum.nu_beta = curr_beta;
		      break;
		    }

		  //std::cerr << i << ' ' << curr_xi << ' ' << curr_beta << '\n';
		}
#endif


#ifdef TESTMAIN_0

#include<iostream>
#include<iomanip>

int main(int argc, char** argv)
{
  std::cout << std::fixed << std::setprecision(8);
  for(double xi=0;xi<=9.0;xi+=0.25)
    {
      std::cout << xi << ' ' << F0(xi) << ' ' << E0(xi) << ' '
		<< E2(xi) << ' ' << E4(xi) << '\n';
      double e0;
      double e2;
      double e4;
      E024(xi, e0, e2, e4);
      std::cout << xi << ' ' << 0.0 << ' ' 
		<< e0 << ' ' << e2 << ' ' << e4 << '\n';
    }
}

#endif

#ifdef TESTMAIN_1

#include<iostream>
#include<iomanip>

int main(int argc, char** argv)
{
  std::cout << std::fixed << std::setprecision(8);
  double beta[] = {-0.5, -0.25, 0, 0.5, 1.0 };
  for(double xi=0;xi<=10.0;xi+=0.01)
    {
      double e0;
      double e2;
      double e4;
      E024(xi, e0, e2, e4);
      std::cout << xi;
      for(unsigned i=0;i<sizeof(beta)/sizeof(*beta);i++)
	std::cout << ' ' << Lambda(xi,beta[i],e0,e2,e4)/(1+beta[i]);
      std::cout << '\n';
    }
}

#endif

#ifdef TESTMAIN

#include<iostream>
#include<iomanip>

int main(int argc, char** argv)
{
  std::cout << std::fixed << std::setprecision(8);
  double beta[] = {-0.5, -0.25, 0, 0.5, 1.0 };
  for(double xi=0;xi<=10.0;xi+=0.01)
    {
      double e0;
      double e2;
      double e4;
      E024(xi, e0, e2, e4);
      std::cout << xi;
      for(unsigned i=0;i<sizeof(beta)/sizeof(*beta);i++)
	std::cout << ' ' << Phi(xi,beta[i],e0,e2,e4);
      std::cout << '\n';
    }
}

#endif
