//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAMath.hpp
  Various math functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSAMATH_HPP
#define VSAMATH_HPP

#define VSA_ASSERT

#include<cmath>
#include<cstdlib>
#include<cfloat>
#include<limits>

#ifdef VSA_ASSERT
#include<cassert>
#endif

#include <VSAssert.hpp>

namespace VERITAS
{
  namespace VSAMath
  {
    template<typename T>
    struct numeric_constants
    {
      static T bessel_zero(int m, unsigned n)       
      { 
	vsassert(m < 4 && n < 20);
	return static_cast<T>(bessel_zeros[std::abs(m)][n]); 
      }
      ///  Constant @f$ \pi @f$.
      static T pi() throw()
      { return static_cast<T>(3.1415926535897932384626433832795029L); }
      ///  Constant @f$ \pi / 2 @f$.
      static T pi_2() throw()
      { return static_cast<T>(1.5707963267948966192313216916397514L); }
      ///  Constant @f$ \pi / 3 @f$.
      static T pi_3() throw()
      { return static_cast<T>(1.0471975511965977461542144610931676L); }
      ///  Constant @f$ \pi / 4 @f$.
      static T pi_4() throw()
      { return static_cast<T>(0.7853981633974483096156608458198757L); }
      ///  Constant @f$ 1 / \pi @f$.
      static T _1_pi() throw()
      { return static_cast<T>(0.3183098861837906715377675267450287L); }
      ///  Constant @f$ 2 / \sqrt(\pi) @f$.
      static T _2_sqrtpi() throw()
      { return static_cast<T>(1.1283791670955125738961589031215452L); }
      ///  Constant @f$ \sqrt(2) @f$.
      static T sqrt2() throw()
      { return static_cast<T>(1.4142135623730950488016887242096981L); }
      ///  Constant @f$ \sqrt(3) @f$.
      static T sqrt3() throw()
      { return static_cast<T>(1.7320508075688772935274463415058723L); }
      ///  Constant @f$ \sqrt(\pi/2) @f$.
      static T sqrtpio2() throw()
      { return static_cast<T>(1.2533141373155002512078826424055226L); }
      ///  Constant @f$ 1 / sqrt(2) @f$.
      static T sqrt1_2() throw()
      { return static_cast<T>(0.7071067811865475244008443621048490L); }
      ///  Constant @f$ \log(\pi) @f$.
      static T lnpi() throw()
      { return static_cast<T>(1.1447298858494001741434273513530587L); }
      ///  Constant Euler's constant @f$ \gamma_E @f$.
      static T gamma_e() throw()
      { return static_cast<T>(0.5772156649015328606065120900824024L); }
      ///  Constant Euler-Mascheroni @f$ e @f$
      static T euler() throw()
      { return static_cast<T>(2.7182818284590452353602874713526625L); }

      static const double bessel_zeros[4][20];
    };
    
    template< typename T >
    const double numeric_constants<T>::bessel_zeros[4][20] = 
  {
    {2.404825557695773,5.520078110286311,8.653727912911011,11.79153443901428,
     14.93091770848779,18.07106396791092,21.21163662987926,24.35247153074930,
     27.49347913204025,30.63460646843198,33.77582021357357,36.91709835366404,
     40.05842576462824,43.19979171317673,46.34118837166181,49.48260989739782,
     52.62405184111500,55.76551075501998,58.90698392608094,62.04846919022717},
    {3.831705970207512,7.015586669815613,10.17346813506272,13.32369193631422,
     16.47063005087763,19.61585851046824,22.76008438059277,25.90367208761838,
     29.04682853491686,32.18967991097440,35.33230755008387,38.47476623477162,
     41.61709421281445,44.75931899765282,47.90146088718545,51.04353518357151,
     54.18555364106132,57.32752543790101,60.46945784534749,63.61135669848123},
    {5.135622301840683,8.417244140399855,11.61984117214906,14.79595178235126,
     17.95981949498783,21.11699705302185,24.27011231357310,27.42057354998456,
     30.56920449551640,33.71651950922270,36.86285651128381,40.00844673347819,
     43.15345377837146,46.29799667723692,49.44216411041687,52.58602350681596,
     55.72962705320114,58.87301577261216,62.01622235921765,65.15927319075780},
    {6.380161895923984,9.761023129981667,13.01520072169843,16.22346616031877,
     19.40941522643501,22.58272959310444,25.74816669929498,28.90835078092176,
     32.06485240709771,35.21867073861011,38.37047243475694,41.52071967040678,
     44.66974311661725,47.81778569153330,50.96502990620518,54.11161556982187,
     57.25765160449901,60.40322413847212,63.54840217856721,66.69324166737268}
  };

    // ------------------------------------------------------------------------
    // Roots of Quadratic and Cubic
    // ------------------------------------------------------------------------

    unsigned realRoots2(const double a, const double b, const double c,
			double roots[2]);
    unsigned realRoots3(const double _a, const double _b, 
			const double _c, const double _d,
			double roots[3]);

    // ------------------------------------------------------------------------
    // Chi Squared Probability
    // ------------------------------------------------------------------------
    double chi2p(double chi2, int ndf);
    double chi2q(double p, int ndf);

    double norm_quantile(double p);

    // ------------------------------------------------------------------------
    // Error Function
    // ------------------------------------------------------------------------
    double erf(double x);

    double erfc(double x);

    double erf_inverse(double x);

    // ------------------------------------------------------------------------
    // Gamma Function
    // ------------------------------------------------------------------------
    double lnGammaLanczos(double z);

    inline double lgamma(double z) 
    { 
      if(z > 0.5)
	return lnGammaLanczos(z); 
      else
	{
	  const double sin_fact
	    = std::abs(std::sin(numeric_constants<double>::pi()*z));
	  if (sin_fact == double(0))
	    std::__throw_domain_error(__N("Argument is nonpositive integer "
					  "in __log_gamma"));
	  return numeric_constants<double>::lnpi()
	    - std::log(sin_fact) - lnGammaLanczos(double(1) - z);
	}
    }

    inline double gamma(double z)
    {
      double v = lgamma(z);
      return exp(v);
    }

    inline void gamma_temme(const double mu,
			    double & gam1, double & gam2, 
			    double & gampl, double & gammi)
    {
      gampl = double(1) / gamma(double(1) + mu);
      gammi = double(1) / gamma(double(1) - mu);

      if (std::abs(mu) < std::numeric_limits<double>::epsilon())
        gam1 = -double(numeric_constants<double>::gamma_e());
      else
        gam1 = (gammi - gampl) / (double(2) * mu);

      gam2 = (gammi + gampl) / (double(2));

      return;
    }
    
    // ------------------------------------------------------------------------
    // Incomplete Gamma Function
    // ------------------------------------------------------------------------
    double gamma(double a, double x);
    double gamma_pseries(double a, double x);
    double gamma_cfe(double a, double x);

    // ------------------------------------------------------------------------
    // Beta Function
    // ------------------------------------------------------------------------
    inline double beta(double x, double y) 
    {   
      return std::exp(lnGammaLanczos(x)+lnGammaLanczos(y)-
		      lnGammaLanczos(x+y));   
    }

    // ------------------------------------------------------------------------
    // Incomplete Beta Function
    // ------------------------------------------------------------------------
    double beta( double a, double b, double x );
    double beta_pseries( double a, double b, double x );
    double beta_cfe1( double a, double b, double x );
    double beta_cfe2( double a, double b, double x );
    
    // ------------------------------------------------------------------------
    // Carlson forms of integrals
    // ------------------------------------------------------------------------

    // Carlson's elliptic integral of the first kind
    double ellipticRFCarlson(double x, double y, double z);

    inline double ellipticRFCarlsonZZero(const double x, const double y)
    {
      const double ERRTOL = 2.7*sqrt(DBL_EPSILON);
      double xm = sqrt(x);
      double ym = sqrt(y);
      while(fabs(xm-ym) > ERRTOL*fabs(xm))
	{
	  const double xm_next = 0.5*(xm+ym);
	  ym = sqrt(xm*ym);
	  xm = xm_next;
	}
      return M_PI/(xm+ym);
    }

    // Carlson's elliptic integral of the second kind
    double ellipticRDCarlson(double x, double y, double z);

    inline double ellipticRGCarlsonZZero(const double x, const double y)
    {
      const double ERRTOL = 2.7*sqrt(DBL_EPSILON);
      double xm = sqrt(x);
      double ym = sqrt(y);
      double fac = 0.5;
      double sum = (xm+ym)/2;
      sum *= sum;
      while(fabs(xm-ym) > ERRTOL*fabs(xm))
	{
	  const double xm_next = 0.5*(xm+ym);
	  ym = sqrt(xm*ym);
	  xm = xm_next;
	  const double dxy = xm-ym;
	  sum -= fac*dxy*dxy;
	  fac *= 2.0;
	}
      return 0.5*sum*M_PI/(xm+ym);
    }

    // ------------------------------------------------------------------------
    // Jacobi elliptic functions
    // ------------------------------------------------------------------------

    void ellipticJacobiFunctions(const double _u, const double k,
				 double &sn, double& cn, double& dn);

    // ------------------------------------------------------------------------
    // Elliptic integrals of the first kind
    // ------------------------------------------------------------------------


    // Legendre elliptic integral of the first kind
    inline double elliptic1Legendre(const double phi, const double k)
    {
      const double s    = sin(phi); 
      const double c    = cos(phi);      
      const double sskk = s*s*k*k;
      return s*ellipticRFCarlson(c*c,1.0-sskk,1.0);
    }

    // Complete elliptic integral of the first kind
    inline double elliptic1Complete(const double k)
    { 
      return ellipticRFCarlsonZZero(1.0-k*k,1.0);
    }

    // ------------------------------------------------------------------------
    // Elliptic integrals of the second kind
    // ------------------------------------------------------------------------

    // Legendre elliptic integral of the second kind
    inline double elliptic2Legendre(const double phi, const double k)
    {
      const double s    = sin(phi);
      const double c    = cos(phi);
      const double cc   = c*c;
      const double sskk = s*s*k*k;
      const double q    = 1.0-sskk;
      return s*(ellipticRFCarlson(cc,q,1.0) -
		sskk*ellipticRDCarlson(cc,q,1.0)/3.0);
    }

    // Complete elliptic integral of the second kind
    inline double elliptic2Complete(const double k)
    {
#if 0
      const double k2 = k*k;
      const double q  = 1.0-k2;
      return (ellipticRFCarlson(0.0,q,1.0)
	      - k2*ellipticRDCarlson(0.0,q,1.0)/3.0);
#endif
      return 2.0*ellipticRGCarlsonZZero(1.0-k*k,1.0);
    }

    // ------------------------------------------------------------------------
    // Bessel and Neumann functions
    // ------------------------------------------------------------------------
    void besselJN(const double nu, const double x,
		  double & jnu, double & nnu, double & jpnu, double & npnu);

    inline double besselJ(const double nu, const double x)
    {
      double jnu, nnu, jpnu, npnu;
      besselJN(nu,x,jnu,nnu,jpnu,npnu);
      return jnu;
    }

    inline double besselJp(const double nu, const double x)
    {
      double jnu, nnu, jpnu, npnu;
      besselJN(nu,x,jnu,nnu,jpnu,npnu);
      return jpnu;
    }

    // ------------------------------------------------------------------------
    // Geometry
    // ------------------------------------------------------------------------
    double areaCircleIntersect(double d, double r1, double r2);

    double areaAnnulusIntersect(double d, double r1, double r2, 
				double r3, double r4);
  }
}

#endif // VSAALGEBRA_HPP
