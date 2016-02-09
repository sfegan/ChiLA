//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAMath.cpp
  Various math functions

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#include<cmath>
#include<cfloat>
#include<limits>
#include<algorithm>
#include<exception>

#include<VSAMath.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAMath;

static inline double sgn(double x)
{
  if(x>=0)return 1.0;
  else return -1.0;
}

unsigned VSAMath::realRoots2(const double a, const double b, const double c,
			     double roots[2])
{
  // Reference: Numerical Recipies

  // a - coefficient of x^2
  // b - coefficient of x^1
  // c - coefficient of x^0

  // formula x = ( -b +/- sqrt(b^2 - 4 a c) ) / 2 a

#ifdef VSA_ASSERT
  assert(a!=0);
#endif

  const double R = b*b - 4.0*a*c;
  
  if(R<0)
    {
      roots[0] = roots[1] = FP_NAN;
      return 0;
    }
  else if(R==0)
    {
      roots[0] = roots[1] = -b/2.0/a;
      return 1;
    }
  else
    {
      double q = -0.5 * ( b + sgn(b)*sqrt(R) );
      roots[0] = q/a;
      roots[1] = c/q;
      if(roots[0]>roots[1])std::swap(roots[0],roots[1]);
      return 2;
    }
  
#ifdef VSA_ASSERT
  assert(0); // never reach here
#endif
}

unsigned VSAMath::realRoots3(const double _a, const double _b, 
			     const double _c, const double _d,
			     double roots[3])
{
  // Reference: Numerical Recipies

  // _a - coefficient of x^3
  // _b - coefficient of x^2
  // _c - coefficient of x^1
  // _d - coefficient of x^0

#ifdef VSA_ASSERT
  assert(_a!=0);
#endif

  // NR uses different variables -- translate for clarity
  const double a = _b/_a;
  const double b = _c/_a;
  const double c = _d/_a;

  const double a2 = a*a;
  const double a3 = a2*a;

  const double Q = (a2 - 3.0*b)/9.0;
  const double R = (2.0*a3 - 9.0*a*b + 27.0*c)/54.0;

  const double R2 = R*R;
  const double Q3 = Q*Q*Q;

  double decider = R2-Q3;
  if(fabs(R2+Q3)>0)decider /= fabs(R2+Q3);

  if(decider < -DBL_EPSILON)
    {
      // Three real roots
      const double theta = acos(R/sqrt(Q3));
      const double factor0 = -2.0*sqrt(Q);
      const double factor1 = a/3.0;
      roots[0] = factor0 * cos(theta/3.0) - factor1;
      roots[1] = factor0 * cos((theta+2.0*M_PI)/3.0) - factor1;
      roots[2] = factor0 * cos((theta-2.0*M_PI)/3.0) - factor1;
      if(roots[0]>roots[1])std::swap(roots[0],roots[1]);
      if(roots[1]>roots[2])
	{
	  std::swap(roots[1],roots[2]);
	  if(roots[0]>roots[1])std::swap(roots[0],roots[1]);
	}
      return 3;
    }
  else if(decider < DBL_EPSILON)
    {
      if(R==0)
	{
	  // One real root
	  roots[0] = roots[1] = roots[2] = -a/3.0;
	  return 1;
	}
      else
	{
	  const double A = -sgn(R)*pow(fabs(R),1.0/3.0);
	  // Two real roots
	  roots[0] = 2.0*A-a/3.0;
	  roots[1] = roots[2] = -A-a/3.0;
	  if(roots[0]>roots[2])std::swap(roots[0],roots[2]);
	  return 2;
	}
    }
  else
    {
      const double A = -sgn(R)*pow(fabs(R) + sqrt(R2-Q3),1.0/3.0);
      const double B = (A==0)?0:Q/A;
      roots[0] = (A+B)-a/3.0;

      if(fabs(A-B)<DBL_EPSILON)
	{
	  //.I think this can never happen 

	  // A-B=0 => R^2-Q^3=0  /OR/  R^2-Q^3=|2R|
	  // The first case is manifestly not true in this case (decider>0)
	  // Is the second possible -- have to multiply it out some time

	  roots[1] = roots[2] = -0.5*(A+B)-a/3.0;
	  if(roots[0]>roots[2])std::swap(roots[0],roots[2]);
	  return 2;
	}
      else
	{
	  roots[1] = roots[2] = FP_NAN;
	  return 1;
	}
    }

#ifdef VSA_ASSERT
  assert(0); // never reach here
#endif
}

// ----------------------------------------------------------------------------
// Gamma Function
//
//  C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
//
// ----------------------------------------------------------------------------
double VSAMath::lnGammaLanczos(double z)
{
  const double xm1 = z - double(1);

  static const double lanczos_cheb_7[9] = {
    0.99999999999980993227684700473478L,
    676.520368121885098567009190444019L,
    -1259.13921672240287047156078755283L,
    771.3234287776530788486528258894L,
    -176.61502916214059906584551354L,
    12.507343278686904814458936853L,
    -0.13857109526572011689554707L,
    9.984369578019570859563e-6L,
    1.50563273514931155834e-7L
  };
  
  static const double LOGROOT2PI = 0.9189385332046727417803297364056176L;
  
  double sum = lanczos_cheb_7[0];
  for(unsigned int k = 1; k < 9; ++k)
    sum += lanczos_cheb_7[k] / (xm1 + k);
  
  const double term1 = (xm1 + double(0.5L))
    * std::log((xm1 + double(7.5L))/numeric_constants<double>::euler());
  const double term2 = LOGROOT2PI + std::log(sum);
  const double result = term1 + (term2 - double(7));
  
  return result;
}

// ----------------------------------------------------------------------------
// Chi Squared Probability
// ----------------------------------------------------------------------------
double VSAMath::chi2p(double chi2, int ndf)
{
  if (ndf <= 0) return 0; // Set CL to zero in case ndf<=0

  if (chi2 <= 0) 
    {
      if (chi2 < 0) return 0;
      else          return 1;
    }

   if (ndf==1) 
     {
       double v = 1.-VSAMath::erf(sqrt(chi2)/sqrt(2.));
       return v;
     }

   // Gaussian approximation for large ndf
   double q = sqrt(2*chi2)-sqrt(double(2*ndf-1));
   if (ndf > 30 && q > 5) 
     {
       double v = 0.5*(1-VSAMath::erf(q/sqrt(2.)));
       return v;
     }

   // Evaluate the incomplete gamma function
   return (1-VSAMath::gamma(0.5*ndf,0.5*chi2));
}

double VSAMath::chi2q(double p, int ndf)
{
   // Evaluate the quantiles of the chi-squared probability
   // distribution function.
   // Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
   // implemented by Anna Kreshuk.
   // Incorporates the suggested changes in AS R85 (vol.40(1), pp.233-5, 1991)
   // Parameters:
   //   p   - the probability value, at which the quantile is computed
   //   ndf - number of degrees of freedom

   const double c[]={0, 0.01, 0.222222, 0.32, 0.4, 1.24, 2.2,
		     4.67, 6.66, 6.73, 13.32, 60.0, 70.0,
		     84.0, 105.0, 120.0, 127.0, 140.0, 175.0,
		     210.0, 252.0, 264.0, 294.0, 346.0, 420.0,
		     462.0, 606.0, 672.0, 707.0, 735.0, 889.0,
		     932.0, 966.0, 1141.0, 1182.0, 1278.0, 1740.0,
		     2520.0, 5040.0};
   double e = 5e-7;
   double aa = 0.6931471806;
   const int maxit = 20;
   double ch, p1, p2, q, t, a, b, x;
   double s1, s2, s3, s4, s5, s6;

   if (ndf <= 0) return 0;

   double g = lgamma(0.5*ndf);

   double xx = 0.5 * ndf;
   double cp = xx - 1;
   if (ndf >= std::log(p)*(-c[5]))
     {
       //starting approximation for ndf less than or equal to 0.32
       if (ndf > c[3]) 
	 {
	   x = norm_quantile(p);
	   //starting approximation using Wilson and Hilferty estimate
	   p1=c[2]/ndf;
	   ch = ndf*std::pow((x*sqrt(p1) + 1 - p1), 3);
	   if (ch > c[6]*ndf + 6)
	     ch = -2 * (std::log(1-p) - cp * std::log(0.5 * ch) + g);
	 } 
       else 
	 {
	   ch = c[4];
	   a = std::log(1-p);
	   do{
	     q = ch;
	     p1 = 1 + ch*(c[7]+ch);
	     p2 = ch*(c[9] + ch*(c[8] + ch));
	     t = -0.5 + (c[7] + 2*ch) / p1 - (c[9] + ch*(c[10] + 3*ch))/p2;
	     ch = ch - (1 - std::exp(a + g + 0.5*ch + cp*aa)*p2/p1)/t;
	   }while (fabs(q/ch - 1) > c[1]);
	 }
     } else {
      ch = std::pow((p * xx * std::exp(g + xx * aa)),(1./xx));
      if (ch < e) return ch;
   }
   //call to algorithm AS 239 and calculation of seven term  Taylor series
   for (int i=0; i<maxit; i++)
     {
       q = ch;
       p1 = 0.5 * ch;
       p2 = p - gamma(xx, p1);
       
       t = p2 * std::exp(xx * aa + g + p1 - cp * std::log(ch));
       b = t / ch;
       a = 0.5 * t - b * cp;
       s1 = 
	 (c[19] + a*(c[17] + a*(c[14] + a*(c[13] + a*(c[12] +c[11]*a)))))/c[24];
       s2 = (c[24] + a*(c[29] + a*(c[32] + a*(c[33] + c[35] * a)))) / c[37];
       s3 = (c[19] + a*(c[25] + a * (c[28] + c[31] * a))) / c[37];
       s4 = (c[20] + 
	     a*(c[27] + c[34]*a) + cp*(c[22] + a*(c[30] + c[36] * a)))/c[38];
       s5 = (c[13] + c[21] * a + cp * (c[18] + c[26] * a)) / c[37];
       s6 = (c[15] + cp * (c[23] + c[16] * cp)) / c[38];
       ch = ch + t*(1 + 0.5*t*s1 - 
		    b*cp*(s1 - b*(s2 - b*(s3 - b*(s4 - b*(s5 - b * s6))))));
       if (fabs(q / ch - 1) > e) break;
   }
   return ch;
}

double VSAMath::norm_quantile(double p)
{
  // Computes quantiles for standard normal distribution N(0, 1)
  // at probability p
  // ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-484.

  vsassert((p>0)&&(p<1));


  const double  a0 = 3.3871328727963666080e0;
  const double  a1 = 1.3314166789178437745e+2;
  const double  a2 = 1.9715909503065514427e+3;
  const double  a3 = 1.3731693765509461125e+4;
  const double  a4 = 4.5921953931549871457e+4;
  const double  a5 = 6.7265770927008700853e+4;
  const double  a6 = 3.3430575583588128105e+4;
  const double  a7 = 2.5090809287301226727e+3;
  const double  b1 = 4.2313330701600911252e+1;
  const double  b2 = 6.8718700749205790830e+2;
  const double  b3 = 5.3941960214247511077e+3;
  const double  b4 = 2.1213794301586595867e+4;
  const double  b5 = 3.9307895800092710610e+4;
  const double  b6 = 2.8729085735721942674e+4;
  const double  b7 = 5.2264952788528545610e+3;
  const double  c0 = 1.42343711074968357734e0;
  const double  c1 = 4.63033784615654529590e0;
  const double  c2 = 5.76949722146069140550e0;
  const double  c3 = 3.64784832476320460504e0;
  const double  c4 = 1.27045825245236838258e0;
  const double  c5 = 2.41780725177450611770e-1;
  const double  c6 = 2.27238449892691845833e-2;
  const double  c7 = 7.74545014278341407640e-4;
  const double  d1 = 2.05319162663775882187e0;
  const double  d2 = 1.67638483018380384940e0;
  const double  d3 = 6.89767334985100004550e-1;
  const double  d4 = 1.48103976427480074590e-1;
  const double  d5 = 1.51986665636164571966e-2;
  const double  d6 = 5.47593808499534494600e-4;
  const double  d7 = 1.05075007164441684324e-9;
  const double  e0 = 6.65790464350110377720e0;
  const double  e1 = 5.46378491116411436990e0;
  const double  e2 = 1.78482653991729133580e0;
  const double  e3 = 2.96560571828504891230e-1;
  const double  e4 = 2.65321895265761230930e-2;
  const double  e5 = 1.24266094738807843860e-3;
  const double  e6 = 2.71155556874348757815e-5;
  const double  e7 = 2.01033439929228813265e-7;
  const double  f1 = 5.99832206555887937690e-1;
  const double  f2 = 1.36929880922735805310e-1;
  const double  f3 = 1.48753612908506148525e-2;
  const double  f4 = 7.86869131145613259100e-4;
  const double  f5 = 1.84631831751005468180e-5;
  const double  f6 = 1.42151175831644588870e-7;
  const double  f7 = 2.04426310338993978564e-15;

  double split1 = 0.425;
  double split2=5.;
  double konst1=0.180625;
  double konst2=1.6;

  double q, r, quantile;
  q=p-0.5;
  if (fabs(q)<split1) {
    r=konst1-q*q;
    quantile = q* (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3)
		     * r + a2) * r + a1) * r + a0) /
      (((((((b7 * r + b6) * r + b5) * r + b4) * r + b3)
	 * r + b2) * r + b1) * r + 1.);
  } else {
    if(q<0) r=p;
    else    r=1-p;
    //error case
    if (r<=0)
      quantile=0;
    else {
      r=sqrt(-std::log(r));
      if (r<=split2) {
	r=r-konst2;
	quantile=(((((((c7 * r + c6) * r + c5) * r + c4) * r + c3)
		    * r + c2) * r + c1) * r + c0) /
	  (((((((d7 * r + d6) * r + d5) * r + d4) * r + d3)
	     * r + d2) * r + d1) * r + 1);
      } else{
	r=r-split2;
	quantile=(((((((e7 * r + e6) * r + e5) * r + e4) * r + e3)
		    * r + e2) * r + e1) * r + e0) /
	  (((((((f7 * r + f6) * r + f5) * r + f4) * r + f3)
	     * r + f2) * r + f1) * r + 1);
      }
      if (q<0) quantile=-quantile;
    }
  }
  return quantile;
}

// ----------------------------------------------------------------------------
// Error Function
//
// ----------------------------------------------------------------------------
double VSAMath::erf(double x)
{
  // Computation of the error function erf(x).
  // Erf(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between 0 and x

  //--- NvE 14-nov-1998 UU-SAP Utrecht

  return (1-VSAMath::erfc(x));
}

double VSAMath::erfc(double x)
{
  // Compute the complementary error function erfc(x).
  // Erfc(x) = (2/sqrt(pi)) Integral(exp(-t^2))dt between x and infinity
  //
  //--- Nve 14-nov-1998 UU-SAP Utrecht

  // The parameters of the Chebyshev fit
  const double a1 = -1.26551223,   a2 = 1.00002368,
    a3 =  0.37409196,   a4 = 0.09678418,
    a5 = -0.18628806,   a6 = 0.27886807,
    a7 = -1.13520398,   a8 = 1.48851587,
    a9 = -0.82215223,  a10 = 0.17087277;
  
  double v = 1; // The return value
  double z = std::abs(x);

  if (z <= 0) return v; // erfc(0)=1

  double t = 1/(1+0.5*z);

  v = t*std::exp((-z*z) +
		 a1+t*(a2+t*(a3+t*(a4+t*
				   (a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));

  if (x < 0) v = 2-v; // erfc(-x)=2-erfc(x)

  return v;
}

// ----------------------------------------------------------------------------
// returns  the inverse error function
// x must be  <-1<x<1
// ----------------------------------------------------------------------------
double VSAMath::erf_inverse(double x)
{
  const int kMaxit    = 50;
  const double kEps   = 1e-14;
  const double kConst = 0.8862269254527579;     // sqrt(pi)/2.0

  if(fabs(x) <= kEps) return kConst*x;

  // Newton iterations
  double erfi, derfi, y0,y1,dy0,dy1;
  if(fabs(x) < 1.0) 
    {
      erfi  = kConst*fabs(x);
      y0    = VSAMath::erf(0.9*erfi);
      derfi = 0.1*erfi;
      for (int iter=0; iter<kMaxit; iter++) 
	{
	  y1  = 1. - VSAMath::erfc(erfi);
	  dy1 = fabs(x) - y1;
	  if(fabs(dy1) < kEps)  
	    {if (x < 0) return -erfi; else return erfi;}
	  dy0    = y1 - y0;
	  derfi *= dy1/dy0;
	  y0     = y1;
	  erfi  += derfi;
	  if(fabs(derfi/erfi) < kEps) 
	    {if (x < 0) return -erfi; else return erfi;}
	}
    }
  return 0; //did not converge
}

// ----------------------------------------------------------------------------
// VSAMath::beta( double aa, double bb, double xx )
//
// Incomplete Beta Function
//
// Returns incomplete beta integral of the arguments, evaluated
// from zero to x.  The function is defined as
//
//                  x
//     -            -
//    | (a+b)      | |  a-1     b-1
//  -----------    |   t   (1-t)   dt.
//   -     -     | |
//  | (a) | (b)   -
//                 0
//
// The domain of definition is 0 <= x <= 1.  In this
// implementation a and b are restricted to positive values.
// The integral from x to 1 may be obtained by the symmetry
// relation
//
//    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).
//
// The integral is evaluated by a continued fraction expansion
// or, when b*x is small, by a power series.
//
// ACCURACY:
//
// Tested at uniformly distributed random points (a,b,x) with a and b
// in "domain" and x between 0 and 1.
//                                        Relative error
// arithmetic   domain     # trials      peak         rms
//    IEEE      0,5         10000       6.9e-15     4.5e-16
//    IEEE      0,85       250000       2.2e-13     1.7e-14
//    IEEE      0,1000      30000       5.3e-12     6.3e-13
//    IEEE      0,10000    250000       9.3e-11     7.1e-12
//    IEEE      0,100000    10000       8.7e-10     4.8e-11
// Outputs smaller than the IEEE gradual underflow threshold
// were excluded from these statistics.
//
// Cephes Math Library, Release 2.8:  June, 2000
// Copyright 1984, 1995, 2000 by Stephen L. Moshier
// ----------------------------------------------------------------------------
double VSAMath::beta( double aa, double bb, double xx )
{
  const double eps = std::numeric_limits<double>::epsilon();
  const double MAXSTIR = 108.116855767857671821730036754;
  const double MAXLOG  = 709.782712893383973096206318587;
  const double MINLOG  = -708.396418532264078748994506896;
  
  double a, b, t, x, xc, w, y;
  int flag;

  if( aa <= 0.0 || bb <= 0.0 ) return( 0.0 );

  if( (xx <= 0.0) || ( xx >= 1.0) )
    {
      if( xx ==0 ) return( 0.0 );
      if( xx ==1 ) return( 1.0 );
      return( 0.0 );
    }

  flag = 0;

  /* - to test if that way is better for large b/ (comment out from
       Cephes version)

     if( (bb * xx) <= 1.0 && xx <= 0.95)
     {
     t = pseries(aa, bb, xx);
     goto done;
     }

  **/
  w = 1.0 - xx;

  // --------------------------------------------------------------------------
  // Reverse a and b if x is greater than the mean.  
  // aa,bb > 1 -> sharp rise at x=aa/(aa+bb)
  // flag = 1 designates that a and b have been reversed.
  // --------------------------------------------------------------------------
  if( xx > (aa/(aa+bb)) )
    {
      flag = 1;
      a = bb;
      b = aa;
      xc = xx;
      x = w;
    }
  else
    {
      a = aa;
      b = bb;
      xc = w;
      x = xx;
    }

  if( flag == 1 && (b * x) <= 1.0 && x <= 0.95)
    {
      t = beta_pseries(a, b, x);
      if( t <= eps ) t = 1.0 - eps;
      else t = 1.0 - t;
      return( t );
    }

  // Choose expansion for better convergence. ---------------------------------
  y = x * (a+b-2.0) - (a-1.0);
  if( y < 0.0 ) w = beta_cfe1( a, b, x );
  else w = beta_cfe2( a, b, x ) / xc;

  // --------------------------------------------------------------------------
  // Multiply w by the factor
  // a      b   _             _     _
  // x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   
  // --------------------------------------------------------------------------
  y = a * std::log(x);
  t = b * std::log(xc);
  if( (a+b) < MAXSTIR && std::abs(y) < MAXLOG && std::abs(t) < MAXLOG )
    {
      t = std::pow(xc,b);
      t *= std::pow(x,a);
      t /= a;
      t *= w;
      t *= gamma(a+b) / (gamma(a) * gamma(b));
    }
  else
    {
      // Resort to logarithms. ------------------------------------------------
      y += t + lgamma(a+b) - lgamma(a) - lgamma(b);
      y += std::log(w/a);
      if( y < MINLOG ) t = 0.0;
      else t = std::exp(y);
    }      

  if( flag == 1 )
    {
      if( t <= eps ) t = 1.0 - eps;
      else t = 1.0 - t;
    }
  return( t );
}


// ----------------------------------------------------------------------------
// VSAMath::beta_pseries( double a, double b, double x )
//
// Power series for incomplete beta integral.  Use when b*x is small
// and x not too close to 1.
// ----------------------------------------------------------------------------

double VSAMath::beta_pseries( double a, double b, double x )
{
  const double eps = std::numeric_limits<double>::epsilon();
  const double MAXSTIR = 108.116855767857671821730036754;
  const double MAXLOG  = 709.782712893383973096206318587;
  const double MINLOG  = -708.396418532264078748994506896;
  double s, t, u, v, n, t1, z, ai;

  ai = 1.0 / a;
  u = (1.0 - b) * x;
  v = u / (a + 1.0);
  t1 = v;
  t = u;
  n = 2.0;
  s = 0.0;
  z = eps * ai;
  while( std::abs(v) > z )
    {
      u = (n - b) * x / n;
      t *= u;
      v = t / (a + n);
      s += v; 
      n += 1.0;
    }
  s += t1;
  s += ai;

  u = a * log(x);
  if( (a+b) < MAXSTIR && std::abs(u) < MAXLOG )
    {
      t = gamma(a+b)/(gamma(a)*gamma(b));
      s = s * t * std::pow(x,a);
    }
  else
    {
      t = lgamma(a+b) - lgamma(a) - lgamma(b) + u + std::log(s);
      if( t < MINLOG ) s = 0.0;
      else s = std::exp(t);
    }
  return(s);
}

// ----------------------------------------------------------------------------
// VSAMath::beta_cfe( double a, double b, double x )
//
// Continued fraction expansion for incomplete beta integral.
// ----------------------------------------------------------------------------
double VSAMath::beta_cfe1( double a, double b, double x )
{
  const double eps    = std::numeric_limits<double>::epsilon();
  const double BIG    = 4.503599627370496e15;
  const double BIGINV = 2.22044604925031308085e-16;

  double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
  double k1, k2, k3, k4, k5, k6, k7, k8;
  double r, t, ans, thresh;
  int n;

  k1 = a;
  k2 = a + b;
  k3 = a;
  k4 = a + 1.0;
  k5 = 1.0;
  k6 = b - 1.0;
  k7 = k4;
  k8 = a + 2.0;

  pkm2 = 0.0;
  qkm2 = 1.0;
  pkm1 = 1.0;
  qkm1 = 1.0;
  ans = 1.0;
  r = 1.0;
  n = 0;
  thresh = 3.0 * eps;
  do
    {	
      xk = -( x * k1 * k2 )/( k3 * k4 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      xk = ( x * k5 * k6 )/( k7 * k8 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      if( qk !=0 )
	r = pk/qk;
      if( r != 0 )
	{
	  t = std::abs( (ans - r)/r );
	  ans = r;
	}
      else t = 1.0;

      if( t < thresh ) break;

      k1 += 1.0;
      k2 += 1.0;
      k3 += 2.0;
      k4 += 2.0;
      k5 += 1.0;
      k6 -= 1.0;
      k7 += 2.0;
      k8 += 2.0;

      if( (std::abs(qk) + std::abs(pk)) > BIG )
	{
	  pkm2 *= BIGINV;
	  pkm1 *= BIGINV;
	  qkm2 *= BIGINV;
	  qkm1 *= BIGINV;
	}
      if( (std::abs(qk) < BIGINV) || (std::abs(pk) < BIGINV) )
	{
	  pkm2 *= BIG;
	  pkm1 *= BIG;
	  qkm2 *= BIG;
	  qkm1 *= BIG;
	}
    }
  while( ++n < 300 );

  return(ans);
}

double VSAMath::beta_cfe2( double a, double b, double x )
{
  const double eps    = std::numeric_limits<double>::epsilon();
  const double BIG    = 4.503599627370496e15;
  const double BIGINV = 2.22044604925031308085e-16;

  double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
  double k1, k2, k3, k4, k5, k6, k7, k8;
  double r, t, ans, z, thresh;
  int n;

  k1 = a;
  k2 = b - 1.0;
  k3 = a;
  k4 = a + 1.0;
  k5 = 1.0;
  k6 = a + b;
  k7 = a + 1.0;;
  k8 = a + 2.0;

  pkm2 = 0.0;
  qkm2 = 1.0;
  pkm1 = 1.0;
  qkm1 = 1.0;
  z = x / (1.0-x);
  ans = 1.0;
  r = 1.0;
  n = 0;
  thresh = 3.0 * eps;
  do
    {	
      xk = -( z * k1 * k2 )/( k3 * k4 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      xk = ( z * k5 * k6 )/( k7 * k8 );
      pk = pkm1 +  pkm2 * xk;
      qk = qkm1 +  qkm2 * xk;
      pkm2 = pkm1;
      pkm1 = pk;
      qkm2 = qkm1;
      qkm1 = qk;

      if( qk != 0 )
	r = pk/qk;
      if( r != 0 )
	{
	  t = std::abs( (ans - r)/r );
	  ans = r;
	}
      else t = 1.0;

      if( t < thresh ) break;

      k1 += 1.0;
      k2 -= 1.0;
      k3 += 2.0;
      k4 += 2.0;
      k5 += 1.0;
      k6 += 1.0;
      k7 += 2.0;
      k8 += 2.0;

      if( (std::abs(qk) + std::abs(pk)) > BIG )
	{
	  pkm2 *= BIGINV;
	  pkm1 *= BIGINV;
	  qkm2 *= BIGINV;
	  qkm1 *= BIGINV;
	}
      if( (std::abs(qk) < BIGINV) || (std::abs(pk) < BIGINV) )
	{
	  pkm2 *= BIG;
	  pkm1 *= BIG;
	  qkm2 *= BIG;
	  qkm1 *= BIG;
	}
    }
  while( ++n < 300 );
  return(ans);
}

// ----------------------------------------------------------------------------
// Incomplete Gamma Function
//
//  Abramowitz and Stegun, Handbook of Mathematical Functions
//
// Alternately use power series or continued fraction expansion.
//
// ----------------------------------------------------------------------------
double VSAMath::gamma(double a, double x)
{
  if (a <= 0 || x <= 0) return 0;

  if (x < (a+1)) return gamma_pseries(a,x);
  else           return gamma_cfe(a,x);
}

// ----------------------------------------------------------------------------
// Computation of the incomplete gamma function P(a,x)
// via its continued fraction representation.
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
double VSAMath::gamma_cfe(double a,double x)
{
  const int itmax    = 100;      // Maximum number of iterations
  const double eps   = std::numeric_limits<double>::epsilon();
  const double fpmin = std::numeric_limits<double>::min();

  if (a <= 0 || x <= 0) return 0;

  double gln = lnGammaLanczos(a);
  double b   = x+1-a;
  double c   = 1/fpmin;
  double d   = 1/b;
  double h   = d;
  double an,del;
  for (int i=1; i<=itmax; i++) 
    {
      an = double(-i)*(double(i)-a);
      b += 2;
      d  = an*d+b;
      if (std::abs(d) < fpmin) d = fpmin;
      c = b+an/c;
      if (std::abs(c) < fpmin) c = fpmin;
      d   = 1/d;
      del = d*c;
      h   = h*del;
      if (std::abs(del-1) < eps) break;
      //if (i==itmax) cout << "*GamCf(a,x)* a too large or itmax too small" << endl;
    }
  double v = std::exp(-x+a*std::log(x)-gln)*h;
  return (1-v);
}

// ----------------------------------------------------------------------------
// Computation of the incomplete gamma function P(a,x)
// via its series representation.
//
//--- Nve 14-nov-1998 UU-SAP Utrecht
double VSAMath::gamma_pseries(double a,double x)
{


  const int itmax  = 100;    // Maximum number of iterations
  const double eps = std::numeric_limits<double>::epsilon();

  if (a <= 0 || x <= 0) return 0;

  double gln = lnGammaLanczos(a);
  double ap  = a;
  double sum = 1/a;
  double del = sum;
  for (int n=1; n<=itmax; n++) 
    {
      ap  += 1;
      del  = del*x/ap;
      sum += del;
      if (std::abs(del) < std::abs(sum*eps)) break;
      //if (n==itmax) cout << "*GamSer(a,x)* a too large or itmax too small" << endl;
    }
  double v = sum*std::exp(-x+a*std::log(x)-gln);
  return v;
}

// ----------------------------------------------------------------------------
// Carlson forms of integrals
// ----------------------------------------------------------------------------

double VSAMath::ellipticRFCarlson(double x, double y, double z)
{
  static const double ERRTOL = 0.0025;
  static const double TINY   = // 1.5e-38
    1.1*std::numeric_limits<double>::min()*5.0; 
  static const double BIG    = // 3.0e37
    0.9*std::numeric_limits<double>::max()*0.2;
  static const double THIRD  = 1.0/3.0;
  static const double C1     = 1.0/24.0;
  static const double C2     = 0.1;
  static const double C3     = 3.0/44.0;
  static const double C4     = 1.0/14.0;

  assert((x>=0)&&(y>=0)&&(z>=0)&&(x<=BIG)&&(y<=BIG)&&(z<=BIG)
	 &&(std::min(std::min(x+y,x+z),y+z)>=TINY));

  double xt = x;
  double yt = y;
  double zt = z;

  double delx;
  double dely;
  double delz;
  double ave;

  do {
    const double sqrtx = sqrt(xt);
    const double sqrty = sqrt(yt);
    const double sqrtz = sqrt(zt);
    const double alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    xt = 0.25*(xt+alamb);
    yt = 0.25*(yt+alamb);
    zt = 0.25*(zt+alamb);
    ave = THIRD*(xt+yt+zt);
    delx = (ave-xt)/ave;
    dely = (ave-yt)/ave;
    delz = (ave-zt)/ave;
  } while((fabs(delx)>ERRTOL)||(fabs(dely)>ERRTOL)||(fabs(delz)>ERRTOL));

  const double e2 = delx*dely - delz*delz;
  const double e3 = delx*dely*delz;
  return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}

// Carlson's elliptic integral of the second kind
double VSAMath::ellipticRDCarlson(double x, double y, double z)
{
  static const double ERRTOL = 0.0015;
  static const double TINY   = //1e-25;
    1.1*pow(std::numeric_limits<double>::max(),-2.0/3.0)*2.0;
  static const double BIG    = //4.5e21;
    0.9*pow(std::numeric_limits<double>::min(),-2.0/3.0)*0.1*ERRTOL;
  static const double C1     = 3.0/14.0;
  static const double C2     = 1.0/6.0;
  static const double C3     = 9.0/22.0;
  static const double C4     = 3.0/26.0;
  static const double C5     = 0.25*C3;
  static const double C6     = 1.5*C4;

  assert((x>=0)&&(y>=0)&&(x+y>=TINY)&&(z>=TINY)&&(x<=BIG)&&(y<=BIG)&&(z<=BIG));
  
  double xt = x;
  double yt = y;
  double zt = z;
  double sum = 0.0;
  double fac = 1.0;

  double delx;
  double dely;
  double delz;
  double ave;

  do {
    const double sqrtx = sqrt(xt);
    const double sqrty = sqrt(yt);
    const double sqrtz = sqrt(zt);
    const double alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    sum += fac/(sqrtz*(zt+alamb));
    fac  = 0.25*fac;
    xt   = 0.25*(xt+alamb);
    yt   = 0.25*(yt+alamb);
    zt   = 0.25*(zt+alamb);
    ave  = 0.2*(xt+yt+3.0*zt);
    delx = (ave-xt)/ave;
    dely = (ave-yt)/ave;
    delz = (ave-zt)/ave;
  } while((fabs(delx)>ERRTOL)||(fabs(dely)>ERRTOL)||(fabs(delz)>ERRTOL));
  
  const double ea = delx*dely;
  const double eb = delz*delz;
  const double ec = ea-eb;
  const double ed = ea-6.0*eb;
  const double ee = ed+ec+ec;

  return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		      + delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}

// ----------------------------------------------------------------------------
// Jacobi elliptic functions
// ----------------------------------------------------------------------------
#include<iostream>
void VSAMath::ellipticJacobiFunctions(const double _u, const double k,
				      double &sn, double& cn, double& dn)
{
  static const double CA = 2.7*sqrt(DBL_EPSILON);

  double emc = 1.0-k*k;

  if(emc) 
    {
      bool bo = (emc < 0.0);
      double u = _u;
      double d = 0;
      if(bo)
	{
	  d = 1.0-emc;
	  emc = -d;
	  d = sqrt(d);
	  u *= d;
	}


      dn = 1.0;

      double a = 1.0;
      double c;
      double em[13];
      double en[13];
      int l;
      for(int i=0;i<13;i++)
	{
	  l = i;
	  em[i] = a;
	  emc = sqrt(emc);
	  en[i] = emc;
	  c = 0.5*(a+emc);
	  if(fabs(a-emc) <= CA*a)break;
	  emc *= a;
	  a = c;
	}
      u *= c;

      sn = sin(u);
      cn = cos(u);

      if(sn) 
	{
	  a = cn/sn;
	  c *= a;
	  for(int i=l;i>=0;i--)
	    {
	      double b = em[i];
	      a *= c;
	      c *= dn;
	      dn = (en[i]+a)/(b+a);
	      a = c/b;
	    }
	  a = 1.0/sqrt(c*c+1.0);
	  sn=(sn >= 0.0 ? a : -a);
	  cn=c*sn;
	}
      if (bo)
	{
	  a = dn;
	  dn = cn;
	  cn = a;
	  sn /= d;
	}
    } 
  else 
    {
      cn = 1.0/cosh(_u);
      dn = cn;
      sn = tanh(_u);
    }
}

// ----------------------------------------------------------------------------
// Bessel and Neumann Functions
// ----------------------------------------------------------------------------
void VSAMath::besselJN(const double nu, const double x,
		       double & Jnu, double & Nnu, double & Jpnu, double & Npnu)
{
  if (x == double(0))
    {
      if (nu == double(0))
	{
	  Jnu = double(1);
	  Jpnu = double(0);
	}
      else if (nu == double(1))
	{
	  Jnu = double(0);
	  Jpnu = double(0.5L);
	}
      else
	{
	  Jnu = double(0);
	  Jpnu = double(0);
	}
      Nnu = -std::numeric_limits<double>::infinity();
      Npnu = std::numeric_limits<double>::infinity();
      return;
    }

  const double eps = std::numeric_limits<double>::epsilon();
  //  When the multiplier is N i.e.
  //  fp_min = N * min()
  //  Then J_0 and N_0 tank at x = 8 * N (J_0 = 0 and N_0 = nan)!
  //const double fp_min = double(20) * std::numeric_limits<double>::min();
  const double fp_min = std::sqrt(std::numeric_limits<double>::min());
  const int max_iter = 15000;
  const double x_min = double(2);

  const int nl = (x < x_min
                    ? static_cast<int>(nu + double(0.5L))
                    : std::max(0, static_cast<int>(nu - x + double(1.5L))));

  const double mu = nu - nl;
  const double mu2 = mu * mu;
  const double xi = double(1) / x;
  const double xi2 = double(2) * xi;
  double w = xi2 / numeric_constants<double>::pi();
  int isign = 1;
  double h = nu * xi;
  if (h < fp_min)
    h = fp_min;
  double b = xi2 * nu;
  double d = double(0);
  double c = h;
  int i;
  for (i = 1; i <= max_iter; ++i)
    {
      b += xi2;
      d = b - d;
      if (std::abs(d) < fp_min)
	d = fp_min;
      c = b - double(1) / c;
      if (std::abs(c) < fp_min)
	c = fp_min;
      d = double(1) / d;
      const double del = c * d;
      h *= del;
      if (d < double(0))
	isign = -isign;
      if (std::abs(del - double(1)) < eps)
	break;
    }
  if (i > max_iter)
    std::__throw_runtime_error(__N("Argument x too large in bessel_jn; "
				 "try asymptotic expansion."));
  double Jnul = isign * fp_min;
  double Jpnul = h * Jnul;
  double Jnul1 = Jnul;
  double Jpnu1 = Jpnul;
  double fact = nu * xi;
  for ( int l = nl; l >= 1; --l )
    {
      const double Jnutemp = fact * Jnul + Jpnul;
      fact -= xi;
      Jpnul = fact * Jnutemp - Jnul;
      Jnul = Jnutemp;
    }
  if (Jnul == double(0))
    Jnul = eps;
  double f= Jpnul / Jnul;
  double Nmu, Nnu1, Npmu, Jmu;
  if (x < x_min)
    {
      const double x2 = x / double(2);
      const double pimu = numeric_constants<double>::pi() * mu;
      double fact = (std::abs(pimu) < eps
		       ? double(1) : pimu / std::sin(pimu));
      double d = -std::log(x2);
      double e = mu * d;
      double fact2 = (std::abs(e) < eps
			? double(1) : std::sinh(e) / e);
      double gam1, gam2, gampl, gammi;
      gamma_temme(mu, gam1, gam2, gampl, gammi);
      double ff = (double(2) / numeric_constants<double>::pi())
	* fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
      e = std::exp(e);
      double p = e / (numeric_constants<double>::pi() * gampl);
      double q = double(1) / (e * numeric_constants<double>::pi() * gammi);
      const double pimu2 = pimu / double(2);
      double fact3 = (std::abs(pimu2) < eps
			? double(1) : std::sin(pimu2) / pimu2 );
      double r = numeric_constants<double>::pi() * pimu2 * fact3 * fact3;
      double c = double(1);
      d = -x2 * x2;
      double sum = ff + r * q;
      double sum1 = p;
      for (i = 1; i <= max_iter; ++i)
	{
	  ff = (i * ff + p + q) / (i * i - mu2);
	  c *= d / double(i);
	  p /= double(i) - mu;
	  q /= double(i) + mu;
	  const double del = c * (ff + r * q);
	  sum += del; 
	  const double del1 = c * p - i * del;
	  sum1 += del1;
	  if ( std::abs(del) < eps * (double(1) + std::abs(sum)) )
	    break;
	}
      if ( i > max_iter )
	std::__throw_runtime_error(__N("Bessel y series failed to converge "
				       "in bessel_jn."));
      Nmu = -sum;
      Nnu1 = -sum1 * xi2;
      Npmu = mu * xi * Nmu - Nnu1;
      Jmu = w / (Npmu - f * Nmu);
    }
  else
    {
      double a = double(0.25L) - mu2;
      double q = double(1);
      double p = -xi / double(2);
      double br = double(2) * x;
      double bi = double(2);
      double fact = a * xi / (p * p + q * q);
      double cr = br + q * fact;
      double ci = bi + p * fact;
      double den = br * br + bi * bi;
      double dr = br / den;
      double di = -bi / den;
      double dlr = cr * dr - ci * di;
      double dli = cr * di + ci * dr;
      double temp = p * dlr - q * dli;
      q = p * dli + q * dlr;
      p = temp;
      int i;
      for (i = 2; i <= max_iter; ++i)
	{
	  a += double(2 * (i - 1));
	  bi += double(2);
	  dr = a * dr + br;
	  di = a * di + bi;
	  if (std::abs(dr) + std::abs(di) < fp_min)
	    dr = fp_min;
	  fact = a / (cr * cr + ci * ci);
	  cr = br + cr * fact;
	  ci = bi - ci * fact;
	  if (std::abs(cr) + std::abs(ci) < fp_min)
	    cr = fp_min;
	  den = dr * dr + di * di;
	  dr /= den;
	  di /= -den;
	  dlr = cr * dr - ci * di;
	  dli = cr * di + ci * dr;
	  temp = p * dlr - q * dli;
	  q = p * dli + q * dlr;
	  p = temp;
	  if (std::abs(dlr - double(1)) + std::abs(dli) < eps)
	    break;
	}
      if (i > max_iter)
	std::__throw_runtime_error(__N("Lentz's method failed "
				   "in bessel_jn."));
      const double gam = (p - f) / q;
      Jmu = std::sqrt(w / ((p - f) * gam + q));
      if (Jmu * Jnul < double(0))
	Jmu = -Jmu;
      Nmu = gam * Jmu;
      Npmu = (p + q / gam) * Nmu;
      Nnu1 = mu * xi * Nmu - Npmu;
    }
  fact = Jmu / Jnul;
  Jnu = fact * Jnul1;
  Jpnu = fact * Jpnu1;
  for (i = 1; i <= nl; ++i)
    {
      const double Nnutemp = (mu + i) * xi2 * Nnu1 - Nmu;
      Nmu = Nnu1;
      Nnu1 = Nnutemp;
    }
  Nnu = Nmu;
  Npnu = nu * xi * Nmu - Nnu1;

  return;
}

// ----------------------------------------------------------------------------
// Geometry
// ----------------------------------------------------------------------------
double VSAMath::areaCircleIntersect(double d, double r1, double r2)
{
  if(r1 == 0 || r2 == 0 || d>=r1+r2 ) return 0;
  else if(r2+d <= r1) return M_PI*std::pow(r2,2);
  else if(r1+d <= r2) return M_PI*std::pow(r1,2);
  else return r1*r1*acos((d*d+r1*r1-r2*r2)/(2*d*r1)) +
	 r2*r2*acos((d*d+r2*r2-r1*r1)/(2*d*r2)) -
	     0.5*sqrt((-d+r1+r2)*(d+r1-r2)*(d-r1+r2)*(d+r1+r2));
}

double VSAMath::areaAnnulusIntersect(double d, double r1, double r2, 
				     double r3, double r4)
{
  return areaCircleIntersect(d,r2,r4) + areaCircleIntersect(d,r1,r3) -
    areaCircleIntersect(d,r1,r4) - areaCircleIntersect(d,r2,r3);
}

#ifdef TESTMAIN
#include<iostream>
int main(int argc, char** argv)
{
  for(double x=0;x<1.0;x+=0.01)
    std::cout << x << ' ' 
	      << VSAMath::elliptic1Complete(x) << ' ' 
	      << VSAMath::elliptic2Complete(x) << '\n'; 
}
#endif
