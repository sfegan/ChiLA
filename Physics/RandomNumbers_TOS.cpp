//-*-mode:c++; mode:font-lock;-*-

// The old supplier (TOS) has been obsoleted by the new generator (TNG)
// classes - SJF 2007-11-03

/**
  * :FILENAME:     RandomNumbers_TOS.cpp
  *
  * :PURPOSE:      C++ implementation of the random number 
  *                generators as a class based on a long period 
  *                (>2e+18) uniform random number generator of 
  *                L'Ecuyer with Bays-Durham shuffle. 
  *
  * :REFERENCES:   Numerical Recipes
  *
  * :AUTHOR:       Name            Vladimir V. Vassiliev
  *                Institution     U of U
  *                E-mail          vvv@physics.utah.edu
  *
  * :DATE:         Saturday, March 11, 2001
  *
  * :NOTES:        Many modifications for cluster file locking SJF 2005/2006
  *
  */

/*! \example example/random.cpp
    This is an example of how to use this class
 */

#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include "RandomNumbers_TOS.hpp"

static const double PI=M_PI;

static const int32_t RAN2_A1=40014;
static const int32_t RAN2_A2=40692;
static const int32_t RAN2_Q1=53668;
static const int32_t RAN2_Q2=52774;
static const int32_t RAN2_R1=12211;
static const int32_t RAN2_R2=3791;
static const int32_t RAN2_M1=2147483563;
static const int32_t RAN2_M2=2147483399;

static const double RNMX=(1.0-FLT_EPSILON);

#ifndef NOTHREADS
pthread_once_t RandomNumbers_TOS::sTDOnceControl = PTHREAD_ONCE_INIT;
pthread_mutex_t RandomNumbers_TOS::sTDMutex = PTHREAD_MUTEX_INITIALIZER;
std::set<std::string> RandomNumbers_TOS::sTDLocks;
#endif

//-------------------------------------------------------------------
// constructor //////////////////////////////////////////////////////
//-------------------------------------------------------------------

RandomNumbers_TOS::RandomNumbers_TOS(int32_t seed):
  ran2_idum1(), ran2_idum2(), ran2_iy(), ran2_iv(),
  fFileName(), fLockFileName(), fLockFileFD(-1), fLock()
{
  if(seed<0)ran2_init(&seed);
  else if(seed>0){ seed=-seed; ran2_init(&seed); }
  else
    {
      seed = -systemRand()%1000000;
      ran2_init(&seed);
    }
}

RandomNumbers_TOS::RandomNumbers_TOS(const char* seeds_filename):
  ran2_idum1(), ran2_idum2(), ran2_iy(), ran2_iv(),
  fFileName(seeds_filename), fLockFileName(), fLockFileFD(-1), fLock()
{
  assert(seeds_filename != 0);
  ran2_lock();                            // lock the lockfile
  ran2_read();                            // initialize seeds
  for(unsigned i=0;i<100;i++)ran2();      // generate some random numbers to
  ran2_write();                           // prevent repetitive crashes
}

RandomNumbers_TOS::RandomNumbers_TOS(const std::string seeds_filename):
  ran2_idum1(), ran2_idum2(), ran2_iy(), ran2_iv(),
  fFileName(seeds_filename), fLockFileName(), fLockFileFD(-1), fLock()
{
  ran2_lock();                            // lock the lockfile
  ran2_read();                            // initialize seeds
  for(unsigned i=0;i<100;i++)ran2();      // generate some random numbers to
  ran2_write();                           // prevent repetitive crashes
}

//-------------------------------------------------------------------
// destructor ///////////////////////////////////////////////////////
//-------------------------------------------------------------------

RandomNumbers_TOS::~RandomNumbers_TOS()
{
  if(!fFileName.empty())
    {
      ran2_write();                           // save seeds
      ran2_unlock();                          // unlock the lockfile
    }
}

//-------------------------------------------------------------------
// accessor functions ///////////////////////////////////////////////
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// public methods ///////////////////////////////////////////////////
//-------------------------------------------------------------------
//-------------------- RandomNumbers_TOS::Uniform ----------------------- 
double RandomNumbers_TOS::Uniform(void)
{
  return ran2();
}
//-------------------- RandomNumbers_TOS::Exponential ------------------- 
double RandomNumbers_TOS::Exponential(void)
/** Returns an exponentially distributed, positive, random deviate 
  * of unit mean, using ran2() as the source of uniform deviates.
  * Waiting times between independent Poisson-random events is
  * exponentially distributed, for example.
  */
{
  double dum;

  do
    dum=ran2();
  while (dum == 0.0);
  return -log(dum);
}
//-------------------- RandomNumbers_TOS::Normal ------------------------ 
double RandomNumbers_TOS::Normal(void)
/** Returns a normally distributed deviate with zero mean and unit variance, 
  * using ran2() as the source of uniform deviates. Algorithm is based on
  * the Box-Muller transformation to get two normal deviates.
  */
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if(iset == 0) {
    do {
      v1=2.0*ran2()-1.0; 
      v2=2.0*ran2()-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } 
  else {
    iset=0;
    return gset;
  }
}
//-------------------- RandomNumbers_TOS::Gamma ------------------------- 
double RandomNumbers_TOS::Gamma(int ia)
/** Returns a deviate distributed as a gamma distribution of integer order ia, 
  * i.e., a waiting time to the iath event in a Poisson process of unit mean, 
  * using ran2() as the source of uniform deviates.
  * pdf=x^(ia-1)*exp(-x)/Gamma(ia) 
  */
{
  int j;
  double am,e,s,v1,v2,x,y;
  if (ia < 1) {
    fprintf(stderr," Error in Gamma: integer order %d is not > 0\n",ia); 
    exit (1);
  }
  if(ia < 6) {
    x=1.0;
    for (j=1;j<=ia;j++) x *= ran2();
    x = -log(x);
  } 
  else {
    do {
      do {
        do { 
          v1=ran2();
          v2=2.0*ran2()-1.0;
        } while (v1*v1+v2*v2 > 1.0);
        y=v2/v1;
        am=ia-1;
        s=sqrt(2.0*am+1.0);
        x=s*y+am; 
      } while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (ran2() > e);
  }
  return x;
} 

//------------------ RandomNumbers_TOS::IncompleteGamma ------------------- 

// Asapted from GSL v1.6 -- original copyright notice below

/* randist/gamma.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* The Gamma distribution of order a>0 is defined by:

   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a1 and a2 with the same scale parameter b, then
   X+Y has gamma distribution of order a1+a2.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

/* The Gamma distribution of order a>0 is defined by:

   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a1 and a2 with the same scale parameter b, then
   X+Y has gamma distribution of order a1+a2.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

double
RandomNumbers_TOS::Gamma(const double a, const double b)
{
  /* assume a > 0 */
  unsigned int na = static_cast<unsigned int>(floor (a));

  if (a == na)
    {
      return b * gsl_ran_gamma_int(na);
    }
  else if (na == 0)
    {
      return b * gsl_gamma_frac(a);
    }
  else
    {
      return b * (gsl_ran_gamma_int(na) + gsl_gamma_frac(a - na)) ;
    }
}

double
RandomNumbers_TOS::GammaByMeanAndStdDev(const double mean, const double stddev)
{
  double a = mean*mean/stddev/stddev;
  double b = stddev*stddev/mean;
  return Gamma(a,b);
}

double
RandomNumbers_TOS::gsl_ran_gamma_int (const unsigned int a)
{
  if (a < 12)
    {
      unsigned int i;
      double prod = 1;

      for (i = 0; i < a; i++)
        {
          prod *= gsl_rng_uniform_pos();
        }

      /* Note: for 12 iterations we are safe against underflow, since
         the smallest positive random number is O(2^-32). This means
         the smallest possible product is 2^(-12*32) = 10^-116 which
         is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gsl_gamma_large((double) a);
    }
}

double
RandomNumbers_TOS::gsl_gamma_large(const double a)
{
  /* Works only if a > 1, and is most efficient if a is large

     This algorithm, reported in Knuth, is attributed to Ahrens.  A
     faster one, we are told, can be found in: J. H. Ahrens and
     U. Dieter, Computing 12 (1974) 223-246.  */

  double sqa, x, y, v;
  sqa = sqrt (2 * a - 1);
  do
    {
      do
        {
          y = tan (M_PI * Uniform());
          x = sqa * y + a - 1;
        }
      while (x <= 0);
      v = Uniform();
    }
  while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

  return x;
}

double
RandomNumbers_TOS::gsl_gamma_frac (const double a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;
  p = M_E / (a + M_E);
  do
    {
      u = Uniform();
      v = gsl_rng_uniform_pos();

      if (u < p)
        {
          x = exp ((1 / a) * log (v));
          q = exp (-x);
        }
      else
        {
          x = 1 - log (v);
          q = exp ((a - 1) * log (x));
        }
    }
  while (Uniform() >= q);

  return x;
}

#if 0
double
RandomNumbers_TOS::gsl_ran_gamma_pdf (const double a, const double b)
{
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (a == 1)
        return 1/b ;
      else
        return 0 ;
    }
  else if (a == 1)
    {
      return exp(-x/b)/b ;
    }
  else 
    {
      double p;
      double lngamma = gammaln(a);
      p = exp ((a - 1) * log (x/b) - x/b - lngamma)/b;
      return p;
    }
}
#endif

double
RandomNumbers_TOS::gsl_rng_uniform_pos()
{
  double x = Uniform();
  while(x <= 0)x=Uniform();
  return x;
}

//-------------------- RandomNumbers_TOS::Poisson ----------------------- 
int RandomNumbers_TOS::Poisson(double xm)
/** Returns random deviate drawn from a Poisson distribution of mean xm, 
  * using ran2() as a source of uniform random deviates. 
  */
{
  static double sq,alxm,g,oldm=(-1.0);
  double em,t,y;

  if (xm < 12.0) {           /* Use direct method */
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    t=1.0;
    do {
      ++em;
      t *= ran2();
    } while (t > g);
  } 
  else {                     /* Use rejection method */
    if (xm != oldm) { 
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
        y=tan(PI*ran2());
        em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (ran2() > t);
  }
  return (int)em;
}
//-------------------- RandomNumbers_TOS::Binomial ---------------------- 
int RandomNumbers_TOS::Binomial(double pp, int n)
/** Returns an integer value that is a random deviate drawn from a 
  * binomial distribution of n trials each of probability pp, using
  * ran2() as a source of uniform random deviates. 
  */
{
  int j;
  static int nold=(-1);
  double am,em,g,angle,p,bnl,sq,t,y;
  static double pold=(-1.0),pc,plog,pclog,en,oldg;

  p=(pp <= 0.5 ? pp : 1.0-pp);
  am=n*p;
  if (n < 25) {              /* Use direct method */
    bnl=0.0;
    for (j=1;j<=n;j++) if (ran2() < p) ++bnl;
  } 
  else if (am < 1.0) {
    g=exp(-am);
    t=1.0;
    for (j=0;j<=n;j++) {
      t *= ran2();
      if (t < g) break;
    }
    bnl=(j <= n ? j : n);
  } 
  else {                     /* Use rejection method */
    if (n != nold) { 
      en=n;
      oldg=gammln(en+1.0);
      nold=n;
    } 
    if (p != pold) { 
      pc=1.0-p;
      plog=log(p);
      pclog=log(pc);
      pold=p;
    }
    sq=sqrt(2.0*am*pc);
    do {
      do {
        angle=PI*ran2();
        y=tan(angle);
        em=sq*y+am;
      } while (em < 0.0 || em >= (en+1.0));
      em=floor(em);
      t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
        -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
    } while (ran2() > t);
    bnl=em;
  }
  if (p != pp) bnl=n-bnl;
  return (int)bnl;
}

//-------------------- RandomNumbers_TOS::InverseCDF ----------------------------- 
double 
RandomNumbers_TOS::InverseCDF(const std::vector< RandomNumbers_TOS::Pair > &inv_cdf)
/**  Returns a continuous random deviate drawn from the inverse CDF
  *  provided in the argument.  Note that this function assumes that
  *  the inverse CDF has equidistant steps in cumulative probability.
  */
{
  double x = ran2();
  unsigned i = (unsigned)ceil(x*(double)(inv_cdf.size()-1));
  return interpolate(x,inv_cdf.begin()+i);
}

//-------------------- RandomNumbers_TOS::GenerateInverseCDF --------------------- 
void 
RandomNumbers_TOS::GenerateInverseCDF(std::vector< RandomNumbers_TOS::Pair > &cdf,
				  unsigned nbins)
/**  Overwrites the CDF provided in the argument with its inverse.
  *  The second argument specifies the number of equidistant bins in
  *  cumulative probability that will be generated.
  */
{
  if(nbins == 0)
    nbins = cdf.size();

  for(std::vector< RandomNumbers_TOS::Pair >::iterator itr = cdf.begin();
      itr != cdf.end(); ++itr)
    *itr = std::make_pair(itr->second,itr->first);

  std::vector< Pair > inv_cdf;
  std::vector< Pair >::const_iterator itr = cdf.begin()+1;

  for(unsigned i = 0; i < nbins; i++) 
    {
      double x;

      if(i==0)
	x = 0.0;
      else if(i == nbins-1)
	x = 1.0;
      else
	x = (double)i/(double)(nbins-1);

      while(x > itr->first && itr+1 != cdf.end())
	itr++;

      inv_cdf.push_back(Pair(x,interpolate(x,itr)));
    }

  cdf = inv_cdf;
}

/* 

   PMTpe(int type=0) generates a PMT like response to a "single"
   photon being converted to electrons at the photo-cathode. This can
   result in single photo-electron (p.e.) production, 2 p.e., 3
   p.e. or afterpulsing which can result in large signals being
   produced.

   On average a single photo-electron will result, i.e. this function,
   when run a sufficiently large number of times, will result in a
   distribution of detected signals whose mean is one. The variance of
   the distribution and its "afterpulsing tail" is dependent on the
   type of tube being simulated.

   At present there are five tubes modelled in this function. They are
   selected by choosing a value for the "type" parameter,

   0 - Hamamatsu R960  - Serial #fa0385 - Inner tube on Granite III camera
   1 - Hamamatsu R960  - Serial #fa0473 - 
   2 - Hamamatsu R960  - Serial #fa0717 - 

   3 - Hamamatsu R1398 - Serial #lb5385 - Outer tube on Granite III camera
   4 - Hamamatsu R1398 - Serial #lb5449 - also in 331 pixel Whipple camera

   Check the individual file (e.g. fa0385.c) for details on the
   distributions.

   These files make up elements of a spline used to digitise the
   inverted, normalised integral rate curves. On the x-axis (the first
   column) is the probability (in log base-e space) of getting a
   certain charge (in p.e.s). On the y-axis (2nd column) is the charge
   (in log p.e.s). Also in the 3rd column is the differential of the
   curve, ie d(ln p.e.)/d(ln prob) at each point. This is used by the
   spline functions.

   So we choose a number at random between 0 and 1. Then we take the
   log and map it into log photo-electrons using the spline, which is
   very fast. Then we take the anti-log and hey-presto we know what
   charge resulted.

   Stephen Fegan - 2001-07-16

*/

  
//-------------------------------------------------------------------
// private methods //////////////////////////////////////////////////
//-------------------------------------------------------------------
//-------------------- RandomNumbers_TOS::ran2 -------------------------- 
double RandomNumbers_TOS::ran2(void)
/**  long period (>2e+18) random number generator of L'Ecuyer with  
  *  Bays-Durham shuffle. Returns a uniform random deviate between 0.0 and 1.0 
  *  (exclusive of the endpoint values). To initialize call ran2_init(&idum)
  *  with any idum (sign of idum makes no difference, idum equal zero is 
  *  the same as idum equal 1). Call to ran2_write() will write state
  *  of generator into RAN2_FILE, and call to ran2_read() will load it
  *  back to continue pseudo-random sequence. RNMX approximates the largest 
  *  floating value that is less than 1 which is defined in <float.h>. 
  */
{
  int    j;
  int32_t   k;
  double am=(1.0/RAN2_M1);
  int32_t   imm1=(RAN2_M1-1);
  int32_t   ndiv=(1+(RAN2_M1-1)/_RANDOMNUMBERS_TOS_NTAB);
  double tmp;

  k=ran2_idum1/RAN2_Q1;
  ran2_idum1=RAN2_A1*(ran2_idum1-k*RAN2_Q1)-RAN2_R1*k; 
  if (ran2_idum1<0) ran2_idum1+=RAN2_M1;
  k=ran2_idum2/RAN2_Q2;
  ran2_idum2=RAN2_A2*(ran2_idum2-k*RAN2_Q2)-RAN2_R2*k;
  if (ran2_idum2<0) ran2_idum2+=RAN2_M2;
  j=ran2_iy/ndiv;
  ran2_iy=ran2_iv[j]-ran2_idum2;
  ran2_iv[j]=ran2_idum1;
  if(ran2_iy<1) ran2_iy+=imm1;
  tmp=(double)(am*ran2_iy);
  if(tmp>RNMX) return RNMX;
  else return tmp;
}

//-------------------- RandomNumbers_TOS::ran2_init --------------------- 
void RandomNumbers_TOS::ran2_init(int32_t *idum)
{
  int32_t k;

  if(*idum <0) *idum=-(*idum);
  if(*idum==0) *idum=1;
  ran2_idum2=*idum;
  for (int j=_RANDOMNUMBERS_TOS_NTAB+7;j>=0;j--) {
    k=(*idum)/RAN2_Q1;
    *idum=RAN2_A1*(*idum-k*RAN2_Q1)-RAN2_R1*k; 
    if (*idum<0) *idum+=RAN2_M1;
    if (j<_RANDOMNUMBERS_TOS_NTAB) ran2_iv[j]= *idum;
  }
  ran2_iy=ran2_iv[0];
  ran2_idum1=*idum;
  return;
}
//-------------------- RandomNumbers_TOS::ran2_read --------------------- 

void RandomNumbers_TOS::ran2_read(void)
{
  unsigned ngood = 0;

#if 1
  std::ifstream fstream(fFileName.c_str());
  if(fstream)
    {
      if(fstream >> ran2_idum1)ngood++;
      if(fstream >> ran2_idum2)ngood++;
      if(fstream >> ran2_iy)ngood++;
      for(unsigned j=0;j<_RANDOMNUMBERS_TOS_NTAB;j++)
	if(fstream >> ran2_iv[j])ngood++;
    }
#else
  int j;
  FILE *fp;

  fp = fopen(fFileName.c_str(),"r");
  if(fp)
    {
      ngood += fscanf(fp,"%d\n",&ran2_idum1);
      ngood += fscanf(fp,"%d\n",&ran2_idum2);
      ngood += fscanf(fp,"%d\n",&ran2_iy);
      for (j=0;j<_RANDOMNUMBERS_TOS_NTAB;j++)
	ngood += fscanf(fp,"%d\n",&ran2_iv[j]);
      fclose(fp);
    }
#endif

  if(ngood != _RANDOMNUMBERS_TOS_NTAB+3)
    {
      std::cerr 
	<< "RandomNumbers_TOS::ran2_read(): could not read RNG state from file"
	<< std::endl
	<< "RandomNumbers_TOS::ran2_read(): " << fFileName << std::endl
	<< "RandomNumbers_TOS::ran2_read(): initializing state from default state"
	<< std::endl;

      int32_t dummy2=-systemRand()%1000000; // initialize seeds
      ran2_init(&dummy2);
      
      unsigned burn = systemRand()%1000000;
      if(burn<500000)burn = systemRand()%1000000;
      while(burn--)ran2();
    }

  return;
}
//-------------------- RandomNumbers_TOS::ran2_write -------------------- 
void RandomNumbers_TOS::ran2_write(void)
{
#if 1
  std::ofstream fstream(fFileName.c_str());
  if(!fstream)
    {
      std::cerr << "ran2_write: could not open " << fFileName << " for writing"
		<< std::endl;
      perror("ran2_write");
      exit(EXIT_FAILURE);
    }
  else
    {
      fstream << ran2_idum1 << std::endl
	      << ran2_idum2 << std::endl
	      << ran2_iy << std::endl;
      for(unsigned j=0;j<_RANDOMNUMBERS_TOS_NTAB;j++)
	fstream << ran2_iv[j] << std::endl;
    }
#else
  int j;
  FILE *fp;

  fp = fopen(fFileName.c_str(),"w");
  if(fp==0)
    {
      std::cerr << "ran2_write: could not open " << fFileName << " for writing"
		<< std::endl;
      perror("ran2_write");
      exit(EXIT_FAILURE);
    }
  else
    {
      fprintf(fp,"%d\n",ran2_idum1);
      fprintf(fp,"%d\n",ran2_idum2);
      fprintf(fp,"%d\n",ran2_iy);
      for (j=0;j<_RANDOMNUMBERS_TOS_NTAB;j++)fprintf(fp,"%d\n",ran2_iv[j]);
      fclose(fp);
    }
#endif
  return;
}

#ifndef NOTHREADS
void RandomNumbers_TOS::initializeThreadOnce(void)
{
  pthread_mutex_init(&sTDMutex,NULL);
}
#endif

void RandomNumbers_TOS::ran2_lock()
{
#ifndef NOTHREADS
  pthread_once(&sTDOnceControl, initializeThreadOnce);
#endif

  std::string file = fFileName;
  std::string::size_type ifind;

  std::string path;
  ifind = file.length();
  for(std::string::size_type ichar = 0; ichar<file.length(); ichar++)
    if(file[ichar]=='/')ifind=ichar;
  if(ifind!=file.length())
    {
      path=file.substr(0,ifind+1);
      file=file.substr(ifind+1);
    }
  
  std::string ext;
  ifind = file.length();
  for(std::string::size_type ichar = 0; ichar<file.length(); ichar++)
    if(file[ichar]=='.')ifind=ichar;
  if(ifind!=file.length())
    {
      ext=file.substr(ifind);
      file=file.substr(0,ifind);
    }

  unsigned instance = 0;
  bool continue_looping = true;

  while(continue_looping)
    {
      char buffer[20];
      sprintf(buffer,"_%u",instance);

      fLockFileName = 
	path+std::string(".lock_")+file+std::string(buffer)+ext;

#ifndef NOTHREADS
      pthread_mutex_lock(&sTDMutex);
      if(sTDLocks.find(fLockFileName) != sTDLocks.end())
	{
	  pthread_mutex_unlock(&sTDMutex);
	  instance++;
	  continue;
	}
      sTDLocks.insert(fLockFileName);
      pthread_mutex_unlock(&sTDMutex);
#endif

      fLockFileFD = open(fLockFileName.c_str(),
			 O_WRONLY|O_CREAT,S_IRUSR|S_IWUSR);

      if(fLockFileFD<0)
	{
	  std::cerr << "ran2_lock: could not open lockfile: " << fLockFileName
		    << std::endl;
	  perror("ran2_lock");
	  exit(EXIT_FAILURE);      
	}

      fLock.l_type   = F_WRLCK;
      fLock.l_whence = SEEK_SET;
      fLock.l_start  = 0;
      fLock.l_len    = 0;
      fLock.l_pid    = 0;
      
      if(fcntl(fLockFileFD, F_SETLK, &fLock)<0)
	{
	  if((errno==EACCES)||(errno==EAGAIN))
	    {
#ifndef NOTHREADS
	      pthread_mutex_lock(&sTDMutex);
	      sTDLocks.erase(fLockFileName);
	      pthread_mutex_unlock(&sTDMutex);
#endif
	      close(fLockFileFD);
	      instance++;
	    }
	  else
	    {
	      std::cerr << "ran2_lock: could not lock lockfile: " 
			<< fLockFileName << std::endl;
	      perror("ran2_lock");
	      exit(EXIT_FAILURE);      
	    }
	}
      else
	{

	  fFileName = path+file+std::string(buffer)+ext;
	  continue_looping=false;
	}
    }
}

void RandomNumbers_TOS::ran2_unlock()
{
  fLock.l_type   = F_UNLCK;
  fLock.l_whence = SEEK_SET;
  fLock.l_start  = 0;
  fLock.l_len    = 0;
  fLock.l_pid    = 0;
  fcntl(fLockFileFD, F_SETLKW, &fLock);

  close(fLockFileFD);
  fLockFileFD=-1;

#ifndef NOTHREADS
  pthread_mutex_lock(&sTDMutex);
  sTDLocks.erase(fLockFileName);
  pthread_mutex_unlock(&sTDMutex);
#endif
}

unsigned int RandomNumbers_TOS::systemRand()
{
  unsigned int sys_rand;

  struct timeval init_tv;
  gettimeofday(&init_tv,0);
  srandom(init_tv.tv_usec^init_tv.tv_sec);
  sys_rand = random();
  
  struct stat stat_buf;
  if(stat("/dev/random", &stat_buf)>=0)
    {
      int fd = open("/dev/random",O_RDONLY);
      read(fd,&sys_rand,sizeof(sys_rand));
      close(fd);
    }

  return sys_rand;
}

std::string RandomNumbers_TOS::defaultFilename()
{
  std::ostringstream stream;
  const char* tmpdir = getenv("TMPDIR");
  if(tmpdir)
    {
      stream << tmpdir;
      if((*tmpdir!='\0')&&(tmpdir[strlen(tmpdir)-1]!='/'))stream << '/';
    }
  else
    {
      stream << "/tmp/";
    } 
  stream << "rng_state_uid" << getuid() << ".dat";
  return stream.str();
}

//-------------------- RandomNumbers_TOS::gammln ------------------------ 
double RandomNumbers_TOS::gammln(double xx)
/* Returns the value ln[Gamma(xx)] for xx > 0. */
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146, -86.50532032941677,
                        24.01409824083091, -1.231739572450155,
                    0.1208650973866179e-2, -0.5395239384953e-5};
  int j;

  if (xx <= 0) {
    fprintf(stderr," Error in gammln: argument %e is not > 0\n",xx); 
    exit (1);
  }

  y=x=xx;
  tmp=x+5.5;
  tmp-=(x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser+=cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

//-------------------- RandomNumbers_TOS::interpolate ---------------------------- 
double RandomNumbers_TOS::interpolate(double x,
				  std::vector<Pair>::const_iterator itr)
/**  Performs a linear interpolation at the coordinate x using two
  *  points along a vector of pairs.
  */ 
{
  return (itr-1)->second + 
    (x-(itr-1)->first)*((itr)->second-(itr-1)->second)/
    ((itr)->first-(itr-1)->first);
}

#ifdef TEST_MAIN_GAMMA
#include<iostream>
#include<sstream>
int main(int argc, char** argv)
{
  double mean=0.5;
  double stddev=0.2;
  unsigned n=10;
  argc--,argv++;
  if(argc)
    {
      std::istringstream str(*argv);
      str >> mean;
      argc--,argv++;
    }
  if(argc)
    {
      std::istringstream str(*argv);
      str >> stddev;
      argc--,argv++;
    }
  if(argc)
    {
      std::istringstream str(*argv);
      str >> n;
      argc--,argv++;
    }
  RandomNumbers_TOS rng("random.seeds");
  for(unsigned i=0;i<n;i++)
    std::cout << rng.GammaByMeanAndStdDev(mean,stddev) << std::endl;
}
#endif

#ifdef TEST_MAIN_LOCK
#include<iostream>
#include<sstream>
int main(int argc, char** argv)
{
  std::ostringstream stream;
  stream << "/tmp/rng_state_uid" << getuid() << ".dat";
  std::string rngstatefile(stream.str());

  RandomNumbers_TOS rng1(rngstatefile.c_str());
  //  sleep(100);
}
#endif
