//-*-mode:c++; mode:font-lock;-*-

#ifndef RANDOMNUMBERS_TOS_HPP
#define RANDOMNUMBERS_TOS_HPP

#include <string>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <stdint.h>
#include <vector>
#include <set>

#include <sys/types.h>
#include <fcntl.h>

#ifndef NOTHREADS
#include <pthread.h>
#endif

#define _RANDOMNUMBERS_TOS_NTAB 32

class RandomNumbers_TOS
{
public:
  typedef std::pair< double, double > Pair;

  RandomNumbers_TOS(int32_t seed);
  RandomNumbers_TOS(const std::string seeds_filename);
  RandomNumbers_TOS(const char* seeds_filename);
  virtual ~RandomNumbers_TOS();
  double Uniform();
  double Exponential();
  double Normal();
  double Gamma( int );
  double Gamma(const double a, const double b);
  double GammaByMeanAndStdDev(const double mean, const double stddev);
  int Poisson( double );
  int Binomial( double , int );
  double InverseCDF(const std::vector< Pair > &inv_cdf);
  void GenerateInverseCDF(std::vector< Pair > &cdf, unsigned nbins = 0);

  static std::string defaultFilename();

private:

#ifndef NOTHREADS
  static void initializeThreadOnce(void);
  static pthread_once_t sTDOnceControl;
  static pthread_mutex_t sTDMutex;
  static std::set<std::string> sTDLocks;
#endif

  RandomNumbers_TOS( const RandomNumbers_TOS & );
  RandomNumbers_TOS operator=( const RandomNumbers_TOS & );

  unsigned int systemRand();

  void ran2_lock();
  void ran2_unlock();

  void ran2_read();
  void ran2_write();
  void ran2_init(int32_t*);

  double ran2();

  double gammln( double );

  double gsl_rng_uniform_pos();
  double gsl_ran_gamma_int(const unsigned int a);
  double gsl_gamma_large(const double a);
  double gsl_gamma_frac(const double a);

  double interpolate(double x,
		     std::vector<Pair>::const_iterator itr); 

  int32_t ran2_idum1;
  int32_t ran2_idum2;
  int32_t ran2_iy;
  int32_t ran2_iv[_RANDOMNUMBERS_TOS_NTAB];

  std::string    fFileName;
  std::string    fLockFileName;
  int            fLockFileFD;
  struct flock   fLock;
};

#endif // ifndef RANDOMNUMBERS_TOS_HPP
