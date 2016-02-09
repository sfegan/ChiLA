//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAAtmosphere.cpp
  Integrated atmospheric thickness profile

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/20/2006
*/

#include<fstream>
#include<cmath>
#include<cerrno>
#include<cstring>
#include<algorithm>

#include<VSLineTokenizer.hpp>

#include<VSAAtmosphere.hpp>

using namespace VERITAS;

VSAAtmosphere::
VSAAtmosphere(const std::string& filename): fGood(false), fAtmo(), fAtmoLn()
{
  VSLineTokenizer tokenizer;
  std::ifstream stream(filename.c_str());
  if(!stream)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": Could not open file \"" 
		<< filename << '"' << std::endl
		<< __PRETTY_FUNCTION__ << ": " << strerror(errno) << std::endl;
      return;
    }

  while(!stream.eof())
    {
      VSTokenList tokens;
      tokenizer.tokenize(stream, tokens);
      if(tokens.size() != 4)continue;

      AtmoSlice slice;
      tokens[0].to(slice.fHeight);
      tokens[1].to(slice.fRho);
      tokens[2].to(slice.fThickness);
      tokens[3].to(slice.fNMinusOne);
      fAtmo.push_back(slice);

      slice.fRho        = log(slice.fRho);
      slice.fThickness  = log(slice.fThickness);
      slice.fNMinusOne  = log(slice.fNMinusOne);
      fAtmoLn.push_back(slice);
    }
  std::sort(fAtmo.begin(), fAtmo.end());
  std::sort(fAtmoLn.begin(), fAtmoLn.end());
  integrate();
  fGood=true;
}

void VSAAtmosphere::integrate()
{
  static const double c = 2.9979246e-04; // c in km/ns

  double delay=0;
  unsigned nslice=fAtmo.size();
  for(unsigned islice=1;islice<nslice;islice++)
    {
      // v = c/n  --  delay is integrated: dT = dH / v - dH / c
      //                                      = dH/c (c/v - 1)
      //                                      = dH/c (n - 1)

      double dH = fAtmo[islice].fHeight - fAtmo[islice-1].fHeight;
      double n_1 = (fAtmo[islice].fNMinusOne+fAtmo[islice-1].fNMinusOne)/2.0;
      delay += dH/c*n_1;
      fAtmo[islice].fDelay = delay;
      fAtmoLn[islice].fDelay = delay;
    }
}

void VSAAtmosphere::
interpolate(const double& height, double& rho, double& thickness, 
	    double& n_minus_one, double& delay) const
{
  unsigned ilo = 0;
  unsigned ihi = fAtmoLn.size()-1;
  if(fAtmoLn.empty())
    { 
      rho         = -1;
      thickness   = -1;
      n_minus_one = -1;
      delay       = -1;
      return;
    }
  else if(height <= fAtmoLn[ilo].fHeight)
    {
      rho         = fAtmo[ilo].fRho;
      thickness   = fAtmo[ilo].fThickness;
      n_minus_one = fAtmo[ilo].fNMinusOne;
      delay       = fAtmo[ilo].fDelay;
      return;
    }
  else if(height >= fAtmoLn[ihi].fHeight)
    {
      rho         = fAtmo[ihi].fRho;
      thickness   = fAtmo[ihi].fThickness;
      n_minus_one = fAtmo[ihi].fNMinusOne;
      delay       = fAtmo[ihi].fDelay;
      return;
    }

  while(ihi-ilo > 1)
    {
      unsigned imid = (ilo+ihi)/2;
#if 0
      std::cerr << "  " << ilo << ' ' << ihi << ' ' << imid << std::endl;
#endif
      if(height>=fAtmoLn[imid].fHeight)ilo=imid;
      else ihi=imid;
    }

  const double hlo = fAtmoLn[ilo].fHeight;
  const double hhi = fAtmoLn[ihi].fHeight;
  const double x = (height-hlo)/(hhi-hlo);

  const double rlo = fAtmoLn[ilo].fRho;
  const double rhi = fAtmoLn[ihi].fRho;
  const double tlo = fAtmoLn[ilo].fThickness;
  const double thi = fAtmoLn[ihi].fThickness;
  const double nlo = fAtmoLn[ilo].fNMinusOne;
  const double nhi = fAtmoLn[ihi].fNMinusOne;
  const double dlo = fAtmoLn[ilo].fDelay;
  const double dhi = fAtmoLn[ihi].fDelay;

#if 0
  std::cerr << height << ' '
	    << ilo << ' ' << ihi << ' '
	    << hlo << ' ' << hhi << ' ' << x << ' '
	    << tlo << ' ' << thi 
	    << std::endl;
#endif

  rho         = exp(rlo + (rhi-rlo)*x);
  thickness   = exp(tlo + (thi-tlo)*x);
  n_minus_one = exp(nlo + (nhi-nlo)*x);
  delay       = dlo + (dhi-dlo)*x;
}

double VSAAtmosphere::rho(const double& height) const
{
  double r, t, n, d;
  interpolate(height,r,t,n,d);
  return r;
}

double VSAAtmosphere::thickness(const double& height) const
{
  double r, t, n, d;
  interpolate(height,r,t,n,d);
  return t;
}

double VSAAtmosphere::nMinusOne(const double& height) const
{
  double r, t, n, d;
  interpolate(height,r,t,n,d);
  return n;
}

double VSAAtmosphere::delay(const double& height) const
{
  double r, t, n, d;
  interpolate(height,r,t,n,d);
  return d;
}

VSAAtmosphere VSAAtmosphere::usStandard()
{
  static double us_std[][4] = {
    { 0.000, 0.12219E-02, 0.10350E+04, 0.28232E-03 },
    { 1.000, 0.11099E-02, 0.91853E+03, 0.25634E-03 },
    { 2.000, 0.10054E-02, 0.81286E+03, 0.23214E-03 },
    { 3.000, 0.90839E-03, 0.71725E+03, 0.20975E-03 },
    { 4.000, 0.81888E-03, 0.63097E+03, 0.18904E-03 },
    { 5.000, 0.73643E-03, 0.55328E+03, 0.16994E-03 },
    { 6.000, 0.66012E-03, 0.48352E+03, 0.15235E-03 },
    { 7.000, 0.59048E-03, 0.42105E+03, 0.13620E-03 },
    { 8.000, 0.52609E-03, 0.36529E+03, 0.12136E-03 },
    { 9.000, 0.46741E-03, 0.31567E+03, 0.10782E-03 },
    { 10.000, 0.41370E-03, 0.27167E+03, 0.95426E-04 },
    { 11.000, 0.36499E-03, 0.23278E+03, 0.84194E-04 },
    { 12.000, 0.31209E-03, 0.19900E+03, 0.71987E-04 },
    { 13.000, 0.26674E-03, 0.17012E+03, 0.61523E-04 },
    { 14.000, 0.22792E-03, 0.14543E+03, 0.52581E-04 },
    { 15.000, 0.19479E-03, 0.12434E+03, 0.44937E-04 },
    { 16.000, 0.16651E-03, 0.10631E+03, 0.38406E-04 },
    { 17.000, 0.14236E-03, 0.90902E+02, 0.32840E-04 },
    { 18.000, 0.12168E-03, 0.77727E+02, 0.28071E-04 },
    { 19.000, 0.10403E-03, 0.66465E+02, 0.23997E-04 },
    { 20.000, 0.88928E-04, 0.56837E+02, 0.20516E-04 },
    { 21.000, 0.75750E-04, 0.48620E+02, 0.17475E-04 },
    { 22.000, 0.64544E-04, 0.41621E+02, 0.14887E-04 },
    { 23.000, 0.55021E-04, 0.35655E+02, 0.12695E-04 },
    { 24.000, 0.46965E-04, 0.30566E+02, 0.10833E-04 },
    { 25.000, 0.40097E-04, 0.26222E+02, 0.92494E-05 },
    { 27.500, 0.27126E-04, 0.17925E+02, 0.62570E-05 },
    { 30.000, 0.18420E-04, 0.12302E+02, 0.42495E-05 },
    { 32.500, 0.12139E-04, 0.85361E+01, 0.28004E-05 },
    { 35.000, 0.84696E-05, 0.59874E+01, 0.19537E-05 },
    { 37.500, 0.59542E-05, 0.42029E+01, 0.13738E-05 },
    { 40.000, 0.39967E-05, 0.29752E+01, 0.92196E-06 },
    { 42.500, 0.27910E-05, 0.21358E+01, 0.64379E-06 },
    { 45.000, 0.19671E-05, 0.15470E+01, 0.45379E-06 },
    { 47.500, 0.14044E-05, 0.11295E+01, 0.32390E-06 },
    { 50.000, 0.10273E-05, 0.82800E+00, 0.23698E-06 },
    { 55.000, 0.56800E-06, 0.44045E+00, 0.13104E-06 },
    { 60.000, 0.30906E-06, 0.22771E+00, 0.71295E-07 },
    { 65.000, 0.16285E-06, 0.11361E+00, 0.37569E-07 },
    { 70.000, 0.82868E-07, 0.54414E-01, 0.19114E-07 },
    { 75.000, 0.40145E-07, 0.24940E-01, 0.92604E-08 },
    { 80.000, 0.18430E-07, 0.10993E-01, 0.42513E-08 },
    { 85.000, 0.82291E-08, 0.46676E-02, 0.18985E-08 },
    { 90.000, 0.34321E-08, 0.19250E-02, 0.79163E-09 },
    { 95.000, 0.14063E-08, 0.78968E-03, 0.32437E-09 },
    { 100.000, 0.57185E-09, 0.32602E-03, 0.13189E-09 },
    { 105.000, 0.24206E-09, 0.13421E-03, 0.55841E-10 },
    { 110.000, 0.10312E-09, 0.52792E-04, 0.23788E-10 },
    { 115.000, 0.46595E-10, 0.17216E-04, 0.10748E-10 },
    { 120.000, 0.24596E-10, 0.00000E+00, 0.56734E-11 } };

  VSAAtmosphere atmo;
  for(unsigned islice=0;islice<sizeof(us_std)/sizeof(*us_std);islice++)
    {
      AtmoSlice slice;
      slice.fHeight     = us_std[islice][0];
      slice.fRho        = us_std[islice][1];
      slice.fThickness  = us_std[islice][2];
      slice.fNMinusOne  = us_std[islice][3];
      atmo.fAtmo.push_back(slice);

      slice.fRho        = log(slice.fRho);
      slice.fThickness  = log(slice.fThickness);
      slice.fNMinusOne  = log(slice.fNMinusOne);
      atmo.fAtmoLn.push_back(slice);
    }
  atmo.integrate();
  atmo.fGood=true;
  return atmo;
}

#ifdef TEST_MAIN

// g++ -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -I. -I../VSUtility -DTEST_MAIN -o test VSAAtmosphere.cpp ../VSUtility/libVSUtility.a

int main(int argc, char** argv)
{
  std::string program(*argv);
  argc--,argv++;

  VSAAtmosphere atmo;

  if(!argc)
    {
      atmo = VSAAtmosphere::usStandard();
    }
  else
    {
      atmo = VSAAtmosphere(*argv);
    }

  for(double h=0; h<150; h+=0.1)
    {
      double r, t, n, d;
      atmo.interpolate(h,r,t,n,d);
      std::cout << h << ' ' << r << ' ' << t << ' ' << n << ' ' << d
		<< std::endl;
    }
}

#endif
