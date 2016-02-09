//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAFunction.cpp

  Various predefined functions.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2008
*/

#include <algorithm>

#include <VSAFunction.hpp>
#include <VSAMath.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAMath;
using namespace VERITAS::VSAFunction;

BicubicInterpolation::BicubicInterpolation():
  m_x1lo(), m_x1hi(), m_nx1(), m_dx1(),
  m_x2lo(), m_x2hi(), m_nx2(), m_dx2(),
  m_y(), m_y1(), m_y2(), m_y12()
{

}	

BicubicInterpolation::
BicubicInterpolation(double x1lo, double x1hi, unsigned nx1,
		     double x2lo, double x2hi, unsigned nx2,
		     const std::vector< std::vector<double> >& y,
		     const std::vector< std::vector<double> >& y1,
		     const std::vector< std::vector<double> >& y2,
		     const std::vector< std::vector<double> >& y12):
  m_x1lo(x1lo), m_x1hi(x1hi), m_nx1(nx1), m_dx1(),
  m_x2lo(x2lo), m_x2hi(x2hi), m_nx2(nx2), m_dx2(),
  m_y(y), m_y1(y1), m_y2(y2), m_y12(y12)
{
  m_dx1 = (x1hi-x1lo)/(double)nx1;
  m_dx2 = (x2hi-x2lo)/(double)nx2;
}			   

double BicubicInterpolation::val(double x1, double x2) const
{
  int ix1 = std::max(0,(int)floor((x1-m_x1lo)/m_dx1));
  int ix2 = std::max(0,(int)floor((x2-m_x2lo)/m_dx2));

  ix1 = std::min(ix1,(int)(m_nx1-2));
  ix2 = std::min(ix2,(int)(m_nx2-2));

  double x1c = m_x1lo + ix1*m_dx1;
  double x2c = m_x2lo + ix2*m_dx2;

  std::vector< std::vector<double> > c;

  bcucof(ix1,ix2,c);

  double t = (x1-x1c)/m_dx1;
  double u = (x2-x2c)/m_dx2;

  double v = 0;
  
  for(int i = 3; i >= 0; i--)
    v=t*v+((c[i][3]*u + c[i][2])*u + c[i][1])*u + c[i][0];

  return v;
}

void BicubicInterpolation::bcucof(unsigned ix, unsigned iy,
				  std::vector< std::vector<double> >& c) const
{
  const unsigned np = 16;
  std::vector<double> x(np);
  std::vector<double> cl(np);

  static const int wt[16][16] =
    {
      { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
      {-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0 },
      { 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
      { 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1 },
      { 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1 },
      {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0 },
      { 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2 },
      {-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2 },
      { 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0 },
      {-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1 },
      { 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1 }
    };

  static int wt2[16][16];

  for(unsigned i = 0; i < 16; i++)
    for(unsigned j = 0; j < 16; j++)
      wt2[i][j] = wt[j][i];

  double dx1dx2 = m_dx1*m_dx2;

  x[0]  = m_y[ix][iy];
  x[1]  = m_y[ix+1][iy];
  x[2]  = m_y[ix+1][iy+1];
  x[3]  = m_y[ix][iy+1];
  x[4]  = m_y1[ix][iy]*m_dx1;
  x[5]  = m_y1[ix+1][iy]*m_dx1;
  x[6]  = m_y1[ix+1][iy+1]*m_dx1;
  x[7]  = m_y1[ix][iy+1]*m_dx1;
  x[8]  = m_y2[ix][iy]*m_dx2;
  x[9]  = m_y2[ix+1][iy]*m_dx2;
  x[10] = m_y2[ix+1][iy+1]*m_dx2;
  x[11] = m_y2[ix][iy+1]*m_dx2;
  x[12] = m_y12[ix][iy]*dx1dx2;
  x[13] = m_y12[ix+1][iy]*dx1dx2;
  x[14] = m_y12[ix+1][iy+1]*dx1dx2;
  x[15] = m_y12[ix][iy+1]*dx1dx2;

  for(unsigned i = 0; i < np; i++)
    {
      double xx = 0;
      for(unsigned k = 0; k < np; k++)
	xx += (double)wt[i][k]*x[k];

      cl[i]=xx;
    }

  c = std::vector< std::vector<double> >(4,std::vector<double>(4));

  unsigned l = 0;
  for(unsigned i = 0; i < 4; i++)
    for(unsigned j = 0; j < 4; j++)
      c[i][j]=cl[l++];
}



// ----------------------------------------------------------------------------
// Spline
// ----------------------------------------------------------------------------
Spline::Spline(BoundaryCondition bc): 
  m_xy(), m_y2(), m_yp1(), m_ypn(), m_bc(bc)
{
  
}

Spline::Spline(const std::vector< std::pair<double,double> >& xy,
	       BoundaryCondition bc):
  m_xy(xy), m_y2(), m_yp1(), m_ypn(), m_bc(bc)
{
  spline();
}

void Spline::spline()
{
  const unsigned n = m_xy.size();
  m_y2.resize(n);
  if(n <= 1) return;
  double p;
  std::vector< double > u(n);

  m_y2[0] = 0;
  u[0] = 0;

  if(m_bc == BC_CLAMPED)
    {
      m_y2[0] = -0.5;
      u[0] = (3.0/(m_xy[1].first-m_xy[0].first))*
	((m_xy[1].second-m_xy[0].second)/(m_xy[1].first-m_xy[0].first)-m_yp1);
    }

  for(unsigned i = 1; i < n-1; i++)
    {
      double sig = 
	(m_xy[i].first - m_xy[i-1].first)/(m_xy[i+1].first - m_xy[i-1].first);
      p = sig*m_y2[i-1]+2.0;
      m_y2[i] = (sig-1.0)/p;
      u[i]=
	(m_xy[i+1].second-m_xy[i].second)/(m_xy[i+1].first - m_xy[i].first) -
	(m_xy[i].second-m_xy[i-1].second)/(m_xy[i].first - m_xy[i-1].first);
      u[i] = (6.0*u[i]/(m_xy[i+1].first - m_xy[i-1].first)-sig*u[i-1])/p;
    }

  double qn = 0;
  double un = 0;

  if(m_bc == BC_CLAMPED)
    {
      qn = 0.5;
      un = (3.0/(m_xy[n-1].first - m_xy[n-2].first))*
	(m_ypn-(m_xy[n-1].second-m_xy[n-2].second)/
	 (m_xy[n-1].first-m_xy[n-2].first));
    }

  std::vector< double >::reverse_iterator kitr = m_y2.rbegin();
  std::vector< double >::reverse_iterator uitr = u.rbegin() + 1;

  m_y2[n-1] = (un-qn*u[n-2])/(qn*m_y2[n-2]+1.0);
  for(kitr = m_y2.rbegin()+1;kitr != m_y2.rend(); ++kitr)
    {
      *kitr = (*kitr)*(*(kitr-1)) + *uitr;
      ++uitr;
    }  
}

double Spline::val(const VSACoord::Coord1D& x) const
{
  return val(x.x());
}

double Spline::val(const double& x) const
{
  const unsigned n = m_xy.size();
  if(n == 0) return 0;
  else if(n == 1) return m_xy[0].second;
  
  vsassert(m_y2.size() == n);

  int k;
  int klo = 0;
  int khi = n-1;

  while(khi-klo > 1)
    {
      k=(khi+klo)>>1;
      if(m_xy[k].first > x) khi=k;
      else klo = k;
    }

  double h = m_xy[khi].first - m_xy[klo].first;
  double a = (m_xy[khi].first-x)/h;
  double b = (x-m_xy[klo].first)/h;
  double y = a*m_xy[klo].second+b*m_xy[khi].second+
    ((a*a*a-a)*m_y2[klo] + (b*b*b-b)*m_y2[khi])*(h*h)/6.0;
  return y;
}

void Spline::setPoint(double x, double y)
{
  m_xy.push_back(std::make_pair(x,y));
}

// ----------------------------------------------------------------------------
// Poly - 1D
// ----------------------------------------------------------------------------
void Poly::operator() (const double& x, VSAAlgebra::VecND& v) const
{
  v.resize(m_n);
  double s = v[0] = 1.0;
  for(unsigned i=1;i<m_n;i++)v[i] = s *= x;
}

double Poly::val(const double& x) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return param()*v;
}

double Poly::val(const double& x, const VSAAlgebra::VecND& a) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return a*v;
}

void Poly::dyda(const double& x, VSAAlgebra::VecND& dyda) const
{
  Poly::dyda(x,param(),dyda);
}

void Poly::dyda(const double& x, const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const
{
  (*this)(x,dyda);
}

void Poly::dydx(const double& x, VSAAlgebra::VecND& dydx) const
{
  dydx.resize(m_n);
  dydx[0] = 0;
  for(unsigned id=1;id<m_n;id++)
    dydx[id] = double(id)*std::pow(x,(double)(id-1));
}

// ----------------------------------------------------------------------------
// Gauss1D
// ----------------------------------------------------------------------------
// double Gauss1D::val(double x) const
// {
//   return m_a*exp(-std::pow((x-m_mu)/m_sigma,2));
// }

// double Gauss1D::val(double x, const VSAAlgebra::VecND& a) const
// {
//   return a(0)*exp(-std::pow((x-a(1))/a(2),2));
// }

// void Gauss1D::dyda(double x, VSAAlgebra::VecND& dyda) const
// {
//   dyda(x,param(),dyda);
// }

// void Gauss1D::dyda(double x, const VSAAlgebra::VecND& a,
// 		   VSAAlgebra::VecND& dyda) const
// {
//   dyda(0) = exp(-std::pow((x-a(1))/a(2),2));
//   dyda(1) = 2*a(0)/std::pow(a(2),2)*(x-a(1))*dyda(0);
//   dyda(2) = 2*a(0)/std::pow(a(2),3)*std::pow(x-a(1),2)*dyda(0);
// }

// ----------------------------------------------------------------------------
// FourierBesselSeries
// ----------------------------------------------------------------------------
FourierBesselSeries::
FourierBesselSeries(unsigned m, unsigned n, double R):
  ParamFn<VSACoord::Coord2D>(), m_mn(), m_R(R), m_norm()
{
  std::vector<std::pair<int,int> > mn;
    
  for(int im = 0; im < (int)m; im++)
    for(int in = 0; in < (int)n; in++)
      {
	mn.push_back(std::make_pair(im,in));
	if(im != 0) mn.push_back(std::make_pair(-im,in));
      }

  set(mn,R);
}

FourierBesselSeries::
FourierBesselSeries(const std::vector<std::pair<int,int> >& mn, double R):
  ParamFn<VSACoord::Coord2D>(mn.size()),
  m_mn(), m_R(R), m_norm()
{
  set(mn,R);
}

void FourierBesselSeries::set(const std::vector<std::pair<int,int> >& mn, 
			      double R)
{
  const double R2 = std::pow(R,2);
  m_mn = mn;
  const unsigned norder = mn.size();
  setNParam(norder);
  m_norm.resize(norder);

  for(unsigned i = 0; i < norder; i++)
    {
      unsigned m = std::abs(m_mn[i].first);
      unsigned n = std::abs(m_mn[i].second);
      double J2 = 
	std::pow(besselJ(m+1,numeric_constants<double>::bessel_zero(m,n)),2);

      m_norm(i) = R2*J2/2.;

      if(m == 0) m_norm(i) *= 2*M_PI;
      else m_norm(i) *= M_PI;
    }
}

void FourierBesselSeries::operator() (const VSACoord::Coord2D& x, 
				      VSAAlgebra::VecND& v) const
{
  v.resize(m_mn.size());
  for(unsigned i=0;i<m_mn.size();i++)
    {
      int m = std::abs(m_mn[i].first);
      int n = m_mn[i].second;

      v[i] = besselJ(m,x.r()*numeric_constants<double>::bessel_zero(m,n)/m_R);
      if(m_mn[i].first > 0) v[i] *= cos(m*x.phi());
      else if(m_mn[i].first < 0) v[i] *= sin(m*x.phi());

      v[i] /= besselJ(m+1,numeric_constants<double>::bessel_zero(m,n));
    }
}

double FourierBesselSeries::val(const VSACoord::Coord2D& x, 
				const VSAAlgebra::VecND& a) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return a * v;
}

double FourierBesselSeries::val(const VSACoord::Coord2D& x) const
{
  VSAAlgebra::VecND v;
  (*this)(x,v);
  return param()*v;
}

void FourierBesselSeries::dyda(const VSACoord::Coord2D& x, 
			       VSAAlgebra::VecND& dyda) const
{
  FourierBesselSeries::dyda(x,param(),dyda);
}

void FourierBesselSeries::dyda(const VSACoord::Coord2D& x, 
			       const VSAAlgebra::VecND& a, 
			       VSAAlgebra::VecND& dyda) const
{
  (*this)(x,dyda);
}

// ----------------------------------------------------------------------------
// FourierBesselSeries2
// ----------------------------------------------------------------------------
FourierBesselSeries2::FourierBesselSeries2(double R):
  ParamFn<VSACoord::Coord2D>(), m_mn(), m_R(R), m_norm(), 
  m_c(), m_dc()
{
  m_norm = 1./(M_PI*std::pow(m_R,2));
}

FourierBesselSeries2::FourierBesselSeries2(unsigned m, unsigned n, double R):
  ParamFn<VSACoord::Coord2D>(), m_mn(), m_R(R), m_norm(),
  m_c(), m_dc()
{
  m_norm = 1./(M_PI*std::pow(m_R,2));
  set(m,n);
}

FourierBesselSeries2::
FourierBesselSeries2(const std::vector<std::pair<unsigned,unsigned> >& mn, 
		     double R):
  ParamFn<VSACoord::Coord2D>(), m_mn(), m_R(R), m_norm(),
  m_c(), m_dc()
{
  m_norm = 1./(M_PI*std::pow(m_R,2));
  set(mn);
}

void FourierBesselSeries2::dyda(const VSACoord::Coord2D& x, 
				VSAAlgebra::VecND& dyda) const
{
  dyda.resize(nparm());
  dyda.set(0.0);

  const unsigned norder = m_mn.size();
  std::vector< double > J(norder);

  getJ(x,J);

  const unsigned nparm = FourierBesselSeries2::nparm();
  double s = m_c;
  for(unsigned i=0;i<nparm;i++) s += param(i)*J[i];
  for(unsigned i=0;i<nparm;i++) dyda[i] = 2*s*(J[i]+m_dc[i]);

  dyda *= m_norm;
}

void FourierBesselSeries2::dyda(const VSACoord::Coord2D& x, 
				const VSAAlgebra::VecND& a, 
				VSAAlgebra::VecND& dyda) const
{
  vsassert(a.size() == nparm());
  dyda.resize(nparm());
  dyda.set(0.0);

  const unsigned norder = m_mn.size();
  std::vector< double > J(norder);

  getJ(x,J);

  double c;
  VSAAlgebra::VecND dc;

  calc(a,c,dc);

  const unsigned nparm = FourierBesselSeries2::nparm();
  double s = c;
  for(unsigned i=0;i<nparm;i++) s += a[i]*J[i];
  for(unsigned i=0;i<nparm;i++) dyda[i] = 2*s*(J[i]+dc[i]);

  dyda *= m_norm;
}

double FourierBesselSeries2::val(const VSACoord::Coord2D& x, 
				 const VSAAlgebra::VecND& a) const
{
  vsassert(a.size() == nparm());
  const unsigned norder = m_mn.size();
  std::vector< double > J(norder);

  getJ(x,J);

  double c;
  VSAAlgebra::VecND dc;

  calc(a,c,dc);

  const unsigned nparm = FourierBesselSeries2::nparm();
  double sum = c;
  for(unsigned i=0;i<nparm;i++) sum += a[i]*J[i];

  return sum*sum*m_norm;
}

double FourierBesselSeries2::val(const VSACoord::Coord2D& x) const
{  
  const unsigned norder = m_mn.size();
  std::vector< double > J(norder);

  getJ(x,J);

  const unsigned nparm = FourierBesselSeries2::nparm();
  double sum = m_c;
  for(unsigned i=0;i<nparm;i++) sum += param(i)*J[i];
  return sum*sum*m_norm;
}

void FourierBesselSeries2::setNParam(unsigned n)
{ 
  ParamFn<VSACoord::Coord2D>::setNParam(n);
  calc(param(),m_c,m_dc);
}

void FourierBesselSeries2::set(unsigned m, unsigned n)
{
  m_mn.clear();
  for(int im = 0; im < (int)m; im++)
    for(int in = 0; in < (int)n; in++)
      {
	unsigned ip = m_mn.size();

	if(im == 0) m_mn.push_back(Param(ip,im,in));
	else
	  {
	    m_mn.push_back(Param(ip,im,in));
	    m_mn.push_back(Param(ip+1,-im,in));
	  }
      }

  setNParam(m_mn.size());
  computeSpline();
}

void FourierBesselSeries2::
set(const std::vector<std::pair<unsigned,unsigned> >& mn)
{
  m_mn.clear();
  for(unsigned i = 0; i < mn.size(); i++)
    {
      int m = (int)mn[i].first;
      int n = (int)mn[i].second;
      unsigned ip = m_mn.size();

      if(m == 0) m_mn.push_back(Param(ip,m,n));
      else
	{
	  m_mn.push_back(Param(ip,m,n));
	  m_mn.push_back(Param(ip+1,-m,n));
	}
    }

  setNParam(m_mn.size());
  computeSpline();
}

void FourierBesselSeries2::computeSpline()
{
  m_bessel_spline.clear();

  for(std::vector< Param >::iterator itr = m_mn.begin();
      itr != m_mn.end(); ++itr)
    {
      int m = std::abs(itr->m);
      int n = itr->n;

      double bzero = numeric_constants<double>::bessel_zero(m,n);
      
//       const unsigned nbinr = 200;
//       const double dr = m_R/(nbinr-1)/2.;
//       const unsigned nbinphi = 200;
//       const double dphi = 2*M_PI/(nbinphi-1)/2.;

//       VSNSpace hist(0-dr,m_R+dr,nbinr,-M_PI-dphi,M_PI+dphi,nbinphi);

//       VSNSpace::Cell c(2);
//       for(c.i[0] = 0; c.i[0] < hist.axis(0).nbin; c.i[0]++)
// 	{
// 	  double r = hist.axis(0).midCoordUnchecked(c.i[0]);
// 	  double J = besselJ(m,r*bzero/m_R)/besselJ(m+1,bzero);

// 	  for(c.i[1] = 0; c.i[1] < hist.axis(1).nbin; c.i[1]++)
// 	    {
// 	      if(m == 0) hist[c] = J;
// 	      else if(itr->m > 0)
// 		{
// 		  double phi = hist.axis(1).midCoordUnchecked(c.i[1]);
// 		  hist[c] = 2*J*cos(m*phi);
// 		}
// 	      else
// 		{
// 		  double phi = hist.axis(1).midCoordUnchecked(c.i[1]);
// 		  hist[c] = 2*J*sin(m*phi);
// 		}
// 	    }
// 	}


      const unsigned npoints = 101;
      std::vector< std::pair< double, double > > xy;
      for(unsigned i = 0; i < npoints; i++)
	{
	  double r = (double)i*m_R/(double)(npoints-1);
	  double J = besselJ(m,r*bzero/m_R)/besselJ(m+1,bzero);
	  xy.push_back(std::make_pair(r,J));
	}       

      Spline spline(xy);      
      m_bessel_spline.push_back(spline);
      //      m_bessel_hist.push_back(hist);
    }
}

void FourierBesselSeries2::getJ(const VSACoord::Coord2D& x,
				std::vector<double>& J) const
{
  const unsigned norder = m_mn.size();
  J.resize(norder);
  for(unsigned i=0;i<norder;i++) 
    {     
      int m = std::abs(m_mn[i].m);
      J[i] = m_bessel_spline[i].val(x.r());
      if(m_mn[i].m > 0) J[i] *= 2*cos(m*x.phi());
      else if(m_mn[i].m < 0) J[i] *= 2*sin(m*x.phi());

//       VSNSpace::Point p(2);
//       p.x[0] = x.r();
//       if(x.phi() < -M_PI) p.x[1] = fmod(x.phi(),2*M_PI) + 2*M_PI;
//       else if(x.phi() > M_PI) p.x[1] = fmod(x.phi(),2*M_PI);
//       else p.x[1] = x.phi();

      //      double J2;
      //      m_bessel_hist[i].interpolateWeight(p,J[i]);
    }
}

void FourierBesselSeries2::calc(const VSAAlgebra::VecND& a, double& c,
				VSAAlgebra::VecND& dc) const
{
  const unsigned nparm = FourierBesselSeries2::nparm();
  double s1 = 0;
  double s2 = 0;
  vsassert(m_mn.size() == nparm);
  for(unsigned i=0;i<nparm;i++)
    {
      double bz = numeric_constants<double>::bessel_zero(m_mn[i].m,m_mn[i].n);

      if(m_mn[i].m == 0)
	{
	  s1 += a(i)/bz;
	  s2 += std::pow(a(i),2);
	}
      else
	s2 += std::pow(a(i),2);
    }
  
  double s3 = sqrt(4*s1*s1 + 1 - s2);
  c = -2*s1 + s3;

  dc.resize(nparm);
  for(unsigned i=0;i<nparm;i++)
    {
      double bz = numeric_constants<double>::bessel_zero(m_mn[i].m,m_mn[i].n);
      if(m_mn[i].m == 0) dc[i] = -2./bz + (4*s1/bz-a[i])/s3;
      else dc[i] = -a[i]/s3;
    }
}

