//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSourceModel.cpp

  Class for generating a CR acceptance model.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.12 $
  \date       07/29/2007

  $Id: VSSourceModel.cpp,v 3.12 2010/10/25 00:42:03 matthew Exp $

*/

#include <fstream>

#include <VSSourceModel.hpp>
#include <VSAMath.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSAQuadrature.hpp>
#include <VSFileUtility.hpp>
#include <VSLineTokenizer.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAFunction;


// ----------------------------------------------------------------------------
// VSSourceModel
// ----------------------------------------------------------------------------
VSSourceModel::VSSourceModel(double bin_size_deg): 
  m_bin_size(bin_size_deg), 
  m_domega(std::pow(bin_size_deg,2)), m_rmax(0.5),
  m_obs_xy(), m_effarea(), m_psf()
{

}

double VSSourceModel::val(const VSAAlgebra::Vec2D& xy) const
{
  const unsigned nobs = m_obs_xy.size();

  double v = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSModelCoord c(xy,iobs);
      v += val(c);
    }
  return v;
}

void VSSourceModel::integrate(const VSAAlgebra::VecND& a,
			      const VSAAlgebra::MatrixND& cov,
			      double& integral, 
			      double& integral_err) const
{
  VSAAlgebra::Vec2D xy(a(1),a(2));
  const double dx = m_bin_size*(unsigned)(m_rmax/m_bin_size);

  const double xlo = xy.x() - dx;
  const double xhi = xy.x() + dx;
  const double ylo = xy.y() - dx;
  const double yhi = xy.y() + dx;

  integral = 0;
  VSAAlgebra::VecND b(a.ndim());

  const unsigned nobs = m_obs_xy.size();
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      for(double x = xlo; x < xhi; x += m_bin_size)
	for(double y = ylo; y < yhi; y += m_bin_size)
	  {
	    VSModelCoord c(x,y,iobs);
	    integral += val(c,a);

	    VSAAlgebra::VecND _dyda;
	    dyda(c,a,_dyda);
	    b += _dyda;	  
	  }
    }

  VSAAlgebra::VecND b2 = cov*b;

  integral_err = sqrt(b*b2);
}

double VSSourceModel::integrate(double R1, double R2) const
{
  const unsigned nobs = m_obs_xy.size();
  VSAAlgebra::Vec2D xy(getXY());

  double z = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSACoord::Polar<VSModelCoord> rphi1(xy,R1,0);
      VSACoord::Polar<VSModelCoord> rphi2(xy,R2,0);

      rphi1.setObs(iobs);
      rphi2.setObs(iobs);

      typedef VSAFunction::ParamMemberFn<VSSourceModel, VSModelCoord > Fn;

      Fn fn(this,&VSSourceModel::val,param());

      try
	{
	  z += 2*M_PI*VSAMath::Quadrature::integrate1D(fn,rphi1,rphi2,0);
	}
      catch(const std::exception& e)
	{
	  std::cerr 
	    << std::string(__PRETTY_FUNCTION__) 
	    << ": Caught exception." << std::endl
	    << "X = " << xy.x() << " Y = " << xy.y() << std::endl
	    << "R1 = " << R1 << " R2 = " << R2 << " IOBS = " << iobs
	    << std::endl
	    << e.what() << std::endl;
	}
    }

  return z/m_domega;
}

double VSSourceModel::integrate(const VSAAlgebra::Vec2D& xy, 
				double R1, double R2) 
  const
{
  const unsigned nobs = m_obs_xy.size();

  double z = 0;
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      VSACoord::Polar<VSModelCoord> rphi1(xy,R1,0);
      VSACoord::Polar<VSModelCoord> rphi2(xy,R2,2*M_PI);

      rphi1.setObs(iobs);
      rphi2.setObs(iobs);

      typedef VSAFunction::ParamMemberFn<VSSourceModel, VSModelCoord > Fn;

      Fn fn(this,&VSSourceModel::val,param());

      try
	{
	  z += VSAMath::Quadrature::integrate(fn,rphi1,rphi2);
	}
      catch(const std::exception& e)
	{
	  std::cerr 
	    << std::string(__PRETTY_FUNCTION__) 
	    << ": Caught exception." << std::endl
	    << "X = " << xy.x() << " Y = " << xy.y() << std::endl
	    << "R1 = " << R1 << " R2 = " << R2 << " IOBS = " << iobs
	    << std::endl
	    << e.what() << std::endl;
	}
    }

  return z/m_domega;
}

void VSSourceModel::clear()
{
  m_effarea.clear();
  m_psf.clear();
  m_obs_xy.clear();
}

void VSSourceModel::setObs(const VSAAlgebra::Vec2D& xy, 
			   const VSNSpace& effarea,
			   const VSNSpace& psf)
{
  m_obs_xy.push_back(xy);
  m_effarea.push_back(effarea);
  m_psf.push_back(psf);

  initialize(psf);
}

// ----------------------------------------------------------------------------
// VSSourceModelGauss
// ----------------------------------------------------------------------------
VSSourceModelGauss::VSSourceModelGauss(double bin_size_deg):
  VSSourceModel(bin_size_deg), m_fn()
{ 
  setNParam(m_fn.nparm());

  setParam(1,0.065);
  fixParam(1,true);
}

double VSSourceModelGauss::val(const VSModelCoord& x) const
{    
  double acceptance = getAcceptance(x,param());    
  return m_fn.val(x)*acceptance*m_domega;
}
    
double VSSourceModelGauss::val(const VSModelCoord& x, 
			       const VSAAlgebra::VecND& a) const
{
  double acceptance = getAcceptance(x,a);  
  return m_fn.val(x,a)*acceptance*m_domega;
}

void VSSourceModelGauss::dyda(const VSModelCoord& x, 
			      VSAAlgebra::VecND& dyda) const
{
  double acceptance = getAcceptance(x,param()); 
  m_fn.dyda(x,dyda);
  dyda*=acceptance*m_domega;
}
    
void VSSourceModelGauss::dyda(const VSModelCoord& x, 
			      const VSAAlgebra::VecND& a,
			      VSAAlgebra::VecND& dyda) const
{
  double acceptance = getAcceptance(x,a);  
  m_fn.dyda(x,a,dyda);
  dyda*=acceptance*m_domega;
}

double VSSourceModelGauss::dyda(const VSModelCoord& x, unsigned ip) const
{
  double acceptance = getAcceptance(x,param());  
  VSAAlgebra::VecND dyda;
  m_fn.dyda(x,dyda);
  dyda*=acceptance;
  return dyda[ip]*m_domega;
}

void VSSourceModelGauss::setSourceParam(const std::vector<std::string>& params)
{
  if(params.size() == 1)
    {      
      double p;
      VSDataConverter::fromString(p, params[0]);

      m_fn.setParam(1,p);
      fixParam(1,true);
    }
  else 
    {
      m_fn.setParam(1,0.065);
      fixParam(1,true);
    }
}

double VSSourceModelGauss::getAcceptance(const VSModelCoord& x, 
					 const VSAAlgebra::VecND& a) const
{
  VSAAlgebra::Vec2D xy(a(2),a(3));
  double dtheta = xy.d(m_obs_xy[x.iobs()]);

  VSNSpace::Point p(1);
  p.x[0] = dtheta;
  double w;

  double lo = m_effarea[x.iobs()].axis(0).lo_bound + 
    m_effarea[x.iobs()].axis(0).bin_size/2.;

  double hi = m_effarea[x.iobs()].axis(0).hi_bound - 
    m_effarea[x.iobs()].axis(0).bin_size/2.;

  if(dtheta <= lo)
    {
      p.x[0] = lo;
      m_effarea[x.iobs()].getWeight(p,w);
    }
  else if(dtheta >= hi)
    {
      p.x[0] = hi;
      m_effarea[x.iobs()].getWeight(p,w);
    }
  else
    m_effarea[x.iobs()].interpolateWeight(p,w);

  return w;
}

// ----------------------------------------------------------------------------
// VSSourceModelPointSource
// ----------------------------------------------------------------------------
VSSourceModelPointSource::VSSourceModelPointSource(double bin_size_deg):
  VSSourceModel(bin_size_deg)
{ 
  setNParam(3);
}

double VSSourceModelPointSource::val(const VSModelCoord& x) const
{    
  double acceptance = getAcceptance(x,param());    
  return param(0)*acceptance*m_domega;
}
    
double VSSourceModelPointSource::val(const VSModelCoord& x, 
				     const VSAAlgebra::VecND& a) const
{
  double acceptance = getAcceptance(x,a);  
  return a(0)*acceptance*m_domega;
}

void VSSourceModelPointSource::dyda(const VSModelCoord& x, 
				    VSAAlgebra::VecND& dyda) const
{
  dyda.resize(3);

  double acceptance = getAcceptance(x,param()); 

  dyda(0) = acceptance*m_domega;
  dyda(1) = 0;
  dyda(2) = 0;
}
    
void VSSourceModelPointSource::dyda(const VSModelCoord& x, 
				    const VSAAlgebra::VecND& a,
				    VSAAlgebra::VecND& dyda) const
{
  dyda.resize(3);

  double acceptance = getAcceptance(x,a);  

  dyda(0) = acceptance*m_domega;
  dyda(1) = 0;
  dyda(2) = 0;
}

double VSSourceModelPointSource::dyda(const VSModelCoord& x, unsigned ip) const
{
  double acceptance = getAcceptance(x,param());  
  VSAAlgebra::VecND dyda;
  dyda*=acceptance;
  return dyda[ip]*m_domega;
}

void VSSourceModelPointSource::
setSourceParam(const std::vector<std::string>& param)
{

}

double VSSourceModelPointSource::getAcceptance(const VSModelCoord& x, 
					       const VSAAlgebra::VecND& a) 
  const
{
  VSAAlgebra::Vec2D src_xy(a(1),a(2));
  VSAAlgebra::Vec2D xy(x.x(),x.y());

  double offset = src_xy.d(m_obs_xy[x.iobs()]);
  double dtheta = src_xy.d(xy);

  VSNSpace::Point p(2);
  p.x[0] = offset;
  p.x[1] = dtheta;
  double w;

  double offset_lo = m_psf[x.iobs()].axis(0).lo_bound + 
    m_psf[x.iobs()].axis(0).bin_size/2.;

  double offset_hi = m_psf[x.iobs()].axis(0).hi_bound - 
    m_psf[x.iobs()].axis(0).bin_size/2.;

  double dtheta_lo = m_psf[x.iobs()].axis(1).lo_bound + 
    m_psf[x.iobs()].axis(1).bin_size/2.;

  double dtheta_hi = m_psf[x.iobs()].axis(1).hi_bound - 
    m_psf[x.iobs()].axis(1).bin_size/2.;

  if(dtheta <= dtheta_lo) p.x[1] = dtheta_lo;
  else if(dtheta >= dtheta_hi) return 0;
  
  if(offset <= offset_lo) p.x[0] = offset_lo;
  else if(offset >= offset_hi) p.x[0] = offset_hi;

  m_psf[x.iobs()].interpolateWeight(p,w);

  //  std::cout << offset << " " << dtheta << " " << w << std::endl;

  return w;
}

// ----------------------------------------------------------------------------
// VSSourceModelExtSource
// ----------------------------------------------------------------------------
double VSSourceModelExtSource::ConvolveFn::val(const VSACoord::Coord2D &a)
{
  double r = a.r();
  double phi = a.phi();
  
  double x = 
    sqrt(r*r + std::pow(m_theta,2) - 2*r*m_theta*cos(phi));
    
  //   std::cout << r << " " << x << " " << phi << std::endl;

  double w = 0;

  VSNSpace::Point p(1);
  p.x[0] = x;
  m_src_model.interpolateWeight(p,w);
  
  VSNSpace::Point p2(2);
  p2.x[0] = 0.5;
  p2.x[1] = r;

  //  double gw = std::exp(-r*r/(2*std::pow(sigma,2)))/(2*M_PI*sigma*sigma);
  double gw = 0;
  m_effarea_domega.interpolateWeight(p2,gw);
  
  return gw*w;
}

VSSourceModelExtSource::
VSSourceModelExtSource(double bin_size_deg,
		       const std::string& ext_model_file):
  VSSourceModel(bin_size_deg), m_model()
{ 
  setNParam(3);

  std::string file = ext_model_file;

  VSFileUtility::expandFilename(file);
  if(!VSFileUtility::isFile(file))
    {
      std::cerr << "Error opening model file " << file << std::endl;
      exit(1);
    }

  std::ifstream datastream(file.c_str());

  std::vector<double> x;
  std::vector<double> y;

  VSLineTokenizer tokenizer;
  VSTokenList tokens;
  std::string line;
  while(getline(datastream,line))
    {
      std::string line_copy = line;
      tokenizer.tokenize(line_copy, tokens);
      if(line_copy.substr(0,1) == "#" || tokens.size() == 0) continue;
      else if(tokens.size() < 2) continue;

      double theta;
      double dpdomega;

      VSDataConverter::fromString(theta,tokens[0].string());
      VSDataConverter::fromString(dpdomega,tokens[1].string());

      x.push_back(theta);
      y.push_back(dpdomega);
    }
  
  VSNSpace::Space space(1);
  
  double delta = *(x.begin()+1) - *x.begin();

  space.axes[0] = VSNSpace::Axis(0.,0.5,delta);
  m_model = VSNSpace(space);

  m_model.clear(0.0);

  for(unsigned i = 0; i < m_model.axis(0).nbin; i++)
    {
      if(i >= y.size()) break;
      
      m_model[i] = y[i];
    }

  double norm = 0;

  const double dtheta = 0.001;
  for(double theta = 0.; theta < 0.5; theta += dtheta)
    {
      VSNSpace::Point p(1);
      p.x[0] = theta;

      double w = 0;
      m_model.interpolateWeight(p,w);

      norm += 2*M_PI*w*theta*dtheta;      
    }

  m_model *= (1./norm);
}

double VSSourceModelExtSource::val(const VSModelCoord& x) const
{    
  double acceptance = getAcceptance(x,param());    
  return param(0)*acceptance*m_domega;
}
    
double VSSourceModelExtSource::val(const VSModelCoord& x, 
				   const VSAAlgebra::VecND& a) const
{
  double acceptance = getAcceptance(x,a);  
  return a(0)*acceptance*m_domega;
}

void VSSourceModelExtSource::dyda(const VSModelCoord& x, 
				    VSAAlgebra::VecND& dyda) const
{
  dyda.resize(3);

  double acceptance = getAcceptance(x,param()); 

  dyda(0) = acceptance*m_domega;
  dyda(1) = 0;
  dyda(2) = 0;
}
    
void VSSourceModelExtSource::dyda(const VSModelCoord& x, 
				    const VSAAlgebra::VecND& a,
				    VSAAlgebra::VecND& dyda) const
{
  dyda.resize(3);

  double acceptance = getAcceptance(x,a);  

  dyda(0) = acceptance*m_domega;
  dyda(1) = 0;
  dyda(2) = 0;
}

double VSSourceModelExtSource::dyda(const VSModelCoord& x, unsigned ip) const
{
  double acceptance = getAcceptance(x,param());  
  VSAAlgebra::VecND dyda;
  dyda*=acceptance;
  return dyda[ip]*m_domega;
}

void VSSourceModelExtSource::
setSourceParam(const std::vector<std::string>& param)
{

}

double VSSourceModelExtSource::getAcceptance(const VSModelCoord& x, 
					       const VSAAlgebra::VecND& a) 
  const
{
  VSAAlgebra::Vec2D src_xy(a(1),a(2));
  VSAAlgebra::Vec2D xy(x.x(),x.y());

  //  double offset = src_xy.d(m_obs_xy[x.iobs()]);
  double dtheta = src_xy.d(xy);

  VSNSpace::Point p(1);
  p.x[0] = dtheta;
  double w;

  double dtheta_lo = m_effarea_domega[x.iobs()].axis(0).lo_bound + 
    m_effarea_domega[x.iobs()].axis(0).bin_size/2.;

  double dtheta_hi = m_effarea_domega[x.iobs()].axis(0).hi_bound - 
    m_effarea_domega[x.iobs()].axis(0).bin_size/2.;

  if(dtheta <= dtheta_lo) p.x[0] = dtheta_lo;
  else if(dtheta >= dtheta_hi) return 0;
  
  m_effarea_domega[x.iobs()].interpolateWeight(p,w);

  return w;
}

void VSSourceModelExtSource::initialize(const VSNSpace& psf)
{
  VSNSpace::Space s(1);
  s.axes[0] = psf.axis(1);

  VSNSpace psf2 = VSNSpace(s);

  VSNSpace::Cell c(1);
  //  c.i[0] = psf.axis(0).indexUnchecked(0.5);

  for(c.i[0] = 0; c.i[0] < psf2.axis(0).nbin; c.i[0]++)
    {
      double theta = psf2.axis(0).midCoordUnchecked(c.i[0]);

      ConvolveFn fn(m_model,psf,theta);

      VSACoord::Polar<VSACoord::Coord2D> x1(0.0,0.0);
      VSACoord::Polar<VSACoord::Coord2D> x2(0.4,M_PI);

      psf2[c] = 2*VSAMath::Quadrature::integrate(fn,x1,x2);

      //      std::cout << theta << " " << psf2[c] << std::endl;
    }


  m_effarea_domega.push_back(psf2);
}

void VSSourceModelExtSource::clear()
{
  m_effarea_domega.clear();
  VSSourceModel::clear();
}
