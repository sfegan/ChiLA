//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSourceInjector.hpp

  Inject a simulated source into the data.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSSourceInjector.cpp,v 3.13 2010/06/21 23:58:22 matthew Exp $

*/

#include <VSAQuadrature.hpp>
#include <VSSourceInjector.hpp>

using namespace VERITAS;
using namespace SEphem;

std::string VSSourceInjector::s_default_src_type = "gauss";
std::string VSSourceInjector::s_default_src_spectrum = "powerlaw,2.5,0";
double VSSourceInjector::s_default_src_rate = 0;
std::pair<double,double> VSSourceInjector::s_default_src_xy = 
std::pair<double,double>(0,0);

VSSourceInjector::
VSSourceInjector(RandomNumbers* rng,
		 double theta,
		 const std::pair< double, double >& ring_radius,
		 const std::string& src_type,
		 const std::string& src_spectrum,
		 double src_rate,
		 const std::pair<double,double>& src_xy):
  m_theta(theta), m_ring_radius(ring_radius), 
  m_src_type(src_type),
  m_src_spectrum(src_spectrum),
  m_src_rate(src_rate),
  m_src_xy(src_xy.first,src_xy.second), 
  m_emin(-1.5),
  m_emax(2.5),
  m_rng(rng),
  m_spectrum_fn()
{
  m_spectrum_fn = VSSpectrumFn::create(m_src_spectrum);
}

VSSourceInjector::~VSSourceInjector()
{
  delete m_spectrum_fn;
}

void VSSourceInjector::initialize(VSAnalysisStage3Data::RunData* data)
{
  const double ebin = 0.005;
  const double offset = (data->obs_xy()-m_src_xy).norm();
  
  VSAFunction::Spline erec_rate;
  for(double erec = m_emin; erec < m_emax; erec += ebin)
    {
      VSAFunction::Spline spline;
      for(double emc = m_emin; emc < m_emax; emc += ebin)
	{
	  double effarea = 
	    std::max(0.,data->irf_calc().getEffectiveArea(offset,m_theta,emc));
	  double dfde = m_spectrum_fn->val(emc);
	  double dpde = data->irf_calc().getKernel(offset,emc,erec);
	  
// 	  std::cout << std::setw(15) << erec 
// 		    << std::setw(15) << emc
// 		    << std::setw(15) << kw
// 		    << std::setw(15) 
// 		    << std::endl;
// 		    << data->egy_kernel().getKernel(offset,emc,erec)
// 		    << std::endl;

	  spline.setPoint(emc,dpde*effarea*dfde*std::pow(10,emc)*log(10));
	}

      spline.spline();
      double s = VSAMath::Quadrature::integrate(spline,m_emin,m_emax,1E-6);
      erec_rate.setPoint(erec,s);
    }

  erec_rate.spline();

  VSAFunction::Spline erec_cdf;
  erec_cdf.setPoint(m_emin,0);

  double rsum = 0;
  std::vector< std::pair<double,double> > xy;
  for(double x = m_emin+ebin; x < m_emax; x+= ebin)
    {
      double s = VSAMath::Quadrature::integrate(erec_rate,x-ebin,x,1E-6);
      if(s<0) s = 0;

      rsum += s;
      xy.push_back(std::make_pair(x,rsum));
    }

  for(std::vector< std::pair<double,double> >::iterator itr = xy.begin();
      itr != xy.end(); ++itr)
    erec_cdf.setPoint(itr->first,itr->second/rsum);

  erec_cdf.setPoint(m_emax,1);
  erec_cdf.spline();


  m_src_rate = VSAMath::Quadrature::integrate(erec_rate,m_emin,m_emax)*60;

  m_erec_rng = 
    RandomNumbers::CDFGenerator<VSAFunction::Spline>(erec_cdf,
 						     m_emin,
 						     m_emax,1000000);
}

void VSSourceInjector::fill(VSAnalysisStage3Data::RunData* data)
{
  unsigned nevent = m_rng->Poisson(m_src_rate*data->livetime_min); 
  if(nevent == 0) return;

  if(m_src_type == "gauss") fillGauss(nevent,data);
  else fillConstant(nevent,data);
}

void VSSourceInjector::fillGauss(unsigned nevent, 
				 VSAnalysisStage3Data::RunData* data)
{
  const double sigma = 0.065;
  for(unsigned ievent = 0; ievent < nevent; ievent++)
    {
      double x =  m_rng->Normal()*sigma + m_src_xy.x();
      double y =  m_rng->Normal()*sigma + m_src_xy.y();
      VSAAlgebra::Vec2D xy(x,y);

      double erec = m_erec_rng.rnd(m_rng);

      data->egy_on_hist.accumulate(erec);
      data->egy_on_np_hist.accumulate(erec);

      data->sky_counts_hist.accumulate(x,y);

      if(xy.d(data->src_xy()) > m_ring_radius.first &&
	 xy.d(data->src_xy()) < m_ring_radius.second)
	data->ring_counts++;

      if(xy.d(data->src_xy()) < m_theta)
	{
	  data->on_counts++;
	}
    }
}

void VSSourceInjector::fillConstant(unsigned nevent, 
				    VSAnalysisStage3Data::RunData* data)
{
  const double rmax = 1.5;

  double obsx = data->obs_xy().x();
  double obsy = data->obs_xy().y();
  
  for(unsigned ievent = 0; ievent < nevent*1.0; ievent++)
    {
      double x = data->obs_xy().x() + (0.5-m_rng->Uniform())*2*rmax;
      double y = data->obs_xy().y() + (0.5-m_rng->Uniform())*2*rmax;
      double r = sqrt(std::pow(x-obsx,2)+std::pow(y-obsy,2));
      while(r > 1.5)
	{
	  x = data->obs_xy().x() + (0.5-m_rng->Uniform())*2*rmax;
	  y = data->obs_xy().y() + (0.5-m_rng->Uniform())*2*rmax;
	  r = sqrt(std::pow(x-obsx,2)+std::pow(y-obsy,2));
	}

      VSAAlgebra::Vec2D xy(x,y);

      data->sky_counts_hist.accumulate(x,y);

      if(xy.d(data->src_xy()) > m_ring_radius.first &&
	 xy.d(data->src_xy()) < m_ring_radius.second)
	data->ring_counts++;

      if(xy.d(data->src_xy()) < m_theta) data->on_counts++;
    }

  double sigma = 1.1;

  for(unsigned ievent = 0; ievent < nevent*3.0; ievent++)
    {
      double x =  obsx + m_rng->Normal()*sigma;
      double y =  obsy + m_rng->Normal()*sigma;
      double r = sqrt(std::pow(x-obsx,2)+std::pow(y-obsy,2));
      while(r > 1.5)
	{
	  x = obsx + m_rng->Normal()*sigma;
	  y = obsy + m_rng->Normal()*sigma;
	  r = sqrt(std::pow(x-obsx,2)+std::pow(y-obsy,2));
	}

      VSAAlgebra::Vec2D xy(x,y);

      data->sky_counts_hist.accumulate(x,y);

      if(xy.d(data->src_xy()) > m_ring_radius.first &&
	 xy.d(data->src_xy()) < m_ring_radius.second)
	data->ring_counts++;

      if(xy.d(data->src_xy()) < m_theta) data->on_counts++;
    }
}

VSSourceInjector::Data 
VSSourceInjector::getData(double ebin, double elo, double ehi) const
{
  Data d;

  d.src_rate = m_src_rate;
  d.integral_flux_hist = VSLimitedErrorsHist<double,double>(ebin,elo,ehi);
  d.diff_flux_hist = VSLimitedErrorsHist<double,double>(ebin,elo,ehi);

  d.integral_flux_hist.fill(0);
  d.diff_flux_hist.fill(0);

  for(VSLimitedErrorsHist<double,double>::iterator itr = 
	d.integral_flux_hist.begin(); itr != d.integral_flux_hist.end(); ++itr)
    {
      double lo = itr->val();
      double hi = itr->val()+ebin;
      double integral_flux = m_spectrum_fn->integralFlux(lo,hi);
      //      double diff_flux = m_spectrum_fn->val(itr->center());
      double de = std::pow(10,hi) - std::pow(10,lo);

      d.integral_flux_hist.accumulate(itr->center(),integral_flux,0);
      d.diff_flux_hist.accumulate(itr->center(),integral_flux/de,0);
    }

  return d;
}

void
VSSourceInjector::fillIntegralHist(VSLimitedErrorsHist<double,double>& h) const
{
  h.clear();
  h.fill(0);

  for(VSLimitedErrorsHist<double,double>::iterator itr = h.begin();
      itr != h.end(); ++itr)
    {
      double elo = itr->val();
      double ehi = elo + h.binSize();

      double integral_flux = m_spectrum_fn->integralFlux(elo,ehi);

      h.accumulate(itr->center(), integral_flux,0);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSSourceInjector::configure(VSOptions& options, 
				 const std::string& profile, 
				 const std::string& opt_prefix)
{
  const std::string category = OPTNAME(opt_prefix,"fake_src");

  options.addCatagory(category,
		      "Options for injecting a source into real data.");

  options.findWithValue(OPTNAME(opt_prefix,"fake_src_xy"),
			s_default_src_xy,
			"Set the position in X/Y of the injected source.",
			category);

  options.findWithValue(OPTNAME(opt_prefix,"fake_src_rate"),
			s_default_src_rate,"Set the fake source injection "
			"rate in g/min.", category);

  options.findWithValue(OPTNAME(opt_prefix,"fake_src_type"),
			s_default_src_type,"",category);

  options.findWithValue(OPTNAME(opt_prefix,"fake_src_spectrum"),
			s_default_src_spectrum,"",category);
}
