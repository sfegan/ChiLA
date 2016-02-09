//-*-mode:c++; mode:font-lock;-*-

/*! \file VSIntegralAnalysis.cpp

  Base class for integral (energy-independent) analysis calculators.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.10 $
  \date       07/29/2007

  $Id: VSIntegralAnalysis.cpp,v 3.10 2010/10/20 04:02:07 matthew Exp $

*/

#include <fstream>
#include <Astro.h>
#include <VSIntegralAnalysis.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSBkgndModel.hpp>
#include <VSOctaveH5Writer.hpp>
#include <VSApertureAnalysis.hpp>
#include <VSIntegralAnalysisData.hpp>
#include <VSReflectedRegionAnalysis.hpp>
#include <VSRingBackgroundAnalysis.hpp>
#include <VSMLMAnalysis.hpp>
#include <VSAStatistics.hpp>

using namespace VERITAS;
using namespace SEphem;

VSIntegralAnalysis::
VSIntegralAnalysis(const std::string bkgnd_model,
		   double bin_size_deg,
		   double theta_cut,
		   const std::pair< double, double >& ring_cut,
		   double offset_max,
		   unsigned max_off_regions,
		   const std::string& spmodel):
  m_bin_size_deg(bin_size_deg),
  m_domega(std::pow(bin_size_deg,2)),
  m_theta_cut(theta_cut),
  m_ring_cut(ring_cut),
  m_offset_max(offset_max),
  m_max_off_regions(max_off_regions),
  m_spmodel(spmodel),
  m_off_xy(),
  m_exclusion_region(),
  m_bkgnd_model(),
  m_spectrum_model()
{
  std::vector< std::string > bkgnd_param;
  VSDataConverter::fromString(bkgnd_param, bkgnd_model);

  vsassert(bkgnd_param.size() >=1);

  std::string model = bkgnd_param[0];
  bkgnd_param.erase(bkgnd_param.begin());

  std::string param = VSDataConverter::toString(bkgnd_param);

  if(model == "bessel2")
    m_bkgnd_model = new VSAcceptanceModelBessel2(param,
						 bin_size_deg, offset_max);
  else if(model == "poly")
    m_bkgnd_model = new VSAcceptanceModelPoly(param,
					      bin_size_deg, offset_max);
  else
    {
      std::cerr << "VSIntegralAnalysis(): Unknown background model: " 
		<< bkgnd_model << std::endl;
      exit(EXIT_FAILURE);
    }

  m_spectrum_model = VSSpectrumFn::create(spmodel);
  m_spectrum_model->setNormalization(1.0);
}

VSIntegralAnalysis::~VSIntegralAnalysis()
{
  delete m_bkgnd_model;
}

void VSIntegralAnalysis::accumulate(VSAnalysisStage3Data::RunData& d,
				    Data& data)
{
  d.irf_calc().getEffectiveAreaOffset(m_spmodel,       
				      data.effarea_nspace,
				      data.effarea_psf_nspace);
}

void VSIntegralAnalysis::accumulate(const VSEventArrayDatum& event,
				    const VSAAlgebra::Vec2D& event_xy,
				    Data& data)
{
  data.sky_counts_hist.accumulate(event_xy.x(),event_xy.y());

  // Accumulate on counts -------------------------------------------------
  double dr_src = data.src_data.src_xy.d(event_xy);
  
  if(dr_src < data.src_data.theta_cut)
    data.src_data.on_counts++;
  
  // Accumulate ring counts -----------------------------------------------
  if(event_xy.d(data.src_data.src_xy) > m_ring_cut.first &&
     event_xy.d(data.src_data.src_xy) < m_ring_cut.second)
    data.src_data.ring_counts++;

  // Accumulate off counts ------------------------------------------------
  for(std::vector< VSAAlgebra::Vec2D >::const_iterator itr =
	data.src_data.off_xy.begin(); itr != data.src_data.off_xy.end(); ++itr)
    {
      double dr_off = itr->d(event_xy);
      if(dr_off < data.src_data.theta_cut) data.src_data.off_counts++;
    }
}

VSAcceptanceData* VSIntegralAnalysis::fitAcceptance(VSAnalysisStage3Data& data)
{
  std::cout << "Entering VSIntegralAnalysis::fitAcceptance()" << std::endl;

  VSAMath::Data<VSModelCoord> fit_data;

  const unsigned nptg = data.m_ptg_xy.size();
  std::vector< VSAAlgebra::Vec2D >& ptg_xy = data.m_ptg_xy;
  std::vector< VSSimple2DHist<double,double> > sky_counts_hists(nptg);

  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      unsigned iptg = 0;

      for(iptg = 0; iptg < nptg; iptg++)
	if(ptg_xy[iptg] == itr->ptg_xy()) break;

      vsassert(iptg != nptg);

      sky_counts_hists[iptg].merge(itr->sky_counts_hist);
    }


  // Construct the data model -------------------------------------------------
  m_bkgnd_model->clear();
  for(unsigned iptg = 0; iptg < nptg; iptg++)
    {
      m_bkgnd_model->setObs(ptg_xy[iptg]);

      double zcount = 0;
      for(VSSimple2DHist<double,double>::iterator itr = 
	    sky_counts_hists[iptg].begin(); itr !=
	    sky_counts_hists[iptg].end(); ++itr)
	{      
	  VSAAlgebra::Vec2D xy(itr->x(),itr->y());

	  if(ptg_xy[iptg].d(xy) > m_offset_max) continue;
	  else if(m_exclusion_region.isExcluded(xy)) continue;
	  zcount+=itr->count();
	}

      // Skip this pointing if all bins are empty -----------------------------
      if(zcount == 0) continue;

      // Create the dataset ---------------------------------------------------
      for(VSSimple2DHist<double,double>::iterator itr = 
	    sky_counts_hists[iptg].begin(); itr !=
	    sky_counts_hists[iptg].end(); ++itr)
	{      
	  VSAAlgebra::Vec2D xy(itr->x(),itr->y());

	  if(ptg_xy[iptg].d(xy) > m_offset_max) continue;
	  else if(m_exclusion_region.isExcluded(xy)) continue;
	
	  double z = itr->count();
	  VSModelCoord c(itr->x(),itr->y(),iptg);
	  VSAMath::DataPoint<VSModelCoord> p(c,z,sqrt(z));
	  fit_data.insert(p);
	}
    }

  VSAcceptanceData* d = m_bkgnd_model->fit(fit_data);
  data.setAcceptanceData(d);

  data.acceptance_data()->fovr2_bkgnd_hist = data.fovr2_counts_hist;
  data.acceptance_data()->fovr2_bkgnd_hist.clear();
  data.acceptance_data()->fovr2_acceptance_hist = data.fovr2_counts_hist;
  data.acceptance_data()->fovr2_acceptance_hist.clear();
  data.acceptance_data()->sky_bkgnd_hist = data.sky_counts_hist;
  data.acceptance_data()->sky_bkgnd_hist.clear();

  std::cout << "VSIntegralAnalysis::fitAcceptance(): Fit Parameters..."
	    << std::endl;

  for(unsigned ip = 0; ip < d->param.ndim(); ip++)
    std::cout << std::setw(10) << ip 
	      << std::setw(15) << d->param(ip) 
	      << std::setw(15) << d->param_err(ip) 
	      << std::endl;

  std::cout << "LNL " << std::setw(15) << d->chi2 << std::endl;

  for(VSSimple2DHist<double,double>::iterator itr = 
	d->fov_acceptance_hist.begin(); itr != d->fov_acceptance_hist.end();
      ++itr)
    {
      VSAAlgebra::Vec2D xy(itr->x(),itr->y());
      double r = sqrt(std::pow(xy.x(),2) + std::pow(xy.y(),2));

      if(r > m_offset_max) continue;

      d->fov_acceptance_hist.setBin(itr->bin(),
				    m_bkgnd_model->getFoVAcceptance(xy));
      d->fov_bkgnd_hist.setBin(itr->bin(),
			       m_domega*m_bkgnd_model->getFoVCounts(xy));
    }

  // Create a sky acceptance hist for every run -------------------------------
//   for(std::vector< VSAnalysisStage3Data::RunData >::iterator itr = 
// 	data.run_data.begin(); itr != data.run_data.end(); ++itr)
//     {      
//       VSSimple2DHist<double,double>& h = itr->sky_acceptance_hist;
//       for(VSSimple2DHist<double,double>::iterator hitr = h.begin(); 
// 	  hitr != h.end(); ++hitr)
// 	{
// 	  if(itr->sky_exposure_hist.countForIndex(hitr->bin())==0) continue;
// 	  VSAAlgebra::Vec2D xy(hitr->x(),hitr->y());
// 	  xy -= itr->obs_xy;
// 	  h.setBin(hitr->bin(),m_bkgnd_model->getAcceptance(xy));
// 	}      
//     }
    
  for(VSLimitedErrorsHist<double,double>::iterator itr = 
	data.fovr2_counts_hist.begin(); itr != 
	data.fovr2_counts_hist.end(); ++itr)
    {
      double R1 = sqrt(itr->val());
      double R2 = sqrt(itr->val()+data.fovr2_counts_hist.binSize());

      double v1 = m_bkgnd_model->integrateFoV(R1,R2);
      double v2 = m_bkgnd_model->getFoVAcceptance(itr->center());

      data.acceptance_data()->fovr2_bkgnd_hist.setBin(itr->bin(),v1,0);
      data.acceptance_data()->fovr2_acceptance_hist.setBin(itr->bin(),v2,0);
    }
  
  // Fill bkgnd distribution for Theta-Squared --------------------------------
  data.acceptance_data()->th2_bkgnd_hist.fill(0);
  for(VSLimitedErrorsHist<double,double>::iterator itr =
	data.acceptance_data()->th2_bkgnd_hist.begin(); 
      itr !=  data.acceptance_data()->th2_bkgnd_hist.end(); ++itr)
    {
      VSAAlgebra::Vec2D xy(0,0);

      double R1 = sqrt(itr->val());
      double R2 = sqrt(itr->val()+ 
		       data.acceptance_data()->th2_bkgnd_hist.binSize());

      double s = m_bkgnd_model->integrate(xy,R1,R2);

      data.acceptance_data()->th2_bkgnd_hist.setBin(itr->bin(),s,0);
    }

  return d;
}

VSIntegralAnalysis::Data 
VSIntegralAnalysis::create(const VSAAlgebra::Vec2D& ptg_xy,
			   const VSAAlgebra::Vec2D& obs_xy)
{
  Data d;
  d.ptg_xy = ptg_xy;
  d.obs_xy = obs_xy;

  d.src_data.src_xy = m_src_xy;
  d.src_data.theta_cut = m_theta_cut;

  double bin_size = m_bin_size_deg;

  double x1 = (lround((obs_xy.x()-m_offset_max-0.2)/bin_size)+0.5)*bin_size;
  double x2 = (lround((obs_xy.x()+m_offset_max+0.2)/bin_size)+0.5)*bin_size;
  double y1 = (lround((obs_xy.y()-m_offset_max-0.2)/bin_size)+0.5)*bin_size;
  double y2 = (lround((obs_xy.y()+m_offset_max+0.2)/bin_size)+0.5)*bin_size;

  d.sky_counts_hist = 
    VSSimple2DHist<double, double>(bin_size,x1,x2,bin_size,y1,y2); 

  VSOffRegion::getRegions(d.src_data.src_xy, d.obs_xy,
			  m_theta_cut, m_max_off_regions,
			  d.src_data.off_xy);

  return d;
}


void VSIntegralAnalysis::
setSourcePosition(const SEphem::SphericalCoords& origin_radec,
		  const SEphem::SphericalCoords& src_radec)
{
  m_origin_radec = origin_radec;
  m_src_radec = src_radec;
  std::pair< Angle, Angle > xy;
  Astro::raDecToXY(m_src_radec,m_origin_radec,xy);  
  m_src_xy.set(xy.first.degPM(),xy.second.degPM());
}

void VSIntegralAnalysis::calcFlux(VSIntegralAnalysisDatum& o)
{
  o.dfde100 = m_spectrum_model->val(-1)*o.dfde;
  o.dfde100_err = fabs(o.dfde100*o.dfde_err/o.dfde);
  o.dfde316 = m_spectrum_model->val(-0.5)*o.dfde;
  o.dfde316_err = fabs(o.dfde316*o.dfde_err/o.dfde);
  o.dfde1000 = m_spectrum_model->val(0.0)*o.dfde;
  o.dfde1000_err = fabs(o.dfde1000*o.dfde_err/o.dfde);
  o.dfde_eth = m_spectrum_model->val(std::log10(o.egy_threshold))*o.dfde;
  o.dfde_eth_err = fabs(o.dfde_eth*o.dfde_err/o.dfde);

  o.flux = m_spectrum_model->integralFlux(m_spectrum_model->normEnergy(),3.)*
    o.dfde;
  o.flux_err = fabs(o.flux*o.dfde_err/o.dfde);
  o.flux100 = m_spectrum_model->integralFlux(-1.,3.)*o.dfde;
  o.flux100_err = fabs(o.flux100*o.dfde_err/o.dfde);
  o.flux316 = m_spectrum_model->integralFlux(-0.5,3.)*o.dfde;
  o.flux316_err = fabs(o.flux316*o.dfde_err/o.dfde);
  o.flux1000 = m_spectrum_model->integralFlux(0.0,3.)*o.dfde;
  o.flux1000_err = fabs(o.flux1000*o.dfde_err/o.dfde);
  o.flux_eth = 
    m_spectrum_model->integralFlux(std::log10(o.egy_threshold),3.)*o.dfde;
  o.flux_eth_err = fabs(o.flux_eth*o.dfde_err/o.dfde);

  o.flux_ul95 = VSAStatistics::heleneUL(o.flux,o.flux_err,0.95);
  o.flux_ul99 = VSAStatistics::heleneUL(o.flux,o.flux_err,0.99);
  o.flux100_ul95 = VSAStatistics::heleneUL(o.flux100,o.flux100_err,0.95);
  o.flux100_ul99 = VSAStatistics::heleneUL(o.flux100,o.flux100_err,0.99);
  o.flux316_ul95 = VSAStatistics::heleneUL(o.flux316,o.flux316_err,0.95);
  o.flux316_ul99 = VSAStatistics::heleneUL(o.flux316,o.flux316_err,0.99);
  o.flux1000_ul95 = VSAStatistics::heleneUL(o.flux1000,o.flux1000_err,0.95);
  o.flux1000_ul99 = VSAStatistics::heleneUL(o.flux1000,o.flux1000_err,0.99);
  o.flux_eth_ul95 = VSAStatistics::heleneUL(o.flux_eth,o.flux_eth_err,0.95);
  o.flux_eth_ul99 = VSAStatistics::heleneUL(o.flux_eth,o.flux_eth_err,0.99);

  o.dfde_ul95 = VSAStatistics::heleneUL(o.dfde,o.dfde_err,0.95);
  o.dfde_ul99 = VSAStatistics::heleneUL(o.dfde,o.dfde_err,0.99);
  o.dfde100_ul95 = VSAStatistics::heleneUL(o.dfde100,o.dfde100_err,0.95);
  o.dfde100_ul99 = VSAStatistics::heleneUL(o.dfde100,o.dfde100_err,0.99);
  o.dfde316_ul95 = VSAStatistics::heleneUL(o.dfde316,o.dfde316_err,0.95);
  o.dfde316_ul99 = VSAStatistics::heleneUL(o.dfde316,o.dfde316_err,0.99);
  o.dfde1000_ul95 = VSAStatistics::heleneUL(o.dfde1000,o.dfde1000_err,0.95);
  o.dfde1000_ul99 = VSAStatistics::heleneUL(o.dfde1000,o.dfde1000_err,0.99);
  o.dfde_eth_ul95 = VSAStatistics::heleneUL(o.dfde_eth,o.dfde_eth_err,0.95);
  o.dfde_eth_ul99 = VSAStatistics::heleneUL(o.dfde_eth,o.dfde_eth_err,0.99);

  o.e2dfde_ul95 = o.dfde_ul95*std::pow(10,2*m_spectrum_model->normEnergy());
  o.e2dfde_ul99 = o.dfde_ul99*std::pow(10,2*m_spectrum_model->normEnergy());
  o.e2dfde100_ul95 = o.dfde100_ul95*std::pow(10,-2.);
  o.e2dfde100_ul99 = o.dfde100_ul99*std::pow(10,-2.);
  o.e2dfde316_ul95 = o.dfde316_ul95*std::pow(10,-1.);
  o.e2dfde316_ul99 = o.dfde316_ul99*std::pow(10,-1.);
  o.e2dfde1000_ul95 = o.dfde1000_ul95;
  o.e2dfde1000_ul99 = o.dfde1000_ul99;
  o.e2dfde_eth_ul95 = 
    o.dfde_eth_ul95*std::pow(10,2.*std::log10(o.egy_threshold));
  o.e2dfde_eth_ul99 = 
    o.dfde_eth_ul99*std::pow(10,2.*std::log10(o.egy_threshold));
}

// ============================================================================
// VSIntegralAnalysisFactory
// ============================================================================
VSIntegralAnalysisFactory::Options::Options():
  method("rbm"),
  bkgnd_model("bessel2"),
  theta_cut(0.15),
  ring_cut(0.4,0.5),
  max_noff_region(10),
  integral_source_spectrum("powerlaw,2.5"),
  mlm_fit_skymap(true),
  mlm_source_model("pointsource"),
  mlm_ext_model(),
  mlm_skymap_fit_radius(0.3)
{

}

std::auto_ptr<VSIntegralAnalysisFactory> VSIntegralAnalysisFactory::s_instance;

VSIntegralAnalysisFactory::Options 
VSIntegralAnalysisFactory::s_default_options = 
		   VSIntegralAnalysisFactory::Options();

VSIntegralAnalysisFactory::VSIntegralAnalysisFactory(const Options& opt):
  m_options(opt)
{

}

VSIntegralAnalysisFactory::~VSIntegralAnalysisFactory()
{

}

VSIntegralAnalysisFactory* VSIntegralAnalysisFactory::getInstance()
{
  if(s_instance.get() == 0)s_instance.reset(new VSIntegralAnalysisFactory());
  return s_instance.get();
}

VSIntegralAnalysis* VSIntegralAnalysisFactory::create(double max_offset_cut, 
					      double bin_width_deg)
{
  if(m_options.method == "rrm")
    {
      return new VSReflectedRegionAnalysis(m_options.bkgnd_model,
					   bin_width_deg,
					   m_options.theta_cut,
					   m_options.ring_cut,
					   max_offset_cut,
					   m_options.max_noff_region,
					   m_options.integral_source_spectrum); 

    }
  else if(m_options.method == "rbm")
    {
      return new VSRingBackgroundAnalysis(m_options.bkgnd_model,
					  bin_width_deg,
					  m_options.theta_cut,
					  m_options.ring_cut,
					  max_offset_cut,
					  m_options.max_noff_region,
					  m_options.integral_source_spectrum);

    }
 else if(m_options.method == "mlm")
    {
      return new VSMLMAnalysis(m_options.bkgnd_model,
			       bin_width_deg,
			       m_options.theta_cut,
			       m_options.ring_cut,
			       max_offset_cut,
			       m_options.integral_source_spectrum,
			       m_options.mlm_fit_skymap,
			       m_options.mlm_source_model,
			       m_options.mlm_ext_model,
			       m_options.mlm_skymap_fit_radius);
    }
  else
    {
      std::cerr 
	<< __PRETTY_FUNCTION__
	<< ": Unrecognized integral analysis method: " 
	<< m_options.method << std::endl
	<< " Valid methods are: rrm, rbm, mlm" 
	<< std::endl;
      exit(EXIT_FAILURE);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSIntegralAnalysisFactory::
configure(VSOptions& options, const std::string& profile, 
	  const std::string& opt_prefix)
{
  const std::string category = "integral";

  options.addCatagory(OPTNAME(opt_prefix,category),
		      "Options related to integral analysis.");

  options.findWithValue(OPTNAME(opt_prefix,"bkgnd_model"),
			s_default_options.bkgnd_model,
			"Set the parametric model for the cosmic-ray "
			"acceptance. "
			"Options are: poly, bessel2.",
			OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"method"),
			s_default_options.method,
			"Set the integral analysis method used for background "
			"modeling and signal extraction.  Available "
			"methods are: \"rrm\" (reflected region), "
			"\"rbm\" (ring background), "
			"\"mlm\" (maximum likelihood)",
			OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"theta_cut"), 
			s_default_options.theta_cut,
			"Signal aperture radius in degrees used to calculate "
			"rates and "
			"significances with the reflected region and ring "
			"background methods.",
			OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"ring_cut"), 
			s_default_options.ring_cut,
			"Define inner and outer ring radius in degrees used "
			"for background estimation by the ring background "
			"method.",OPTNAME(opt_prefix,category));
  
  options.findWithValue(OPTNAME(opt_prefix,"max_noff_region"), 
			s_default_options.max_noff_region,
			"Maximum number of off regions that will be used "
			"when estimating background with the reflected "
			"region method.",OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"integral_source_spectrum"),
			s_default_options.integral_source_spectrum,
			"Set the spectral parameterization of "
			"the integral source model.  Fluxes "
			"will be calculated from the integral counts using "
			"this spectrum.  This option also determines "
			"the psf model used by the \"mlm\" method.  "
			"Spectral parameterization should "
			"be specified as a spectral type followed by one or "
			"more parameters (e.g. powerlaw,2.5). "
			"Available spectral types are: "
			"powerlaw,powerlawexp",OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"mlm_source_model"),
			s_default_options.mlm_source_model,
			"Set the morphological parameterization "
			"for the MLM gamma-ray source model. "
			"Options are: gauss,pointsource,ext",
			OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"mlm_ext_model"),
			s_default_options.mlm_ext_model,
			"Set the file containing the definition of the "
			"extended source model.",
			OPTNAME(opt_prefix,category));

  options.findWithValue(OPTNAME(opt_prefix,"mlm_skymap_fit_radius"),
			s_default_options.mlm_skymap_fit_radius,
			"Set the radius in degrees out to which the source "
			"model should be fit when generating sky maps.  This "
			"option should generally be set to a value slightly "
			"larger than the extent of the source model.",
			OPTNAME(opt_prefix,category));

  options.findBoolValue(OPTNAME(opt_prefix,"mlm_fit_skymap"),
			s_default_options.mlm_fit_skymap, true,
			"Generate sky maps when running with mlm method.",
			OPTNAME(opt_prefix,category));

}
