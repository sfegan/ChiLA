//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMLMAnalysis.cpp

  Integral analysis with maximum likelihood method.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.7 $
  \date       07/29/2007

  $Id: VSMLMAnalysis.cpp,v 3.7 2010/10/20 03:33:16 matthew Exp $

*/

#include <Astro.h>
#include <VSAAlgebra.hpp>
#include <VSAStatistics.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSMLMAnalysis.hpp>
#include <VSAFunction.hpp>
#include <VSSourceModel.hpp>
#include <VSBkgndModel.hpp>

using namespace VERITAS;

VSMLMAnalysis::VSMLMAnalysis(const std::string& bkgnd_model,
			     double bin_size_deg,
			     double theta_cut, 
			     const std::pair< double, double >& ring_cut,
			     double theta_max,
			     const std::string& source_spectrum,
			     bool fit_skymap,
			     const std::string& source_model,
			     const std::string& ext_model,
			     double fit_radius):
  VSIntegralAnalysis(bkgnd_model,bin_size_deg,theta_cut, ring_cut,
		     theta_max,5,source_spectrum),
  m_domega(bin_size_deg*bin_size_deg),
  m_fit_skymap(fit_skymap), m_fit_radius(fit_radius), 
  m_source_model(), m_data_model(), m_fitter()
{
  std::vector< std::string > p;
  VSDataConverter::fromString(p, source_model);

  vsassert(p.size() >= 1);
  std::string model = p[0];

  p.erase(p.begin());

  if(model == "gauss")
    m_source_model = new VSSourceModelGauss(bin_size_deg);
  else if(model == "pointsource")
    m_source_model = new VSSourceModelPointSource(bin_size_deg);
  else if(model == "ext")
    m_source_model = new VSSourceModelExtSource(bin_size_deg,ext_model);
  else
    {
      std::cerr << "Error. Invalid source model type: " << model
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  m_source_model->setSourceParam(p);
}

VSMLMAnalysis::~VSMLMAnalysis()
{
  delete m_fitter;
}

void VSMLMAnalysis::analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			    const VSAnalysisStage3Data& data,
			    const VSAcceptanceData& acceptance,
			    VSIntegralAnalysisData& o)
{
  o.sky_livetime_hist = data.sky_livetime_hist;
  o.sky_significance_hist = data.sky_counts_hist;
  o.sky_on_counts_hist = data.sky_counts_hist;
  o.sky_excess_hist = data.sky_counts_hist;
  o.sky_excess_rate_hist = data.sky_counts_hist;
  o.sky_source_hist = data.sky_counts_hist;
  o.sky_flux_hist = data.sky_counts_hist;
  o.sky_flux_ul95_hist = data.sky_counts_hist;
  o.sky_bkgnd_hist = data.sky_counts_hist;
  o.sky_bkgnd_rate_hist = data.sky_counts_hist;
  o.sky_bkgnd_density_rate_hist = data.sky_counts_hist;

  o.sky_significance_hist.clear();
  o.sky_on_counts_hist.clear();
  o.sky_excess_hist.clear();
  o.sky_excess_rate_hist.clear();
  o.sky_source_hist.clear();
  o.sky_flux_hist.clear();
  o.sky_flux_ul95_hist.clear();
  o.sky_bkgnd_hist.clear();
  o.sky_bkgnd_rate_hist.clear();
  o.sky_bkgnd_density_rate_hist.clear();


  if(!loadFitData(data,acceptance)) return;

  // Fit the source model -----------------------------------------------------
  std::cout << "VSMLMAnalysis::analyze(): Fitting source model" << std::endl;
  analyze(d,data,acceptance,o.src_data);

  // Generate sky map of the best fit source model ----------------------------
  m_data_model->setParam(m_source_fit_data.param);
  for(VSSimple2DHist<double,double>::iterator itr = 
	o.sky_source_hist.begin(); itr != 
	o.sky_source_hist.end(); ++itr)
    {
      if(o.sky_livetime_hist.countForIndex(itr->bin()) == 0)
	continue;

      VSAAlgebra::Vec2D xy(itr->x(),itr->y());
      o.sky_source_hist.setBin(itr->bin(),m_data_model->val(xy));
    }

  // Perform fits on map ------------------------------------------------------
  VSAcceptanceModelFixed* bkgnd_fn = 
    new VSAcceptanceModelFixed(m_bkgnd_model,data.sky_counts_hist);
  m_data_model->setBkgndFn(bkgnd_fn);
  m_data_model->setSourceFn(m_source_model);

  // Instantiate the fitter ---------------------------------------------------
  m_fitter = 
    VSAMath::NLFitterFactory::createLnL<VSDataModel,VSModelCoord>
    (m_data_model);

  std::cout << "VSMLMAnalysis::analyze(): Generating sky maps" << std::endl;

  // Loop over sky map bins ---------------------------------------------------
  const unsigned nbin = data.sky_counts_hist.nBins();
  for(unsigned ibin = 0; ibin < nbin; ibin++)
    {
      if(m_fit_skymap && ibin % 1000 == 0) 
	std::cout << "VSMLMAnalysis::analyze(): BIN " 
		  << std::setw(12) << ibin 
		  << std::setw(12) << nbin << std::endl;

      if(o.sky_livetime_hist.countForIndex(ibin) == 0)
	continue;

      double x = 0, y = 0;
      data.sky_counts_hist.center(ibin,x,y);
      VSAAlgebra::Vec2D xy(x,y);	  
      //double livetime = o.sky_livetime_hist.countForIndex(ibin);
      double livetime = 0;

      const unsigned nptg = data.m_ptg_xy.size();
      for(unsigned iptg = 0; iptg < nptg; iptg++)
	{
	  if(xy.d(data.m_ptg_xy[iptg]) < 3*m_offset_max)
	    livetime += m_livetime[iptg];
	}

      if(m_fit_skymap && ibin) 
	{
	  fit(m_fit_data,xy,m_fit_radius,
	      false,m_bkgnd_fit_data,m_source_fit_data);
	  double lambda = m_bkgnd_fit_data.chi2-m_source_fit_data.chi2;
	  double sigma = sqrt(lambda);

	  if(m_source_fit_data.source_param(0) < 0) sigma *= -1;

	  o.sky_significance_hist.setBin(ibin,sigma);
	  o.significance_hist.accumulate(sigma);	  
	  if(!m_exclusion_region.isExcluded(xy))
	    o.significance_excluded_hist.accumulate(sigma);

	  double sum_counts, sum_counts_err;
	  m_source_model->integrate(m_source_fit_data.source_param,
				    m_source_fit_data.source_cov,
				    sum_counts,sum_counts_err);


	  double flux = m_source_fit_data.source_param(0);
	  double flux_err = sqrt(m_source_fit_data.source_cov(0,0));
	  double flux_ul95 = VSAStatistics::heleneUL(flux,flux_err,0.95);

	  o.sky_excess_hist.setBin(ibin,sum_counts);
	  o.sky_excess_rate_hist.setBin(ibin,sum_counts/livetime);
	  o.sky_flux_hist.setBin(ibin,flux);
	  o.sky_flux_ul95_hist.setBin(ibin,flux_ul95);
	}

      double bkgnd = m_bkgnd_model->val(xy);

      o.sky_bkgnd_hist.setBin(ibin,bkgnd);
      o.sky_bkgnd_rate_hist.setBin(ibin,bkgnd/livetime);
      o.sky_bkgnd_density_rate_hist.setBin(ibin,bkgnd/(livetime*m_domega));
    }

  delete m_fitter;
  m_fitter = NULL;
}

void VSMLMAnalysis::analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			    const VSAnalysisStage3Data& data,
			    const VSAcceptanceData& acceptance,
			    VSIntegralAnalysisDatum& o)
{
  o.method = "mlm";

  m_bkgnd_model->setAcceptanceParam(acceptance.acceptance_param);

  if(!loadFitData(data,acceptance)) return;

  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
   	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      o.livetime_min += itr->livetime_min;
      o.elaptime_min += itr->elaptime_min;      
    }

  // ==========================================================================
  // Obtain effective area for integral flux calculation
  // ==========================================================================
  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      double offset = itr->src_xy().d(itr->obs_xy());
      VSNSpace v = itr->irf_calc().getEffectiveArea(offset,m_fit_radius);
      v *= itr->livetime_min;

      if(o.effarea.size() == 0) o.effarea = v;
      else o.effarea += v;

      VSNSpace vpsf = itr->irf_calc().getPSF(m_spectrum_model,offset);
      vpsf *= itr->livetime_min;

      if(o.psf.size() == 0) o.psf = vpsf;
      else o.psf += vpsf;
    }

  o.effarea *= (1./o.livetime_min);
  o.psf *= (1./o.livetime_min);

  // Fit the model ------------------------------------------------------------
  std::cout << "VSMLMAnalysis::analyze(): Fitting source model" << std::endl;
  fitSource(data,o);

  o.drde = o.effarea;
  o.edrde = o.effarea;

  VSNSpace::Cell c(1);
  for(c.i[0] = 0; c.i[0] < o.effarea.axis(0).nbin; c.i[0]++)
    {
      double loge = o.effarea.axis(0).midCoordUnchecked(c.i[0]);
      o.drde[c] = o.effarea[c]*m_spectrum_model->val(loge);
      o.edrde[c] = o.drde[c]*std::pow(10,loge);
    }

  for(c.i[0] = 1; c.i[0] < o.drde.axis(0).nbin; c.i[0]++)
    {
      
      VSNSpace::Cell c2 = c;
      c2.i[0]++;

      if(o.drde[c2] < o.drde[c])
	{
	  VSAMath::Data<double> data;
	  
	  unsigned ilo = (unsigned)std::max(0,(int)c.i[0]-2);
	  unsigned ihi = 
	    (unsigned)std::min((int)o.drde.axis(0).nbin-1,(int)c.i[0]+2);

	  for(c2.i[0] = ilo; c2.i[0] <= ihi; c2.i[0]++)
	    {
	      double loge = o.drde.axis(0).midCoordUnchecked(c2.i[0]);
	      data.insert(VSAMath::DataPoint<double>(loge,o.drde[c2],1));

	      if(o.drde[c2] == 0) break;
	    }

	  VSAAlgebra::MatrixND cov;
	  VSAAlgebra::VecND param;

	  try
	    {
	      VSAMath::PolyFit::fit(2,data,param,&cov);
	    }
	  catch(const std::string& s)
	    {
	      std::cerr << s << std::endl;
	    }

	  o.egy_threshold = std::pow(10,-0.5*param[1]/param[2]);

	  break;
	}

    }

  calcFlux(o);

  o.drde *= o.dfde_ul95;
  o.edrde *= o.dfde_ul95;

  // Fill bkgnd distribution for Theta-Squared --------------------------------
  o.th2_bkgnd_counts_model_hist.fill(0);
  o.th2_source_counts_model_hist.fill(0);
  o.th2_source_model_hist.fill(0);
  o.th2_counts_model_hist.fill(0);

  double s_signal_sum = 0;

  for(VSLimitedErrorsHist<double,double>::iterator itr =
	o.th2_counts_model_hist.begin(); 
      itr != o.th2_counts_model_hist.end(); ++itr)
    {
      m_data_model->setParam(m_source_fit_data.param);
      VSAAlgebra::Vec2D xy = m_data_model->getSourceModel()->getXY();

      double R1 = sqrt(itr->val());
      double R2 = sqrt(itr->val()+o.th2_counts_model_hist.binSize());

      double s_signal = m_data_model->getSourceModel()->integrate(R1,R2);
      double s_bkgnd = m_data_model->getBkgndModel()->integrate(xy,R1,R2);

      s_signal_sum += s_signal;

      o.th2_bkgnd_counts_model_hist.setBin(itr->bin(),s_bkgnd,0);
      o.th2_source_counts_model_hist.setBin(itr->bin(),s_signal,0);
      o.th2_source_model_hist.setBin(itr->bin(),s_signal,0);
      o.th2_counts_model_hist.setBin(itr->bin(),s_signal+s_bkgnd,0);
    }

  o.th2_source_model_hist *= 
    std::pow(s_signal_sum*o.th2_counts_model_hist.binSize(),-1);
}

bool VSMLMAnalysis::loadFitData(const VSAnalysisStage3Data& data,
				const VSAcceptanceData& acceptance)
{
  m_fit_data.clear();

  const unsigned nptg = data.m_ptg_xy.size();

  std::map< unsigned, VSNSpace > ptg_effarea;
  std::map< unsigned, VSNSpace > ptg_psf;
  std::map< unsigned, VSSimple2DHist<double,double> > sky_counts_hists;
  const std::vector< VSAAlgebra::Vec2D >& ptg_xy = data.m_ptg_xy;

  m_livetime.clear();
  m_livetime.resize(nptg);

  // --------------------------------------------------------------------------
  // Create dataset and generate effarea/psf models for each pointing
  // --------------------------------------------------------------------------
  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      unsigned iptg = 0;
      for(iptg = 0; iptg < nptg; iptg++)
	if(ptg_xy[iptg] == itr->ptg_xy()) break;

      vsassert(iptg != nptg);

      VSNSpace effarea = itr->int_effarea_nspace;
      VSNSpace psf = itr->int_effarea_psf_nspace;

      effarea *= itr->livetime_min*60;
      psf *= itr->livetime_min*60;

      sky_counts_hists[iptg].merge(itr->sky_counts_hist);
      m_livetime[iptg] += itr->livetime_min;

      if(ptg_effarea[iptg].size() == 0)
	{
	  ptg_effarea[iptg] = effarea;
	  ptg_psf[iptg] = psf;
	}
      else
	{
	  ptg_effarea[iptg] += effarea;
	  ptg_psf[iptg] += psf;
	}      
    }

  // Construct the data model -------------------------------------------------
  m_bkgnd_model->clear();
  m_source_model->clear();

  std::vector< double > ncount(nptg);
  double zcount = 0;
  for(unsigned iptg = 0; iptg < nptg; iptg++)
    {
      m_source_model->setObs(ptg_xy[iptg],
			     ptg_effarea[iptg],ptg_psf[iptg]);
      m_bkgnd_model->setObs(ptg_xy[iptg]);

      double zsum = 0;
      for(VSSimple2DHist<double,double>::iterator hitr = 
	    sky_counts_hists[iptg].begin(); hitr !=
	    sky_counts_hists[iptg].end(); ++hitr)
	{
	  VSAAlgebra::Vec2D xy(hitr->x(),hitr->y());
	  if(ptg_xy[iptg].d(xy) > m_offset_max) continue;
	  zsum += hitr->count();
	}

      ncount[iptg] = zsum;
      if(zsum == 0) continue;

      // Create the dataset ---------------------------------------------------
      for(VSSimple2DHist<double,double>::iterator hitr = 
	    sky_counts_hists[iptg].begin(); hitr !=
	    sky_counts_hists[iptg].end(); ++hitr)
	{
	  VSAAlgebra::Vec2D xy(hitr->x(),hitr->y());
	  if(ptg_xy[iptg].d(xy) > m_offset_max) continue;

	  double z = hitr->count();
	 
	  zcount += z;
	  VSModelCoord c(hitr->x(),hitr->y(),iptg);
	  VSAMath::DataPoint<VSModelCoord> p(c,z,sqrt(z));
	  m_fit_data.insert(p);
	}
    }

  // Construct signal + bkgnd model -------------------------------------------
  const unsigned np_acceptance = acceptance.acceptance_param.ndim();
  for(unsigned ip = 0; ip < np_acceptance; ip++)
    m_bkgnd_model->setParam(ip,acceptance.acceptance_param(ip));

  for(unsigned iptg = 0; iptg < nptg; iptg++)
    m_bkgnd_model->setParam(np_acceptance+iptg,ncount[iptg]);

  m_data_model = new VSDataModel(m_source_model,m_bkgnd_model);

  if(zcount == 0) return false;
  else return true;
}

void VSMLMAnalysis::fitSource(const VSAnalysisStage3Data& data,
			      VSIntegralAnalysisDatum& o)
{
  o.src_xy = m_src_xy;
  o.src_dec_J2000 = m_src_radec.latitudeRad();
  o.src_ra_J2000 = m_src_radec.longitudeRad();  

  m_fitter = 
    VSAMath::NLFitterFactory::createLnL<VSDataModel,VSModelCoord>
    (m_data_model);

  fit(m_fit_data,o.src_xy,0,true,m_bkgnd_fit_data,m_source_fit_data);

  delete m_fitter;
  m_fitter = NULL;

  for(unsigned ip = 0; ip < m_source_fit_data.param.ndim(); ip++)
    std::cout << std::setw(15) << m_source_fit_data.param(ip)
	      << std::setw(15) << sqrt(m_source_fit_data.cov(ip,ip))
	      << std::endl;

  double sigma = sqrt(m_bkgnd_fit_data.chi2-m_source_fit_data.chi2);

  if(m_source_fit_data.source_param(0) < 0) sigma *= -1;

  o.significance = sigma;
  o.sigma_sqrthr = sigma/sqrt(o.elaptime_min/60.);

  double excess, excess_err;
  m_source_model->integrate(m_source_fit_data.source_param,
			    m_source_fit_data.source_cov,
			    excess,excess_err);

  // Calculate background per unit solid angle --------------------------------
  //  m_data_model->getSourceModel()->setAmplitude(0);

  m_bkgnd_model->setParam(m_source_fit_data.bkgnd_param);
  double bkgnd, bkgnd_err;
  m_bkgnd_model->val(o.src_xy,m_source_fit_data.bkgnd_param,
		     m_source_fit_data.bkgnd_cov,bkgnd,bkgnd_err);

  o.bkgnd_density = bkgnd/m_domega;
  o.bkgnd_density_err = bkgnd_err/m_domega;
  o.bkgnd_density_rate = o.bkgnd_density/o.livetime_min;  
  o.bkgnd_density_rate_err = o.bkgnd_density_err/o.livetime_min;
  
  // Calculate background within aperture -------------------------------------
  o.bkgnd = m_bkgnd_model->integrate(o.src_xy,0,m_theta_cut);
  o.bkgnd_rate = o.bkgnd/o.livetime_min;

  o.excess = excess;
  o.excess_err = excess_err;
  o.excess_rate = excess/o.livetime_min;
  o.excess_rate_err = excess_err/o.livetime_min;
  o.excess_rate_ul95 = VSAStatistics::heleneUL(o.excess_rate,
					       o.excess_rate_err,0.95);
  o.excess_rate_ul99 = VSAStatistics::heleneUL(o.excess_rate,
					       o.excess_rate_err,0.99);

  o.dfde = m_source_fit_data.source_param(0);
  o.dfde_err = sqrt(m_source_fit_data.source_cov(0,0));
}

void VSMLMAnalysis::excludeData(const VSAMath::Data<VSModelCoord>& data,
				const VSAAlgebra::Vec2D& xy, double R,
				VSAMath::Data<VSModelCoord>& d)
{
  for(VSAMath::Data<VSModelCoord>::const_iterator itr = data.begin();
      itr != data.end(); ++itr)
    {
      double r = sqrt(std::pow(xy.x()-itr->x.x(),2)+
 		      std::pow(xy.y()-itr->x.y(),2));
      
      VSAAlgebra::Vec2D xy(itr->x.x(),itr->x.y());

      if(r > R && R != 0) continue;
      else if(m_exclusion_region.isExcluded(xy) && 
	      r > m_fit_radius) continue;
      else d.insert(*itr);
    }
}

void VSMLMAnalysis::excludeData(const VSAMath::Data<VSModelCoord>& data,
				VSAMath::Data<VSModelCoord>& d)
{
  for(VSAMath::Data<VSModelCoord>::const_iterator itr = data.begin();
      itr != data.end(); ++itr)
    {
      VSAAlgebra::Vec2D xy(itr->x.x(),itr->x.y());
      if(m_exclusion_region.isExcluded(xy)) continue;
      else d.insert(*itr);
    }
}
			       
void VSMLMAnalysis::fit(const VSAMath::Data<VSModelCoord>& data,
			const VSAAlgebra::Vec2D& xy,
			double R, bool fit_bkgnd,
			FitResults& bkgnd_fit_data,
			FitResults& source_fit_data)
			   
{
  VSAMath::Data<VSModelCoord> d;
  excludeData(data,xy,R,d);

  if(d.size() == 0) return;

  // --------------------------------------------------------------------------
  // Fit Null Hypothesis 
  // --------------------------------------------------------------------------
  m_fitter->fn().setData(d);
  const unsigned nbkgnd_param = m_data_model->nBkgndParm();
  const unsigned nsource_param = m_data_model->nSourceParm();

  m_fitter->setScaling(nbkgnd_param,1E6);
  m_fitter->setMinTolerance(0.01);
  m_fitter->setMaxTolerance(0.1);

  m_data_model->fixBkgndParam(false);
  m_data_model->fixSourceParam(true);

  m_data_model->getSourceModel()->setXY(xy);
  m_data_model->getSourceModel()->setAmplitude(0);

  if(fit_bkgnd)
    {
      m_fitter->initialize(m_data_model->param(),m_data_model->fixed());
      m_fitter->fit();

      bkgnd_fit_data.chi2  = m_fitter->chi2();
      bkgnd_fit_data.param = m_fitter->param();
      bkgnd_fit_data.cov   = m_fitter->cov();
    }
  else
    {
      bkgnd_fit_data.chi2 = m_fitter->chi2(m_data_model->param());
      bkgnd_fit_data.param = m_data_model->param();
    }

  // --------------------------------------------------------------------------
  // Fit Source Hypothesis 
  // --------------------------------------------------------------------------

  // Hold background parameters -----------------------------------------------
  m_data_model->fixBkgndParam(true);
  m_data_model->setParam(bkgnd_fit_data.param);

//   const unsigned nparm = m_source_model_param.size();
//   for(unsigned ip = 0; ip < nparm; ip++)
//     {
//       m_data_model->getSourceModel()->
// 	setSourceParam(ip,m_source_model_param[ip].first,
// 		       m_source_model_param[ip].second);
//     }

  m_data_model->setParam(nbkgnd_param,0.0);
  m_data_model->fixParam(nbkgnd_param,false);
  m_fitter->initialize(m_data_model->param(),m_data_model->fixed());
  m_fitter->fit();

  m_data_model->setParam(m_fitter->param());

  // Free background parameters -----------------------------------------------
  if(fit_bkgnd)
    {
      m_data_model->fixBkgndParam(false);
      m_fitter->initialize(m_data_model->param(),m_data_model->fixed());
      m_fitter->fit();
    }

  source_fit_data.chi2  = m_fitter->chi2();
  source_fit_data.param = m_fitter->param();
  source_fit_data.cov   = m_fitter->cov();

  source_fit_data.bkgnd_param = source_fit_data.param.subVector(0,nbkgnd_param);
  source_fit_data.bkgnd_cov =    
    source_fit_data.cov.subMatrix(0,0,nbkgnd_param,nbkgnd_param);
  source_fit_data.source_param =
    source_fit_data.param.subVector(nbkgnd_param,nsource_param);
  source_fit_data.source_cov =    
    source_fit_data.cov.subMatrix(nbkgnd_param,nbkgnd_param,
				  nsource_param,nsource_param);
}

