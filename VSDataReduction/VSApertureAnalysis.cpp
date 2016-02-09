//-*-mode:c++; mode:font-lock;-*-

/*! \file VSApertureAnalysis.cpp

  Base call for aperture-based integral analysis calcultors (ring
  background and reflected region).

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.9 $
  \date       07/29/2007

  $Id: VSApertureAnalysis.cpp,v 3.9 2010/10/20 03:30:40 matthew Exp $

*/

#include <fstream>
#include <Astro.h>
#include <VSAStatistics.hpp>
#include <VSApertureAnalysis.hpp>
#include <VSOctaveH5Writer.hpp>

using namespace VERITAS;

VSApertureAnalysis::
VSApertureAnalysis(const std::string bkgnd_model,
		   double bin_size_deg,
		   double theta_cut,
		   const std::pair< double, double >& ring_cut,
		   double theta_max,
		   unsigned max_nregion,
		   const std::string& spmodel):
  VSIntegralAnalysis(bkgnd_model,bin_size_deg,theta_cut, ring_cut,
		     theta_max,max_nregion,spmodel),
  m_sky_on_counts_hist(),
  m_sky_off_counts_hist(),
  m_sky_domega_hist(),
  m_sky_livetime_hist()
{

}

VSApertureAnalysis::~VSApertureAnalysis()
{

}

void VSApertureAnalysis::analyze(const std::vector<VSIntegralAnalysis::Data>& d,
				 const VSAnalysisStage3Data& data,
				 const VSAcceptanceData& acceptance,
				 VSIntegralAnalysisData& o)
{
  const unsigned nptg = m_sky_on_counts_hist.size();
  for(unsigned iptg = 0; iptg < nptg; iptg++)
    {
      o.sky_on_counts_hist.merge(m_sky_on_counts_hist[iptg]);
      o.sky_off_counts_hist.merge(m_sky_off_counts_hist[iptg]);
      o.sky_livetime_hist.merge(m_sky_livetime_hist[iptg]);
    }

  o.sky_alpha_hist = o.sky_on_counts_hist;
  o.sky_significance_hist = o.sky_on_counts_hist;
  o.sky_bkgnd_hist = o.sky_on_counts_hist;
  o.sky_bkgnd_rate_hist = o.sky_on_counts_hist;
  o.sky_bkgnd_density_rate_hist = o.sky_on_counts_hist;
  o.sky_excess_hist = o.sky_on_counts_hist;
  o.sky_excess_rate_hist = o.sky_on_counts_hist;
  o.sky_flux_hist = o.sky_on_counts_hist;
  o.sky_flux_ul95_hist = o.sky_on_counts_hist;

  o.sky_alpha_hist.clear();
  o.sky_significance_hist.clear();
  o.sky_bkgnd_hist.clear();
  o.sky_bkgnd_rate_hist.clear();
  o.sky_bkgnd_density_rate_hist.clear();
  o.sky_excess_hist.clear();
  o.sky_excess_rate_hist.clear();
  o.sky_flux_hist.clear();
  o.sky_flux_ul95_hist.clear();



  for(VSSimple2DHist<double, double>::iterator hitr = 
	o.sky_on_counts_hist.begin(); 
      hitr != o.sky_on_counts_hist.end(); ++hitr)
    {  
      std::vector< double > on_counts;
      std::vector< double > off_counts;
      std::vector< double > alpha;

      int bin = hitr->bin();
      double x = hitr->x();
      double y = hitr->y();
      VSAAlgebra::Vec2D xy(hitr->x(),hitr->y());

      bool skip = true;
      for(unsigned iptg = 0; iptg < nptg; iptg++)
	if(xy.d(m_ptg_xy[iptg]) < m_offset_max)
	  {
	    skip = false;
	    break;
	  }

      if(skip) continue;

      double excess_sum = 0;
      double excess_var = 0;
      double bkgnd_sum = 0;
      double atot = 0;
      double livetime = 0;
      double effarea_sum = 0;

      const unsigned nptg = m_sky_on_counts_hist.size();
      for(unsigned iptg = 0; iptg < nptg; iptg++)
	{
	  if(xy.d(m_ptg_xy[iptg]) < 3*m_offset_max)
	    livetime += m_livetime[iptg];

	  double on = m_sky_on_counts_hist[iptg].countForVal(x,y);
	  double off = m_sky_off_counts_hist[iptg].countForVal(x,y);
	  double a = m_sky_alpha_hist[iptg].countForVal(x,y);
	  double domega = m_sky_domega_hist[iptg].countForVal(x,y);
	  double effarea = m_sky_effarea_hist[iptg].countForVal(x,y);

	  if(domega == 0) continue;

	  excess_sum += on-off*a;
	  excess_var += on + a*a*off;
	  effarea_sum += effarea;

	  atot += a;
	  //	  livetime += m_sky_livetime_hist[iptg].countForVal(x,y);

	  on_counts.push_back(on);
	  off_counts.push_back(off);
	  alpha.push_back(a);

	  o.sky_excess_hist.accumulateBin(bin,on-off*a);	  
	  o.sky_bkgnd_hist.accumulateBin(bin,off*a*m_domega/domega);
	  bkgnd_sum += off*a*m_domega/domega;
	}

      double sig = 
	VSAStatistics::limaSignificance(on_counts,off_counts,alpha);
      o.sky_significance_hist.setBin(bin,sig);
      o.significance_hist.accumulate(sig);
      if(!m_exclusion_region.isExcluded(xy))
	o.significance_excluded_hist.accumulate(sig);
      
      o.sky_bkgnd_rate_hist.setBin(bin,bkgnd_sum/livetime);
      o.sky_bkgnd_density_rate_hist.
	setBin(bin,bkgnd_sum/(livetime*m_domega));
      o.sky_excess_rate_hist.setBin(bin,excess_sum/livetime);
      o.sky_alpha_hist.setBin(bin,atot/alpha.size());
      o.sky_flux_hist.setBin(bin,excess_sum/effarea_sum);
      
      double excess_ul95 =
	VSAStatistics::heleneUL(excess_sum,sqrt(excess_var),0.95);
      
      o.sky_flux_ul95_hist.setBin(bin,excess_ul95/effarea_sum);
    }
}


void VSApertureAnalysis::analyze(const std::vector<VSIntegralAnalysis::Data>& d,
				 const VSAnalysisStage3Data& data,
				 const VSAcceptanceData& acceptance,
				 VSIntegralAnalysisDatum& o)
{
  double domega = M_PI*std::pow(m_theta_cut,2);
  
  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
   	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      o.livetime_min += itr->livetime_min;
      o.elaptime_min += itr->elaptime_min;
    }

  o.bkgnd_density = o.bkgnd/domega;
  o.bkgnd_density_err = o.bkgnd_err/domega;

  VSSpectrumFn* spectrum_fn = VSSpectrumFn::create(m_spmodel);
  spectrum_fn->setNormalization(1);

  // ==========================================================================
  // Obtain effective area for integral flux calculation
  // ==========================================================================
  double effarea_sum = 0;
  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      double offset = itr->src_xy().d(itr->obs_xy());
      
      double effarea = 
	itr->irf_calc().getEffectiveAreaAperture(m_spmodel,
						 offset,m_theta_cut);

      effarea_sum += effarea*itr->livetime_min*60;

      VSNSpace v = itr->irf_calc().getEffectiveArea(offset,m_theta_cut);
      v *= itr->livetime_min;

      if(o.effarea.size() == 0) o.effarea = v;
      else o.effarea += v;

      VSNSpace vpsf = itr->irf_calc().getPSF(spectrum_fn,offset);
      vpsf *= itr->livetime_min;

      if(o.psf.size() == 0) o.psf = vpsf;
      else o.psf += vpsf;
    }

  o.effarea *= (1./o.livetime_min);
  o.psf *= (1./o.livetime_min);

  delete spectrum_fn;

  // Calculate Rates ----------------------------------------------------------
  o.on_rate      = o.on_counts/o.livetime_min;
  o.on_rate_err  = o.on_counts_err/o.livetime_min;
  o.off_rate     = o.off_counts/o.livetime_min;
  o.off_rate_err = o.off_counts_err/o.livetime_min;
  o.bkgnd_rate      = o.bkgnd/o.livetime_min;
  o.bkgnd_rate_err  = o.bkgnd_err/o.livetime_min;
  o.bkgnd_density_rate = o.bkgnd_density/o.livetime_min;
  o.bkgnd_density_rate_err = o.bkgnd_density_err/o.livetime_min;
  o.excess_rate     = o.excess/o.livetime_min;
  o.excess_rate_err = o.excess_err/o.livetime_min;
  o.excess_rate_ul95 = 
    VSAStatistics::heleneUL(o.excess_rate,
			    o.excess_rate_err,0.95);
  o.excess_rate_ul99 = 
    VSAStatistics::heleneUL(o.excess_rate,
			    o.excess_rate_err,0.99);

  o.sigma_sqrthr = o.significance/sqrt(o.elaptime_min/60.);

  // Calculate flux -----------------------------------------------------------
  o.dfde = o.excess/effarea_sum;
  o.dfde_err = o.excess_err/effarea_sum;
  o.dfde_ul95 = VSAStatistics::heleneUL(o.dfde,o.dfde_err,0.95);
  o.dfde_ul99 = VSAStatistics::heleneUL(o.dfde,o.dfde_err,0.99);
}
