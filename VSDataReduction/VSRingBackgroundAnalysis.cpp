//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRingBackgroundAnalysis.cpp

  Integral analysis with ring background algorithm.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.8 $
  \date       07/29/2007

  $Id: VSRingBackgroundAnalysis.cpp,v 3.8 2010/10/20 03:10:06 matthew Exp $

*/

#include <Astro.h>
#include <VSAStatistics.hpp>
#include <VSAAlgebra.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSRingBackgroundAnalysis.hpp>

using namespace VERITAS;


VSRingBackgroundAnalysis::
VSRingBackgroundAnalysis(const std::string bkgnd_model,
			 double bin_size_deg,
			 double theta_cut, 
			 const std::pair< double, double >& ring_cut,
			 double max_offset_cut,
			 unsigned max_nregion,
			 const std::string& spmodel):
  VSApertureAnalysis(bkgnd_model,bin_size_deg,theta_cut, ring_cut,
		     max_offset_cut, max_nregion, spmodel),
  m_xbinner(bin_size_deg,0), m_ybinner(bin_size_deg,0)
{
  int xlo = m_xbinner.valToBin(-ringCut().second)-1;
  int ylo = m_ybinner.valToBin(-ringCut().second)-1;
  int xhi = m_xbinner.valToBin(ringCut().second)+1;
  int yhi = m_ybinner.valToBin(ringCut().second)+1;
  
  for(int xbin = xlo; xbin <= xhi; xbin++)
    {
      for(int ybin = ylo; ybin <= yhi; ybin++)
	{
	  double bin_x = bin_size_deg*(double)xbin;
	  double bin_y = bin_size_deg*(double)ybin;
	  double dr = sqrt(bin_x*bin_x + bin_y*bin_y);

	  if(dr > ringCut().first && dr < ringCut().second)
	    m_bkgnd_offsets.push_back(std::make_pair(xbin,ybin));
	  else if(dr < m_theta_cut)
	    m_signal_offsets.push_back(std::make_pair(xbin,ybin));       
	}
    }
}

VSRingBackgroundAnalysis::~VSRingBackgroundAnalysis()
{

}

void VSRingBackgroundAnalysis::
analyze(const std::vector<VSIntegralAnalysis::Data>& d,
	const VSAnalysisStage3Data& data,
	const VSAcceptanceData& acceptance,
	VSIntegralAnalysisData& o)
{
  o.src_data.method = "rbm";

  m_ptg_xy = data.m_ptg_xy;
  m_livetime.resize(m_ptg_xy.size());

  m_bkgnd_model->setAcceptanceParam(acceptance.acceptance_param);

  const unsigned nrun = data.nrun();
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      unsigned iptg = 0;
      const unsigned nptg = m_ptg_xy.size();
      for(iptg = 0; iptg < nptg; iptg++)
	if(m_ptg_xy[iptg] == data.run_data(irun).ptg_xy()) break;

      vsassert(iptg != nptg);

      accumulate(iptg,data.run_data(irun));
    }

  // Perform analysis on sky maps ---------------------------------------------
  VSApertureAnalysis::analyze(d,data,acceptance,o);

  // Perform analysis on source location --------------------------------------
  analyze(d,data,acceptance,o.src_data);
}

void VSRingBackgroundAnalysis::
analyze(const std::vector<VSIntegralAnalysis::Data>& d,
	const VSAnalysisStage3Data& data,
	const VSAcceptanceData& acceptance,
	VSIntegralAnalysisDatum& o)
{
  o.src_xy = m_src_xy;
  o.src_dec_J2000 = m_src_radec.latitudeRad();
  o.src_ra_J2000 = m_src_radec.longitudeRad();

  std::vector< double > on_counts;
  std::vector< double > off_counts;
  std::vector< double > alpha;

  m_bkgnd_model->setAcceptanceParam(acceptance.acceptance_param);

  for(std::vector< VSAnalysisStage3Data::RunData >::const_iterator itr =
   	data.run_data().begin(); itr != data.run_data().end(); ++itr)
    {
      double acceptance_signal = 
	m_bkgnd_model->integrate(itr->src_xy(),0,m_theta_cut);
      double acceptance_bkgnd = 
	m_bkgnd_model->integrate(itr->src_xy(),ringCut().first,
				 ringCut().second);
      double a = acceptance_signal/acceptance_bkgnd;

      double on = itr->on_counts;
      double off = itr->ring_counts;

      on_counts.push_back(on);
      off_counts.push_back(off);
      alpha.push_back(a);

      double bkgnd_var = std::pow(a,2)*off;
      double excess_var = on + std::pow(a,2)*off;

      o.on_counts += itr->on_counts;
      o.on_counts_err = sqrt(o.on_counts);
      o.off_counts += itr->ring_counts;
      o.off_counts_err = sqrt(o.off_counts);
      o.bkgnd += a*off;
      o.bkgnd_err = sqrt(std::pow(o.bkgnd_err,2) + bkgnd_var);

      o.excess += on - a*off;
      o.excess_err = sqrt(std::pow(o.excess_err,2) + excess_var);
    }

  o.significance = 
    VSAStatistics::limaSignificance(on_counts,off_counts,alpha);

  VSApertureAnalysis::analyze(d,data,acceptance,o);
}

void VSRingBackgroundAnalysis::
accumulate(unsigned iptg,
	   const VSAnalysisStage3Data::RunData& data)
{
  // Calculate results for source position ------------------------------------
  m_acceptance_hist = data.sky_counts_hist;
  m_acceptance_hist.clear();

  for(VSSimple2DHist<double, double>::iterator itr = m_acceptance_hist.begin();
      itr != m_acceptance_hist.end(); ++itr)
    {
      VSAAlgebra::Vec2D xy(itr->x(),itr->y());
      xy -= data.ptg_xy();
      m_acceptance_hist.setBin(itr->bin(),m_bkgnd_model->getFoVAcceptance(xy));
    }

  if(m_sky_on_counts_hist.size() <= iptg)
    {
      double xlo = data.sky_counts_hist.xLoLimit();
      double xhi = data.sky_counts_hist.xHiLimit();
      double xbin = data.sky_counts_hist.xBinSize();
      double ylo = data.sky_counts_hist.yLoLimit();
      double yhi = data.sky_counts_hist.yHiLimit();
      double ybin = data.sky_counts_hist.yBinSize();

      VSSimple2DHist<double,double> h(xbin,xlo,xhi,ybin,ylo,yhi);

      m_sky_on_counts_hist.resize(iptg+1,h);
      m_sky_off_counts_hist.resize(iptg+1,h);
      m_sky_domega_hist.resize(iptg+1,h);
      m_sky_livetime_hist.resize(iptg+1,h);
      m_sky_alpha_hist.resize(iptg+1,h);
      m_sky_effarea_hist.resize(iptg+1,h);
    } 

  m_livetime[iptg] += data.livetime_min;

  VSSimple2DHist<double,double>& on_hist = m_sky_on_counts_hist[iptg];
  VSSimple2DHist<double,double>& off_hist = m_sky_off_counts_hist[iptg];
  VSSimple2DHist<double,double>& domega_hist = m_sky_domega_hist[iptg];
  VSSimple2DHist<double,double>& livetime_hist = m_sky_livetime_hist[iptg];
  VSSimple2DHist<double,double>& alpha_hist = m_sky_alpha_hist[iptg];
  VSSimple2DHist<double,double>& effarea_hist = m_sky_effarea_hist[iptg];

  // Calculate results from 2D maps -------------------------------------------
  for(VSSimple2DHist<double, double>::iterator itr = 
	data.sky_counts_hist.begin(); 
      itr != data.sky_counts_hist.end(); ++itr)
    {
      VSAAlgebra::Vec2D xy(itr->x(),itr->y());

      if(xy.d(data.ptg_xy()) < m_offset_max)  
	{     
	  livetime_hist.accumulateBin(itr->bin(),data.livetime_min);
	}

      // Accumulate counts in the signal histogram ----------------------------
      double effarea = 0;
      std::vector< unsigned > on_bins;
      for(std::vector< std::pair< int, int > >::iterator offset_itr =
	    m_signal_offsets.begin(); offset_itr != m_signal_offsets.end(); 
	  ++offset_itr)
	{  
	  int xbin = data.sky_counts_hist.indexToXBin(itr->bin()) + 
	    offset_itr->first;
	  int ybin = data.sky_counts_hist.indexToYBin(itr->bin()) +
	    offset_itr->second;

	  VSAAlgebra::Vec2D xy2(data.sky_counts_hist.xBinToCenter(xbin),
			       data.sky_counts_hist.yBinToCenter(ybin));

	  double dr_obs = data.ptg_xy().d(xy2);
	  if(dr_obs > m_offset_max) continue;

	  int bin = data.sky_counts_hist.xyBinToIndex(xbin,ybin);
	  if(bin < (int)on_hist.nBins())
	    {
	      on_hist.accumulateBin(itr->bin(),
				    data.sky_counts_hist.countForIndex(bin));
	      on_bins.push_back(bin);

	      VSNSpace::Point p(2);

	      p.x[0] = dr_obs;
	      p.x[1] = xy.d(xy2);

	      double w = 0;
	      data.int_effarea_psf_nspace.interpolateWeight(p,w);
	      effarea += w*m_domega;
	    }
	}

      // Accumulate solid angle of signal aperture ----------------------------
      domega_hist.setBin(itr->bin(),m_domega*on_bins.size());


      // Accumulate counts in the background histogram ------------------------
      std::vector< unsigned > off_bins;      
      for(std::vector< std::pair< int, int > >::iterator offset_itr =
	    m_bkgnd_offsets.begin(); offset_itr != m_bkgnd_offsets.end();
	  ++offset_itr)
	{
	  int xbin = data.sky_counts_hist.indexToXBin(itr->bin()) + 
	    offset_itr->first;
	  int ybin = data.sky_counts_hist.indexToYBin(itr->bin()) +
	    offset_itr->second;

	  VSAAlgebra::Vec2D xy2(data.sky_counts_hist.xBinToCenter(xbin),
			       data.sky_counts_hist.yBinToCenter(ybin));

	  double dr_obs = data.ptg_xy().d(xy2);
	  
	  if(dr_obs > m_offset_max || m_exclusion_region.isExcluded(xy2))
	    continue;

	  int bin = data.sky_counts_hist.xyBinToIndex(xbin,ybin);

	  if(bin < (int)off_hist.nBins())
	    {
	      off_hist.accumulateBin(itr->bin(),
				     data.sky_counts_hist.countForIndex(bin));
	      off_bins.push_back(bin);
	    }
	}

      // Calculate the acceptance correction for this bin ---------------------
      double alpha = 
 	integrateAcceptance(m_acceptance_hist,on_bins)/
 	integrateAcceptance(m_acceptance_hist,off_bins);

      alpha_hist.setBin(itr->bin(), alpha);
      effarea_hist.setBin(itr->bin(), effarea*data.livetime_min*60);
    }
}

double VSRingBackgroundAnalysis::
integrateAcceptance(const VSSimple2DHist<double,double>& h,
		    const std::vector<unsigned>& bins)
{
  double sum = 0;
  for(unsigned ibin = 0; ibin < bins.size(); ibin++)
    {
      sum += h.countForIndex(bins[ibin]);
      //      std::cout << h.countForIndex(bins[ibin]) << std::endl;
    }
  return sum;
}

