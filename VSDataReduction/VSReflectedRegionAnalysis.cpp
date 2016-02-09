//-*-mode:c++; mode:font-lock;-*-

/*! \file VSReflectedRegionAnalysis.cpp

  Integral analysis with reflected region algorithm.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.5 $
  \date       07/29/2007

  $Id: VSReflectedRegionAnalysis.cpp,v 3.5 2010/10/20 03:10:26 matthew Exp $

*/

#include <VSReflectedRegionAnalysis.hpp>

using namespace VERITAS;


VSReflectedRegionAnalysis::
VSReflectedRegionAnalysis(const std::string bkgnd_model,
			  double bin_size_deg,
			  double theta_cut, 
			  const std::pair< double, double >& ring_cut,
			  double theta_max,
			  unsigned max_nregion,
			  const std::string& spmodel):
  VSApertureAnalysis(bkgnd_model,bin_size_deg,theta_cut,ring_cut,
		     theta_max,max_nregion,spmodel)
{
  
}

VSReflectedRegionAnalysis::~VSReflectedRegionAnalysis()
{

}

void VSReflectedRegionAnalysis::
analyze(const std::vector<VSIntegralAnalysis::Data>& d,
	const VSAnalysisStage3Data& data,
	const VSAcceptanceData& acceptance,
	VSIntegralAnalysisData& o)
{
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

  VSApertureAnalysis::analyze(d,data,acceptance,o);
  analyze(d,data,acceptance,o.src_data);
}

void VSReflectedRegionAnalysis::
analyze(const std::vector<VSIntegralAnalysis::Data>& d,
	const VSAnalysisStage3Data& data,
	const VSAcceptanceData& acceptance,
	VSIntegralAnalysisDatum& o)
{
  // o->alpha = 1./(double)data.off_xy().size();
  // o->on_counts = data.on_counts;
  // o->off_counts = data.off_counts;

  VSApertureAnalysis::analyze(d,data,acceptance,o);
}

void VSReflectedRegionAnalysis::
accumulate(unsigned iptg,
	   const VSAnalysisStage3Data::RunData& data)
{
  // Calculate results for source position ------------------------------------
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

  VSSimple2DHist<double,double>& on_hist = m_sky_on_counts_hist[iptg];
  VSSimple2DHist<double,double>& off_hist = m_sky_off_counts_hist[iptg];
  VSSimple2DHist<double,double>& domega_hist = m_sky_domega_hist[iptg];
  VSSimple2DHist<double,double>& livetime_hist = m_sky_livetime_hist[iptg];
  VSSimple2DHist<double,double>& alpha_hist = m_sky_alpha_hist[iptg];

  // Calculate results from 2D maps -------------------------------------------
  for(VSSimple2DHist<double, double>::iterator itr = 
	data.sky_counts_hist.begin(); 
      itr != data.sky_counts_hist.end(); ++itr)
    {
      VSAAlgebra::Vec2D xy(itr->x(),itr->y());

      if(xy.d(data.ptg_xy()) < 2*m_theta_cut || 
	 xy.d(data.ptg_xy()) > m_offset_max)
	continue;

      livetime_hist.accumulateBin(itr->bin(),data.livetime_min);

      std::vector< unsigned > on_bins;
      getSignalBins(data.sky_counts_hist,xy,data.ptg_xy(),on_bins);

      const unsigned non_bins = on_bins.size();
      for(unsigned ibin = 0; ibin < non_bins; ibin++)
	on_hist.accumulateBin(itr->bin(),
			      data.sky_counts_hist.
			      countForIndex(on_bins[ibin]));

      // Accumulate solid angle of signal aperture ----------------------------
      domega_hist.setBin(itr->bin(),m_domega*on_bins.size());

      std::vector< unsigned > off_bins;

      std::vector< VSAAlgebra::Vec2D > off_coords;
      VSOffRegion::getRegions(xy,data.ptg_xy(),m_theta_cut, 
			      maxOffRegions(), off_coords);

      for(std::vector< VSAAlgebra::Vec2D >::iterator itr2 = off_coords.begin();
	  itr2 != off_coords.end(); ++itr2)
	getBkgndBins(data.sky_counts_hist,*itr2,data.ptg_xy(),off_bins);

      double alpha = (double)on_bins.size()/(double)off_bins.size();
      const unsigned noff_bins = off_bins.size();
      for(unsigned ibin = 0; ibin < noff_bins; ibin++)
	off_hist.accumulateBin(itr->bin(),
			       data.sky_counts_hist.
			       countForIndex(off_bins[ibin]));

      alpha_hist.accumulateBin(itr->bin(), alpha);
    }
}

void VSReflectedRegionAnalysis::
getBkgndBins(const VSSimple2DHist<double,double>& h,
	     const VSAAlgebra::Vec2D& coord,
	     const VSAAlgebra::Vec2D& obs_xy,
	     std::vector< unsigned > &bins)
{    
  double rmax = m_theta_cut+h.xBinSize();

  int xlo = std::max((int)0,(int)h.xValToBin(coord.x()-rmax));
  int xhi = std::min((int)h.xValToBin(coord.x()+rmax),(int)h.nXBins()-1);
  int ylo = std::max((int)0,(int)h.yValToBin(coord.y()-rmax));
  int yhi = std::min((int)h.yValToBin(coord.y()+rmax),(int)h.nYBins()-1);

  for(int xbin = xlo; xbin <= xhi; xbin++)
    {
      for(int ybin = ylo; ybin <= yhi; ybin++)
	{
	  int bin = h.xyBinToIndex(xbin,ybin);   
	  VSAAlgebra::Vec2D xy(h.xBinToCenter(xbin),h.yBinToCenter(ybin));
	  double dbin = coord.d(xy);
	  double dobs = obs_xy.d(xy);

	  if(dobs < m_offset_max && !m_exclusion_region.isExcluded(xy) &&
	     dbin < m_theta_cut)
	    bins.push_back(bin);
	}
    }
}

void VSReflectedRegionAnalysis::
getSignalBins(const VSSimple2DHist<double,double>& h,
	      const VSAAlgebra::Vec2D& coord,
	      const VSAAlgebra::Vec2D& obs_xy,
	      std::vector< unsigned > &bins)
{    
  double rmax = m_theta_cut+h.xBinSize();

  int xlo = std::max((int)0,(int)h.xValToBin(coord.x()-rmax));
  int xhi = std::min((int)h.xValToBin(coord.x()+rmax),(int)h.nXBins()-1);
  int ylo = std::max((int)0,(int)h.yValToBin(coord.y()-rmax));
  int yhi = std::min((int)h.yValToBin(coord.y()+rmax),(int)h.nYBins()-1);

  for(int xbin = xlo; xbin <= xhi; xbin++)
    {
      for(int ybin = ylo; ybin <= yhi; ybin++)
	{
	  int bin = h.xyBinToIndex(xbin,ybin); 
	  VSAAlgebra::Vec2D xy(h.xBinToCenter(xbin),h.yBinToCenter(ybin));
	  double dobs = obs_xy.d(xy);
	  double dbin = coord.d(xy);
	  if(dobs < m_offset_max && dbin < m_theta_cut) bins.push_back(bin);
	}
    }
}

