//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMuonAnalysis.hpp

  Try to fit a ring to an image, cut on the image parameters, do the muon
  analysis and return the muon paramaters

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/14/2007

  $Id: VSMuonAnalysis.hpp,v 3.1 2007/04/13 23:32:41 sfegan Exp $

*/

#define MUON_TEST_NON_UNIFORMITY

#ifndef VSMUONANALYSIS_HPP
#define VSMUONANALYSIS_HPP

#include<vector>

#include<VSRingFitter.hpp>
#include<VSMuonAnalysisData.hpp>

namespace VERITAS
{

  class VSMuonAnalysis
  {
  public:
    VSMuonAnalysis(const VSAReconstruction::ScopeInfo& scope,
		   double raw_ring_width,
		   unsigned nimage_cut,
		   double radius_min_cut, double radius_max_cut,
		   double rms_max_cut, 
		   double ring_edge_dist_max_cut,
		   double centroid_radius_ratio_max_cut
#ifdef MUON_TEST_NON_UNIFORMITY
		   , double nu_beta = 0
#endif		   
		   ):
      m_fitter(new VSRingFitter(scope)), m_scope(scope),
      m_raw_ring_width(raw_ring_width),
      m_nimage_cut(nimage_cut), 
      m_radius_min_cut(radius_min_cut), m_radius_max_cut(radius_max_cut),
      m_rms_max_cut(rms_max_cut), 
      m_ring_edge_dist_max_cut(ring_edge_dist_max_cut),
      m_centroid_radius_ratio_max_cut(centroid_radius_ratio_max_cut)
#ifdef MUON_TEST_NON_UNIFORMITY
      , m_nu_beta(nu_beta>-0.5?nu_beta:-0.5)
#endif
    {
      // nothing to see here
    }

    ~VSMuonAnalysis()
    {
      delete m_fitter;
    }

    class RawDataFetcher
    {
    public:
      virtual ~RawDataFetcher();
      virtual bool fetchRawData(unsigned npixel, 
				double* raw_data, double* raw_total_data) = 0;
    };

    bool analyze(VSMuonAnalysisDatum& datum,
		 unsigned ievent,
		 const VSAReconstruction::ScopeImage& image,
		 RawDataFetcher* raw_image_fetcher = 0) const;

  private:
    const VSRingFitter*                  m_fitter;
    const VSAReconstruction::ScopeInfo&  m_scope;
    const double                         m_raw_ring_width;
    const unsigned                       m_nimage_cut;
    const double                         m_radius_min_cut;
    const double                         m_radius_max_cut;
    const double                         m_rms_max_cut;
    const double                         m_ring_edge_dist_max_cut;
    const double                         m_centroid_radius_ratio_max_cut;
#ifdef MUON_TEST_NON_UNIFORMITY
    const double                         m_nu_beta;
#endif
  };

}

#endif // VSMUONANALYSIS_HPP
