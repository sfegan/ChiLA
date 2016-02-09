//-*-mode:c++; mode:font-lock;-*-

/*! \file VSIntegralAnalysis.hpp

  Base class for integral (energy-independent) analysis calculators.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.8 $
  \date       07/29/2007

  $Id: VSIntegralAnalysis.hpp,v 3.8 2010/10/20 04:02:07 matthew Exp $

*/

#ifndef VSINTEGRALANALYSIS_HPP
#define VSINTEGRALANALYSIS_HPP

#include <SphericalCoords.h>
#include <VSAAlgebra.hpp>
#include <VSEventData.hpp>
#include <VSIntegralAnalysisData.hpp>
#include <VSDataModel.hpp>
#include <VSBkgndModel.hpp>
#include <VSAnalysisStage3Data.hpp>

namespace VERITAS
{

  //! Base class for integral (energy-independent) stage3 analysis.
  //! The two analyze() methods accept a VSIntegralAnalysisData or
  //! VSIntegralAnalysisDatum structure which is filled with the
  //! results of the analysis.
  class VSIntegralAnalysis 
  {
  public:

    struct SourceData
    {
      SourceData():
	on_counts(), off_counts(), ring_counts(), src_xy(),
	off_xy(), iptg(), theta_cut()
      { }

      unsigned                             on_counts;
      unsigned                             off_counts;
      unsigned                             ring_counts;
      VSAAlgebra::Vec2D                    src_xy;
      std::vector< VSAAlgebra::Vec2D >     off_xy;
      unsigned                             iptg;
      double                               theta_cut;
    };

    struct Data
    {
      Data(): livetime_min(), elaptime_min(), obs_xy(), ptg_xy() { }

      double                         livetime_min;
      double                         elaptime_min;
      VSAAlgebra::Vec2D              obs_xy;
      VSAAlgebra::Vec2D              ptg_xy;

      VSNSpace                       effarea_nspace;
      VSNSpace                       effarea_psf_nspace;

      VSSimple2DHist<double,double>  sky_counts_hist;
      SourceData                     src_data;
    };

    VSIntegralAnalysis(const std::string bkgnd_model,
		       double bin_size_deg,
		       double theta_cut,
		       const std::pair< double, double >& ring_cut,
		       double offset_max,
		       unsigned max_nregion,
		       const std::string& spmodel);
  
    virtual ~VSIntegralAnalysis();

    VSAcceptanceData* fitAcceptance(VSAnalysisStage3Data& data);

    //! Generate results for every position in a two-dimensional field
    //! (sky maps).
    //! @param d Input vector of integral analysis data structure
    //! containing histograms of event parameters.
    //! @param o Output data structure with results.
    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisData& o) = 0;

    //! Generate results for a single source position.
    //! @param d Input vector of integral analysis data structure
    //! containing histograms of event parameters.
    //! @param o Output data structure with results.
    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisDatum& o) = 0;
    
    //! Method for filling parameters of each event into
    //! the VSIntegralAnalysis::Data structure.
    virtual void accumulate(const VSEventArrayDatum& event,
			    const VSAAlgebra::Vec2D& event_xy,
			    Data& data);

    virtual void accumulate(VSAnalysisStage3Data::RunData& d,
			    Data& data);

    void calcFlux(VSIntegralAnalysisDatum& o);

    //! Create a new datum.
    Data create(const VSAAlgebra::Vec2D& ptg_xy,
		const VSAAlgebra::Vec2D& obs_xy);

    //! Set the source coordinate.
    void setSourcePosition(const SEphem::SphericalCoords& origin_radec,
			   const SEphem::SphericalCoords& src_radec);

    //! Define a region in RA/DEC to be excluded from the background
    //! modeling analysis.
    void setExclusionRegion(const VSExclusionRegion& exclusion_region)
    {
      m_exclusion_region = exclusion_region;
    }
          
    void setSpectrumModel(VSSpectrumFn* sp)
    {
      delete m_spectrum_model;
      m_spectrum_model = sp->clone();
    }

    double thetaCut() const { return m_theta_cut; }
    const std::pair< double, double >& ringCut() const { return m_ring_cut; }
    unsigned maxOffRegions() const { return m_max_off_regions; }
    const std::string& spmodel() const { return m_spmodel; }

  protected:

    double                               m_bin_size_deg;
    double                               m_domega;
    double                               m_theta_cut;
    std::pair< double, double >          m_ring_cut;
    double                               m_offset_max;
    unsigned                             m_max_off_regions;
    std::string                          m_spmodel;

    SEphem::SphericalCoords              m_origin_radec;
    SEphem::SphericalCoords              m_src_radec;
    VSAAlgebra::Vec2D                    m_src_xy;
    std::vector< VSAAlgebra::Vec2D >     m_off_xy;

    VSExclusionRegion                    m_exclusion_region;
    VSAcceptanceModel*                   m_bkgnd_model;
    VSSpectrumFn*                        m_spectrum_model;
  };  


  //! Factory class for creating instances of integral analysis
  //! calculator (VSIntegralAnalysis).
  class VSIntegralAnalysisFactory
  {
  public:
    struct Options
    {
      Options(); 

      std::string                 method;
      std::string                 bkgnd_model;
      double                      theta_cut;
      std::pair< double, double > ring_cut;
      unsigned                    max_noff_region;
      std::string                 integral_source_spectrum; 
      bool                        mlm_fit_skymap;
      std::string                 mlm_source_model;
      std::string                 mlm_ext_model;
      double                      mlm_skymap_fit_radius;
    };

    virtual ~VSIntegralAnalysisFactory();

    static VSIntegralAnalysisFactory* getInstance();

    VSIntegralAnalysis* create(double max_offset_cut, double bin_width_deg);

    void setOptions(const Options& opt) { m_options = opt; }

    const Options& options() const { return m_options; }

    //! Configure the factory to create a given type of integral
    //! analysis calculator with the options defined on the
    //! command-line.
    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="s3_");

  private:
    VSIntegralAnalysisFactory(const Options& opt = s_default_options);

    Options                     m_options;

    static Options              s_default_options;

    static std::auto_ptr<VSIntegralAnalysisFactory> s_instance;
  };

};

#endif // VSINTEGRALANALYSIS_HPP
