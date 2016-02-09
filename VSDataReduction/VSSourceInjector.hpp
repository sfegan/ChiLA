//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSourceInjector.hpp

  Perform coordinate transforms on simulation events.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSSourceInjector.hpp,v 3.5 2009/09/20 01:46:25 matthew Exp $

*/

#ifndef VSSOURCEINJECTOR_HPP
#define VSSOURCEINJECTOR_HPP

#include <VSAAlgebra.hpp>
#include <SphericalCoords.h>
#include <VSOptions.hpp>
#include <VSSimple2DHist.hpp>
#include <RandomNumbers.hpp>
#include <VSAnalysisStage3Data.hpp>

namespace VERITAS
{

  class VSSourceInjector
  {
  public:
    class Data
    {
    public:

      Data(): src_rate(), integral_flux_hist(), diff_flux_hist() { }

      double                             src_rate;

      VSLimitedErrorsHist<double,double> integral_flux_hist;
      VSLimitedErrorsHist<double,double> diff_flux_hist;

      void save(VSOctaveH5Writer* writer) const
      {
	writer->writeCompositeHere(*this);
	integral_flux_hist.save(writer->writeStruct("integral_flux_hist"));
	diff_flux_hist.save(writer->writeStruct("diff_flux_hist"));
      }
      
      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,Data,src_rate);
      }
    };

    VSSourceInjector(RandomNumbers* rng,
		     double theta,
		     const std::pair< double, double >& ring_radius,
		     const std::string& src_type = s_default_src_type,
		     const std::string& src_spectrum = s_default_src_spectrum,
		     double src_rate = s_default_src_rate,
		     const std::pair<double,double>& src_xy =
		     s_default_src_xy);
    virtual ~VSSourceInjector();

    virtual void initialize(VSAnalysisStage3Data::RunData* data);
    virtual void fill(VSAnalysisStage3Data::RunData* data);
    virtual void fillGauss(unsigned nevent, 
			   VSAnalysisStage3Data::RunData* data);
    virtual void fillConstant(unsigned nevent, 
			      VSAnalysisStage3Data::RunData* data);

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

    Data getData(double ebin, double elo, double ehi) const;

    void fillIntegralHist(VSLimitedErrorsHist<double,double>& h) const;

  private:
    double                      m_theta;
    std::pair< double, double > m_ring_radius;

    std::string                 m_src_type;
    std::string                 m_src_spectrum;
    double                      m_src_rate;
    VSAAlgebra::Vec2D           m_src_xy;
    double                      m_emin;
    double                      m_emax;

    RandomNumbers*              m_rng;
    VSSpectrumFn*               m_spectrum_fn;

    RandomNumbers::CDFGenerator<VSAFunction::Spline> m_emc_rng;
    RandomNumbers::CDFGenerator<VSAFunction::Spline> m_erec_rng;
    VSNSpace m_krn;


    // Default options --------------------------------------------------------
    static std::string                 s_default_src_type;
    static std::string                 s_default_src_spectrum;
    static double                      s_default_src_rate;
    static unsigned                    s_default_rng_seed;
    static std::pair< double, double > s_default_src_xy;
  };
 
} // namespace VERITAS


#endif // VSSOURCEINJECTOR_HPP
