//-*-mode:c++; mode:font-lock;-*-

/*! \file VSInstrumentResponseCalc.hpp

  General purpose class for generating detector response functions.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSInstrumentResponseCalc.hpp,v 3.3 2010/10/20 04:50:57 matthew Exp $

*/

#ifndef VSINSTRUMENTRESPONSECALC_HPP
#define VSINSTRUMENTRESPONSECALC_HPP

#include <VSNSpace.hpp>
#include <VSSimEnergyWeightCalc.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSLTLibraryCalc.hpp>

namespace VERITAS
{
  //! Parameterized angular resolution function.
  class VSPSFFn : public VSAFunction::ParamFn<double>
  {
  public:
    VSPSFFn();

    // Function Evalulation -------------------------------------------------
    double val(const double& x) const;
    double val(const double& x, const VSAAlgebra::VecND& a) const;
    
    void dyda(const double& x, VSAAlgebra::VecND& dyda) const; 
    void dyda(const double& x, const VSAAlgebra::VecND& a, 
	      VSAAlgebra::VecND& dyda) const; 
    
    double integrate(double xlo, double xhi);
    double integrate(const VSAAlgebra::VecND& a, double xlo, double xhi);

    // Virtual Constructor --------------------------------------------------
    VSPSFFn* clone() const { return new VSPSFFn(*this); }
  };

  //! Parameterized energy resolution function.
  class VSEnergyResponseFn : public VSAFunction::ParamFn<VSACoord::CoordND>
  {
  public:
    VSEnergyResponseFn(double dx = 1, double norm = 1);

    void setNorm(double norm) { m_norm = norm; }

    // Function Evalulation -------------------------------------------------
    double val(const VSACoord::CoordND& x) const;
    double val(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a) const;
    double val(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a, 
	       double norm) const;
      
    void dyda(const VSACoord::CoordND& x, VSAAlgebra::VecND& dyda) const; 
    void dyda(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a, 
	      VSAAlgebra::VecND& dyda) const; 
      
    // Virtual Constructor --------------------------------------------------
    VSEnergyResponseFn* clone() const 
    { return new VSEnergyResponseFn(*this); }
      
    double norm() const { return m_norm; }
      
  private:
    VSAFunction::Gauss1D<VSACoord::Coord1D> m_fn;
      
    double m_norm;
    double m_dx;
  };


  class VSEnergyKernelEffareaFn 
  {
  public:

    VSEnergyKernelEffareaFn(const VSNSpace&  effarea,
			    const VSNSpace&  krn_sigma1,
			    const VSNSpace&  krn_bias1,
			    const VSNSpace&  krn_sigma2,
			    const VSNSpace&  krn_bias2,
			    const VSNSpace&  krn_alpha,
			    VSSpectrumFn* sp_fn,
			    double offset,
			    double theta_cut);

    double val(const VSACoord::Coord2D& x) const;
      
  private:

    VSNSpace                m_effarea;
    VSNSpace                m_krn_sigma1;
    VSNSpace                m_krn_bias1;
    VSNSpace                m_krn_sigma2;
    VSNSpace                m_krn_bias2;
    VSNSpace                m_krn_alpha;
    VSEnergyResponseFn      m_dpdloge_fn;
    VSSpectrumFn*           m_sp_fn;
    double                  m_offset;
    double                  m_theta_cut;
  };

  //! Helper class for generating effective area, PSF, and energy
  //! response detector models from simulation lookup tables.
  class VSInstrumentResponseCalc : public VSLTLibraryCalc
  {
  public:

    struct Options
    {
      Options(): irf_lookup_file("") { }

      std::string              irf_lookup_file;
    };

    VSInstrumentResponseCalc(const Options& opt = s_default_options);

    // VSInstrumentResponseCalc(const VSNSpace& effarea,
    // 			     const VSNSpace& psf_sigma1,
    // 			     const VSNSpace& psf_sigma2,
    // 			     const VSNSpace& psf_alpha,
    // 			     double theta_max,
    // 			     const Options& opt = s_default_options);
    

    virtual ~VSInstrumentResponseCalc();

    // ========================================================================
    // Effective Area Functions
    // ========================================================================

    //! Returns total point-source effective area as a function of
    //! gamma-ray energy and offset angle.
    //!
    //! @param spectrum_fn Spectral energy distribution
    //! @param effarea Returned 1D nspace with total effective area
    //! vs. offset angle
    //! @param effarea Returned 2D nspace of total effective area
    //! vs. offset angle and reconstructed direction angle
    void getEffectiveAreaOffset(const std::string& spectrum_fn,
				VSNSpace& effarea, 
				VSNSpace& effarea_psf) const;

    void getEffectiveAreaOffset(const VSSpectrumFn* sp_fn,
				VSNSpace& effarea, 
				VSNSpace& effarea_psf) const;

    //! Returns an nspace containing point-source effective area
    //! within a circular aperture as a function of gamma-ray energy.
    void getEffectiveAreaAperture(const std::string& spectrum_fn,
				  double theta_cut, VSNSpace& effarea) const;

    //! Returns the point-source effective area within a circular
    //! aperture as a function of gamma-ray energy and offset angle.
    //!
    //! @param spectrum_fn Spectral energy distribution of the source
    //! model.
    double getEffectiveAreaAperture(const std::string& spectrum_fn,
				    double offset, double theta_cut) const;


    //! Returns the point-source effective area as a function of MC
    //! energy.
    VSNSpace getEffectiveArea(double ebin, double elo, double ehi,
			      double offset, double theta_cut) const;

    //! Returns the point-source effective area as a function of MC
    //! energy with default energy binning.
    VSNSpace getEffectiveArea(double offset, double theta_cut) const;

    double getEffectiveArea(double offset, double theta_cut, double egy) const;

    // ========================================================================
    // PSF Functions
    // ========================================================================

    VSNSpace getPSF(const VSSpectrumFn* spectrum_fn, double offset) const;

    // ========================================================================
    // Energy Kernel Functions
    // ========================================================================

    VSNSpace getKernel(double emc_ebin, double emc_elo, double emc_ehi,
		       double erec_ebin, double erec_elo, double erec_ehi,
		       double offset, double theta_cut) const;

    VSNSpace getKernel(double ebin, double elo, double ehi,
		       double offset, double theta_cut) const;

    VSNSpace getKernel(double ebin, double elo, double ehi,
		       double offset) const;

    double getKernel(double offset, double emc, double erec) const;

    // Load/Save --------------------------------------------------------------
    virtual void clear() { }

    virtual bool load(const std::string& lookup_file,
		      const std::vector<unsigned>& nchan, 
		      double zn_deg, double az_deg, double ped_rms);

    bool load(const std::vector<unsigned>& nchan, double zn_deg, 
	      double az_deg, double ped_rms);

    bool load(const VSTime& date,
	      const std::vector<unsigned>& nchan, double zn_deg, 
	      double az_deg, double ped_rms);

    void save(VSScaledParameterLibraryWriter& writer,
	      double zenith_deg,
	      double azimuth_deg,
	      double ped_dev) const;

    void load(VSScaledParameterLibraryReader* reader,
	      double zenith_deg,
	      double azimuth_deg,
	      double ped_dev);    

    bool load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="s3_");

  private:

    Options                        m_options;

    VSNSpace m_effarea_nspace;
    VSNSpace m_psf_sigma1_nspace;
    VSNSpace m_psf_sigma2_nspace;
    VSNSpace m_psf_alpha_nspace;
    VSNSpace m_krn_sigma1_nspace;
    VSNSpace m_krn_bias1_nspace;
    VSNSpace m_krn_sigma2_nspace;
    VSNSpace m_krn_bias2_nspace;
    VSNSpace m_krn_alpha_nspace;

    double   m_theta_max;

    static Options                 s_default_options;
  };
}

#endif // VSINSTRUMENTRESPONSECALC_HPP
