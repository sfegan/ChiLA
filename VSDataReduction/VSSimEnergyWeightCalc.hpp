//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimEnergyWeightCalc.hpp

  Class that reweights simulated energies.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Timothy C. Arlen            \n
              UCLA                        \n
              arlen@astro.ucla.edu        \n

  \version    1.0
  \date       03/09/2008

  $Id: VSSimEnergyWeightCalc.hpp,v 3.11 2010/10/20 17:51:19 matthew Exp $

*/

#ifndef VSSIMENERGYWEIGHTCALC_HPP
#define VSSIMENERGYWEIGHTCALC_HPP

#include<VSOptions.hpp>
#include<VSSimulationData.hpp>
#include<VSAFunction.hpp>

namespace VERITAS
{

  //! Spectrum base class.
  class VSSpectrumFn : public VSAFunction::ParamFn<double>
  {
  public:
    VSSpectrumFn(unsigned nparm, const std::string& name);
    virtual ~VSSpectrumFn() {}

    // Evaluation -------------------------------------------------------------

    //! Return the differential flux (dF/dE) at the given energy.  
    //! @param log10_egy_tev Energy at which the differential flux
    //! should be evaluated in log10(E/TeV).
    //! @param a Vector of free parameters of the spectral model.
    virtual double val(const double& log10_egy_tev, 
		       const VSAAlgebra::VecND& a) const = 0;
    virtual double val(const double& log10_egy_tev) const
    {
      return val(log10_egy_tev,param());
    }

    virtual void dyda(const double& log10_egy_tev, 
		      VSAAlgebra::VecND& dyda) const = 0;
    virtual void dyda(const double& log10_egy_tev, 
		      const VSAAlgebra::VecND& a,
		      VSAAlgebra::VecND& dyda) const = 0;

    virtual double diffFlux(double log10_egy_tev,
			    const VSAAlgebra::VecND& a) const
    {
      return val(log10_egy_tev,a);
    }

    virtual double diffFlux(double log10_egy_tev) const
    {
      return val(log10_egy_tev,param());
    }

    virtual double dfdlog10e(const double& log10_egy_tev,
			     const VSAAlgebra::VecND& a) const
    {
      return val(log10_egy_tev,a)*
	std::pow(10,log10_egy_tev)*std::log(10);
    }

    //! Return the integral of the spectrum between two energies.
    virtual double integralFlux(double log10_elo_tev,
				double log10_ehi_tev,
				const VSAAlgebra::VecND& a) const = 0;

    virtual double integralFlux(double log10_elo_tev,
				double log10_ehi_tev) const
    {
      return integralFlux(log10_elo_tev,log10_ehi_tev,param());
    }

    virtual void integralFluxDerivative(double log10_elo_tev,
					double log10_ehi_tev,
					VSAAlgebra::VecND& dyda) const = 0;

    const std::string& name() const { return m_name; }

    double normEnergy() const { return m_log10_norm_egy; }

    // Setters ----------------------------------------------------------------
    virtual void setNormalization(double norm) = 0;
    virtual void setNormEnergy(double egy) { m_log10_norm_egy = egy; }

    // Factory Methods --------------------------------------------------------
    static VSSpectrumFn* create(const std::string& sp);
    static VSSpectrumFn* create(const std::vector<std::string>& sp);
    static VSSpectrumFn* create(const std::string& sp_type,
				const VSAAlgebra::VecND& param);

    // Virtual Constructor ----------------------------------------------------
    virtual VSSpectrumFn* clone() const = 0;

  protected:

    std::string m_name;
    double      m_log10_norm_egy;
  };

  //! Two parameter spectral model corresponding to a power-law
  //! parameterized by a spectral index and normalization (dF/dE = F_0
  //! E^(-Lambda)).
  class VSSpectrumFnPowerLaw : public VSSpectrumFn
  {
  public:
    VSSpectrumFnPowerLaw(const std::vector<std::string>& args);

    VSSpectrumFnPowerLaw(double index = 2.5, 
			 double flux_constant = 3.2E-7);
    virtual ~VSSpectrumFnPowerLaw() {}

    // Evaluation -------------------------------------------------------------
    virtual double val(const double& log10_egy_tev, 
		       const VSAAlgebra::VecND& a) const;

    virtual void dyda(const double& log10_egy_tev, 
		      VSAAlgebra::VecND& dyda) const;
    virtual void dyda(const double& log10_egy_tev, 
		      const VSAAlgebra::VecND& a,
		      VSAAlgebra::VecND& dyda) const;

    virtual double integralFlux(double log10_elo_tev,
				double log10_ehi_tev,
				const VSAAlgebra::VecND& a) const;

    virtual void integralFluxDerivative(double log10_elo_tev,
					double log10_ehi_tev,
					VSAAlgebra::VecND& dyda) const;

    // Setters ----------------------------------------------------------------
    virtual void setNormalization(double norm)
    {
      setParam(0,norm);
      m_flux_constant = norm;
    }

    // Virtual Constructor ----------------------------------------------------
    virtual VSSpectrumFnPowerLaw* clone() const 
    { return new VSSpectrumFnPowerLaw(*this); }

  private:

    double      m_index;
    double      m_flux_constant;
    double      m_norm_egy;
  };

  //! Three parameter spectral model corresponding to power-law with
  //! an exponential cut-off (dF/dE).
  class VSSpectrumFnPowerLawExp : public VSSpectrumFn
  {
  public:
    VSSpectrumFnPowerLawExp(const std::vector<std::string>& args);
    VSSpectrumFnPowerLawExp(double index = 2.5, 
			    double log10_ecut_tev = 0, 
			    double flux_constant = 3.2E-7);
    virtual ~VSSpectrumFnPowerLawExp() {}

    // Evaluation -------------------------------------------------------------
    virtual double val(const double& log10_egy_tev, 
		       const VSAAlgebra::VecND& a) const;

    virtual void dyda(const double& log10_egy_tev, 
		      VSAAlgebra::VecND& dyda) const;
    virtual void dyda(const double& log10_egy_tev, 
		      const VSAAlgebra::VecND& a,
		      VSAAlgebra::VecND& dyda) const;

    virtual double integralFlux(double log10_elo_tev,
				double log10_ehi_tev,
				const VSAAlgebra::VecND& a) const;

    virtual void integralFluxDerivative(double log10_elo_tev,
					double log10_ehi_tev,
					VSAAlgebra::VecND& dyda) const;    

    // Setters ----------------------------------------------------------------
    virtual void setNormalization(double norm)
    {
      setParam(0,norm);
      m_flux_constant = norm;
    }

    // Virtual Constructor ----------------------------------------------------
    virtual VSSpectrumFnPowerLawExp* clone() const 
    { return new VSSpectrumFnPowerLawExp(*this); }

  private:
    double m_index;
    double m_log10_ecut_tev;
    double m_flux_constant;
  };

  // ==========================================================================
  // VSSpectrumFnPowerLawExpCutoff
  // ==========================================================================
  class VSSpectrumFnPowerLawExpCutoff : public VSSpectrumFn
  {
  public:
    VSSpectrumFnPowerLawExpCutoff(double index = -2.5, 
				  double energy_cutoff = 0.2, 
				  double alpha = 1.0,
				  double flux_constant = 3.2E-7);
    virtual ~VSSpectrumFnPowerLawExpCutoff() {}

    // Evaluation -------------------------------------------------------------
    virtual double val(const double& log10_egy_tev, 
		       const VSAAlgebra::VecND& a) const;

    virtual void dyda(const double& log10_egy_tev, 
		      VSAAlgebra::VecND& dyda) const;
    virtual void dyda(const double& log10_egy_tev, 
		      const VSAAlgebra::VecND& a,
		      VSAAlgebra::VecND& dyda) const;

    virtual double integralFlux(double log10_elo_tev,
				double log10_ehi_tev,
				const VSAAlgebra::VecND& a) const;

    virtual void integralFluxDerivative(double log10_elo_tev,
					double log10_ehi_tev,
					VSAAlgebra::VecND& dyda) const;   

    // Setters ----------------------------------------------------------------
    virtual void setNormalization(double norm) 
    {
      setParam(0,norm);
      m_flux_constant = norm;
    }

    // Virtual Constructor ----------------------------------------------------
    virtual VSSpectrumFnPowerLawExpCutoff* clone() const 
    { return new VSSpectrumFnPowerLawExpCutoff(*this); }

  private:
    double m_index;
    double m_energy_cutoff;
    double m_alpha;
    double m_flux_constant;
  };


  // ==========================================================================
  // VSSimEnergyWeightCalc
  // ==========================================================================
  class VSSimEnergyWeightCalc
  {
  public:
    struct Options
    {
      Options():
	sim_energy_weight()
      { } 

      std::vector< std::string > sim_energy_weight;
    };

    VSSimEnergyWeightCalc(const Options& opt = s_default_options);
    ~VSSimEnergyWeightCalc();

    void calcWeighting(const VSHeaderSimulationDatum& sim_header);

    double calcWeight(double egy_tev);
    VSSpectrumFn* getSpectrum();
    
    static void configure(VSOptions& options,
			  const std::string& profile="", 
			  const std::string& opt_prefix="");

    static Options& getDefaultOptions() { return s_default_options; }

  private:
    VSSpectrumFn*                  m_spectrum;
    std::map<double,double>        m_energy_weight;
    double                         m_energy_norm;
    bool                           m_apply_weighting;

    Options                        m_options;

    static Options                 s_default_options;
  };
}


#endif // VSSIMENERGYWEIGHTCALC_HPP
