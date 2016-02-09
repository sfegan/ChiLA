//-*-mode:c++; mode:font-lock;-*-

/*! \file VSPSFCalc.hpp

  Model gamma-ray psf from simulations.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSPSFCalc.hpp,v 3.4 2010/05/31 00:30:11 matthew Exp $

*/

#ifndef VSPSFCALC_HPP
#define VSPSFCALC_HPP

#include <VSNSpace.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSSimulationData.hpp>
#include <VSAFunction.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSBinnedPSFFn
  // ==========================================================================
  class VSBinnedPSFFn : public VSAFunction::ParamFn<double>
  {
  public:
    VSBinnedPSFFn(double dthsq);
    VSBinnedPSFFn(double dthsq, double counts);

    // Function Evalulation -------------------------------------------------
    double val(const double& x) const;
    double val(const double& x, const VSAAlgebra::VecND& a) const;
    
    void dyda(const double& x, VSAAlgebra::VecND& dyda) const; 
    void dyda(const double& x, const VSAAlgebra::VecND& a, 
	      VSAAlgebra::VecND& dyda) const; 
    
    double integrate(const VSAAlgebra::VecND& a, double xlo, double xhi);

    void setCounts(double counts) { m_counts = counts; }

    // Virtual Constructor --------------------------------------------------
    VSBinnedPSFFn* clone() const { return new VSBinnedPSFFn(*this); }

  private:
    double m_counts;
    double m_dthsq;
  };

  // ==========================================================================
  // VSPSFEnergyFn
  // ==========================================================================
  class VSPSFEnergyFn : public VSAFunction::ParamFn<VSACoord::CoordND>
  {
  public:
    VSPSFEnergyFn(double dthsq, unsigned negy, double dloge, double logemin);

    void setCounts(unsigned iegy, double counts);

    // Function Evalulation -------------------------------------------------
    double val(const VSACoord::CoordND& x) const;
    double val(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a) const;
    
    void dyda(const VSACoord::CoordND& x, VSAAlgebra::VecND& dyda) const; 
    void dyda(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a, 
	      VSAAlgebra::VecND& dyda) const; 
    
    // Virtual Constructor --------------------------------------------------
    VSPSFEnergyFn* clone() const { return new VSPSFEnergyFn(*this); }

  private:
    double                        m_dthsq;
    std::vector<VSBinnedPSFFn>    m_fn;
    VSBinCalcLinear<double>       m_ebin;
  };


  // ==========================================================================
  // VSPSFCalc
  // ==========================================================================
  class VSPSFCalc
  {
  public:
    struct Data
      {
      public:

	struct EnergyData
	{
	public:

	  EnergyData(double log10_egy):
	    log10_egy_tev(log10_egy),
	    thetasq_hist(0.0002,0.,0.08),
	    thetasq_fit_hist(0.0002,0.,0.08)
	  { }

	  double                      log10_egy_tev;

	  VSLimitedErrorsHist<double> thetasq_hist;
	  VSLimitedErrorsHist<double> thetasq_fit_hist;
	  VSLimitedErrorsHist<double> thetasq_fit2_hist;

	  void save(VSOctaveH5WriterStruct* writer) const;
	};


	Data(double offset_deg, double ebin, double elo, double ehi): 
	  offset_deg(offset_deg),
	  egymc_hist(ebin,elo,ehi),
	  egymc_thetasq_hist(ebin,elo,ehi,0.0002,0.,0.08),
	  sigma1_hist(ebin,elo,ehi),
	  sigma2_hist(ebin,elo,ehi),
	  alpha_hist(ebin,elo,ehi),
	  sigma1_fit_hist(ebin,elo,ehi),
	  sigma2_fit_hist(ebin,elo,ehi),
	  alpha_fit_hist(ebin,elo,ehi),
	  th68_fit_hist(ebin,elo,ehi),
	  m_data()
	{ }

	double                             offset_deg;
	VSLimitedErrorsHist<double>        egymc_hist;
	VSSimple2DHist<double,double>      egymc_thetasq_hist;
	VSLimitedErrorsHist<double>        sigma1_hist;
	VSLimitedErrorsHist<double>        sigma2_hist;
	VSLimitedErrorsHist<double>        alpha_hist;

	VSLimitedErrorsHist<double>        sigma1_fit_hist;
	VSLimitedErrorsHist<double>        sigma2_fit_hist;
	VSLimitedErrorsHist<double>        alpha_fit_hist;
	VSLimitedErrorsHist<double>        th68_fit_hist;

	std::vector<EnergyData>            m_data;

	void save(VSOctaveH5WriterStruct* writer) const;
      };


    VSPSFCalc();
    ~VSPSFCalc();

    void loadSimInfo(VSSimInfoData* sim_info);
    void loadHeader(const VSHeaderSimulationDatum& sim_header);

    void setEnergyBinning(double ebin, double elo, double ehi);

    void accumulate(double energy_tev, double thetasq);

    void fit();
    void fit(Data *data);

    void save(VSOctaveH5WriterStruct* writer) const;

    const VSNSpace& getSigma1() const { return m_sigma1_offset_nspace; }
    const VSNSpace& getSigma2() const { return m_sigma2_offset_nspace; }
    const VSNSpace& getAlpha() const { return m_alpha_offset_nspace; }

  private:

    VSSimple2DHist<double,double>   m_sigma1_offset_hist;
    VSSimple2DHist<double,double>   m_sigma2_offset_hist;
    VSSimple2DHist<double,double>   m_alpha_offset_hist;
    VSSimple2DHist<double,double>   m_th68_offset_hist;

    VSNSpace                        m_sigma1_offset_nspace;
    VSNSpace                        m_sigma2_offset_nspace;
    VSNSpace                        m_alpha_offset_nspace;

    double                          m_log10_egybin;
    double                          m_log10_egylo;
    double                          m_log10_egyhi;
    unsigned                        m_negy;

    std::vector<Data*>              m_data;
    Data*                           m_data_ptr;
  };

} // namespace VERITAS

#endif // VSPSFCALC_HPP
