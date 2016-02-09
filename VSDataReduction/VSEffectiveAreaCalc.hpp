//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEffectiveAreaCalc.hpp

  Calculate effective areas

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSEffectiveAreaCalc.hpp,v 3.9 2010/06/25 01:53:30 matthew Exp $

*/

#ifndef VSEFFECTIVEAREACALC_HPP
#define VSEFFECTIVEAREACALC_HPP

#include <VSNSpace.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSSimulationData.hpp>
#include <VSAFunction.hpp>
#include <VSInstrumentResponseCalc.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSEffectiveAreaCalc
  // ==========================================================================
  class VSEffectiveAreaCalc
  {
  public:

    class PolyFn : public VSAFunction::ParamFn<double>
    {
    public:
      PolyFn(unsigned n1, unsigned n2, double x0);

      void operator() (const double& x, VSAAlgebra::VecND& v) const;
      double val(const double& x) const;
      double val(const double& x, const VSAAlgebra::VecND& a) const;
      void dyda(const double& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const double& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& dyda) const;

      // Virtual Constructor --------------------------------------------------
      PolyFn* clone() const 
      { return new PolyFn(*this); }

      double            m_x0;
      unsigned          m_n;
      VSAFunction::Poly m_p1;
      VSAFunction::Poly m_p2;
    };

    struct Data
      {
      public:
	Data(double offset_deg, double ebin, double elo, double ehi): 
	  offset_deg(offset_deg),
	  egymc_selected_hist(ebin,elo,ehi), 
	  egymc_total_hist(ebin,elo,ehi), 
	  egymc_area_hist(ebin,elo,ehi),
	  egymc_effarea_hist(ebin,elo,ehi), 
	  egymc_log_effarea_hist(ebin,elo,ehi), 
	  egymc_effarea_fit_hist(ebin,elo,ehi),
	  //	  egymc_effarea_fit2_hist(ebin,elo,ehi),
	  egymc_log_effarea_fit_hist(ebin,elo,ehi),
	  effarea_fit_param(), 
	  effarea_fit_cov(),
	  log_effarea_fit_param(), 
	  log_effarea_fit_cov()
	{ }

	double                             offset_deg;	
	VSLimitedErrorsHist<double,double> egymc_selected_hist;
	VSLimitedErrorsHist<double,double> egymc_total_hist;
	VSLimitedErrorsHist<double,double> egymc_area_hist;
	VSLimitedErrorsHist<double,double> egymc_effarea_hist;
	VSLimitedErrorsHist<double,double> egymc_log_effarea_hist;
	VSLimitedErrorsHist<double,double> egymc_effarea_fit_hist;
	//	VSLimitedErrorsHist<double,double> egymc_effarea_fit2_hist;
	VSLimitedErrorsHist<double,double> egymc_log_effarea_fit_hist;
	VSAAlgebra::VecND                  effarea_fit_param;
	VSAAlgebra::MatrixND               effarea_fit_cov;
	VSAAlgebra::VecND                  log_effarea_fit_param;
	VSAAlgebra::MatrixND               log_effarea_fit_cov;

	void save(VSOctaveH5WriterStruct* writer) const;
      };


    VSEffectiveAreaCalc();
    ~VSEffectiveAreaCalc();

    void loadSimInfo(VSSimInfoData* sim_info);
    void loadHeader(const VSHeaderSimulationDatum& sim_header);

    void setEnergyBinning(double ebin, double elo, double ehi);

    void accumulate(double energy_tev);

    void calcEffarea();
    void calcEffarea(Data* data);
    void fit();
    void fit(Data* data);
    void save(VSOctaveH5WriterStruct* writer) const;

    const VSNSpace& getEffarea() const { return m_effarea_offset_nspace; }

  private:

    VSSimple2DHist<double,double>   m_effarea_offset_hist;
    VSNSpace                        m_effarea_offset_nspace;
    std::vector< VSLimitedErrorsHist<double,double> > m_effarea_param_hist;

    double                          m_log10_egybin;
    double                          m_log10_egylo;
    double                          m_log10_egyhi;
    unsigned                        m_negy;

    std::vector<Data*>              m_data;
    Data*                           m_data_ptr;
  };

  // ==========================================================================
  // VSEnergyKernelFn
  // ==========================================================================
  class VSPolyFn : public VSAFunction::ParamFn<double>
  {
  public:
    VSPolyFn(unsigned n1, unsigned n2, double x0);

    void operator() (const double& x, VSAAlgebra::VecND& v) const;
    double val(const double& x) const;
    double val(const double& x, const VSAAlgebra::VecND& a) const;
    void dyda(const double& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const double& x, const VSAAlgebra::VecND& a, 
	      VSAAlgebra::VecND& dyda) const;

    // Virtual Constructor --------------------------------------------------
    VSPolyFn* clone() const 
    { return new VSPolyFn(*this); }

    double            m_x0;
    unsigned          m_n;
    VSAFunction::Poly m_p1;
    VSAFunction::Poly m_p2;
  };

  class VSEnergyKernelFn : public VSAFunction::ParamFn<VSACoord::CoordND>
  {
  public:


    VSEnergyKernelFn(double dx, unsigned negy,
		     double dloge, double logemin);

    ~VSEnergyKernelFn();

    void setNorm(unsigned iegy, double counts);

    // Function Evalulation -------------------------------------------------
    double val(const VSACoord::CoordND& x) const;
    double val(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a) const;
    double val(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a, 
	       double counts) const;

    void dyda(const VSACoord::CoordND& x, VSAAlgebra::VecND& dyda) const; 
    void dyda(const VSACoord::CoordND& x, const VSAAlgebra::VecND& a, 
	      VSAAlgebra::VecND& dyda) const; 

    // Virtual Constructor --------------------------------------------------
    VSEnergyKernelFn* clone() const { return new VSEnergyKernelFn(*this); }

    
    // Copy-Constructor -----------------------------------------------------
    VSEnergyKernelFn(const VSEnergyKernelFn& o);

    void coeff(double x, 
	       const VSAAlgebra::VecND& a,
	       VSAAlgebra::VecND& afn) const;

    void coeff_var(double x, 
		   const VSAAlgebra::VecND& a,
		   const VSAAlgebra::MatrixND& cov,
		   VSAAlgebra::VecND& var) const;

    double norm(unsigned iegy) const { return m_fn[iegy].norm(); }

  private:
    
    // Assignment Operator --------------------------------------------------
    VSEnergyKernelFn& operator= (const VSEnergyKernelFn& o);

    std::vector< VSAFunction::BaseFn<double>* >  m_pol;
    std::vector<VSEnergyResponseFn>              m_fn;
    VSBinCalcLinear<double>                      m_ebin;
  };

  // ==========================================================================
  // VSEnergyKernelCalc
  // ==========================================================================
  class VSEnergyKernelCalc
  {
  public:
    struct Data
      {
      public:
	Data(double offset_deg, double ebin, double elo, double ehi): 
	  offset_deg(offset_deg),
	  egymc_egy_hist(ebin,elo,ehi,0.02,-1.5,2.5),
	  egymc_egyerr_hist(ebin,elo,ehi,0.02,-1.5,2.5),
	  egymc_biasdev_hist(ebin,elo,ehi),
	  egymc_mse_hist(ebin,elo,ehi),
	  egymc_bias_hist(ebin,elo,ehi),
	  egymc_dev_hist(ebin,elo,ehi),
	  egymc_res_hist(ebin,elo,ehi),
	  kernel_hist(ebin,elo,ehi,0.02,-1.5,2.5),
	  kernel_fit_hist(ebin,elo,ehi,0.02,-1.5,2.5),
	  kernel_fit_counts_hist(ebin,elo,ehi,0.02,-1.5,2.5),
	  egy_hist(),
	  egy_fit_hist(),
	  egy_fit2_hist(),
	  param_hist(5,VSLimitedErrorsHist<double,double>(ebin,elo,ehi)),
	  param_fit_hist(5,VSLimitedErrorsHist<double,double>(ebin,elo,ehi))
	{ }

	double                              offset_deg;	
	VSSimple2DHist<double,double>       egymc_egy_hist;
	VSSimple2DHist<double, double>      egymc_egyerr_hist;
	VSLimitedErrorsHist<double,double>  egymc_biasdev_hist;
	VSLimitedErrorsHist<double,double>  egymc_mse_hist;
	VSLimitedErrorsHist<double,double>  egymc_bias_hist;
	VSLimitedErrorsHist<double,double>  egymc_dev_hist;
	VSLimitedErrorsHist<double,double>  egymc_res_hist;
	VSSimple2DHist<double,double>       kernel_hist;
	VSSimple2DHist<double,double>       kernel_fit_hist;
	VSSimple2DHist<double,double>       kernel_fit_counts_hist;
	VSNSpace                            kernel_nspace;

	std::vector< VSLimitedErrorsHist<double,double> > egy_hist;
	std::vector< VSLimitedErrorsHist<double,double> > egy_fit_hist;
	std::vector< VSLimitedErrorsHist<double,double> > egy_fit2_hist;

	std::vector<VSLimitedErrorsHist<double,double> >  param_hist;
	std::vector<VSLimitedErrorsHist<double,double> >  param_fit_hist;

	void save(VSOctaveH5WriterStruct* writer) const;
      };

    VSEnergyKernelCalc();
    ~VSEnergyKernelCalc();

    void loadSimInfo(VSSimInfoData* sim_info);
    void loadHeader(const VSHeaderSimulationDatum& sim_header);

    void setEnergyBinning(double ebin, double elo, double ehi);

    void accumulate(double emc_tev, double erec_tev);

    void calcKernel();
    void calcKernel(Data* data);

    void fit();
    void fit(Data* data);

    void save(VSOctaveH5WriterStruct* writer) const;

    const VSNSpace& getSigma1() const { return m_sigma1_offset_nspace; }
    const VSNSpace& getBias1() const { return m_bias1_offset_nspace; }
    const VSNSpace& getSigma2() const { return m_sigma2_offset_nspace; }
    const VSNSpace& getBias2() const { return m_bias2_offset_nspace; }
    const VSNSpace& getAlpha() const { return m_alpha_offset_nspace; }

  private:

    VSSimple2DHist<double,double>   m_sigma1_offset_hist;
    VSSimple2DHist<double,double>   m_bias1_offset_hist;
    VSSimple2DHist<double,double>   m_sigma2_offset_hist;
    VSSimple2DHist<double,double>   m_bias2_offset_hist;
    VSSimple2DHist<double,double>   m_alpha_offset_hist;

    VSNSpace                        m_sigma1_offset_nspace;
    VSNSpace                        m_bias1_offset_nspace;
    VSNSpace                        m_sigma2_offset_nspace;
    VSNSpace                        m_bias2_offset_nspace;
    VSNSpace                        m_alpha_offset_nspace;

    double                          m_log10_egybin;
    double                          m_log10_egylo;
    double                          m_log10_egyhi;
    unsigned                        m_negy;

    std::vector<Data*>              m_data;
    Data*                           m_data_ptr;
  };

  

} // namespace VERITAS

#endif // VSEFFECTIVEAREACALC_HPP
