//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSpectrumCalc.hpp

  Various classes for spectral reconstruction.

  \author     Tim Arlen                   \n
              UCLA                        \n
              arlen@astro.ucla.edu        \n

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       10/03/2008

  $Id: VSSpectrumCalc.hpp,v 3.17 2010/07/16 21:01:28 matthew Exp $

*/

#ifndef VSSPECTRUMCALC_HPP
#define VSSPECTRUMCALC_HPP

#include <VSSimpleErrorsHist.hpp>
#include <VSSimpleGraph.hpp>
#include <VSNSpace.hpp>
#include <VSAAlgebra.hpp>
#include <VSOptions.hpp>
#include <VSAnalysisStage3Data.hpp>

namespace VERITAS
{
  class VSSpectrumFitData
  {
  public:
    VSSpectrumFitData();
    ~VSSpectrumFitData() { }

    std::string                        model;
    double                             log10_enorm;

    // Spectral fit parameters and errors -------------------------------------
    VSAAlgebra::VecND                  param;
    VSAAlgebra::VecND                  param_err;
    VSAAlgebra::MatrixND               param_cov;

    // Differential Fluxes ----------------------------------------------------
    double                             dfde;
    double                             dfde_err;
    double                             e2dfde;
    double                             e2dfde_err;
    double                             dfde100;
    double                             dfde100_err;
    double                             dfde200;
    double                             dfde200_err;
    double                             dfde316;
    double                             dfde316_err;
    double                             dfde1000;
    double                             dfde1000_err;

    // Integral fluxes --------------------------------------------------------
    double                             flux100;
    double                             flux200;
    double                             flux316;
    double                             flux1000;

    // Goodness of fit --------------------------------------------------------
    double                             chi2_pearson;
    double                             chi2_pearson_pval;
    double                             chi2_ml;
    double                             chi2_ml_pval;
    unsigned                           ndf;

    VSNSpace                           kernel_nspace;
    VSLimitedErrorsHist<double,double> on_hist;
    VSLimitedErrorsHist<double,double> off_hist;
    VSLimitedErrorsHist<double,double> mu_on_hist;
    VSLimitedErrorsHist<double,double> mu_off_hist;
    VSLimitedErrorsHist<double,double> mu_bkgnd_hist;

    VSSimple2DHist<double,double>               lhist;
    std::vector< VSSimpleGraph<double,double> > lcont68;
    std::vector< VSSimpleGraph<double,double> > lcont90;

    VSSimpleGraph<double,double>                dfde_butt68;
    VSSimpleGraph<double,double>                e2dfde_butt68;

    // Load/Save --------------------------------------------------------------
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

  };

  class VSSpectrumData
  {
  public:
    VSSpectrumData();
    ~VSSpectrumData() { }

    double                             chi2;
    double                             lnl;
    double                             lambda;

    VSLimitedErrorsHist<double,double> gcc_hist;

    VSAAlgebra::MatrixND               krn;
    VSNSpace                           krn_nspace;

    VSLimitedErrorsHist<double,double> egy_on_hist;
    VSLimitedErrorsHist<double,double> egy_off_hist;
    VSLimitedErrorsHist<double,double> egy_excess_hist;
    VSLimitedErrorsHist<double,double> egy_nu_hist;
    VSLimitedErrorsHist<double,double> egy_resid_hist;
    VSLimitedErrorsHist<double,double> resid_hist;
    VSLimitedErrorsHist<double,double> dfde_hist;
    VSLimitedErrorsHist<double,double> edfde_hist;
    VSLimitedErrorsHist<double,double> e2dfde_hist;
    VSLimitedErrorsHist<double,double> flux_hist;
    VSSimple2DHist<double,double>      flux_cov_hist;
    VSSimple2DHist<double,double>      mu_cov_hist;
    VSSimple2DHist<double,double>      nu_cov_hist;
    VSSimple2DHist<double,double>      k_hist;
    VSSimple2DHist<double,double>      ksq_hist;
    VSSimple2DHist<double,double>      ksq_u_hist;
    VSLimitedHist<double,double>       ksq_w_hist;
    VSLimitedHist<double,double>       ksq_sqw_hist;
    VSLimitedErrorsHist<double,double> ksq_c_hist;
    VSLimitedErrorsHist<double,double> ksq_nc_hist;
    VSSimple2DHist<double,double>      ksqr_hist;
    VSSimple2DHist<double,double>      ksqr_u_hist;
    VSLimitedErrorsHist<double,double> ksqr_c_hist;
    VSLimitedErrorsHist<double,double> ksqr_f_hist;
    VSLimitedErrorsHist<double,double> flux_bias_hist;
    VSLimitedErrorsHist<double,double> flux_rbias_hist;

    VSLimitedErrorsHist<double,double> lambda_gcc_hist;
    VSLimitedErrorsHist<double,double> lambda_chi2_hist;
    VSLimitedErrorsHist<double,double> lambda_chi2b_hist;
    VSLimitedErrorsHist<double,double> lambda_chi2_ml_hist;
    VSLimitedErrorsHist<double,double> lambda_msbias_hist;
    VSLimitedErrorsHist<double,double> lambda_mvar_hist;
    VSLimitedErrorsHist<double,double> lambda_mse_hist;
    VSSimpleGraph<double,double>       chi2_mse_graph;
    VSSimpleGraph<double,double>       chi2_ml_mse_graph;
    VSLimitedErrorsHist<double,double> lambda_wmse_hist;
    VSSimpleGraph<double,double>       chi2_wmse_graph;
    VSSimpleGraph<double,double>       chi2_ml_wmse_graph;

    // Load/Save --------------------------------------------------------------
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;
  };

  // ==========================================================================
  //! VSSpectrumFit: Class for parametric spectrum modeling (forward folding).
  // ==========================================================================
  class VSSpectrumFit
  {
  public:

    class Data
    {
    public:
      
      Data(): on(), off(), alpha(), krn() { }

      VSAAlgebra::VecND on;
      VSAAlgebra::VecND off;
      double            alpha;
      VSNSpace          krn;
    };

    //! Functor for evaluating spectrum likelihood.  Background is evaluated
    //! for each reconstructed energy bin independently.
    class LnLFn 
    {
    public:

      LnLFn();
      LnLFn(const std::vector<Data>& data,
	    VSSpectrumFn* sp_fn);

      ~LnLFn();

      // Accessors ------------------------------------------------------------
      unsigned nparm() const { return m_sp_fn->nparm(); }
      VSAAlgebra::VecND param() const { return m_sp_fn->param(); }
      double param(unsigned ip) const { return m_sp_fn->param(ip); }
      std::vector<bool> fixed() const { return m_sp_fn->fixed(); }
      bool fixed(unsigned ip) const { return m_sp_fn->fixed(ip); }

      double operator() (const VSAAlgebra::VecND& a) const { return val(a); }

      double val() const;
      double val(const VSAAlgebra::VecND& a) const;      

      void gof(const VSAAlgebra::VecND& a,
	       double& chi2_pearson, double& chi2_ml,
	       unsigned& ndf) const;     

      void mu(const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& mu_on,
	      VSAAlgebra::VecND& mu_off,
	      VSAAlgebra::VecND& mu_bkgnd);
 
      void val(const VSAAlgebra::VecND& a, 
	       double& chi2,
	       VSAAlgebra::MatrixND& beta,
	       VSAAlgebra::MatrixND& alpha) const;

      void dyda(const VSAAlgebra::VecND& a, VSAAlgebra::VecND& dyda) const
      {

      }

      unsigned ndata() const 
      {
	unsigned ndata = 0;
	for(std::vector<Data>::const_iterator itr = m_data.begin();
	    itr != m_data.end(); ++itr)
	  ndata += itr->on.ndim() + itr->off.ndim();
	return ndata;
      }

      // Copy-Constructor -----------------------------------------------------
      LnLFn(const LnLFn& o)
      {
	m_data    = o.m_data;
	m_lny     = o.m_lny;
	m_negy    = o.m_negy;
	m_sp_fn   = o.m_sp_fn->clone();
      }

      // Assignment Operator --------------------------------------------------
      LnLFn& operator=(const LnLFn& o)
      {
	m_data  = o.m_data;
	m_lny   = o.m_lny;
	m_negy  = o.m_negy;
	delete m_sp_fn;
	m_sp_fn = o.m_sp_fn->clone();

	return *this;
      }

      
      static void mu_bkgnd(const VSAAlgebra::VecND& on,
			   const VSAAlgebra::VecND& off,
			   double alpha,
			   const VSAAlgebra::VecND& mus,
			   VSAAlgebra::VecND& mub);
    private:

      void mu_signal(const VSNSpace& krn,
		     VSAAlgebra::VecND& mus) const;


      

      double val(const VSAAlgebra::VecND& a,
		 const Data& data, unsigned iegylo, unsigned iegyhi) const;
      void val(const VSAAlgebra::VecND& a, 
	       double& chi2,
	       VSAAlgebra::MatrixND& beta,
	       VSAAlgebra::MatrixND& alpha,
	       const Data& data) const;


      std::vector<Data> m_data;
      double            m_lny;
      unsigned          m_negy;
      VSSpectrumFn*     m_sp_fn;

    };

    VSSpectrumFit(double theta_cut, const std::string& sp_model,
		  double log10_norm_egy = 0);
    virtual ~VSSpectrumFit();

    void fit(const VSAnalysisStage3Data& stage3_data,
	     VSSpectrumFitData& data);

  private:

    double        m_theta_cut;
    VSSpectrumFn* m_sp_fn;
  };

  /////////////////////////////////////////////////////////////////////////////
  //! Abstract base class for non-parametric spectrum calculator.
  //! Instances can be generated with VSSpectrumCalcFactory.
  /////////////////////////////////////////////////////////////////////////////
  class VSSpectrumUnfolding
  {
  public:
    
    class Data
    {
    public:
      
      Data(): 
	flux(), flux_cov(), bias(), bias_cov(), nu(), gcc(),
	ndf(), rnorm(), chi2(), chi2b(), chi2_ml(),
	lambda(), msbias(), mvar(), mse(), wmse(), mean_gcc()
      { }
	
      VSAAlgebra::VecND     flux;
      VSAAlgebra::MatrixND  flux_cov;
      VSAAlgebra::VecND     bias;
      VSAAlgebra::MatrixND  bias_cov;
      VSAAlgebra::VecND     nu;
      VSAAlgebra::VecND     resid;
      VSAAlgebra::MatrixND  c;
      VSAAlgebra::MatrixND  kc;
      VSAAlgebra::MatrixND  k;
      VSAAlgebra::MatrixND  ksq;
      VSAAlgebra::VecND     ksq_w;  // Eigenvalues
      VSAAlgebra::VecND     ksq_sqw;  // Square Root of Eigenvalues
      VSAAlgebra::VecND     ksq_c;  // Fourier Coefficients
      VSAAlgebra::VecND     ksq_nc; // Fourier Coefficients
      VSAAlgebra::VecND     ksq_nc_err; // Fourier Coefficients Error
      VSAAlgebra::MatrixND  ksq_u;      // Eigenvectors
      VSAAlgebra::MatrixND  ksqr;
      VSAAlgebra::VecND     ksqr_c;     // Fourier Coefficients
      VSAAlgebra::MatrixND  ksqr_u;     // Eigenvectors
      VSAAlgebra::VecND     ksqr_f;     // Filter Factors
      VSAAlgebra::VecND     gcc;    // Global correlation coefficient

      unsigned  ksq_ndf;
      double    ksqr_ndf;
      unsigned  ndf;
      double    rnorm;
      double    chi2;
      double    chi2b;
      double    chi2_ml;
      double    dchi2_eff;
      double    lambda;
      double    msbias;
      double    mvar;
      double    mse;
      double    wmse;
      double    mean_gcc;

    };

    VSSpectrumUnfolding(double theta_cut, const std::string& sp_model,
			double ebin, double emin, double emax);
    virtual ~VSSpectrumUnfolding();

    //! Perform non-parametric spectral reconstruction on the data
    //! contained in stage3_data.  Reconstruced spectrum is written to
    //! VSSpectrumData.
    virtual void reconstruct(const VSAnalysisStage3Data& stage3_data,
			     VSSpectrumData& data) = 0;
    
    void finalize(VSSpectrumData& data);

  protected:

    double                         m_ebin;
    double                         m_emin;
    double                         m_emax;

  private:

    double                         m_theta_cut;
    unsigned                       m_max_off_regions;
 
  };

  /////////////////////////////////////////////////////////////////////////////
  //! Class for accumulating spectrum data and performing both parametric
  //! and non-parametric spectral reconstruction.
  //! Instances can be generated with VSSpectrumCalcFactory.
  /////////////////////////////////////////////////////////////////////////////
  class VSSpectrumCalc
  {
  public:
    VSSpectrumCalc(double theta_cut, const std::string& sp_model,
		   double log10_enorm);
    virtual ~VSSpectrumCalc();

    //! Perform non-parametric spectral reconstruction on the data
    //! contained in stage3_data.  Reconstruced spectrum is contained in
    //! VSSpectrumData.
    void reconstruct(const VSAnalysisStage3Data& stage3_data,
		     VSSpectrumData& data);

    void fit(const VSAnalysisStage3Data& stage3_data,
	     VSSpectrumFitData& data);
    
    void finalize(VSSpectrumData& data);

    //    void setRun(const VSAnalysisStage3Data::RunData& data);

    void accumulate(const VSEventArrayDatum& event,
		    const VSAAlgebra::Vec2D& event_xy,
		    VSAnalysisStage3Data::RunData& data);
      
    static VSAAlgebra::VecND 
    toVecND(const VSLimitedErrorsHist<double, double>& h);
    static VSAAlgebra::VecND toVecND(const VSLimitedHist<double, double>& h);
    static VSAAlgebra::VecND toVecND(const VSNSpace& nspace);
    static VSAAlgebra::MatrixND toMatrixND(const VSNSpace& nspace);

    
    unsigned maxOffRegions() const { return m_max_off_regions; }
    double thetaCut() const { return m_theta_cut; }

  private:

    VSSpectrumUnfolding*           m_sp_unfold;
    VSSpectrumFit*                 m_sp_fit;
    double                         m_theta_cut;
    unsigned                       m_max_off_regions;
  };


  //! Non-parametric spectrum calculation by correction factor method.
  class VSSpectrumUnfoldingCFM : public VSSpectrumUnfolding
  {
  public:
    VSSpectrumUnfoldingCFM(double theta_cut,
			   const std::string& sp_model,
			   double ebin, double emin, double emax,
			   double log10_enorm,
			   double ul_threshold,
			   unsigned niter);

    virtual ~VSSpectrumUnfoldingCFM() {}

    virtual void reconstruct(const VSAnalysisStage3Data& stage3_data,
			     VSSpectrumData& data);
    

  private:
    double         m_theta_cut;
    double         m_ul_threshold;
    unsigned       m_niter;
  };


  /////////////////////////////////////////////////////////////////////////////
  //! Non-parametric spectrum calculation with chi2 functional and
  //! regularization.
  /////////////////////////////////////////////////////////////////////////////
  class VSSpectrumUnfoldingChi2 : public VSSpectrumUnfolding
  {
  public:
    VSSpectrumUnfoldingChi2(double theta_cut,
			    const std::string& sp_model,
			    double ebin, double emin, double emax,
			    double log10_enorm,
			    double ul_threshold,
			    const std::string& reg_method,
			    double lambda, double index);
    virtual ~VSSpectrumUnfoldingChi2() {}

    virtual void reconstruct(const VSAnalysisStage3Data& stage3_data,
			     VSSpectrumData& data);

    void reconstruct(const VSNSpace& kernel,
		     const VSAAlgebra::VecND& excess,
		     const VSAAlgebra::VecND& var,
		     const VSAAlgebra::VecND& on,
		     const VSAAlgebra::VecND& off,
		     double alpha,
		     VSSpectrumData& data);

    static VSAAlgebra::MatrixND calcSmoothingMatrix(const std::string& method,
						    const VSNSpace& krn);

  private:

    void calcSpectrum(const VSAAlgebra::MatrixND& k,
		      const VSAAlgebra::MatrixND& s,
		      const VSAAlgebra::MatrixND& vi,
		      const VSAAlgebra::MatrixND& v,
		      const VSAAlgebra::VecND& excess,
		      double lambda,
		      Data& data);

    double      m_theta_cut;
    double      m_ul_threshold;
    std::string m_reg_method;
    double      m_lambda;
    double      m_scaling_index;

  };

#if 0
  /////////////////////////////////////////////////////////////////////////////
  //! Non-parametric spectrum calculation with chi2 functional and
  //! regularization.
  /////////////////////////////////////////////////////////////////////////////
  class VSSpectrumUnfoldingML : public VSSpectrumUnfolding
  {
  public:

    class LnLFn
    {
    public:
      LnLFn();
      LnLFn(const std::vector<VSSpectrumFit::Data>& data,
	    VSSpectrumFn* sp_fn,
	    const VSAAlgebra::MatrixND& s,
	    const VSAAlgebra::MatrixND& sm,
	    double lambda);

      ~LnLFn() { }

      // Accessors ------------------------------------------------------------
      unsigned nparm() const { return m_fn.nparm(); }
      VSAAlgebra::VecND param() const { return m_fn.param(); }
      double param(unsigned ip) const { return m_fn.param(ip); }
      std::vector<bool> fixed() const { return m_fn.fixed(); }
      bool fixed(unsigned ip) const { return m_fn.fixed(ip); }

      double val(const VSAAlgebra::VecND& a) const;     

      void dyda(const VSAAlgebra::VecND& a, VSAAlgebra::VecND& dyda) const
      {

      }

      void val(const VSAAlgebra::VecND& a, 
	       double& chi2,
	       VSAAlgebra::MatrixND& beta,
	       VSAAlgebra::MatrixND& alpha) const;

      unsigned ndata() const 
      {
	return m_fn.ndata();
      }

    private:
      VSSpectrumFit::LnLFn m_fn;
      VSAAlgebra::MatrixND m_s;
      VSAAlgebra::MatrixND m_sm;
      double               m_lambda;
    };

    class SpectrumFn : public VSSpectrumFn
    {
    public:
      SpectrumFn(unsigned nparm, double ebin, double emin);
      
      virtual double val(const double& log10_egy_tev, 
		       const VSAAlgebra::VecND& a) const;
      virtual double val(const double& log10_egy_tev) const;

      virtual void dyda(const double& log10_egy_tev, 
			VSAAlgebra::VecND& dyda) const;
      virtual void dyda(const double& log10_egy_tev, 
			const VSAAlgebra::VecND& a,
			VSAAlgebra::VecND& dyda) const;
      
      virtual double integralFlux(double log10_elo_tev,
				  double log10_ehi_tev,
				  const VSAAlgebra::VecND& a) const;
      virtual double integralFlux(double log10_elo_tev,
				  double log10_ehi_tev) const;

      virtual void integralFluxDerivative(double log10_elo_tev,
					  double log10_ehi_tev,
					  VSAAlgebra::VecND& dyda) const;

      // Virtual Constructor --------------------------------------------------
      virtual SpectrumFn* clone() const { return new SpectrumFn(*this); }
      
    private:

      VSBinCalcLinear<double> m_binner;
      double                  m_ebin;
    };


    VSSpectrumUnfoldingML(double theta_cut,
			  const std::string& sp_model,
			  double ebin, double emin, double emax,
			  double log10_enorm,
			  double ul_threshold,
			  const std::string& reg_method,
			  double lambda);
    virtual ~VSSpectrumUnfoldingML() {}

    virtual void reconstruct(const VSAnalysisStage3Data& stage3_data,
			     VSSpectrumData& data);

    void reconstruct(const VSNSpace& kernel,
		     const VSAAlgebra::VecND& excess,
		     const VSAAlgebra::VecND& var,
		     const VSAAlgebra::VecND& on,
		     const VSAAlgebra::VecND& off,
		     double alpha,
		     VSSpectrumData& data);

    static VSAAlgebra::MatrixND calcSmoothingMatrix(const std::string& method,
						    double bin_width,
						    unsigned nrow,
						    unsigned ncol);

  private:

    void calcSpectrum(const VSAAlgebra::MatrixND& k,
		      const VSAAlgebra::MatrixND& s,
		      const VSAAlgebra::VecND& ep,
		      const VSAAlgebra::VecND& var,		 
		      double lambda,
		      Data& data);

    double      m_theta_cut;
    double      m_ul_threshold;
    std::string m_reg_method;
    double      m_lambda;

  };

#endif

  /////////////////////////////////////////////////////////////////////////////
  //! Singleton factory class for creating instances of 
  //! spectrum calculator (VSSpectrumCalc).  A pointer to the factory is
  //! obtained by calling getInstance().
  /////////////////////////////////////////////////////////////////////////////
  class VSSpectrumCalcFactory
  {
  public:
    struct Options
    {
      Options(); 

      std::string                 np_method;
      std::string                 reg_matrix;
      double                      lambda;
      double                      scaling_index;
      unsigned                    cfm_niter;
      std::string                 spectrum_model;
      double                      spectrum_log10_enorm;
      std::pair< double, double > spectrum_bin_scheme;
      double                      spectrum_bin_width;
      double                      theta_cut;
      double                      upper_limit_threshold;
    };


    virtual ~VSSpectrumCalcFactory();

    static VSSpectrumCalcFactory* getInstance();

    VSSpectrumCalc* create();
    VSSpectrumUnfolding* createSpectrumUnfolding();

    double egyBinWidth() const { return m_egy_bin_width; }
    double egyMin() const { return m_egy_min; }
    double egyMax() const { return m_egy_max; }
    double thetaCut() const { return m_options.theta_cut; }

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="s3_");

  private:
    VSSpectrumCalcFactory(const Options& opt = s_default_options);

    double                      m_egy_bin_width;
    double                      m_egy_min;
    double                      m_egy_max;

    Options                     m_options;

    static Options              s_default_options;

    static std::auto_ptr<VSSpectrumCalcFactory> s_instance;


  };



}

#endif // VSSPECTRUMCALC_HPP
