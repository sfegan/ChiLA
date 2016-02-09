//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLTLibraryVisitor.hpp

  Generate a library of lookup tables.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSLTLibraryVisitor.hpp,v 3.11 2010/07/12 23:07:49 matthew Exp $

*/

#ifndef VSLTLIBRARYVISITOR_HPP
#define VSLTLIBRARYVISITOR_HPP

#include <VSAFunction.hpp>
#include <VSNSpace.hpp>
#include <VSGenSpace.hpp>
#include <VSEventDataVisitor.hpp>
#include <VSEffectiveAreaCalc.hpp>
#include <VSPSFCalc.hpp>
#include <VSSimpleErrorsHist.hpp>

namespace VERITAS
{
  //! Class responsible for generating energy lookup tables.  After
  //! data from a set of simulation files is accumulated, lookup
  //! tables are generated with the call to createLibrary().
  class VSEnergyLTLibraryVisitor : public VSEventDataVisitor
  {
  public:
    
    class Fn : public VSAFunction::ParamFn<VSAAlgebra::Vec2D>
    {
    public:

      Fn();
      void operator() (const VSAAlgebra::Vec2D& x, VSAAlgebra::VecND& v) const;

      double val(const VSAAlgebra::Vec2D& x, const VSAAlgebra::VecND& a) const;

    };

    struct Options
    {
      Options(): 
	min_count(5), use_size_tables(true),
	poly_max_order(8), poly_rchi2_max_change(1.5),	
	median(true), no_scope_tables(true), space("N_R_disp"),
	no_reconstruction_cut(false),
	no_reconstructed_impact(true),
	theta_cut(0.3),
	dist_cut(1.3),
	no_extrapolation(false),
	no_diagnostics(true)
      { }

      unsigned   min_count;
      bool       use_size_tables;
      unsigned   poly_max_order;
      double     poly_rchi2_max_change;
      bool       median;
      bool       no_scope_tables;
      std::string space;
      bool       no_reconstruction_cut;
      bool       no_reconstructed_impact;
      double     theta_cut;
      double     dist_cut;
      bool       no_extrapolation;
      bool       no_diagnostics;
    };

    struct SizeFitData
    {
      SizeFitData(): 
	ndata(), R(), disp(),
	exp_fit_param(), exp_fit_nparam(), exp_fit_chi2(), exp_fit_rchi2(),
	dx(), dy(), ds(), dp(), 
	size_exp_hist(), size_rms_hist(),
	size_fit_exp_hist(), size_fit_rms_hist()	
      { } 

      unsigned            ndata;
      double              R;
      double              disp;

      VSAAlgebra::VecND   exp_fit_param;
      unsigned            exp_fit_nparam;
      double              exp_fit_chi2;
      double              exp_fit_rchi2;

      std::vector<double> dx;
      std::vector<double> dy;
      std::vector<double> ds;
      std::vector<double> dp;

      VSLimitedErrorsHist<double,double> size_exp_hist;
      VSLimitedErrorsHist<double,double> size_rms_hist;
      VSLimitedErrorsHist<double,double> size_fit_exp_hist;
      VSLimitedErrorsHist<double,double> size_fit_rms_hist;
      VSLimitedErrorsHist<double,double> exp_hist;
      VSLimitedErrorsHist<double,double> rms_hist;

      void save(VSOctaveH5WriterStruct* writer) const;
    };

    struct Spaces
    {
      Spaces(const VSNSpace::Space& energy_def, 
	     const VSNSpace::Space& size_def, 
	     unsigned _min_count, const std::string& _comment_base);

      void accumulate(const VSNSpace::Point& p_energy, double log10_energy,
		      const VSNSpace::Point& p_size, double log10_size,
		      double weight = 1.0);

      void calculateMomentTables(const VSNSpace& n,
				 VSNSpace& exp, VSNSpace& exp_var, 
				 VSNSpace& rms, VSNSpace& rms_var);

      void calculateMedianTables(const VSNSpace& n,
				 VSNSpace& exp, VSNSpace& exp_var, 
				 VSNSpace& rms, VSNSpace& rms_var);

      void fit(const VSAMath::Data<double>& xys,
	       double& fit_chi2,
	       double& fit_rchi2,
	       VSAAlgebra::VecND& fit_param,
	       VSAAlgebra::MatrixND& fit_cov,
	       double poly_rchi2_max_change,
	       unsigned poly_max_order);

      void extrapolate(VSNSpace& exp, VSNSpace& exp_var);
      void extrapolate2(VSNSpace& exp, VSNSpace& exp_var);
      void smooth(VSNSpace& exp, VSNSpace& exp_var,
		  double dx1, double dx2, unsigned norder);      

      void calculateMomentEnergyTables(VSNSpace& exp, VSNSpace& exp_var, 
				       VSNSpace& rms, VSNSpace& rms_var);
      void calculateMedianEnergyTables(VSNSpace& exp, VSNSpace& exp_var, 
				       VSNSpace& rms, VSNSpace& rms_var);
      void calculateMomentSizeTables(VSNSpace& exp, VSNSpace& exp_var, 
				     VSNSpace& rms, VSNSpace& rms_var);
      void calculateMedianSizeTables(VSNSpace& exp, VSNSpace& exp_var, 
				     VSNSpace& rms, VSNSpace& rms_var);
     
      void size2Energy(const VSNSpace& s_exp, const VSNSpace& s_exp_var, 
		       const VSNSpace& s_rms, const VSNSpace& s_rms_var, 
		       const VSNSpace& s_n, 
		       VSNSpace& e_exp, VSNSpace& e_exp_var, 
		       VSNSpace& e_rms, VSNSpace& e_rms_var, 
		       std::vector<SizeFitData>& size_fit_data,
		       double poly_rchi2_max_change = 1.5,
		       unsigned poly_max_order = 7);

      unsigned min_count;
      std::string comment_base;
      
      VSNSpace::Space                    energy_space;
      VSNSpace                           n_energy;
      VSGenSpace<MeanWeight>             log10_energy_mean;
      VSGenSpace<StandardDevWeight>      log10_energy_rms;
      VSGenSpace<OneSidedIntervalWeight> log10_energy_med;
      VSGenSpace<TwoSidedIntervalWeight> log10_energy_i68;
      
      VSNSpace                           n_size;
      VSGenSpace<MeanWeight>             log10_size_mean;
      VSGenSpace<StandardDevWeight>      log10_size_rms;
      VSGenSpace<OneSidedIntervalWeight> log10_size_med;
      VSGenSpace<TwoSidedIntervalWeight> log10_size_i68;
 
    };

    struct LTData
    {
      LTData():
	m_size_n(), 
	m_size_exp(), m_size_exp_var(), 
	m_size_rms(), m_size_rms_var(),
	m_energy_n(), 
	m_energy_mask(), 
	m_energy_exp(), m_energy_exp_var(), 
	m_energy_rms(), m_energy_rms_var(),
	m_size_fit_data()
      { }

      VSNSpace                      m_size_n;
      VSNSpace                      m_size_exp;
      VSNSpace                      m_size_exp_var;
      VSNSpace                      m_size_rms;
      VSNSpace                      m_size_rms_var;
      VSNSpace                      m_energy_n;
      VSNSpace                      m_energy_mask;
      VSNSpace                      m_energy_exp;
      VSNSpace                      m_energy_exp_var;
      VSNSpace                      m_energy_rms;
      VSNSpace                      m_energy_rms_var;

      std::vector<SizeFitData>      m_size_fit_data;

      std::vector< VSLimitedErrorsHist<double,double> > m_energy_fit_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_energy_fit_rms_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_energy_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_energy_rms_hist;

      void save(VSOctaveH5WriterStruct* writer) const;
    };

    struct Data
    {
      Data(const VSNSpace::Space& energy_space, 
	   const VSNSpace::Space& size_space, 
	   unsigned min_count,
	   double zn, double az, double ped);
      ~Data();

      double                        zenith_deg;
      double                        azimuth_deg;
      double                        ped_dev;

      SEphem::SphericalCoords       m_azel;

      LTData                        m_array_lt;
      std::vector<LTData>           m_scope_lt;

      Spaces                        m_sp_all;
      std::vector<Spaces*>          m_sp_tel;
    };

    VSEnergyLTLibraryVisitor(const Options& opt = s_default_options);
    virtual ~VSEnergyLTLibraryVisitor();

    virtual void getSimSet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("energy_tev");
      if(m_options.no_reconstructed_impact)
	{
	  s.insert("core_east_m");
	  s.insert("core_north_m");
	  s.insert("core_elevation_asl_m");
	  s.insert("primary_azimuth_deg");
	  s.insert("primary_zenith_deg");
	}
    }

    virtual void getScopeSet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("used_in_reconstruction");
      s.insert("R");
      s.insert("N");
      s.insert("G");
      s.insert("theta1");
      s.insert("fp_disp");
      s.insert("fp_dist");
      s.insert("lambdad");
      m_cuts->getScopeParamSet(s);
    }

    virtual void getArraySet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("event_num");
      s.insert("scope");
      s.insert("theta1");
      m_cuts->getArrayParamSet(s);
    }

    virtual void visitRun(const VSAnalysisStage1Data& stage1,
			  const VSTargetTable::Observation& obs,
			  const VSArrayMergedCalibrationData& cal);
  
    virtual void leaveRun();
  
    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();
  
    virtual void visitScope(const VSEventScopeDatum& scope, unsigned iscope);
    virtual void leaveScope();

    virtual void visitSimEvent(const VSArraySimulationDatum& sim);
    virtual void leaveSimEvent();
  
    virtual void visitSimHeader(const VSHeaderSimulationDatum& header);
    virtual void leaveSimHeader();
    
    void createLibrary();
    void createLibrary(Spaces& sp, LTData& data);

    void save(VSOctaveH5WriterStruct* writer) const;

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");
  private:

    VSCutsEvaluator*       m_cuts;
    VSSimEnergyWeightCalc* m_egywt_calc;

    unsigned               m_ndim;
    VSNSpace::Space        m_energy_space;
    VSNSpace::Space        m_size_space;

    double                 m_zn_deg;
    double                 m_az_deg;
    double                 m_ped_dev;
    double                 m_offset_deg;

    double                 m_log10_egybin;
    double                 m_log10_egylo;
    double                 m_log10_egyhi;

    std::vector< Data* >   m_data;
    Data*                  m_data_ptr;
    double                 m_wt;
    bool                   m_is_selected;

    std::vector<VSAAlgebra::Vec3D> m_scope_pos;
    VSAAlgebra::Vec3D              m_sim_core;
    VSAAlgebra::Vec3D              m_sim_azel;

    VSArraySimulationDatum  m_sim;
    VSEventArrayDatum       m_evt;
    VSHeaderSimulationDatum m_sim_header;

    Options                 m_options;

    static Options          s_default_options;
  };


  //! Class responsible for generating the mean scaled parameter
  //! lookup tables (width, length, disp).  After data from a set of
  //! simulation files is accumulated, lookup tables are generated
  //! with the call to createLibrary().
  class VSScaledParameterLibraryVisitor : public VSEventDataVisitor
  {
  public:
   
    struct Options
    {
      Options(): 
	min_count(5), theta_cut(0.30), median(true), no_scope_tables(true), 
	no_diagnostics(true), ndim(2)
      { }

      unsigned   min_count;
      double     theta_cut;
      bool       median;
      bool       no_scope_tables;
      bool       no_diagnostics;
      unsigned   ndim;
    };

    struct Spaces
    {
      Spaces(const VSNSpace::Space& def, unsigned _min_count,
	     const std::string& _comment_base);

      void accumulate(const VSNSpace::Point& p, 
		      double width, double length, double disp, 
		      double weight = 1.0);

      void calculateMomentWidthTables(VSNSpace& exp, VSNSpace& exp_var,
				      VSNSpace& rms, VSNSpace& rms_var);
      void calculateMomentLengthTables(VSNSpace& exp, VSNSpace& exp_var,
				       VSNSpace& rms, VSNSpace& rms_var);
      void calculateMomentDispTables(VSNSpace& exp, VSNSpace& exp_var,
				     VSNSpace& rms, VSNSpace& rms_var);
      void calculateMedianWidthTables(VSNSpace& exp, VSNSpace& exp_var,
				      VSNSpace& rms, VSNSpace& rms_var);
      void calculateMedianLengthTables(VSNSpace& exp, VSNSpace& exp_var,
				       VSNSpace& rms, VSNSpace& rms_var);
      void calculateMedianDispTables(VSNSpace& exp, VSNSpace& exp_var,
				     VSNSpace& rms, VSNSpace& rms_var);

      unsigned min_count;
      std::string comment_base;

      VSNSpace n;

      VSGenSpace<MeanWeight>             width_mean;
      VSGenSpace<StandardDevWeight>      width_rms;
      VSGenSpace<MeanWeight>             length_mean;
      VSGenSpace<StandardDevWeight>      length_rms;
      VSGenSpace<MeanWeight>             disp_mean;
      VSGenSpace<StandardDevWeight>      disp_rms;

      VSGenSpace<OneSidedIntervalWeight> width_med;
      VSGenSpace<TwoSidedIntervalWeight> width_i68;
      VSGenSpace<OneSidedIntervalWeight> length_med;
      VSGenSpace<TwoSidedIntervalWeight> length_i68;
      VSGenSpace<OneSidedIntervalWeight> disp_med;
      VSGenSpace<TwoSidedIntervalWeight> disp_i68;
    };

    struct LTData
    {
      LTData(); 

      VSNSpace                      m_width_mask;
      VSNSpace                      m_width_exp;
      VSNSpace                      m_width_exp_var;
      VSNSpace                      m_width_rms;
      VSNSpace                      m_width_rms_var;
      VSNSpace                      m_width_fit_exp;
      VSNSpace                      m_width_fit_rms;

      VSNSpace                      m_length_mask;
      VSNSpace                      m_length_exp;
      VSNSpace                      m_length_exp_var;
      VSNSpace                      m_length_rms;
      VSNSpace                      m_length_rms_var;
      VSNSpace                      m_length_fit_exp;
      VSNSpace                      m_length_fit_rms;

      VSNSpace                      m_disp_mask;
      VSNSpace                      m_disp_exp;
      VSNSpace                      m_disp_exp_var;
      VSNSpace                      m_disp_rms;
      VSNSpace                      m_disp_rms_var;
      VSNSpace                      m_disp_fit_exp;
      VSNSpace                      m_disp_fit_rms;

      std::vector< VSLimitedErrorsHist<double,double> > m_width_fit_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_width_fit_rms_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_width_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_width_rms_hist;

      std::vector< VSLimitedErrorsHist<double,double> > m_length_fit_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_length_fit_rms_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_length_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_length_rms_hist;

      std::vector< VSLimitedErrorsHist<double,double> > m_disp_fit_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_disp_fit_rms_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_disp_exp_hist;
      std::vector< VSLimitedErrorsHist<double,double> > m_disp_rms_hist;

      void save(VSOctaveH5WriterStruct* writer) const;
    };


    struct Data
    {
      Data(const VSNSpace::Space& space, unsigned min_count,
	   double zn, double az, double ped);
      ~Data();

      double                        zenith_deg;
      double                        azimuth_deg;
      double                        ped_dev;

      SEphem::SphericalCoords       m_azel;

      LTData                        m_array_lt;
      std::vector<LTData>           m_scope_lt;

      Spaces                        m_sp_all;
      std::vector<Spaces*>          m_sp_tel;
    };

    VSScaledParameterLibraryVisitor(const Options& opt = s_default_options);
    virtual ~VSScaledParameterLibraryVisitor();

    virtual void visitRun(const VSAnalysisStage1Data& stage1,
			  const VSTargetTable::Observation& obs,
			  const VSArrayMergedCalibrationData& cal);
  
    virtual void leaveRun();
  
    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();
  
    virtual void visitScope(const VSEventScopeDatum& scope, unsigned iscope);
    virtual void leaveScope();

    virtual void visitSimEvent(const VSArraySimulationDatum& sim);
    virtual void leaveSimEvent();
  
    virtual void visitSimHeader(const VSHeaderSimulationDatum& header);
    virtual void leaveSimHeader();

    //! Generate the set of lookup tables using the accumulated
    //! simulation data.
    void createLibrary();

    void createLibrary(Spaces& sp, LTData& data);

    void fit(LTData& data);
    void fit(VSNSpace& exp, VSNSpace& var, VSNSpace& fit, unsigned npoly);

    void smooth(VSNSpace& exp, VSNSpace& var, VSNSpace& fit);

    void save(VSOctaveH5WriterStruct* writer) const;

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

  private:
    VSCutsEvaluator*       m_cuts;
    VSSimEnergyWeightCalc* m_egywt_calc;

    VSNSpace::Space        m_space;

    double                 m_zn_deg;
    double                 m_az_deg;
    double                 m_ped_dev;
    double                 m_offset_deg;

    std::vector< Data* >   m_data;
    Data*                  m_data_ptr;
    double                 m_wt;
    bool                   m_is_selected;

    VSArraySimulationDatum  m_sim;
    VSEventArrayDatum       m_evt;
    VSHeaderSimulationDatum m_sim_header;

    Options                 m_options;

    static Options          s_default_options;
  };


  // ==========================================================================
  // VSEffectiveAreaLibraryVisitor
  // ==========================================================================
  class VSEffectiveAreaLibraryVisitor : public VSEventDataVisitor
  {
  public:

    struct Options
    {
      Options(): 
	theta_cut(0.3),
	kernel_theta_cut(0.15),
	no_diagnostics(true)
      { }

      double     theta_cut;
      double     kernel_theta_cut;
      bool       no_diagnostics;
    };    

    struct Data
    {
      Data(): 
	zenith_deg(), azimuth_deg(), ped_dev(), m_azel(),
	m_effarea_calc()
      { }
      Data(double zn, double az, double ped): 
	zenith_deg(zn), azimuth_deg(az), ped_dev(ped), m_azel(),
	m_effarea_calc()
      { 
	m_azel = SEphem::SphericalCoords::makeDeg(zn,az);
      }

      ~Data();

      double                        zenith_deg;
      double                        azimuth_deg;
      double                        ped_dev;

      SEphem::SphericalCoords       m_azel;

      VSEffectiveAreaCalc           m_effarea_calc;
      VSEffectiveAreaCalc           m_effarea_trigger_calc;
      VSEnergyKernelCalc            m_kernel_calc;
      VSPSFCalc                     m_psf_calc;

      VSScaledParameterCalc*        m_sp_calc;
      VSEnergyCalcLT*               m_egy_calc;

      VSSimple2DHist<double,double> m_effarea_offset_nspace;
    };

    VSEffectiveAreaLibraryVisitor(const Options& opt = s_default_options);
    virtual ~VSEffectiveAreaLibraryVisitor();

    virtual void getSimSet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("energy_tev");
    }

    virtual void getScopeSet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("used_in_reconstruction");
      s.insert("R");
      s.insert("N");
      s.insert("G");
      s.insert("theta1");
      s.insert("intrinsic_width");
      s.insert("intrinsic_length");
      s.insert("sc_width");
      s.insert("sc_length");
      s.insert("sc_disp");
      s.insert("fp_disp");
      s.insert("fp_dist");
      s.insert("lambdad");
      s.insert("lt_log10_energy");
      s.insert("lt_log10_energy_err");
      m_cuts->getScopeParamSet(s);
    }

    virtual void getArraySet(VSEventDataReader::MemberSubset& s) const
    { 
      s.insert("event_num");
      s.insert("scope");
      s.insert("theta1");
      s.insert("mlt_log10_energy");
      s.insert("mlt_log10_energy_chi2");
      s.insert("used_in_reconstruction_mask");
      m_cuts->getArrayParamSet(s);
    }

    virtual void visitRun(const VSAnalysisStage1Data& stage1,
			  const VSTargetTable::Observation& obs,
			  const VSArrayMergedCalibrationData& cal);
  
    virtual void leaveRun();
  
    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();
  
    virtual void visitSimEvent(const VSArraySimulationDatum& sim);
    virtual void leaveSimEvent();
  
    virtual void visitSimHeader(const VSHeaderSimulationDatum& header);
    virtual void leaveSimHeader();

    void createLibrary();

    void save(VSOctaveH5WriterStruct* writer) const;

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

  private:

    VSCutsEvaluator*       m_cuts;

    double                 m_zn_deg;
    double                 m_az_deg;
    double                 m_ped_dev;
    double                 m_offset_deg;

    double                 m_log10_egybin;
    double                 m_log10_egylo;
    double                 m_log10_egyhi;

    std::vector< Data* >   m_data;
    Data*                  m_data_ptr;

    VSArraySimulationDatum  m_sim;
    VSEventArrayDatum       m_evt;
    VSHeaderSimulationDatum m_sim_header;

    Options                 m_options;

    static Options          s_default_options;
  };

}

#endif // VSLTLIBRARYVISITOR_HPP
