//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage3Data.hpp

  Stage 1 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/12/2006

  $Id: VSAnalysisStage3Data.hpp,v 3.20 2010/10/20 03:59:11 matthew Exp $

*/

#ifndef VSANALYSISSTAGE3DATA_HPP
#define VSANALYSISSTAGE3DATA_HPP

#include <VSAAlgebra.hpp>

#include <VSPointing.hpp>

#include <VSSimple2DHist.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSNSpace.hpp>
#include <VSInstrumentResponseCalc.hpp>
#include <VSScaledParameterCalc.hpp>
#include <VSEnergyCalcLT.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSOnOffHist1D
  // ==========================================================================
  template< typename T, typename count_type, 
	    template <typename, typename, typename BIN=VSBinCalcLinear<T> > 
	    class HIST >
  class VSOnOffHist1D
  {
  public:
    VSOnOffHist1D(): 
      on_hist(), off_hist(), excess_hist(), total_hist()
    { }

    VSOnOffHist1D(const std::string& _name, T bin_size, T lo, T hi): 
      name(_name),
      on_hist(bin_size,lo,hi), 
      off_hist(bin_size,lo,hi), 
      excess_hist(bin_size,lo,hi), 
      total_hist(bin_size,lo,hi)
    { }

    ~VSOnOffHist1D() {}

    void accumulate(T x, double wt = 1.0) { total_hist.accumulate(x,wt); }
    void accumulateOn(T x, double wt = 1.0) { on_hist.accumulate(x,wt); }
    void accumulateOff(T x, double wt = 1.0) { off_hist.accumulate(x,wt); }
    void calcExcess(double alpha);

    bool load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    VSOnOffHist1D& operator+= (const VSOnOffHist1D& o);

    std::string         name;

    HIST<T,count_type>  on_hist;
    HIST<T,count_type>  off_hist;
    HIST<T,count_type>  excess_hist;
    HIST<T,count_type>  total_hist;
  };

  template< typename T, typename count_type, 
	    template <typename, typename, typename BIN=VSBinCalcLinear<T> > 
	    class HIST >
  void VSOnOffHist1D<T,count_type,HIST>::calcExcess(double alpha)
  {
    off_hist *= alpha;
    excess_hist = on_hist - off_hist;
  }

  template< typename T, typename count_type, 
	    template <typename, typename, typename BIN=VSBinCalcLinear<T> > 
	    class HIST >
  bool VSOnOffHist1D<T,count_type,HIST>::load(VSOctaveH5ReaderStruct* reader)
  {
    if(reader == NULL) return false;
    VSOctaveH5ReaderStruct *s = reader->readStruct(name);
    if(s == NULL) return false;

    on_hist.load(s->readStruct("on_hist"));
    off_hist.load(s->readStruct("off_hist"));
    excess_hist.load(s->readStruct("excess_hist"));
    total_hist.load(s->readStruct("total_hist"));
    
    delete s;
    return true;
  }
  
  template< typename T, typename count_type, 
	    template <typename, typename, typename BIN=VSBinCalcLinear<T> > 
	    class HIST >
  void VSOnOffHist1D<T,count_type,HIST>::save(VSOctaveH5WriterStruct* writer) 
    const
  {
    VSOctaveH5WriterStruct *s = writer->writeStruct(name);
    vsassert(s);
    on_hist.save(s->writeStruct("on_hist"));
    off_hist.save(s->writeStruct("off_hist"));
    excess_hist.save(s->writeStruct("excess_hist"));
    total_hist.save(s->writeStruct("total_hist"));
    delete s;
  }

  template< typename T, typename count_type, 
	    template <typename, typename, typename BIN=VSBinCalcLinear<T> > 
	    class HIST >
  VSOnOffHist1D<T,count_type,HIST>&
  VSOnOffHist1D<T,count_type,HIST>::operator+=
  (const VSOnOffHist1D<T,count_type,HIST>& o)
  {
    on_hist += o.on_hist;
    off_hist += o.off_hist;
    excess_hist += o.excess_hist;
    total_hist += o.total_hist;
    return *this;
  }

  typedef VSOnOffHist1D< double, double, VSLimitedErrorsHist > Hist1D;

  // ==========================================================================
  // VSAnalysisStage3ScopeDatum
  //
  // Scope based distributions go here
  // ==========================================================================
  class VSAnalysisStage3OnOffHistDatum
  {
  public:
//     class Scope
//     {
//     public:
//       Scope();
//       ~Scope() {}

//       void accumulate(const VSEventScopeDatum& scope_data, bool is_on_event,
// 		      bool is_off_event);
//       void calcExcess(double alpha);
//       void merge(Scope* results_data);

//       void load(VSOctaveH5ReaderStruct* reader);
//       void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

//       static void _compose(VSOctaveH5CompositeDefinition& c)
//       { }

//       // Parameter Distributions ----------------------------------------------
//       Hist1D              fp_dist;
//       Hist1D              fp_disp;
//       Hist1D              fp_width;
//       Hist1D              fp_length;
//       Hist1D              log10_N;
//       Hist1D              G;
//       Hist1D              delta2l;
//       Hist1D              delta2m;
//       Hist1D              log10_lambdac;
//       Hist1D              log10_lambdad;
//       Hist1D              sc_disp;
//       Hist1D              sc_width;
//       Hist1D              sc_length;
//     };

    VSAnalysisStage3OnOffHistDatum();

    // Parameter Distributions ------------------------------------------------
//     Hist1D                               thetasq;
//     Hist1D                               core_R;
    Hist1D                               msc_width;
    Hist1D                               msc_length;
    Hist1D                               msc_disp;
    Hist1D                               log10_N2;

    void accumulate(const VSEventArrayDatum& event, double wt = 1.0);
    void accumulateOn(const VSEventArrayDatum& event, double wt = 1.0);
    void accumulateOff(const VSEventArrayDatum& event, double wt = 1.0);
    void calcExcess(double alpha);

    VSAnalysisStage3OnOffHistDatum& operator+= 
    (const VSAnalysisStage3OnOffHistDatum& o);

    bool load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;
  };

  class VSAnalysisStage3DataBase
  {
  public:

    VSAnalysisStage3DataBase();
    VSAnalysisStage3DataBase(const SEphem::SphericalCoords& origin_radec);
    VSAnalysisStage3DataBase(const SEphem::SphericalCoords& origin_radec,
			     double offset_max, double bin_size,
			     double log10_egy_binsize, double log10_egy_min,
			     double log10_egy_max);

    // Position of coordinate origin ------------------------------------------
    double                               origin_ra_rad;
    double                               origin_dec_rad;
    std::string                          origin_ra_hms;
    std::string                          origin_dec_dms;

    // Various distributions in sky and camera coordinates --------------------
    VSSimple2DHist<double,double>        sky_counts_hist;
    VSSimple2DHist<double,double>        sky_livetime_hist;
    VSSimple2DHist<double,double>        cam_counts_hist;
    VSSimple2DHist<double,double>        fov_counts_hist;
    
    VSLimitedErrorsHist<double, double>  fovr2_counts_hist;
    VSLimitedErrorsHist<double, double>  th2_on_counts_hist;
    VSLimitedErrorsHist<double, double>  th2_off_counts_hist;

    // Histograms for Energy Spectrum -----------------------------------------
    VSLimitedErrorsHist<double, double>  egy_on_hist;
    VSLimitedErrorsHist<double, double>  egy_off_hist;
    VSLimitedErrorsHist<double, double>  egy_ring_hist;   
 
    VSSimple2DHist<double,double>        egy_th2_on_hist;
    VSSimple2DHist<double,double>        egy_th2_off_hist;

    VSLimitedErrorsHist<double, double>  egy_on_np_hist;
    VSLimitedErrorsHist<double, double>  egy_off_np_hist;
    VSLimitedErrorsHist<double, double>  egy_ring_np_hist;    
    
    VSSimple2DHist<double,double>        egy_th2_on_np_hist;
    VSSimple2DHist<double,double>        egy_th2_off_np_hist;

    VSNSpace                             egy_excess_counts;
    VSNSpace                             egy_excess_counts_var;

    // Parameter Histograms ---------------------------------------------------
    VSLimitedErrorsHist<double,double>   log10_N2_hist;
    VSLimitedErrorsHist<double,double>   msc_width_hist;
    VSLimitedErrorsHist<double,double>   msc_length_hist;
    VSLimitedErrorsHist<double,double>   msc_disp_hist;

    VSLimitedErrorsHist<unsigned,double> triggered_hist;
    VSLimitedErrorsHist<unsigned,double> has_image_hist;
    VSLimitedErrorsHist<unsigned,double> used_in_reconstruction_hist;

    VSAnalysisStage3OnOffHistDatum       hist_datum;

    // Effarea/PSF/Energy Kernel  ---------------------------------------------
    VSNSpace                             int_effarea_nspace;
    VSNSpace                             int_effarea_aperture_nspace;
    VSNSpace                             int_effarea_psf_nspace;
    
    VSNSpace                             effarea_nspace;

    const SEphem::SphericalCoords& origin_radec() const 
    { return m_origin_radec; }

    // Misc Methods ---------------------------------------------------------
    void accumulate(const VSEventArrayDatum& event, double wt = 1.0);
    void accumulateOn(const VSEventArrayDatum& event, double wt = 1.0);
    void accumulateOff(const VSEventArrayDatum& event, double wt = 1.0);

    // Load/Save --------------------------------------------------------------
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSAnalysisStage3DataBase,origin_ra_rad);
      H5_ADDMEMBER(c,VSAnalysisStage3DataBase,origin_dec_rad);
      H5_ADDMEMBER(c,VSAnalysisStage3DataBase,origin_ra_hms);
      H5_ADDMEMBER(c,VSAnalysisStage3DataBase,origin_dec_dms);
    }

  protected:
    SEphem::SphericalCoords              m_origin_radec;

  };

  // ==========================================================================
  // VSAcceptanceData
  // ==========================================================================
  class VSAcceptanceData
  {
  public:
    VSAcceptanceData(double offset_max = 1.8,
		     const std::string& model = "",
		     const VSAAlgebra::VecND& param = VSAAlgebra::VecND(),
		     const VSAAlgebra::MatrixND& cov = VSAAlgebra::MatrixND());
    virtual ~VSAcceptanceData() { }


    virtual VSAcceptanceData* clone();

    // Load/Save --------------------------------------------------------------
    virtual bool load(VSOctaveH5ReaderStruct* reader);
    virtual void save(VSOctaveH5WriterStruct* writer) const;

    std::string                                  model;

    double                                       chi2;
    VSAAlgebra::VecND                            param;
    VSAAlgebra::VecND                            param_err;
    VSAAlgebra::MatrixND                         param_cov;
    
    VSAAlgebra::VecND                            acceptance_param;
    VSAAlgebra::VecND                            acceptance_param_err;
    VSAAlgebra::MatrixND                         acceptance_param_cov;

    VSSimple2DHist<double,double>                fov_acceptance_hist;
    VSSimple2DHist<double,double>                fov_bkgnd_hist;
    VSLimitedErrorsHist<double, double>          fovr2_acceptance_hist;
    VSLimitedErrorsHist<double, double>          fovr2_bkgnd_hist;
    VSSimple2DHist<double,double>                sky_bkgnd_hist;
    VSLimitedErrorsHist<double, double>          th2_bkgnd_hist;
    VSSimple2DHist<double,double>                lnl_hist;
    VSSimple2DHist<double,double>                lnl_diff_hist;
  };

  // ==========================================================================
  // VSAnalysisStage3Data
  //
  // This is the base data class used by other stage3 calculators to
  // generate sky maps and/or create spectra.
  // ==========================================================================
  class VSAnalysisStage3Data : public VSAnalysisStage3DataBase
  {
  public:

    static const unsigned EC_TRIGGERED          = 0x00000001;
    static const unsigned EC_RECONSTRUCTED      = 0x00000002;
    static const unsigned EC_SELECTED           = 0x00000004;
    static const unsigned EC_ON                 = 0x00000008;
    static const unsigned EC_OFF                = 0x00000010;

    // ========================================================================
    // VSAnalysisStage3Data::RunData
    // ========================================================================
    class RunData : public VSAnalysisStage3DataBase
    {
    public:
      RunData();
      RunData(double offset_max, double bin_size,
	      double log10_egy_binsize, double log10_egy_min,
	      double log10_egy_max,
	      const SEphem::SphericalCoords& src_radec,
	      const SEphem::SphericalCoords& origin_radec,
	      const VSTargetTable::Observation& obs);
      ~RunData();

      // ======================================================================
      // Public Data
      // ======================================================================
      std::string                   src_name;

      double                        src_ra_rad;
      double                        src_dec_rad;
      std::string                   src_ra_hms;
      std::string                   src_dec_dms;

      double                        obs_ra_rad;
      double                        obs_dec_rad;
      std::string                   obs_ra_hms;
      std::string                   obs_dec_dms;

      double                        ptg_ra_rad;
      double                        ptg_dec_rad;
      std::string                   ptg_ra_hms;
      std::string                   ptg_dec_dms;

      unsigned                      run_number;
      std::string                   mode;
      double                        wobble_theta_deg;
      double                        wobble_phi_deg;

      double                        zn_mean_deg;
      double                        zn_rms_deg;
      double                        az_mean_deg;
      double                        az_rms_deg;
      double                        pangle_mean_deg;
      double                        pangle_rms_deg;
      double                        mean_scaled_dev;

      VSTime                        lo_event_time;
      std::string                   lo_event_time_string;
      VSTime                        hi_event_time;
      std::string                   hi_event_time_string;

      int64_t                       livetime_ticks;
      int64_t                       elaptime_ticks;
      double                        livetime_min;
      double                        elaptime_min;
      unsigned                      ring_counts;
      unsigned                      on_counts;
      unsigned                      off_counts;

      // Accessor Methods -----------------------------------------------------
      const VSAAlgebra::Vec2D& src_xy() const { return m_src_xy; }
      const VSAAlgebra::Vec2D& obs_xy() const { return m_obs_xy; }
      const VSAAlgebra::Vec2D& ptg_xy() const { return m_ptg_xy; }

      const unsigned noff() const { return m_off_xy.size(); }
      const VSAAlgebra::Vec2D& off_xy(unsigned i) const { return m_off_xy[i]; }
      const std::vector<VSAAlgebra::Vec2D>& off_xy() const { return m_off_xy; }

      const unsigned sp_noff() const { return m_sp_off_xy.size(); }
      const VSAAlgebra::Vec2D& sp_off_xy(unsigned i) const 
      { return m_sp_off_xy[i]; }
      const std::vector<VSAAlgebra::Vec2D>& sp_off_xy() const 
      { return m_sp_off_xy; }

      const SEphem::SphericalCoords& obs_radec() const { return m_obs_radec; }

      const VSInstrumentResponseCalc& irf_calc() const { return m_irf_calc; }
      const VSScaledParameterCalc* sp_calc() const { return m_sp_calc; }
      const VSEnergyCalcLT* egy_calc() const { return m_egy_calc; }

      // Set Methods ----------------------------------------------------------
      std::vector<VSAAlgebra::Vec2D>& off_xy() { return m_off_xy; }
      std::vector<VSAAlgebra::Vec2D>& sp_off_xy() { return m_sp_off_xy; }

      VSInstrumentResponseCalc& irf_calc() { return m_irf_calc; }
      VSScaledParameterCalc* sp_calc() { return m_sp_calc; }
      VSEnergyCalcLT* egy_calc() { return m_egy_calc; }

      void setPointing(unsigned ipointing, 
		       const SEphem::SphericalCoords& radec);
		
      void calcExcess()
      {
	hist_datum.calcExcess(1./(double)noff());
      }
		
      // Load/Save ------------------------------------------------------------
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

      // Copy Constructor / Assignment Operator -------------------------------
      RunData(const RunData& o);
      RunData& operator=(const RunData& o);

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,src_name);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,src_ra_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,src_dec_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,src_ra_hms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,src_dec_dms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,obs_ra_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,obs_dec_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,obs_ra_hms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,obs_dec_dms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,ptg_ra_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,ptg_dec_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,ptg_ra_hms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,ptg_dec_dms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,run_number);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,mode);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,wobble_theta_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,wobble_phi_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,zn_mean_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,zn_rms_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,az_mean_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,az_rms_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,pangle_mean_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,pangle_rms_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,mean_scaled_dev);
	H5_ADDSIMPLECOMPOSITE(c,VSAnalysisStage3Data::RunData,lo_event_time);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,lo_event_time_string);
	H5_ADDSIMPLECOMPOSITE(c,VSAnalysisStage3Data::RunData,hi_event_time);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,hi_event_time_string);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,livetime_ticks);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,elaptime_ticks);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,livetime_min);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,elaptime_min);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,ring_counts);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,on_counts);
	H5_ADDMEMBER(c,VSAnalysisStage3Data::RunData,off_counts);
      }

    private:
      void copy(const RunData& o);

      // Miscellaneous coordinates --------------------------------------------
      VSAAlgebra::Vec2D                    m_src_xy;
      VSAAlgebra::Vec2D                    m_obs_xy;
      VSAAlgebra::Vec2D                    m_ptg_xy;
      std::vector< VSAAlgebra::Vec2D >     m_off_xy;
      std::vector< VSAAlgebra::Vec2D >     m_sp_off_xy;
      SEphem::SphericalCoords              m_src_radec;
      SEphem::SphericalCoords              m_obs_radec;
      SEphem::SphericalCoords              m_ptg_radec;      

      // Effarea, PSF, and Energy Kernel Calculator ---------------------------
      VSInstrumentResponseCalc             m_irf_calc;

      // Scaled Parameter and Energy Calculators ------------------------------
      VSScaledParameterCalc*               m_sp_calc;
      VSEnergyCalcLT*                      m_egy_calc;
    };
    
    VSAnalysisStage3Data();
    VSAnalysisStage3Data(const std::string& src_name,
			 const SEphem::SphericalCoords& src_radec,
			 const SEphem::SphericalCoords& origin_radec);
    VSAnalysisStage3Data(const std::string& src_name,
			 const SEphem::SphericalCoords& src_radec,
			 const std::vector< SEphem::SphericalCoords>&
			 ptg_radec,
			 const SEphem::SphericalCoords& origin_radec);
    ~VSAnalysisStage3Data();

    // Accessors --------------------------------------------------------------
    unsigned nrun() const { return m_run_data.size(); }
    const std::vector<RunData>& run_data() const { return m_run_data; }
    const RunData& run_data(unsigned irun) const { return m_run_data[irun]; }
    VSAcceptanceData* acceptance_data() const { return m_acceptance_data; }
    

    // Setters ----------------------------------------------------------------
    std::vector<RunData>& run_data() { return m_run_data; }
    RunData& run_data(unsigned irun) { return m_run_data[irun]; }

    void addData(const std::vector< RunData* >& run_data);
    void addData(const RunData& data);
    void setAcceptanceData(VSAcceptanceData* data)
    {
      delete m_acceptance_data;
      m_acceptance_data = data;
    }

    void clear();

    // Load/Save --------------------------------------------------------------
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer, bool write_hists = true) const;

    // Copy Constructor / Assignment Operator ---------------------------------
    VSAnalysisStage3Data(const VSAnalysisStage3Data& o);
    VSAnalysisStage3Data& operator=(const VSAnalysisStage3Data& o);
    
    // Data -------------------------------------------------------------------
    std::string                          src_name;
    double                               src_ra_rad;
    double                               src_dec_rad;
    std::string                          src_ra_hms;
    std::string                          src_dec_dms;

    double                               pangle_mean_deg;
    double                               pangle_rms_deg;

    int64_t                              livetime_ticks;
    int64_t                              elaptime_ticks;
    double                               livetime_min;
    double                               elaptime_min;

    SEphem::SphericalCoords              src_radec;
    std::vector<SEphem::SphericalCoords> ptg_radec;
    std::vector<VSAAlgebra::Vec2D>       m_ptg_xy;

    
    static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,VSAnalysisStage3Data,src_name);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,src_ra_rad);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,src_dec_rad);	
	H5_ADDMEMBER(c,VSAnalysisStage3Data,src_ra_hms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,src_dec_dms);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,pangle_mean_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,pangle_rms_deg);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,livetime_ticks);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,elaptime_ticks);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,livetime_min);
	H5_ADDMEMBER(c,VSAnalysisStage3Data,elaptime_min);
      }

  private:
    void copy(const VSAnalysisStage3Data& o);

    VSAcceptanceData*                    m_acceptance_data;
    std::vector<RunData>                 m_run_data;

  };

  class VSAnalysisStage3SimData
  {
  public:
    VSAnalysisStage3SimData() { }
    ~VSAnalysisStage3SimData() { }


  };

}

#endif // VSANALYSISSTAGE3DATA_HPP
