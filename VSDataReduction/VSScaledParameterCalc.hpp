//-*-mode:c++; mode:font-lock;-*-

/*! \file VSScaledParameterCalc.hpp

  Class which calculates mean scaled parameters

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/09/2008

  $Id: VSScaledParameterCalc.hpp,v 3.10 2010/06/22 00:00:32 matthew Exp $

*/

#ifndef VSSCALEDPARAMETERCALC_HPP
#define VSSCALEDPARAMETERCALC_HPP

#include<vector>
#include<memory>
#include<cmath>
#include<iostream>

#include<VSOptions.hpp>
#include<VSEventData.hpp>
#include<VSNSpace.hpp>
#include<VSLTLibraryCalc.hpp>

namespace VERITAS
{

  //! Class responsible for calculating mean scaled event parameters
  //! (width,length,disp).  The relative weighting of the scaled
  //! parameters is controlled by the sp_weight_power option.
  class VSScaledParameterCalc : public VSLTLibraryCalc
  {
  public:

    struct Options
    {
      Options():
	sp_lookup_file(""), 
	sp_weight_power(1.0)
	//	sp_quality_cuts(0,std::numeric_limits<double>::infinity(),
	//			std::numeric_limits<double>::infinity(),
	//			std::numeric_limits<double>::infinity()) 
      { }

      typedef quad<double,double,double,double> EnergyQualityCuts;

      std::string              sp_lookup_file;
      double                   sp_weight_power;
      //      double                   sp_rms_cut;
      //      EnergyQualityCuts        sp_quality_cuts;
    };

    struct State
    {
      State(): sum_width_w(), sum_width(), sum_length_w(), sum_length(),
	       sum_disp_w(), sum_disp() { }

      void reset()
      {
	sum_width_w   = 0;
	sum_width     = 0;
	sum_length_w  = 0;
	sum_length    = 0;
	sum_disp_w    = 0;
	sum_disp      = 0;
      }
      
      double sum_width_w;
      double sum_width;
      double sum_length_w;
      double sum_length;
      double sum_disp_w;
      double sum_disp;
    };

    struct SCParameterPair
    {
      SCParameterPair(): exp(), rms(), mask() { }
      SCParameterPair(VSNSpace* _exp, VSNSpace* _rms, VSNSpace* _mask):
	exp(_exp), rms(_rms), mask(_mask)
      { }

      VSNSpace* exp;
      VSNSpace* rms;
      VSNSpace* mask;
    };

    VSScaledParameterCalc(const Options& opt = s_default_options);

    VSScaledParameterCalc(const std::vector<unsigned>& nchan,
			  double zn_deg, double az_deg, double ped_rms,
			  std::ostream* stream = 0,
			  const Options& opt = s_default_options);

    virtual ~VSScaledParameterCalc();

    bool getScopeSP(State& state, unsigned iscope, 
		    double R, double N, 
		    double fp_width_deg, double fp_length_deg, 
		    double fp_disp_deg, 
		    double fp_dist_deg,
		    double& scope_sc_width, double& scope_sc_length, 
		    double& scope_sc_disp);

    bool getArraySP(State& state, 
		    double& array_msc_width, double& array_msc_length, 
		    double& array_msc_disp);
    
    bool calcSP(VSEventArrayDatum& event_data)
    {
      if(!m_has_lt) return false;

      State state;
      unsigned nscope = event_data.scope.size();
      for(unsigned iscope=0;iscope<nscope; iscope++)
	{
	  VSEventScopeDatum* scope_data = event_data.scope[iscope];

	  bool used_in_reconstruction = 
	    ((0x1<<iscope)&event_data.used_in_reconstruction_mask);

	  if(scope_data && used_in_reconstruction)
	    getScopeSP(state, iscope,
		       scope_data->R, scope_data->N,
		       scope_data->intrinsic_width, 
		       scope_data->intrinsic_length, 
		       scope_data->fp_disp,
		       scope_data->fp_dist,
		       scope_data->sc_width, scope_data->sc_length, 
		       scope_data->sc_disp);
	}
      return getArraySP(state, 
			event_data.msc_width, event_data.msc_length, 
			event_data.msc_disp);
    }
    
    // Load/Save --------------------------------------------------------------
    virtual void clear();

    virtual bool load(const std::string& lookup_file,
		      const std::vector<unsigned>& nchan, 
		      double zn_deg, double az_deg, double ped_rms);

    bool load(const std::vector<unsigned>& nchan, double zn_deg, 
	      double az_deg, double ped_rms);

    bool load(const VSTime& date,
	      const std::vector<unsigned>& nchan, double zn_deg, 
	      double az_deg, double ped_rms);
    bool load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    // Copy Constructor -------------------------------------------------------
    VSScaledParameterCalc(const VSScaledParameterCalc& o);

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="s2_");

  private:

    VSScaledParameterCalc& operator=(const VSScaledParameterCalc& o);
    
    inline bool setScaledParameter(const SCParameterPair& sc_tables,
				   const VSNSpace::Point& p, double data,
				   double& sc_val, double& msc_w, 
				   double& msc_sum);
      
    Options                        m_options;
    

    std::vector<SCParameterPair>   m_lt_width;
    std::vector<SCParameterPair>   m_lt_length;
    std::vector<SCParameterPair>   m_lt_disp;

    static Options                 s_default_options;
  };

  // ==========================================================================
  // Scaled parameter helper
  // ==========================================================================

  inline bool VSScaledParameterCalc::
  setScaledParameter(const SCParameterPair& sc_tables,
		     const VSNSpace::Point& p, double data,
		     double& sc_val, double& sum_w, double& sum)
  {
    if(sc_tables.exp)
      {
	double expected = 0;
	double rms = 0;
#if 0
	if((sc_tables.exp->getWeight(p, expected))&&
	   (sc_tables.rms->getWeight(p, rms))&&(rms>0))
#endif
	if((sc_tables.exp->getWeight(p, expected))&&
	   (sc_tables.rms->getWeight(p, rms))&&(rms>0))
	  {
	    sc_val = (data-expected)/rms;

	    if(std::isfinite(m_options.sp_weight_power))
	      {
		double w = 1.0;
		if(m_options.sp_weight_power > 0)
		  w = pow(rms,-m_options.sp_weight_power);
		sum_w += w;
		sum += w*sc_val;
	      }
	    else if((sum_w == 0)||(rms < sum_w))
	      {
		sum_w = rms;
		sum = rms*sc_val;
	      }
	  }
	else
	  {
	    sc_val = std::numeric_limits<double>::infinity();
	    return false;
	  }
      }
    else 
      {
	sc_val = 0;
	return true; // NB!
      }

    return true;
  }
  
}

#endif // VSSCALEDPARAMETERCALC_HPP
