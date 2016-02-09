//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEnergyCalcLT.hpp

  Class which calculates reconstructed energy.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/09/2008

  $Id: VSEnergyCalcLT.hpp,v 3.11 2010/10/20 20:16:30 matthew Exp $

*/

#ifndef VSENERGYCALCLT_HPP
#define VSENERGYCALCLT_HPP

#include<vector>
#include<memory>
#include<cmath>
#include<ostream>

#include<VSOptions.hpp>
#include<VSEventData.hpp>
#include<VSNSpace.hpp>
#include<VSLTLibraryCalc.hpp>

namespace VERITAS
{

  //! Class that performs energy reconstruction on invidual stage2
  //! events.
  class VSEnergyCalcLT : public VSLTLibraryCalc
  {
  public:

    struct Options
    {
      Options():
	energy_lookup_file(""), 
	energy_weight_power(1.0), energy_rms_cut(0.0),
	energy_quality_cuts(1.3,250,
			    std::numeric_limits<double>::infinity(),
			    std::numeric_limits<double>::infinity()) 
      { }

      typedef quad<double,double,double,double> EnergyQualityCuts;

      std::string              energy_lookup_file;
      double                   energy_weight_power;
      double                   energy_rms_cut;
      EnergyQualityCuts        energy_quality_cuts;
    };

    enum EnergyEstimator
      {
	EE_N_R,
	EE_N_R_DISP,
	EE_LAMBDAD_R,
	EE_LAMBDAD_G,
	EE_LAMBDAD_ETA,
	EE_LAMBDAD_G_ETA
      };

    struct State
    {
    public:
      State(): 
	sum_e(), sum_w(), sum_1_sigma2(), sum_e_sigma2(), sum_e2_sigma2() { }

      void reset()
      {
	sum_e         = 0.0;
	sum_w         = 0.0;
	sum_1_sigma2  = 0.0;
	sum_e_sigma2  = 0.0;
	sum_e2_sigma2 = 0.0;
      }

      double                   sum_e;
      double                   sum_w;
      double                   sum_1_sigma2;
      double                   sum_e_sigma2;
      double                   sum_e2_sigma2;
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

    VSEnergyCalcLT(const Options& opt = s_default_options);

    VSEnergyCalcLT(const std::vector<unsigned>& nchan,
		   double zn_deg, double az_deg, double ped_rms,
		   std::ostream* stream = 0,
		   const Options& opt = s_default_options);
    
    ~VSEnergyCalcLT();

    //! Get the telescope energy estimator.
    //! @param iscope           Telescope index
    //! @param R                Telescope impact distance
    //! @param N                Telescope image size
    //! @param fp_disp_deg      Telescope image displacement in degrees
    //! @param fp_dist_deg      Telescope image distance in degrees
    //! @param scope_log10_e    Estimated energy for this telescope
    //! @param scope_log10_e_rms Energy estimate dispersion for this
    //! telescope
    bool getScopeEnergy(State& state, 
			unsigned iscope, 
			double R, double N, double fp_disp_deg, 
			double fp_dist_deg,
			double& scope_log10_e, double& scope_log10_e_rms);

    //! Get the telescope energy estimator.
    bool getScopeEnergy(State& state, 
			unsigned iscope, VSEventScopeDatum* sd,
			double& scope_log10_e, double& scope_log10_e_rms);

    //! Get the array energy estimator.
    bool getArrayEnergy(State& state, 
			double& array_log10_e, double& array_log10_e_chi2);
    
    static void point(VSNSpace::Point& p, VSEventScopeDatum* sd, 
		      EnergyEstimator ee);

    //! Calculate the array and telescope energies for this event.
    //! @param event_data stage2 event data structure
    bool calcEnergy(VSEventArrayDatum& event_data);
    
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
    VSEnergyCalcLT(const VSEnergyCalcLT& o);

    static void configure(VSOptions& options,
			  const std::string& profile="", 
			  const std::string& opt_prefix="s2_");

  private:

    VSEnergyCalcLT& operator=(const VSEnergyCalcLT& o);

    Options                        m_options;

    std::vector<SCParameterPair>   m_lt;
    EnergyEstimator                m_egy_estimator;

    static Options                 s_default_options;
  };

  typedef VSEnergyCalcLT VSEnergyCalc;
  
}

#endif // VSENERGYCALCLT_HPP
