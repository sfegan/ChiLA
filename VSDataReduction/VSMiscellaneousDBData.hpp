//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMiscellaneousDBData.hpp

  Various data from the run that comes from the database but is not
  strictly necessary for the analaysis. (Happy birthday Eoin)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       02/24/2007 

  $Id: VSMiscellaneousDBData.hpp,v 3.6 2008/04/14 22:29:59 sfegan Exp $

*/

#ifndef VSMISCELLANEOUSDBDATA_HPP
#define VSMISCELLANEOUSDBDATA_HPP

#include <vector>

#include <VSCentralizedDBAccess.hpp>
#include <VSOctaveIO.hpp>

namespace VERITAS
{

  class VSMiscellaneousDBData
  {
  public:
    VSMiscellaneousDBData(): scope(), fir(), target_table()
    { /* nothing to see here */ }
    
    // ========================================================================
    // CLASS DEFINITIONS
    // ========================================================================

    // Tracking target --------------------------------------------------------

    typedef VSCentralizedDBAccess::TrackingTargetDatum TrackingTargetDatum;

    // Correction parameters --------------------------------------------------

    typedef VSCentralizedDBAccess::CorrectionParametersDatum 
    CorrectionParametersDatum;

    // VPM data ---------------------------------------------------------------

    typedef VSCentralizedDBAccess::VPMCentroids VPMCentroids;
    typedef VSCentralizedDBAccess::VPMLEDs VPMLEDs;

    // HV Status data ---------------------------------------------------------

    struct OneHVMeasChan
    {
      OneHVMeasChan(): voltage(), current() { /* nothing to see here */ }
      float                            voltage;
      float                            current;      
    };

    struct OneHVMeasScope
    {
      OneHVMeasScope(unsigned nchan): timestamp(), chan(nchan) 
      { /* nothing to see here */ }
      VATime                           timestamp;
      std::vector<OneHVMeasChan>       chan;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,OneHVMeasScope,timestamp);
      }
    };

    // L1 Rate data -----------------------------------------------------------

    struct OneL1RateMeasScope
    {
      OneL1RateMeasScope(unsigned nchan): timestamp(), rate(nchan) 
      { /* nothing to see here */ }
      VATime                           timestamp;
      std::vector<float>               rate;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDSIMPLECOMPOSITE(c,OneL1RateMeasScope,timestamp);
      }
    };

    // L3 Telescope data ------------------------------------------------------

    typedef VSCentralizedDBAccess::L3ScopeDatum L3ScopeDatum;

    // Per telescope data -----------------------------------------------------

    class Scope
    {
    public:
      Scope(): 
	has_scope(), tracking_targets(), hv_status(), l1_rate(),
	l3_scope(), correction_parameters(), vpm_stars(), vpm_leds()
      { /* nothing to see here */ }

      bool getCorrectionParameters(const VSTime& time,
				   CorrectionParametersDatum& cpd) const;

      bool getNextCorrectionParameters(const VSTime& time,
				       CorrectionParametersDatum& cpd) const;

      bool                                   has_scope;
      std::vector<TrackingTargetDatum>       tracking_targets;
      std::vector<OneHVMeasScope>            hv_status;
      std::vector<OneL1RateMeasScope>        l1_rate;
      std::vector<L3ScopeDatum>              l3_scope;
      std::vector<CorrectionParametersDatum> correction_parameters;
      std::vector<VPMCentroids>              vpm_stars;
      std::vector<VPMLEDs>                   vpm_leds;
    };

    // FIR data ---------------------------------------------------------------

    typedef VSCentralizedDBAccess::FIRDatum FIRDatum;
    typedef std::vector<FIRDatum> FIRDataSet;

    // Target table entries ---------------------------------------------------

    typedef VSCentralizedDBAccess::TargetTableCoord TargetTableCoord;

    // ========================================================================
    // DATA
    // ========================================================================

    std::vector<Scope>                       scope;
    std::vector<FIRDataSet>                  fir;
    std::vector<TargetTableCoord>            target_table;

    // ========================================================================
    // MEMBER FUNCTIONS
    // ========================================================================

    void clear();
    void fillFromDB(const VSTime& start_time, const VSTime& stop_time,
		    const std::vector<unsigned>& telescopes);
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;
  };

}

#endif // VSMISCELLANEOUSDBDATA_HPP
