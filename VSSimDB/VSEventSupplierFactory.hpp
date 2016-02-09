//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierFactory.hpp

  Factory for creating event suppliers

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       19/07/2006
  \note
*/

// g++ -I../Physics -I/usr/include/mysql -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -o test test.cpp -L. -L../VSUtility -L/usr/lib/mysql -L../Physics -I. -I../VSUtility -lVSSimDB -lVSUtility -lPhysics -lmysqlclient -lhdf5

#ifndef VSEVENTSUPPLIERFACTORY_HPP
#define VSEVENTSUPPLIERFACTORY_HPP

#include <vector>
#include <stdint.h>

#include <RandomNumbers.hpp>
#include <VSSimDB.hpp>

#include "VSEventSupplier.hpp"

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierFactory
  {
  public:
    VSEventSupplierFactory(VSSimDB* sim_db = 0, RandomNumbers* rng = 0,
                           const std::string& dir = "");
    ~VSEventSupplierFactory();

    // Empty (pedestal) supplier ----------------------------------------------
    VSEventSupplier* newEmptySupplier(unsigned nevents,
				      unsigned nscope = 4, 
				      unsigned nchan = 500);

    // Laser supplier ---------------------------------------------------------
    VSEventSupplier* newLaserSupplier(unsigned nevents,
				      double mean_pe_per_pixel = 100.0,
				      double mean_pe_per_pixel_dev = 10.0,
				      double mean_pulse_fwhm_ns = 5.0,
				      unsigned nscope = 4, 
				      unsigned nchan = 500);

    // One table supplier -----------------------------------------------------
    VSEventSupplier* newOneTableSupplier(const std::string& table_name);

    VSEventSupplier* newOneTableSupplier(const std::string& table_name,
					 uint32_t nevent_to_supply, 
					 double ievent_to_start_fraction);

    // All events supplier ----------------------------------------------------
    VSEventSupplier* 
    newAllEventsSupplier(bool randomize_if_available = true);


    // Supply events weighted to high energies --------------------------------
    VSEventSupplier* newHEWeightedSupplier(uint32_t nevent_to_supply,
					   bool randomize_if_available = true);

    // Supply events with a monoenergetic distribution ------------------------
    VSEventSupplier* newMonoenergeticSupplier(double target_energy_gev,
					      uint32_t nevent_to_supply = 0,
					      bool randomize_if_available = 
					      true);

    // Supply events with a power-law distribution ----------------------------
    VSEventSupplier* newPLSpectrumSupplier(double spectral_index,
					   uint32_t nevent_to_supply = 0,
					   bool randomize_if_available = true);
    
    // Legacy methods (for backwards compatibility) ---------------------------
    VSEventSupplier* 
    newPLSpectrumClosestZenithSupplier(double spectral_index,
				       double min_energy_gev,
				       double max_energy_gev,
				       double zn_rad,
				       uint32_t nevent_to_supply = 0,
				       bool randomize_if_available = true);

    VSEventSupplier* 
    newPLSpectrumClosestPointingSupplier(double spectral_index,
					 double min_energy_gev,
					 double max_energy_gev,
					 double zn_rad,
					 double az_rad,
					 uint32_t nevent_to_supply = 0,
					 bool randomize_if_available = true);

    VSEventSupplier* 
    newAllTablesClosestZenithSupplier(double min_energy_gev,
				      double max_energy_gev,
				      double zn_rad, 
				      uint32_t nevent_to_supply,
				      bool randomize_if_available = true);

    VSEventSupplier* 
    newAllTablesClosestPointingSupplier(double min_energy_gev,
					double max_energy_gev,
					double zn_rad, 
					double az_rad,
					uint32_t nevent_to_supply,
					bool randomize_if_available = true);

    // Methods to configure the factory ---------------------------------------
    void setCheckHDFEventCount(bool check_hdf_event_count)
    {
      fCheckHDFEventCount = check_hdf_event_count;
    }
    
    void setEnergyRange(double min_energy_gev, double max_energy_gev)
    {
      fMinEnergyGeV = min_energy_gev;
      fMaxEnergyGeV = max_energy_gev;
    }

    void setTargetZenith(double target_zenith_rad)
    {
      fTargetZenithRad = target_zenith_rad;
      fHasTargetZenith = true;
    }

    void setTargetAzimuth(double target_azimuth_rad)
    {
      fTargetAzimuthRad = target_azimuth_rad;
      fHasTargetAzimuth = true;
    }

  private:
    
    void selectByEnergy(std::vector<VSSimDBTableParam>& tables);
    void selectByZenith(std::vector<VSSimDBTableParam>& tables);
    void selectByAzimuth(std::vector<VSSimDBTableParam>& tables);

    double angleDistance(double target_angle, 
			 double angle_range_min, double angle_range_max);

    double findClosestEnergy(double target_energy_gev,
			     const std::vector<VSSimDBTableParam>& tables);
    double findMinZenith(double target_angle,
			 const std::vector<VSSimDBTableParam>& table_params);
    double findMinAzimuth(double target_angle,
			  double zn_angle,
			  const std::vector<VSSimDBTableParam>& table_params);

    void findPLWeighting(double spectral_index,
			 uint32_t nevent_to_supply,
			 bool randomize_if_available,
			 std::vector<VSSimDBTableParam> &table_params,
			 std::vector< std::pair< std::string, uint32_t > >& 
			 tables);
    
    void findHEWeighting(uint32_t nevent_to_supply,
			 bool randomize_if_available,
			 std::vector<VSSimDBTableParam> &table_params,
			 std::vector< std::pair< std::string, uint32_t > >& 
			 tables);


    VSSimDB*             fSimDB;
    RandomNumbers*       fRNG;

    std::string          fOverrideDataDirectory;
    uint32_t             fMaxBufferSize;
    bool                 fCheckHDFEventCount;
    double               fMinEnergyGeV;
    double               fMaxEnergyGeV;
    double               fTargetZenithRad;
    double               fTargetAzimuthRad;
    bool                 fHasTargetZenith;
    bool                 fHasTargetAzimuth;
  };

} // namespace VERITAS

#endif // VSEVENTSUPPLIERFACTORY_HPP
