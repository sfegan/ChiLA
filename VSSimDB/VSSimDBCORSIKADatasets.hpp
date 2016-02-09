//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimDBWavelengthTable.hpp
  Database access to wavelength/value tables

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       08/19/2005
*/

#ifndef VSWAVELENGTHTABLE_HPP
#define VSWAVELENGTHTABLE_HPP

#include <map>
#include <vector>
#include <set>
#include <string>

#include "VSDatabase.hpp"

//! VERITAS namespace
namespace VERITAS
{
  typedef std::map<unsigned, float>              VSSimDBWavelengthData;
  typedef std::map<unsigned, float>              VSSimDBAltitudeData;
  typedef std::map<unsigned,VSSimDBAltitudeData> VSSimDBWavelengthAltitudeData;
  
  class VSSimDBWavelengthDataset
  {
  public:
    VSSimDBWavelengthDataset(): comment(), data() 
    { /* nothing to see here */ }

    std::string comment;
    VSSimDBWavelengthData data;

    bool writeToCORSIKA(const std::string& filename) const;
    static VSSimDBWavelengthDataset* 
    createFromCORSIKA(const std::string& filename);
  };

  class VSSimDBWavelengthAltitudeDataset
  {
  public:
    VSSimDBWavelengthAltitudeDataset(): comment(), data() 
    { /* nothing to see here */ }

    std::string comment;
    VSSimDBWavelengthAltitudeData data;

    bool writeToCORSIKA(const std::string& filename) const;
    static VSSimDBWavelengthAltitudeDataset* 
    createFromCORSIKA(const std::string& filename);
  };

  class VSSimDBModtranProfileDatum
  {
  public:
    VSSimDBModtranProfileDatum():
      altitude(), rho(), thick(), n_minus_one()
    { /* nothing to see here */ }

    float altitude;
    float rho;
    float thick;
    float n_minus_one;
    bool operator< (const VSSimDBModtranProfileDatum& o) const 
    { return altitude < o.altitude; }
  };

  typedef std::vector<VSSimDBModtranProfileDatum>    VSSimDBModtranProfileData;

  class VSSimDBModtranProfileDataset
  {
  public:
    VSSimDBModtranProfileDataset(): comment(), data() 
    { /* nothing to see here */ }

    std::string comment;
    VSSimDBModtranProfileData data;

    bool writeToCORSIKA(const std::string& filename) const;
    static VSSimDBModtranProfileDataset* 
    createFromCORSIKA(const std::string& filename);
  };

  class VSSimDBCORSIKADatasets
  {
  public:
    VSSimDBCORSIKADatasets(VSDatabase* db): fDB(db) { }
    virtual ~VSSimDBCORSIKADatasets();

    // WAVELENGTH DATA

    void storeWavelengthData(const std::string& table_name,
			     const VSSimDBWavelengthData& data);
    VSSimDBWavelengthData* 
    retrieveWavelengthData(const std::string& table_name);

    void storeWavelengthDataset(const std::string& table_name,
				const VSSimDBWavelengthDataset& dataset);
    VSSimDBWavelengthDataset* 
    retrieveWavelengthDataset(const std::string& table_name);

    // WAVELENGTH-ALTITUDE DATA

    void storeWavelengthAltitudeData(const std::string& table_name,
				    const VSSimDBWavelengthAltitudeData& data);
    VSSimDBWavelengthAltitudeData* 
    retrieveWavelengthAltitudeData(const std::string& table_name);    

    void storeWavelengthAltitudeDataset(const std::string& table_name,
			      const VSSimDBWavelengthAltitudeDataset& dataset);
    VSSimDBWavelengthAltitudeDataset* 
    retrieveWavelengthAltitudeDataset(const std::string& table_name);    

    // MODTRAN ATMOSPHERIC PROFILE

    void storeModtranProfileData(const std::string& table_name,
				 const VSSimDBModtranProfileData& data);
    VSSimDBModtranProfileData* 
    retrieveModtranProfileData(const std::string& table_name);    

    void storeModtranProfileDataset(const std::string& table_name,
				const VSSimDBModtranProfileDataset& dataset);
    VSSimDBModtranProfileDataset* 
    retrieveModtranProfileDataset(const std::string& table_name);    

  private:
    VSDatabase* fDB;
  };


}

#endif // VSWAVELENGTHTABLE_HPP
