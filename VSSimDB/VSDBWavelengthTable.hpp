//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBWavelengthTable.hpp
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
  typedef std::map<unsigned, float>              VSDBWavelengthData;
  typedef std::map<unsigned, float>              VSDBAltitudeData;
  typedef std::map<unsigned, VSDBAltitudeData>   VSDBWavelengthAltitudeData;
  
  class VSDBWavelengthDataset
  {
  public:
    std::string comment;
    VSDBWavelengthData data;

    bool writeToCORSIKA(const std::string& filename) const;
    static VSDBWavelengthDataset* 
    createFromCORSIKA(const std::string& filename);
  };

  class VSDBWavelengthAltitudeDataset
  {
  public:
    std::string comment;
    VSDBWavelengthAltitudeData data;

    bool writeToCORSIKA(const std::string& filename) const;
    VSDBWavelengthAltitudeDataset* 
    createFromCORSIKA(const std::string& filename);
  };

  class VSDBModtranProfileDatum
  {
  public:
    float altitude;
    float rho;
    float thick;
    float n_minus_one;
  };

  typedef std::vector<VSDBModtranProfileDatum>    VSDBModtranProfileData;

  class VSDBModtranProfileDataset
  {
  public:
    std::string comment;
    VSDBModtranProfileData data;

    bool writeToCORSIKA(const std::string& filename) const;
    static VSDBModtranProfileDataset* 
    createFromCORSIKA(const std::string& filename);
  };

  class VSDBCORSIKADatasets
  {
  public:
    VSDBCORSIKADatasets(VSDatabase* db): fDB(db) { }
    virtual ~VSDBCORSIKADatasets();

    // WAVELENGTH DATA

    void storeWavelengthData(const std::string& table_name,
			     const VSDBWavelengthData& data);
    VSDBWavelengthData retrieveWavelengthData(const std::string& table_name);

    void storeWavelengthDataset(const std::string& table_name,
				const VSDBWavelengthDataset& dataset);
    VSDBWavelengthDataset 
    retrieveWavelengthDataset(const std::string& table_name);

    // WAVELENGTH-ALTITUDE DATA

    void storeWavelengthAltitudeData(const std::string& table_name,
				     const VSDBWavelengthAltitudeData& data);
    VSDBWavelengthAltitudeData 
    retrieveWavelengthAltitudeData(const std::string& table_name);    

    void storeWavelengthAltitudeDataset(const std::string& table_name,
				const VSDBWavelengthAltitudeDataset& dataset);
    VSDBWavelengthAltitudeDataset 
    retrieveWavelengthAltitudeDataset(const std::string& table_name);    

    // MODTRAN ATMOSPHERIC PROFILE

    void storeModtranProfileData(const std::string& table_name,
				 const VSDBModtranProfileData& data);
    VSDBModtranProfileData 
    retrieveModtranProfileData(const std::string& table_name);    

    void storeModtranProfileDataset(const std::string& table_name,
				const VSDBModtranProfileDataset& dataset);
    VSDBModtranProfileDataset 
    retrieveModtranProfileDataset(const std::string& table_name);    

  private:
    VSDatabase* fDB;
  };


}

#endif // VSWAVELENGTHTABLE_HPP
