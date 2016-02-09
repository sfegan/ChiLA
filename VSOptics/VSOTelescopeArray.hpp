//-*-mode:c++; mode:font-lock;-*-

/*! \file OSTelescopeArray.hpp
  Array class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz             \n
           UCLA                        \n
	   nicewicz@physics.ucla.edu   \n

  \date    11/30/2004
  \version 0.2
  \note
*/

#ifndef VSOTELESCOPEARRAY_HPP
#define VSOTELESCOPEARRAY_HPP

#include <string>

#include <Vec3D.hpp>

#include <VSDatabase.hpp>

#include "VSOArrayParameters.hpp"
#include "VSOTelescope.hpp"

namespace VERITAS
{

  class VSOTelescopeArray
  {
  public:
    VSOTelescopeArray();
#if 0
    VSOTelescopeArray(Database *db);
#endif
    virtual ~VSOTelescopeArray();
    
   // ************************************************************************
    // Read and Write to the database
    // ************************************************************************
    bool readFromDatabase(VSDatabase* db, uint32_t optics_id);
    void writeToDatabase(VSDatabase* db, uint32_t optics_id) const;
    static void createTelescopeArrayTable(VSDatabase* db);
    
    // ************************************************************************
    // Create a new array randomly using parameters
    // ************************************************************************
    void generateFromArrayParameters(const VSOArrayParameters& param, 
				     RandomNumbers& rng);
    
    // ************************************************************************
    // Other stuff
    // ************************************************************************
    bool pointTelescopesAzEl(const double az_rad, const double el_rad);
    bool pointTelescopes(const Physics::Vec3D& v);

    // ************************************************************************
    // Accessors
    // ************************************************************************
    double latitude() const { return fLatitude; }
    double longitude() const { return fLatitude; }
    double altitude() const { return fAltitude; }
    double spacing() const { return fSpacing; }
    double arrayParity() const { return fArrayParity; }
    
    unsigned numTelescopes() const { return fTelescopes.size(); }
    unsigned numTelescopeHexSites() const { return fTelescopesByHexID.size(); }
    
    inline const VSOTelescope* telescope(unsigned id) const;
    inline const VSOTelescope* telescopeByHexID(unsigned hexID) const;
    
    inline VSOTelescope* telescope(unsigned id);
    inline VSOTelescope* telescopeByHexID(unsigned hexID);
    
    // ************************************************************************
    // Dump
    // ************************************************************************
    void dump(std::ostream& stream, unsigned l=0) const;
    void dumpShort(const std::string& filename) const;
    void dumpShort(std::ostream& stream) const;
    void writeToShortDump(const std::string& filename) const { dumpShort(filename); }
    void writeToShortDump(std::ostream& stream) const { dumpShort(stream); }
    bool readFromShortDump(std::istream& stream);
    bool readFromShortDump(const std::string& filename);
    
  private:
    double                     fLatitude;
    double                     fLongitude;
    double                     fAltitude;
    double                     fSpacing;
    bool                       fArrayParity;
    
    std::vector<VSOTelescope*>  fTelescopes;
    std::vector<VSOTelescope*>  fTelescopesByHexID;
  };

  inline const VSOTelescope* VSOTelescopeArray::telescope(unsigned id) const
  {
    if(id>=fTelescopes.size())return 0;
    else return fTelescopes[id];
  }
  
  inline const VSOTelescope* VSOTelescopeArray::telescopeByHexID(unsigned hexID) const
  {
    if((hexID<1)||(hexID>fTelescopes.size()))return 0;
    else return fTelescopes[hexID-1];
  }
  
  inline VSOTelescope* VSOTelescopeArray::telescope(unsigned id)
  {
    if(id>=fTelescopes.size())return 0;
    else return fTelescopes[id];
  }
  
  inline VSOTelescope* VSOTelescopeArray::telescopeByHexID(unsigned hexID)
  {
    if((hexID<1)||(hexID>fTelescopes.size()))return 0;
    else return fTelescopes[hexID-1];
  }

}

#endif // VSOTELESCOPEARRAY_HPP
