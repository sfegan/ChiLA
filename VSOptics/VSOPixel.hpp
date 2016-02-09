//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOPixel.hpp
  Pixel class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz
           UCLA
	   nicewicz@physics.ucla.edu

  \date    11/30/2004
  \version 0.2
  \note
*/

#ifndef VSOPIXEL_HPP
#define VSOPIXEL_HPP

#include <iostream>

#include <Vec3D.hpp>

#include <VSDatabase.hpp>

namespace VERITAS 
{

  class VSOTelescope;
  
  //! Class for an individual pixel in the camera
  class VSOPixel
  {
  public:
    // ************************************************************************
    // Constructor and Destructor
    // ************************************************************************
    VSOPixel();     //!<default constructor
    VSOPixel(const VSOTelescope* T, unsigned PID, unsigned PHID, bool REM,
	     const Physics::Vec3D& P);
    virtual ~VSOPixel();
        
    // ************************************************************************
    // Write to and Read from the database
    // ************************************************************************
    static VSDBStatement* createInsertQuery(VSDatabase* db);
    static VSDBStatement* createSelectQuery(VSDatabase* db);
    void writeToDatabase(VSDBStatement* stmt, uint32_t optics_id) const;
    static VSOPixel* createFromDatabaseRow(VSDBStatement* stmt, 
					   const VSOTelescope* T);
    static void createPixelTable(VSDatabase* db);
    
    // ************************************************************************
    // Dump
    // ************************************************************************
    void dumpShort(std::ostream& stream) const;
    void dump(std::ostream& stream, unsigned l=0) const;
    static VSOPixel* createFromShortDump(std::istream& stream,
					 const VSOTelescope* T);
    
    // ************************************************************************
    // Accessors
    // ************************************************************************
    const VSOTelescope* telescope() const { return fTelescope; }
    unsigned            id() const { return fID; }
    unsigned            hexID() const { return fHexID; }
    const Physics::Vec3D& pos() const { return fPos; }

    Physics::Vec3D      incomingSkyVectorAtZenith(double plate_scale=1.0) 
      const;
    
  private:
    const VSOTelescope* fTelescope;         //!< Telescope 
    unsigned            fID;                //!< Hex ID
    unsigned            fHexID;             //!< Hex ID
    bool                fRemoved;           //!< Pixel removed from camera
    Physics::Vec3D      fPos;               //!< Position
  };

}
#endif // VSOPIXEL_HPP
