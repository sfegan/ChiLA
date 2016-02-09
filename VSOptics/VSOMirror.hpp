//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOMirror.hpp
  Mirror class header file

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

#ifndef VSOMIRROR_HPP
#define VSOMIRROR_HPP

#include <iostream>

#include <Vec3D.hpp>
#include <Particle.hpp>

#include <VSDatabase.hpp>

namespace VERITAS
{
  class VSOTelescope;
  
  //! VSOMirror class for telescope, stores information about individual mirror
  class VSOMirror 
  {
  public:
    // ************************************************************************
    // Constructor and Destructor
    // ************************************************************************
    VSOMirror();
    VSOMirror(const VSOTelescope* T, unsigned ID, unsigned MHID, bool REM,
	      const Physics::Vec3D& P, const Physics::Vec3D& A, 
	      double FL, double SS, double DF);
    virtual ~VSOMirror();
    
    // ************************************************************************
    // Write to and Read from the database
    // ************************************************************************
    static VSDBStatement* createInsertQuery(VSDatabase* db);
    static VSDBStatement* createSelectQuery(VSDatabase* db);
    void writeToDatabase(VSDBStatement* stmt, uint32_t optics_id) const;
    static VSOMirror* createFromDatabaseRow(VSDBStatement* stmt, 
					    const VSOTelescope* T);
    static void createMirrorTable(VSDatabase* db);
    
    // ************************************************************************
    // Coordinate transformations
    // ************************************************************************
    void reflectorToMirror(Physics::Vec3D& v) const;    //!< Transform vector from reflector to mirror
    void mirrorToReflector(Physics::Vec3D& v) const ;   //!< Transform vector from mirror to reflector
    void reflectorToMirror(Physics::Particle& p) const; //!< Transform particle reflector to mirror
    void mirrorToReflector(Physics::Particle& p) const; //!< Transform particle mirror to reflector
    
    // ************************************************************************
    // Dump
    // ************************************************************************
    void dumpShort(std::ostream& stream) const;
    void dump(std::ostream& stream, unsigned l=0) const;
    static VSOMirror* createFromShortDump(std::istream& stream,
					  const VSOTelescope* T);

    // ************************************************************************
    // Accessors
    // ************************************************************************
    const VSOTelescope* telescope() const { return fTelescope; }
    unsigned           id() const { return fID; }
    unsigned           hexID() const { return fHexID; }
    bool               removed() const { return fRemoved; }
    const Physics::Vec3D& pos() const { return fPos; }
    const Physics::Vec3D& align() const { return fAlign; }
    double             focalLength() const { return fFocalLength; }
    double             spotSize() const { return fSpotSize; }
    double             degradingFactor() const { return fDegradingFactor; }
    
  private:
    const VSOTelescope* fTelescope;         //!< Telescope
    unsigned            fID;                //!< Sequential ID (starting at 0)
    unsigned            fHexID;             //!< Hex index on reflector
    bool                fRemoved;           //!< Mirror is removed from scope
    Physics::Vec3D      fPos;               //!< Position
    Physics::Vec3D      fAlign;             //!< Alignment angles
    double              fFocalLength;       //!< Focal length
    double              fSpotSize;          //!< Spot size
    double              fDegradingFactor;   //!< Degrading factor of mirror
    
    Physics::Vec3D      fRotationVector;
    void calculateRotationVector();
  };
}

#endif // VSOMIRROR_H
