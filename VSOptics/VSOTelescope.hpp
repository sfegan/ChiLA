//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOTelescope.hpp
  Reflector class header file

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

#ifndef VSOTELESCOPE_HPP
#define VSOTELESCOPE_HPP

#include <iostream>
#include <string>
#include <vector>

#include <Vec3D.hpp>
#include <Particle.hpp>
#include <RandomNumbers.hpp>
#include <xytohex.h>

#include <VSDatabase.hpp>
#include <VSSimDBTables.hpp>

#include "VSOArrayParameters.hpp"
#include "VSOMirror.hpp"
#include "VSOPixel.hpp"

namespace VERITAS
{

  /*! \class VSOTelescope
    \brief VSOTelescope class;
  */
  
  class VSOTelescope
  {
  public:
    VSOTelescope();
    VSOTelescope(unsigned TID, unsigned THID, const Physics::Vec3D&P, 
		 double DY, double AX, double AY, double EL, double AZ,
		 const Physics::Vec3D& T, double CR, double A, 
		 double FSP, double FS, double RR, const Physics::Vec3D& AP,
		 unsigned HRN, double RIP, bool MP, 
		 const Physics::Vec3D& FPT, double CD, double FOV, double D,
		 double PS, double CSP, const Physics::Vec3D& FPR, 
		 double CIP, bool PP 
#if 0
		 , bool SEC, double RI, double RI1, double RI2, double C10, 
		 double C11, double C12, double C13, double C14, double C20,
		 double C21, double C22,double C23, double C24
#endif
		 );
    VSOTelescope(const VSOTelescope& o);
    virtual ~VSOTelescope();
    
    const VSOTelescope& operator =(const VSOTelescope& o);
    
    // ************************************************************************
    // Create a new telescope randomly using parameters
    // ************************************************************************

    void populateMirrorsAndPixelsRandom(const VSOArrayParameters&,
					RandomNumbers& rng);
    
    // ************************************************************************
    // Write to and Read from the database
    // ************************************************************************
    
    static VSDBStatement* createInsertQuery(VSDatabase* db);
    static VSDBStatement* createSelectQuery(VSDatabase* db);
    static VSOTelescope* createFromDatabaseRow(VSDBStatement* stmt);
    void populateMirrorsAndPixelsFromDatabase(VSDatabase* db,
					      uint32_t optics_id);
    void writeToDatabase(VSDatabase* db, VSDBStatement* stmt, 
			 uint32_t optics_id) const;
    static void createTelescopeTable(VSDatabase* db);

    // ************************************************************************
    // Telescope repointing
    // ************************************************************************
    
    //! Points telescope along axis 
    bool pointTelescope(const Physics::Vec3D& v); 
    bool pointTelescopeAzEl(const double az_rad, const double el_rad);

    //! Transform particle global to reflector
    void globalToReflector(Physics::Vec3D& v) const; 
    //! Transform particle reflector to global
    void reflectorToGlobal(Physics::Vec3D& v) const; 
    //! Transform particle global to reflector
    void globalToReflector(Physics::Particle& p) const; 
    //! Transform particle reflector to global
    void reflectorToGlobal(Physics::Particle& p) const; 
    //! Transform from focal plane to reflector
    void focalPlaneToReflector(Physics::Particle& p) const; 
    //! Transform from reflector to focal plane
    void reflectorToFocalPlane(Physics::Particle& p) const; 
    
    // ************************************************************************
    // Dump and Reload
    // ************************************************************************
    
    void dumpShort(std::ostream& stream) const;
    void dump(std::ostream& stream, unsigned l=0) const;
    static VSOTelescope* createFromShortDump(std::istream& stream);

    // ************************************************************************
    // Accessors
    // ************************************************************************

    unsigned       id() const { return fID; }
    unsigned       hexID() const { return fTelescopeHexID; }
    double         altitude() const { return fPos.z; }
    const Physics::Vec3D& pos() const { return fPos; }
    const Physics::Vec3D& position() const { return fPos; }
    double         deltaY() const { return fDeltaY; }
    double         alphaX() const { return fAlphaX; } 
    double         alphaY() const { return fAlphaY; }
    double         azimuth() const { return fAzimuth; }
    double         elevation() const { return fElevation; }
    Physics::Vec3D opticalAxis() const;
    
    const Physics::Vec3D& translation() const { return fTranslation; }
    double         curvatureRadius() const { return fCurvatureRadius; }
    double         aperture() const { return fAperture; }
    double         facetSpacing() const { return  fFacetSpacing; }
    double         facetSize() const { return  fFacetSize; }
    double         reflectorRotation() const { return fReflectorRotation; }
    const Physics::Vec3D& alignmentPoint() const { return fAlignmentPoint; }
    unsigned       mirrorHexRings() const { return fHexagonRingsN; }
    double         reflectorIP() const { return fReflectorIP; }
    bool           mirrorParity() const { return fMirrorParity; }
    
    const Physics::Vec3D& focalPlanePosition() const { return fFPTranslation; }
    double         cameraDiameter() const { return fCameraDiameter; }
    double         fov() const { return fFieldOfView; }
    double         cathodeDiameter() const { return fCathodeDiameter; }
    double         pixelSpacing() const { return fPixelSpacing; }
    double         concentratorSurvivalProb() const { return fConcSurvProb; }
    const Physics::Vec3D& focalPlaneRotion() const { return fFPRotation; }
    double         cameraIP() const { return fCameraIP; }
    bool           pixelParity() const { return fPixelParity; }
    
#if 0
    bool           hasSecondary() const { return fHasSecondary; }
    double         refractiveIndex() const { return fRefractiveIndex; }
    double         refractiveIndex1() const { return fRefractiveIndex1; }
    double         refractiveIndex2() const { return fRefractiveIndex2; }
    double         ce1Parameter0() const { return fCE1Parameter0; }
    double         ce1Parameter2() const { return fCE1Parameter2; }
    double         ce1Parameter3() const { return fCE1Parameter3; }
    double         ce1Parameter4() const { return fCE1Parameter4; }
    double         ce1Parameter5() const { return fCE1Parameter5; }
    double         ce2Parameter0() const { return fCE2Parameter0; }
    double         ce2Parameter2() const { return fCE2Parameter2; }
    double         ce2Parameter3() const { return fCE2Parameter3; }
    double         ce2Parameter4() const { return fCE2Parameter4; }
    double         ce2Parameter5() const { return fCE2Parameter5; }
#endif

    unsigned       numMirrors() const { return fMirrors.size(); }
    unsigned       numMirrorHexSites() const { return fMirrorsByHexID.size(); }
    
    unsigned       numPixels() const { return fPixels.size(); }
    unsigned       numPixelHexSites() const { return fPixelsByHexID.size(); }
    
    inline const VSOMirror* mirror(unsigned id) const;
    inline const VSOMirror* mirrorByHexID(unsigned hexID) const;
    
    inline const VSOPixel* pixel(unsigned id) const;
    inline const VSOPixel* pixelByHexID(unsigned hexID) const;
    
 private:
    // ************************************************************************
    // Telescope Parameters
    // ************************************************************************
    unsigned       fID;                //!< Sequential ID (starting at 0)
    unsigned       fTelescopeHexID;    //!< Hex ID on the honeycomb (from 1)
    Physics::Vec3D fPos;               //!< Position of the telescope
    double         fDeltaY;            //!< Nonorthogonality of rotation planes
    double         fAlphaX;            //!< Alignment angles(?)
    double         fAlphaY;            //!< Alignment angles(?)
    double         fElevation;         //!< Elevation of the telescope -- angle to the x-y plane
    double         fAzimuth;           //!< Azimuthal position -- angle made in the x-y plane with the y axis
    
    // ************************************************************************
    // Reflector Parameters
    // ************************************************************************
    Physics::Vec3D fTranslation;       //!< Vector in reflector r.f. to the intersection of the rotation axes
    double         fCurvatureRadius;   //!< Radius of curvature of reflector
    double         fAperture;          //!< Diameter of projection of reflector onto a plane
    double         fFacetSpacing;      //!< Size of mirror
    double         fFacetSize;         //!< Spacing between mirrors
    double         fReflectorRotation; //!< Rotation about the axis of reflector
    Physics::Vec3D fAlignmentPoint;    //!< Alignment of telescope
    unsigned       fHexagonRingsN;     //!< Number of hexagon rings of mirrors
    double         fReflectorIP;       //!< Diameter of sphere embedding reflector 
    bool           fMirrorParity;      //!< Parity for counting mirrors; 1 => counting toward +x-axis
    //!  in reflector r.f.; -1 => counting toward -x-axis
    
    // ************************************************************************
    // Camera Parameters
    // ************************************************************************
    Physics::Vec3D fFPTranslation;     
    double         fCameraDiameter;    //!< Dependent on the field of view
    double         fFieldOfView;       //!< Field of view
    double         fCathodeDiameter;   //!< Diameter of photocathode
    double         fPixelSpacing;      //!< Spacing of pixels
    double         fConcSurvProb;      //!< Pixel properties; probability of survival
    Physics::Vec3D fFPRotation; 
    double         fCameraIP;
    bool           fPixelParity;       //!< Parity for counting pixels; 1 => counting toward +x-axis
    //!  in reflector r.f.; -1 => counting toward -x-axis

#if 0    
    // ************************************************************************
    // Secondary Optics Parameters
    // ************************************************************************
    bool           fHasSecondary;
    double         fRefractiveIndex;   //!<secondary optics parameters
    double         fRefractiveIndex1;  //!<secondary optics parameters
    double         fRefractiveIndex2;  //!<secondary optics parameters
    double         fCE1Parameter0;
    double         fCE1Parameter2;
    double         fCE1Parameter3;
    double         fCE1Parameter4;
    double         fCE1Parameter5;
    double         fCE2Parameter0;
    double         fCE2Parameter2;
    double         fCE2Parameter3;
    double         fCE2Parameter4;
    double         fCE2Parameter5;
#endif

    // ************************************************************************
    // Mirrors and Pixel data
    // ************************************************************************
    
    std::vector<VSOMirror*> fMirrors;
    std::vector<VSOMirror*> fMirrorsByHexID;
    std::vector<VSOPixel*>  fPixels;
    std::vector<VSOPixel*>  fPixelsByHexID;
    
    // ************************************************************************
    // Precalculated rotation vector -- not stored in DB
    // ************************************************************************
    
    Physics::Vec3D fRotationVector;  //!<the pre-calculated rotation vector
    
    void calculateRotationVector();
  };
  
  inline const VSOMirror* VSOTelescope::mirror(unsigned id) const
  {
    if(id>=fMirrors.size())return 0;
    else return fMirrors[id];
  }
  
  inline const VSOMirror* VSOTelescope::mirrorByHexID(unsigned hexID) const
  {
    if((hexID<1)||(hexID>fMirrorsByHexID.size()))return 0;
    else return fMirrorsByHexID[hexID-1];
  }
  
  inline const VSOPixel* VSOTelescope::pixel(unsigned id) const
  {
    if(id>fPixels.size())return 0;
    else return fPixels[id];
  }
  
  inline const VSOPixel* VSOTelescope::pixelByHexID(unsigned hexID) const
  {
    if((hexID<1)||(hexID>fPixelsByHexID.size()))return 0;
    else return fPixelsByHexID[hexID-1];
  }
}

#endif // VSOTELESCOPE_HPP
