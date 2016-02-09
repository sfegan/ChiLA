//-*-mode:c++; mode:font-lock;-*-

#include <sstream>
#include <algorithm>

#include <VSDataConverter.hpp>
#include <VSSimDBTables.hpp>

#include "VSOTelescope.hpp"

using namespace Physics;
using namespace VERITAS;

VSOTelescope::VSOTelescope():
  fID(), fTelescopeHexID(), fPos(), 
  fDeltaY(), fAlphaX(), fAlphaY(),
  fElevation(), fAzimuth(), 
  fTranslation(), fCurvatureRadius(), fAperture(), 
  fFacetSpacing(), fFacetSize(), 
  fReflectorRotation(), fAlignmentPoint(), fHexagonRingsN(),
  fReflectorIP(), fMirrorParity(), fFPTranslation(), fCameraDiameter(),
  fFieldOfView(), fCathodeDiameter(), fPixelSpacing(), fConcSurvProb(),
  fFPRotation(), fCameraIP(), fPixelParity(),
#if 0
  fHasSecondary(false),
  fRefractiveIndex(), fRefractiveIndex1(), fRefractiveIndex2(), 
  fCE1Parameter0(), fCE1Parameter2(), fCE1Parameter3(), fCE1Parameter4(), 
  fCE1Parameter5(), fCE2Parameter0(), fCE2Parameter2(), fCE2Parameter3(),
  fCE2Parameter4(), fCE2Parameter5(), 
#endif
  fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(), fRotationVector()
{
  calculateRotationVector();
}

VSOTelescope::
VSOTelescope(unsigned TID, unsigned THID, const Vec3D&P, 
	     double DY, double AX, double AY, double EL, double AZ,
	     const Vec3D& T, double CR, double A, double FSP, double FS, 
	     double RR, const Vec3D& AP, unsigned HRN, double RIP, bool MP, 
	     const Vec3D& FPT, double CD, double FOV, double D, double PS, 
	     double CSP, const Vec3D& FPR, double CIP, bool PP
#if 0
	     , bool SEC, double RI, double RI1, double RI2, double C10, 
	     double C12, double C13, double C14, double C15, double C20,
	     double C22, double C23,double C24, double C25
#endif
	     ):
  fID(TID), fTelescopeHexID(THID), fPos(P),
  fDeltaY(DY), fAlphaX(AX), fAlphaY(AY), fElevation(EL), fAzimuth(AZ), 
  fTranslation(T), fCurvatureRadius(CR), fAperture(A), fFacetSpacing(FSP), 
  fFacetSize(FS), fReflectorRotation(RR), fAlignmentPoint(AP), 
  fHexagonRingsN(HRN), fReflectorIP(RIP), fMirrorParity(MP), 
  fFPTranslation(FPT),  fCameraDiameter(CD), fFieldOfView(FOV), 
  fCathodeDiameter(D), fPixelSpacing(PS), fConcSurvProb(CSP),
  fFPRotation(FPR), fCameraIP(CIP), fPixelParity(PP),
#if 0
  fHasSecondary(SEC),
  fRefractiveIndex(RI), fRefractiveIndex1(RI1), fRefractiveIndex2(RI2),
  fCE1Parameter0(C10), fCE1Parameter2(C12), fCE1Parameter3(C13), 
  fCE1Parameter4(C14), fCE1Parameter5(C15), fCE2Parameter0(C20), 
  fCE2Parameter2(C22), fCE2Parameter3(C23), fCE2Parameter4(C24),
  fCE2Parameter5(C25), 
#endif
  fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(), fRotationVector()
{
  calculateRotationVector();
}

VSOTelescope::VSOTelescope(const VSOTelescope& o):
  fID(o.fID), fTelescopeHexID(o.fTelescopeHexID),
  fPos(o.fPos), fDeltaY(o.fDeltaY), fAlphaX(o.fAlphaX), fAlphaY(o.fAlphaY),
  fElevation(o.fElevation), fAzimuth(o.fAzimuth), fTranslation(o.fTranslation),
  fCurvatureRadius(o.fCurvatureRadius), fAperture(o.fAperture), 
  fFacetSpacing(o.fFacetSpacing), fFacetSize(o.fFacetSize), 
  fReflectorRotation(o.fReflectorRotation),
  fAlignmentPoint(o.fAlignmentPoint), fHexagonRingsN(o.fHexagonRingsN),
  fReflectorIP(o.fReflectorIP), fMirrorParity(o.fMirrorParity),
  fFPTranslation(o.fFPTranslation), fCameraDiameter(o.fCameraDiameter),
  fFieldOfView(o.fFieldOfView), fCathodeDiameter(o.fCathodeDiameter),
  fPixelSpacing(o.fPixelSpacing), fConcSurvProb(o.fConcSurvProb),
  fFPRotation(o.fFPRotation), fCameraIP(o.fCameraIP),
  fPixelParity(o.fPixelParity), 
#if 0
  fHasSecondary(o.fHasSecondary),
  fRefractiveIndex(o.fRefractiveIndex),
  fRefractiveIndex1(o.fRefractiveIndex1), 
  fRefractiveIndex2(o.fRefractiveIndex2),
  fCE1Parameter0(o.fCE1Parameter0), fCE1Parameter2(o.fCE1Parameter2),
  fCE1Parameter3(o.fCE1Parameter3), fCE1Parameter4(o.fCE1Parameter4),
  fCE1Parameter5(o.fCE1Parameter5), fCE2Parameter0(o.fCE2Parameter0),
  fCE2Parameter2(o.fCE2Parameter2), fCE2Parameter3(o.fCE2Parameter3),
  fCE2Parameter4(o.fCE2Parameter4), fCE2Parameter5(o.fCE2Parameter5),
#endif
  fMirrors(), fMirrorsByHexID(), fPixels(), fPixelsByHexID(),
  fRotationVector()
{
  fMirrors.resize(o.fMirrors.size());
  fMirrorsByHexID.resize(o.fMirrorsByHexID.size());
  for(std::vector<VSOMirror*>::const_iterator i=o.fMirrors.begin(); 
      i!=o.fMirrors.end(); i++)
    {
      VSOMirror* mirror = new VSOMirror(**i);
      fMirrors[(*i)->id()]=mirror;
      fMirrorsByHexID[(*i)->hexID()]=mirror;
    }

  fPixels.resize(o.fPixels.size());
  fPixelsByHexID.resize(o.fPixelsByHexID.size());
  for(std::vector<VSOPixel*>::const_iterator i=o.fPixels.begin(); 
      i!=o.fPixels.end(); i++)
    {
      VSOPixel* pixel = new VSOPixel(**i);
      fPixels[(*i)->id()]=pixel;
      fPixelsByHexID[(*i)->hexID()]=pixel;
    }

  calculateRotationVector();
}

VSOTelescope::~VSOTelescope()
{
  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
}

// Stroustrup third edition sec 11.3.4 recommends default copy assignment
const VSOTelescope& VSOTelescope::operator =(const VSOTelescope& o)
{
  fID                = o.fID;
  fTelescopeHexID    = o.fTelescopeHexID;
  fPos               = o.fPos;
  fDeltaY            = o.fDeltaY;
  fAlphaX            = o.fAlphaX;
  fAlphaY            = o.fAlphaY;
  fElevation         = o.fElevation;
  fAzimuth           = o.fAzimuth;
  fTranslation       = o.fTranslation;
  fCurvatureRadius   = o.fCurvatureRadius;
  fAperture          = o.fAperture;
  fFacetSpacing      = o.fFacetSpacing;
  fFacetSize         = o.fFacetSize;
  fReflectorRotation = o.fReflectorRotation;
  fAlignmentPoint    = o.fAlignmentPoint;
  fHexagonRingsN     = o.fHexagonRingsN;
  fReflectorIP       = o.fReflectorIP;
  fMirrorParity      = o.fMirrorParity;
  fFPTranslation     = o.fFPTranslation;
  fCameraDiameter    = o.fCameraDiameter;
  fFieldOfView       = o.fFieldOfView;
  fCathodeDiameter   = o.fCathodeDiameter;
  fPixelSpacing      = o.fPixelSpacing;
  fConcSurvProb      = o.fConcSurvProb;
  fFPRotation        = o.fFPRotation;
  fCameraIP          = o.fCameraIP;
  fPixelParity       = o.fPixelParity;
#if 0
  fHasSecondary      = o.fHasSecondary,
  fRefractiveIndex   = o.fRefractiveIndex;
  fRefractiveIndex1  = o.fRefractiveIndex1; 
  fRefractiveIndex2  = o.fRefractiveIndex2;
  fCE1Parameter0     = o.fCE1Parameter0; 
  fCE1Parameter2     = o.fCE1Parameter2;
  fCE1Parameter3     = o.fCE1Parameter3; 
  fCE1Parameter4     = o.fCE1Parameter4;
  fCE1Parameter5     = o.fCE1Parameter5;
  fCE2Parameter0     = o.fCE2Parameter0;
  fCE2Parameter2     = o.fCE2Parameter2;
  fCE2Parameter3     = o.fCE2Parameter3;
  fCE2Parameter4     = o.fCE2Parameter4;
  fCE2Parameter5     = o.fCE2Parameter5;
#endif

  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;

  fMirrors.clear();
  fMirrorsByHexID.clear();

  fMirrors.resize(o.fMirrors.size());
  fMirrorsByHexID.resize(o.fMirrorsByHexID.size());

  for(std::vector<VSOMirror*>::const_iterator i=o.fMirrors.begin(); 
      i!=o.fMirrors.end(); i++)
    {
      VSOMirror* mirror = new VSOMirror(**i);
      fMirrors[(*i)->id()]=mirror;
      fMirrorsByHexID[(*i)->hexID()]=mirror;
    }

  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;

  fPixels.clear();
  fPixelsByHexID.clear();

  fPixels.resize(o.fPixels.size());
  fPixelsByHexID.resize(o.fPixelsByHexID.size());

  for(std::vector<VSOPixel*>::const_iterator i=o.fPixels.begin(); 
      i!=o.fPixels.end(); i++)
    {
      VSOPixel* pixel = new VSOPixel(**i);
      fPixels[(*i)->id()]=pixel;
      fPixelsByHexID[(*i)->hexID()]=pixel;
    }

  calculateRotationVector();

  return *this;
}

// ****************************************************************************
// Accessor
// ****************************************************************************

Vec3D VSOTelescope::opticalAxis() const
{
  return Vec3D(cos(fElevation)*sin(fAzimuth),
	       cos(fElevation)*cos(fAzimuth),
	       sin(fElevation));
}

// ****************************************************************************
// Repoint the telescope along a vector
// ****************************************************************************

bool VSOTelescope::pointTelescope(const Vec3D& v)
{
  if(v.Norm2()==0)return false;
  fElevation = atan2(v.z,sqrt(v.x*v.x + v.y*v.y));
  fAzimuth = fmod(atan2(v.x,v.y)+2.0*M_PI, 2.0*M_PI);
  calculateRotationVector();
  return true;
}

bool VSOTelescope::pointTelescopeAzEl(const double az_rad, const double el_rad)
{
  fElevation = fmod(fmod(el_rad,2.0*M_PI)+2.0*M_PI, 2.0*M_PI);
  fAzimuth = fmod(fmod(az_rad,2.0*M_PI)+2.0*M_PI, 2.0*M_PI);
  calculateRotationVector();
  return true;
}

// ****************************************************************************
// Calculate rotation vector and map between Global and Telescope coordinates
// ****************************************************************************

void VSOTelescope::calculateRotationVector()
{
  // Rotation vector maps from Reflector to Global
  fRotationVector = 
    Vec3D(1,0,0)*fElevation &
    Vec3D(0,1,0)*fDeltaY &
    Vec3D(0,0,-1)*fAzimuth &
    Vec3D(0,1,0)*fAlphaX &
    Vec3D(1,0,0)*fAlphaY;
}

void VSOTelescope::globalToReflector(Physics::Vec3D& v) const
{
  // First: Translate from center of array to drive axes intersection
  v -= fPos;
  // Second: Rotate coordinate system to reflector orientation
  v.Rotate(-fRotationVector);
  // Third: Translate from intersection of drive axes to reflector
  v += fTranslation;
}

void VSOTelescope::reflectorToGlobal(Physics::Vec3D& v) const
{
  // First: Translate from reflector to intersection of drive axes
  v -= fTranslation;
  // Second: Rotate coordinate system to ground based
  v.Rotate(fRotationVector);
  // Third: Translate from drive axes intersection to center of array
  v += fPos;
}

void VSOTelescope::globalToReflector(Particle& p) const
{
  // First: Translate from center of array to drive axes intersection
  p.TranslateOrigin(Vec4D(0,fPos));
  // Second: Rotate coordinate system to reflector orientation
  p.Rotate(-fRotationVector);
  // Third: Translate from intersection of drive axes to reflector
  p.TranslateOrigin(Vec4D(0,-fTranslation));
}

void VSOTelescope::reflectorToGlobal(Particle& p) const
{
  // First: Translate from reflector to intersection of drive axes
  p.TranslateOrigin(Vec4D(0,fTranslation));
  // Second: Rotate coordinate system to ground based
  p.Rotate(fRotationVector);
  // Third: Translate from drive axes intersection to center of array
  p.TranslateOrigin(Vec4D(0,-fPos));
}

void VSOTelescope::focalPlaneToReflector(Physics::Particle& p) const
{
  // First: Rotate coordinate system
  p.Rotate(fFPRotation);
  // Second: Translate from center of Focal Plane
  p.TranslateOrigin(Vec4D(0,-fFPTranslation));
}

void VSOTelescope::reflectorToFocalPlane(Physics::Particle& p) const
{
  // First: Translate to center of Focal Plane
  p.TranslateOrigin(Vec4D(0,fFPTranslation));
  // Second: Rotate coordinate system
  p.Rotate(-fFPRotation);
}

// ****************************************************************************
// Fill the mirrors and pixels tables up with randomly generated objects
// ****************************************************************************

void VSOTelescope::
populateMirrorsAndPixelsRandom(const VSOArrayParameters& param, 
			       RandomNumbers& rng)
{
  // **************************************************************************
  // Clear the MIRRORs table and repopulate it with randomly generated mirrors
  // **************************************************************************
  double reflector_r2 = fAperture*fAperture/4.0;
  double reflector_c2 = fCurvatureRadius*fCurvatureRadius;
  Vec3D reflector_center(0, fCurvatureRadius, 0);

  std::set<unsigned> mirrors_missing;
  Physics::tokenize(param.MirrorMissingList,mirrors_missing);

  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  fMirrors.clear();
  fMirrorsByHexID.clear();  

  int num_hex_mirror_sites = 3*fHexagonRingsN*(fHexagonRingsN+1)+1;

  unsigned id = 0;
  for(int i=0; i<num_hex_mirror_sites; i++)
    {
      int hexid = i+1;

      if(mirrors_missing.find(hexid) != mirrors_missing.end())
	{
	  fMirrorsByHexID.push_back(0);
	  continue; // skip mirror if on the mirring list
	}

      Vec3D nominal_position;

      // Compute the mirror's nominal position
      nh_to_xy(&hexid, &nominal_position.x, &nominal_position.z); 
      if(fMirrorParity)nominal_position.x=-nominal_position.x;
      nominal_position.x *= fFacetSpacing;
      nominal_position.z *= fFacetSpacing;

      double X2 = nominal_position.x*nominal_position.x;
      double Z2 = nominal_position.z*nominal_position.z;

      if( (reflector_r2 - X2 - Z2) <= 0 )
	{
	  fMirrorsByHexID.push_back(0);
	  continue; // skip mirror if projected hex position too far out
	}

      nominal_position.y = fCurvatureRadius-sqrt(reflector_c2 - X2 - Z2);

      // Add Gaussian normal/tangenetial position error
      Vec3D reflector_normal(reflector_center-nominal_position);
      reflector_normal /= reflector_normal.Norm();

      Vec3D position(nominal_position);
      position += reflector_normal*(rng.Normal()*param.MirrorPosNormalDisp);
      position -= reflector_center;
      position.ScatterDirection(param.MirrorPosTangentDisp/fCurvatureRadius,rng);
      position += reflector_center;
      
      // Rotate by global rotation angle
      position.Rotate(Vec3D(0,fReflectorRotation,0));

      // Get the (perturbed) alignment angle of mirror
      Vec3D alignment;

      if(param.ReflectorAlignMode == 0)
	{
	  // Standard DC alignment to a fixed point in space (with scatter)
	  alignment = (fAlignmentPoint-position);
	  alignment /= alignment.Norm();

	  double align_disp = 
	    param.MirrorAlignTangentDisp/fAlignmentPoint.Norm();
	  alignment.ScatterDirection(align_disp,rng);
	}
      else
	{
	  // Alignment to a point in the focal plane
	  
	  double stheta = sin(param.ReflectorFPAlignTheta);
	  double ctheta = cos(param.ReflectorFPAlignTheta);
	  double ttheta = stheta/ctheta;

	  double sphi = sin(param.ReflectorFPAlignPhi);
	  double cphi = cos(param.ReflectorFPAlignPhi);

	  double y_fp = fFPTranslation.y; // Distance to focal plane
	  Vec3D r_fp(y_fp*ttheta*sphi ,y_fp, -y_fp*ttheta*cphi);

	  Vec3D e_in(-stheta*sphi, ctheta, stheta*cphi);
	  Vec3D e_out = r_fp-position;
	  e_out /= e_out.Norm();

	  alignment = e_in;

	  Vec3D e_rot = e_in^e_out;
	  double strot2 = e_rot.Norm();
	  if(strot2 != 0)
	    {
	      double ctrot2 = e_in*e_out;
	      double trot2 = atan2(strot2,ctrot2);
	      e_rot *= 0.5*trot2/strot2;
	      alignment.Rotate(e_rot);
	    }
	}

      double focal_length = 
	param.MirrorFLength + rng.Normal()*param.MirrorFLengthDisp;

      double spot_size;
      if(param.MirrorSpotSizeDisp > 0)
	spot_size = rng.GammaByMeanAndStdDev(param.MirrorSpotSize,
					     param.MirrorSpotSizeDisp);
      else spot_size = param.MirrorSpotSize;

      // If param.MirrorSpotSizePhotonFraction not set -- assume FWHM
      if((param.MirrorSpotSizePhotonFraction>0.0)&&
	 (param.MirrorSpotSizePhotonFraction<1.0))
	spot_size /= 
	  2.0*sqrt(log(1.0/(1.0-param.MirrorSpotSizePhotonFraction)));
      else spot_size /= 2.0*sqrt(log(2.0));
      
      VSOMirror* mirror = 
	new VSOMirror(this, id, hexid, false, position, alignment, 
		     focal_length, spot_size, param.MirrorDegradingFactor);
      
      fMirrors.push_back(mirror);
      fMirrorsByHexID.push_back(mirror);

      id++;
    }
  
  // **************************************************************************
  // Clear the PIXELSs table and repopulate it with randomly generated pixels
  // **************************************************************************

  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
  fPixels.clear();
  fPixelsByHexID.clear();  

  std::set<unsigned> pixels_missing;
  Physics::tokenize(param.PixelMissingList,pixels_missing);

  unsigned num_hex_pixel_rings = 
    unsigned(floor((fCameraDiameter/2.0)/fPixelSpacing))+2;
  unsigned num_hex_pixel_sites = 
    3*num_hex_pixel_rings*(num_hex_pixel_rings+1)+1;

  id = 0;
  for(unsigned i=0; i<num_hex_pixel_sites; i++)
    {
      int hexid = i+1;

      if(pixels_missing.find(hexid) != pixels_missing.end())
	{
	  fPixelsByHexID.push_back(0);
	  continue; // skip pixel if on the missing list
	}

      Vec3D nominal_position;

      // Compute the pixel's nominal position
      nh_to_xy(&hexid, &nominal_position.x, &nominal_position.z); 
      if(fPixelParity)nominal_position.x=-nominal_position.x;
      nominal_position.x *= fPixelSpacing;
      nominal_position.z *= fPixelSpacing;
      nominal_position.y = 0;

      if(nominal_position.Norm() > fCameraDiameter/2.0)
	{
	  fPixelsByHexID.push_back(0);
	  continue; // skip pixel if position too far out
	}

      nominal_position.Rotate(fFPRotation);

      VSOPixel* pixel = new VSOPixel(this, id, hexid, false, nominal_position);
      
      fPixels.push_back(pixel);
      fPixelsByHexID.push_back(pixel);

      id++;
    }
}

// ****************************************************************************
// DATABASE functions
// ****************************************************************************

VSDBStatement* VSOTelescope::createInsertQuery(VSDatabase* db)
{
#if 0
  return db->createInsertQuery(VSIMDB_TABLE_NAME_TELESCOPE,52);
#else
  return db->createInsertQuery(VSIMDB_TABLE_NAME_TELESCOPE,38);
#endif
}

VSDBStatement* VSOTelescope::createSelectQuery(VSDatabase* db)
{
  return db->createSelectQuery(VSIMDB_TABLE_NAME_TELESCOPE,
			       "OpticsID=?");
}

VSOTelescope* VSOTelescope::createFromDatabaseRow(VSDBStatement* stmt)
{
  VSOTelescope* telescope = new VSOTelescope;
  unsigned OID;

  stmt->bindToResult(OID);

  stmt->bindToResult(telescope->fID);
  stmt->bindToResult(telescope->fTelescopeHexID);
  stmt->bindToResult(telescope->fPos.x);
  stmt->bindToResult(telescope->fPos.y);
  stmt->bindToResult(telescope->fPos.z);

  stmt->bindToResult(telescope->fDeltaY);
  stmt->bindToResult(telescope->fAlphaX);
  stmt->bindToResult(telescope->fAlphaY);
  stmt->bindToResult(telescope->fElevation);
  stmt->bindToResult(telescope->fAzimuth);

  stmt->bindToResult(telescope->fTranslation.x);
  stmt->bindToResult(telescope->fTranslation.y);
  stmt->bindToResult(telescope->fTranslation.z);
  stmt->bindToResult(telescope->fCurvatureRadius);
  stmt->bindToResult(telescope->fAperture);

  stmt->bindToResult(telescope->fFacetSpacing);
  stmt->bindToResult(telescope->fFacetSize);
  stmt->bindToResult(telescope->fReflectorRotation);
  stmt->bindToResult(telescope->fAlignmentPoint.x);
  stmt->bindToResult(telescope->fAlignmentPoint.y);

  stmt->bindToResult(telescope->fAlignmentPoint.z);
  stmt->bindToResult(telescope->fHexagonRingsN);
  stmt->bindToResult(telescope->fReflectorIP);
  stmt->bindToResult(telescope->fMirrorParity);
  stmt->bindToResult(telescope->fFPTranslation.x);

  stmt->bindToResult(telescope->fFPTranslation.y);
  stmt->bindToResult(telescope->fFPTranslation.z);
  stmt->bindToResult(telescope->fCameraDiameter);
  stmt->bindToResult(telescope->fFieldOfView);
  stmt->bindToResult(telescope->fCathodeDiameter);

  stmt->bindToResult(telescope->fPixelSpacing);
  stmt->bindToResult(telescope->fConcSurvProb);
  stmt->bindToResult(telescope->fFPRotation.x);
  stmt->bindToResult(telescope->fFPRotation.y);
  stmt->bindToResult(telescope->fFPRotation.z);

  stmt->bindToResult(telescope->fCameraIP);
  stmt->bindToResult(telescope->fPixelParity);
#if 0
  stmt->bindToResult(telescope->fHasSecondary);
  stmt->bindToResult(telescope->fRefractiveIndex);
  stmt->bindToResult(telescope->fRefractiveIndex1);

  stmt->bindToResult(telescope->fRefractiveIndex2);
  stmt->bindToResult(telescope->fCE1Parameter0);
  stmt->bindToResult(telescope->fCE1Parameter2);
  stmt->bindToResult(telescope->fCE1Parameter3);
  stmt->bindToResult(telescope->fCE1Parameter4);

  stmt->bindToResult(telescope->fCE1Parameter5);
  stmt->bindToResult(telescope->fCE2Parameter0);
  stmt->bindToResult(telescope->fCE2Parameter2);
  stmt->bindToResult(telescope->fCE2Parameter3);
  stmt->bindToResult(telescope->fCE2Parameter4);

  stmt->bindToResult(telescope->fCE2Parameter5);
#endif

  int row = stmt->retrieveNextRow();
  if(row == -1)
    {
      std::cerr << stmt->getErrorMessage();
      assert(0);
    }
  else if(row == 0)
    {
      delete telescope;
      return 0;
    }

  return telescope;
}

void VSOTelescope::populateMirrorsAndPixelsFromDatabase(VSDatabase* db,
							uint32_t optics_id)
{
  VSDBStatement* stmt;

  // **************************************************************************
  // Clear the MIRRORs table and then re-populate it from the DB
  // **************************************************************************

  for(std::vector<VSOMirror*>::iterator i=fMirrors.begin();
      i!=fMirrors.end(); i++)delete *i;
  fMirrors.clear();
  fMirrorsByHexID.clear();  

  int num_hex_mirror_sites = 3*fHexagonRingsN*(fHexagonRingsN+1)+1;
  fMirrorsByHexID.resize(num_hex_mirror_sites);

  stmt = VSOMirror::createSelectQuery(db);
  assert(stmt);
  stmt->bindToParam(optics_id);
  stmt->bindToParam(fID);
  assert(stmt->execute() >= 0);

  while(VSOMirror* mirror = VSOMirror::createFromDatabaseRow(stmt,this))
    {
      if(mirror->id() >= fMirrors.size())
	fMirrors.resize(mirror->id()+1);
      fMirrors[mirror->id()]=mirror;
      
      if(mirror->hexID() > fMirrorsByHexID.size())
	fMirrorsByHexID.resize(mirror->hexID());
      fMirrorsByHexID[mirror->hexID()-1]=mirror;
    }

  delete stmt;

  // **************************************************************************
  // Clear the PIXELs table and then re-populate it from the DB
  // **************************************************************************

  for(std::vector<VSOPixel*>::iterator i=fPixels.begin();
      i!=fPixels.end(); i++)delete *i;
  fPixels.clear();
  fPixelsByHexID.clear();  

  unsigned num_hex_pixel_rings = 
    unsigned(floor((fCameraDiameter/2.0)/fPixelSpacing))+2;
  unsigned num_hex_pixel_sites = 
    3*num_hex_pixel_rings*(num_hex_pixel_rings+1)+1;

  fPixelsByHexID.resize(num_hex_pixel_sites);

  stmt = VSOPixel::createSelectQuery(db);
  assert(stmt);
  stmt->bindToParam(optics_id);
  stmt->bindToParam(fID);
  assert(stmt->execute() >= 0);
  
  while(VSOPixel* pixel = VSOPixel::createFromDatabaseRow(stmt,this))
    {
      if(pixel->id() >= fPixels.size())
	fPixels.resize(pixel->id()+1);
      fPixels[pixel->id()]=pixel;

      if(pixel->hexID() > fPixelsByHexID.size())
	fPixelsByHexID.resize(pixel->hexID());
      fPixelsByHexID[pixel->hexID()-1]=pixel;
    }

  delete stmt;

  // Recalculate rotation vector
  calculateRotationVector();

  return;
}

void VSOTelescope::writeToDatabase(VSDatabase* db, VSDBStatement* stmt, 
				   uint32_t optics_id) const
{
  stmt->bindToParam(optics_id);

  stmt->bindToParam(fID);
  stmt->bindToParam(fTelescopeHexID);
  stmt->bindToParam(fPos.x);
  stmt->bindToParam(fPos.y);
  stmt->bindToParam(fPos.z);

  stmt->bindToParam(fDeltaY);
  stmt->bindToParam(fAlphaX);
  stmt->bindToParam(fAlphaY);
  stmt->bindToParam(fElevation);
  stmt->bindToParam(fAzimuth);

  stmt->bindToParam(fTranslation.x);
  stmt->bindToParam(fTranslation.y);
  stmt->bindToParam(fTranslation.z);
  stmt->bindToParam(fCurvatureRadius);
  stmt->bindToParam(fAperture);

  stmt->bindToParam(fFacetSpacing);
  stmt->bindToParam(fFacetSize);
  stmt->bindToParam(fReflectorRotation);
  stmt->bindToParam(fAlignmentPoint.x);
  stmt->bindToParam(fAlignmentPoint.y);

  stmt->bindToParam(fAlignmentPoint.z);
  stmt->bindToParam(fHexagonRingsN);
  stmt->bindToParam(fReflectorIP);
  stmt->bindToParam(fMirrorParity);
  stmt->bindToParam(fFPTranslation.x);

  stmt->bindToParam(fFPTranslation.y);
  stmt->bindToParam(fFPTranslation.z);
  stmt->bindToParam(fCameraDiameter);
  stmt->bindToParam(fFieldOfView);
  stmt->bindToParam(fCathodeDiameter);

  stmt->bindToParam(fPixelSpacing);
  stmt->bindToParam(fConcSurvProb);
  stmt->bindToParam(fFPRotation.x);
  stmt->bindToParam(fFPRotation.y);
  stmt->bindToParam(fFPRotation.z);

  stmt->bindToParam(fCameraIP);
  stmt->bindToParam(fPixelParity);
#if 0
  stmt->bindToParam(fHasSecondary);
  stmt->bindToParam(fRefractiveIndex);
  stmt->bindToParam(fRefractiveIndex1);

  stmt->bindToParam(fRefractiveIndex2);
  stmt->bindToParam(fCE1Parameter0);
  stmt->bindToParam(fCE1Parameter2);
  stmt->bindToParam(fCE1Parameter3);
  stmt->bindToParam(fCE1Parameter4);

  stmt->bindToParam(fCE1Parameter5);
  stmt->bindToParam(fCE2Parameter0);
  stmt->bindToParam(fCE2Parameter2);
  stmt->bindToParam(fCE2Parameter3);
  stmt->bindToParam(fCE2Parameter4);

  stmt->bindToParam(fCE2Parameter5);
#endif
  assert(stmt->execute() == 1);

  VSDBStatement* sub_stmt;

  sub_stmt = VSOMirror::createInsertQuery(db);
  assert(sub_stmt);
  for(std::vector<VSOMirror*>::const_iterator i = fMirrors.begin(); 
      i != fMirrors.end(); i++)(*i)->writeToDatabase(sub_stmt, optics_id);
  delete(sub_stmt);

  sub_stmt = VSOPixel::createInsertQuery(db);
  assert(sub_stmt);
  for(std::vector<VSOPixel*>::const_iterator i = fPixels.begin(); 
      i != fPixels.end(); i++)(*i)->writeToDatabase(sub_stmt, optics_id);
  delete(sub_stmt);
}

void VSOTelescope::createTelescopeTable(VSDatabase* db)
{
  VSOTelescope temp; // stupid to create this --- C++ should have a "typeof" operator (g++ does but it
                    // is not portable .. and doesn't work in a static function anyway)
  unsigned OID;

  db->createTable(VSIMDB_TABLE_NAME_TELESCOPE,

	 db->sqlSpecOf("OpticsID",OID,true,"NOT NULL")+

	 db->sqlSpecOf("TelescopeID",temp.fID,false,"NOT NULL")+
	 db->sqlSpecOf("TelescopeHexID",temp.fTelescopeHexID,false,"NOT NULL")+
	 db->sqlSpecOf("PosX",temp.fPos.x,false,"NOT NULL")+
	 db->sqlSpecOf("PosY",temp.fPos.y,false,"NOT NULL")+
	 db->sqlSpecOf("PosZ",temp.fPos.z,false,"NOT NULL")+

	 db->sqlSpecOf("DeltaY",temp.fDeltaY,false,"NOT NULL")+
	 db->sqlSpecOf("AlphaX",temp.fAlphaX,false,"NOT NULL")+
	 db->sqlSpecOf("AlphaY",temp.fAlphaY,false,"NOT NULL")+
	 db->sqlSpecOf("Elevation",temp.fElevation,false,"NOT NULL")+
	 db->sqlSpecOf("Azimuth",temp.fAzimuth,false,"NOT NULL")+

	 db->sqlSpecOf("TranslationX",temp.fTranslation.x,false,"NOT NULL")+
	 db->sqlSpecOf("TranslationY",temp.fTranslation.y,false,"NOT NULL")+
	 db->sqlSpecOf("TranslationZ",temp.fTranslation.z,false,"NOT NULL")+
	 db->sqlSpecOf("CurvatureRadius",temp.fCurvatureRadius,false,"NOT NULL")+
	 db->sqlSpecOf("Aperture",temp.fAperture,false,"NOT NULL")+
	
	 db->sqlSpecOf("FacetSpacing",temp.fFacetSpacing,false,"NOT NULL")+
	 db->sqlSpecOf("FacetSize",temp.fFacetSize,false,"NOT NULL")+
	 db->sqlSpecOf("ReflectorRotation",temp.fReflectorRotation,false,"NOT NULL")+
	 db->sqlSpecOf("AlignmentPointX",temp.fAlignmentPoint.x,false,"NOT NULL")+
	 db->sqlSpecOf("AlignmentPointY",temp.fAlignmentPoint.y,false,"NOT NULL")+

	 db->sqlSpecOf("AlignmentPointZ",temp.fAlignmentPoint.z,false,"NOT NULL")+
	 db->sqlSpecOf("HexagonRingsN",temp.fHexagonRingsN,false,"NOT NULL")+
	 db->sqlSpecOf("ReflectorIP",temp.fReflectorIP,false,"NOT NULL")+
	 db->sqlSpecOf("MirrorParity",temp.fMirrorParity,false,"NOT NULL")+
 	 db->sqlSpecOf("FPTranslationX",temp.fFPTranslation.x,false,"NOT NULL")+

	 db->sqlSpecOf("FPTranslationY",temp.fFPTranslation.y,false,"NOT NULL")+
	 db->sqlSpecOf("FPTranslationZ",temp.fFPTranslation.z,false,"NOT NULL")+
	 db->sqlSpecOf("CameraDiameter",temp.fCameraDiameter,false,"NOT NULL")+
	 db->sqlSpecOf("FieldOfView",temp.fFieldOfView,false,"NOT NULL")+
	 db->sqlSpecOf("PixelDiameter",temp.fCathodeDiameter,false,"NOT NULL")+

	 db->sqlSpecOf("PixelSpacing",temp.fPixelSpacing,false,"NOT NULL")+
	 db->sqlSpecOf("ConcSurvProb",temp.fConcSurvProb,false,"NOT NULL")+
	 db->sqlSpecOf("FPRotationX",temp.fFPRotation.x,false,"NOT NULL")+
	 db->sqlSpecOf("FPRotationY",temp.fFPRotation.y,false,"NOT NULL")+
	 db->sqlSpecOf("FPRotationZ",temp.fFPRotation.z,false,"NOT NULL")+

	 db->sqlSpecOf("CameraIP",temp.fCameraIP,false,"NOT NULL")+
	 db->sqlSpecOf("PixelParity",temp.fPixelParity,false,"NOT NULL")+
#if 0
	 db->sqlSpecOf("HasSecondary",temp.fHasSecondary,false,"NOT NULL")+
	 db->sqlSpecOf("RefractiveIndex",temp.fRefractiveIndex,false,"NOT NULL")+
	 db->sqlSpecOf("RefractiveIndex1",temp.fRefractiveIndex1,false,"NOT NULL")+

	 db->sqlSpecOf("RefractiveIndex2",temp.fRefractiveIndex2,false,"NOT NULL")+
	 db->sqlSpecOf("CE1Parameter0",temp.fCE1Parameter0,false,"NOT NULL")+
	 db->sqlSpecOf("CE1Parameter2",temp.fCE1Parameter2,false,"NOT NULL")+
	 db->sqlSpecOf("CE1Parameter3",temp.fCE1Parameter3,false,"NOT NULL")+
	 db->sqlSpecOf("CE1Parameter4",temp.fCE1Parameter4,false,"NOT NULL")+

	 db->sqlSpecOf("CE1Parameter5",temp.fCE1Parameter5,false,"NOT NULL")+
	 db->sqlSpecOf("CE2Parameter0",temp.fCE1Parameter0,false,"NOT NULL")+
	 db->sqlSpecOf("CE2Parameter2",temp.fCE1Parameter2,false,"NOT NULL")+
	 db->sqlSpecOf("CE2Parameter3",temp.fCE1Parameter3,false,"NOT NULL")+
	 db->sqlSpecOf("CE2Parameter4",temp.fCE1Parameter4,false,"NOT NULL")+
	
	 db->sqlSpecOf("CE2Parameter5",temp.fCE1Parameter5,false,"NOT NULL")+
#endif

		  ", PRIMARY KEY ( OpticsID, TelescopeID )",

		  VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
}

void VSOTelescope::dumpShort(std::ostream& stream) const
{
  stream
    << "TELESCOPE "
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fTelescopeHexID) << ' '
    << VSDataConverter::toString(fMirrors.size()) << ' '
    << VSDataConverter::toString(fPixels.size()) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << ' '
    << VSDataConverter::toString(fDeltaY) << ' '
    << VSDataConverter::toString(fAlphaX) << ' '
    << VSDataConverter::toString(fAlphaY) << ' '

    << VSDataConverter::toString(fElevation) << ' '
    << VSDataConverter::toString(fAzimuth) << ' '
    << VSDataConverter::toString(fTranslation.x) << ' '
    << VSDataConverter::toString(fTranslation.y) << ' '
    << VSDataConverter::toString(fTranslation.z) << ' '

    << VSDataConverter::toString(fCurvatureRadius) << ' '
    << VSDataConverter::toString(fAperture) << ' '
    << VSDataConverter::toString(fFacetSpacing) << ' '
    << VSDataConverter::toString(fFacetSize) << ' '
    << VSDataConverter::toString(fReflectorRotation) << ' '

    << VSDataConverter::toString(fAlignmentPoint.x) << ' '
    << VSDataConverter::toString(fAlignmentPoint.y) << ' '
    << VSDataConverter::toString(fAlignmentPoint.z) << ' '
    << VSDataConverter::toString(fHexagonRingsN) << ' '
    << VSDataConverter::toString(fReflectorIP) << ' '

    << VSDataConverter::toString(fMirrorParity) << ' '
    << VSDataConverter::toString(fFPTranslation.x) << ' '
    << VSDataConverter::toString(fFPTranslation.y) << ' '
    << VSDataConverter::toString(fFPTranslation.z) << ' '
    << VSDataConverter::toString(fCameraDiameter) << ' '

    << VSDataConverter::toString(fFieldOfView) << ' '
    << VSDataConverter::toString(fCathodeDiameter) << ' '
    << VSDataConverter::toString(fPixelSpacing) << ' '
    << VSDataConverter::toString(fConcSurvProb) << ' '
    << VSDataConverter::toString(fFPRotation.x) << ' '

    << VSDataConverter::toString(fFPRotation.y) << ' '
    << VSDataConverter::toString(fFPRotation.z) << ' '
    << VSDataConverter::toString(fCameraIP) << ' '
    << VSDataConverter::toString(fPixelParity) << ' '
#if 0
    << VSDataConverter::toString(fHasSecondary) << ' '

    << VSDataConverter::toString(fRefractiveIndex) << ' '
    << VSDataConverter::toString(fRefractiveIndex1) << ' '
    << VSDataConverter::toString(fRefractiveIndex2) << ' '
    << VSDataConverter::toString(fCE1Parameter0) << ' '
    << VSDataConverter::toString(fCE1Parameter2) << ' '

    << VSDataConverter::toString(fCE1Parameter3) << ' '
    << VSDataConverter::toString(fCE1Parameter4) << ' '
    << VSDataConverter::toString(fCE1Parameter5) << ' '
    << VSDataConverter::toString(fCE1Parameter0) << ' '
    << VSDataConverter::toString(fCE1Parameter2) << ' '

    << VSDataConverter::toString(fCE1Parameter3) << ' '
    << VSDataConverter::toString(fCE1Parameter4) << ' '
    << VSDataConverter::toString(fCE1Parameter5) 
#endif
    << std::endl;

  for(std::vector<VSOMirror*> ::const_iterator i = fMirrors.begin();
      i!=fMirrors.end(); i++)
    (*i)->dumpShort(stream);

  for(std::vector<VSOPixel*> ::const_iterator i = fPixels.begin();
      i!=fPixels.end(); i++)
    (*i)->dumpShort(stream);
}

void VSOTelescope::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fID, "", 30, l) << std::endl
    << FDAV("Telescope Hex ID", fTelescopeHexID, "", 30, l) << std::endl
    << FDAV("Num Mirrors", fMirrors.size(), "", 30, l) << std::endl
    << FDAV("Num Pixels", fPixels.size(), "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "cm", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "cm", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "cm", 30, l) << std::endl
    << FDAV("Delta Y", fDeltaY, "rad", 30, l) << std::endl
    << FDAV("Alpha X", fAlphaX, "rad", 30, l) << std::endl
    << FDAV("Alpha Y", fAlphaY, "rad", 30, l) << std::endl

    << FDAV("Elevation", fElevation, "rad", 30, l) << std::endl
    << FDAV("Azimuth", fAzimuth, "rad", 30, l) << std::endl
    << FDAV("Translation X", fTranslation.x, "cm", 30, l) << std::endl
    << FDAV("Translation Y", fTranslation.y, "cm", 30, l) << std::endl
    << FDAV("Translation Z", fTranslation.z, "cm", 30, l) << std::endl

    << FDAV("Curvature Radius", fCurvatureRadius, "cm", 30, l) << std::endl
    << FDAV("Aperture", fAperture, "cm", 30, l) << std::endl
    << FDAV("Facet Spacing", fFacetSpacing, "cm", 30, l) << std::endl
    << FDAV("Facet Size", fFacetSize, "cm", 30, l) << std::endl
    << FDAV("Reflector Rotation", fReflectorRotation, "rad", 30, l) << std::endl 

    << FDAV("Alignment Point X", fAlignmentPoint.x, "cm", 30, l) << std::endl
    << FDAV("Alignment Point Y", fAlignmentPoint.y, "cm", 30, l) << std::endl
    << FDAV("Alignment Point Z", fAlignmentPoint.z, "cm", 30, l) << std::endl
    << FDAV("Num Mirror Hexagon Rings", fHexagonRingsN, "", 30, l) << std::endl
    << FDAV("Reflector IP", fReflectorIP, "cm", 30, l) << std::endl

    << FDAV("Mirror Parity", fMirrorParity, "", 30, l) << std::endl
    << FDAV("FP Translation X", fFPTranslation.x, "cm", 30, l) << std::endl
    << FDAV("FP Translation Y", fFPTranslation.y, "cm", 30, l) << std::endl
    << FDAV("FP Translation Z", fFPTranslation.z, "cm", 30, l) << std::endl
    << FDAV("CameraDiameter", fCameraDiameter, "cm", 30, l) << std::endl

    << FDAV("Field Of View", fFieldOfView, "deg", 30, l) << std::endl
    << FDAV("Pixel Diameter", fCathodeDiameter, "cm", 30, l) << std::endl
    << FDAV("Pixel Spacing", fPixelSpacing, "cm", 30, l) << std::endl
    << FDAV("Conc Surv Prob", fConcSurvProb, "", 30, l) << std::endl
    << FDAV("FP Rotation X", fFPRotation.x, "rad", 30, l) << std::endl

    << FDAV("FP Rotation Y", fFPRotation.y, "rad", 30, l) << std::endl
    << FDAV("FP Rotation Z", fFPRotation.z, "rad", 30, l) << std::endl
    << FDAV("Camera IP", fCameraIP, "cm", 30, l) << std::endl
    << FDAV("Pixel Parity", fPixelParity, "", 30, l) << std::endl;
#if 0
    << FDAV("Has Secondary", fHasSecondary, "", 30, l) << std::endl

    << FDAV("Refractive Index", fRefractiveIndex, "", 30, l) << std::endl
    << FDAV("Refractive Index 1", fRefractiveIndex1, "", 30, l) << std::endl
    << FDAV("Refractive Index 2", fRefractiveIndex2, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 0", fCE1Parameter0, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 2", fCE1Parameter2, "", 30, l) << std::endl

    << FDAV("CE1 Parameter 3", fCE1Parameter3, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 4", fCE1Parameter4, "", 30, l) << std::endl
    << FDAV("CE1 Parameter 5", fCE1Parameter5, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 0", fCE1Parameter0, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 2", fCE1Parameter2, "", 30, l) << std::endl

    << FDAV("CE2 Parameter 3", fCE1Parameter3, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 4", fCE1Parameter4, "", 30, l) << std::endl
    << FDAV("CE2 Parameter 5", fCE1Parameter5, "", 30, l) << std::endl;
#endif

  for(std::vector<VSOMirror*> ::const_iterator i = fMirrors.begin();
      i!=fMirrors.end(); i++)
    {
      stream << std::endl;
      (*i)->dump(stream,l+1);
    }

  for(std::vector<VSOPixel*> ::const_iterator i = fPixels.begin();
      i!=fPixels.end(); i++)
    {
      stream << std::endl;
      (*i)->dump(stream,l+1);
    }
}

VSOTelescope* VSOTelescope::createFromShortDump(std::istream& stream)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOTelescope* telescope = new VSOTelescope;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("TELESCOPE"));

  unsigned mirrors_size;
  unsigned pixels_size;

  linestream
    >> telescope->fID
    >> telescope->fTelescopeHexID
    >> mirrors_size
    >> pixels_size
    >> telescope->fPos.x
    
    >> telescope->fPos.y
    >> telescope->fPos.z
    >> telescope->fDeltaY
    >> telescope->fAlphaX
    >> telescope->fAlphaY

    >> telescope->fElevation
    >> telescope->fAzimuth
    >> telescope->fTranslation.x
    >> telescope->fTranslation.y
    >> telescope->fTranslation.z

    >> telescope->fCurvatureRadius
    >> telescope->fAperture
    >> telescope->fFacetSpacing
    >> telescope->fFacetSize
    >> telescope->fReflectorRotation

    >> telescope->fAlignmentPoint.x
    >> telescope->fAlignmentPoint.y
    >> telescope->fAlignmentPoint.z
    >> telescope->fHexagonRingsN
    >> telescope->fReflectorIP

    >> telescope->fMirrorParity
    >> telescope->fFPTranslation.x
    >> telescope->fFPTranslation.y
    >> telescope->fFPTranslation.z
    >> telescope->fCameraDiameter

    >> telescope->fFieldOfView
    >> telescope->fCathodeDiameter
    >> telescope->fPixelSpacing
    >> telescope->fConcSurvProb
    >> telescope->fFPRotation.x

    >> telescope->fFPRotation.y
    >> telescope->fFPRotation.z
    >> telescope->fCameraIP
    >> telescope->fPixelParity
#if 0
    >> telescope->fHasSecondary

    >> telescope->fRefractiveIndex
    >> telescope->fRefractiveIndex1
    >> telescope->fRefractiveIndex2
    >> telescope->fCE1Parameter0
    >> telescope->fCE1Parameter2

    >> telescope->fCE1Parameter3
    >> telescope->fCE1Parameter4
    >> telescope->fCE1Parameter5
    >> telescope->fCE1Parameter0
    >> telescope->fCE1Parameter2

    >> telescope->fCE1Parameter3
    >> telescope->fCE1Parameter4
    >> telescope->fCE1Parameter5
#endif
    ;  
  
  if(!linestream)
    {
      delete telescope;
      return 0;
    }

  for(unsigned i=0; i<mirrors_size; i++)
    {
      VSOMirror* mirror = VSOMirror::createFromShortDump(stream, telescope);
      if(mirror==0)
	{
	  delete telescope;
	  return 0;
	}
      
      if(mirror->id() >= telescope->fMirrors.size())
	telescope->fMirrors.resize(mirror->id()+1);
      telescope->fMirrors[mirror->id()]=mirror;
      
      if(mirror->hexID() > telescope->fMirrorsByHexID.size())
	telescope->fMirrorsByHexID.resize(mirror->hexID());
      telescope->fMirrorsByHexID[mirror->hexID()-1]=mirror;
    }

  for(unsigned i=0; i<pixels_size; i++)
    {
      VSOPixel* pixel = VSOPixel::createFromShortDump(stream, telescope);
      if(pixel==0)
	{
	  delete telescope;
	  return 0;
	}
      
      if(pixel->id() >= telescope->fPixels.size())
	telescope->fPixels.resize(pixel->id()+1);
      telescope->fPixels[pixel->id()]=pixel;
      
      if(pixel->hexID() > telescope->fPixelsByHexID.size())
	telescope->fPixelsByHexID.resize(pixel->hexID());
      telescope->fPixelsByHexID[pixel->hexID()-1]=pixel;
    }

  // Recalculate rotation vector
  telescope->calculateRotationVector();
  
  return telescope;
}

