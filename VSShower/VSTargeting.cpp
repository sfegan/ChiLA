//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTargeting.cpp

  Array and primary tracking algorithms

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/13/2005
*/

#include <cmath>
#include <cassert>
#include <sstream>

#include <Vec3D.hpp>

#include <VSTargeting.hpp>

using namespace VERITAS;
using namespace Physics;

// ----------------------------------------------------------------------------
// VSTargeting
// ----------------------------------------------------------------------------

VSTargeting::~VSTargeting()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// VSTargetingFixed
// ----------------------------------------------------------------------------

VSTargetingFixed::VSTargetingFixed(const VSDBParameterSet& param):
  fFixedZnRad(), fFixedAzRad()
{
  VSDBParameterSet::const_iterator iparam;
  
  iparam = param.find("zenith_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fFixedZnRad, iparam->second);

  iparam = param.find("azimuth_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fFixedAzRad, iparam->second);
}

VSTargetingFixed::~VSTargetingFixed()
{
  // nothing to see here
}

void VSTargetingFixed::setTarget(float& target_zn_rad, float& target_az_rad)
{
  target_zn_rad = fFixedZnRad;
  target_az_rad = fFixedAzRad;
};

void VSTargetingFixed::writeToDatabase(VSDBParameterTable* param_db,
				       const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  VSDataConverter::toString(param["zenith_rad"],fFixedZnRad);
  VSDataConverter::toString(param["azimuth_rad"],fFixedAzRad);
  param_db->storeParameterSet(collection, param);
}

// ----------------------------------------------------------------------------
// VSTargetingTracking
// ----------------------------------------------------------------------------

VSTargetingTracking::VSTargetingTracking(const VSDBParameterSet& param)
{
  // nothing to see here
}

VSTargetingTracking::~VSTargetingTracking()
{
  // nothing to see here
}

void VSTargetingTracking::setTarget(float& target_zn_rad, float& target_az_rad)
{
  // nothing to see here
}

void VSTargetingTracking::writeToDatabase(VSDBParameterTable* param_db,
					  const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  param_db->storeParameterSet(collection, param);
}

// ----------------------------------------------------------------------------
// VSTargetingOffsetTracking
// ----------------------------------------------------------------------------

VSTargetingOffsetTracking::
VSTargetingOffsetTracking(const VSDBParameterSet& param)
{
  VSDBParameterSet::const_iterator iparam;
  
  iparam = param.find("offset_theta_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fOffsetThetaRad, iparam->second);

  iparam = param.find("offset_phi_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fOffsetPhiRad, iparam->second);
}

VSTargetingOffsetTracking::~VSTargetingOffsetTracking()
{
  // nothing to see here
}

void VSTargetingOffsetTracking::
setTarget(float& target_zn_rad, float& target_az_rad)
{
  Vec3D axis(0,0,1);
  axis.Rotate(Vec3D(fOffsetThetaRad,0,0));
  axis.Rotate(Vec3D(0,0,fOffsetPhiRad));
  axis.Rotate(Vec3D(-target_zn_rad,0,0));
  axis.Rotate(Vec3D(0,0,-target_az_rad));
  target_zn_rad=atan2(sqrt(axis.x*axis.x+axis.y*axis.y), axis.z);
  target_az_rad=atan2(axis.x,axis.y);
}

void VSTargetingOffsetTracking::writeToDatabase(VSDBParameterTable* param_db,
						const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  VSDataConverter::toString(param["offset_theta_rad"],fOffsetThetaRad);
  VSDataConverter::toString(param["offset_phi_rad"],fOffsetPhiRad);
  param_db->storeParameterSet(collection, param);
}

// ----------------------------------------------------------------------------
// VSTargetingWobbleTracking
// ----------------------------------------------------------------------------

VSTargetingWobbleTracking::
VSTargetingWobbleTracking(RandomNumbers& rng, const VSDBParameterSet& param)
  : fRNG(rng), fWobbleThetaRad()
{
  VSDBParameterSet::const_iterator iparam;
  
  iparam = param.find("offset_theta_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fWobbleThetaRad, iparam->second);
}

VSTargetingWobbleTracking::~VSTargetingWobbleTracking()
{
  // nothing to see here
}

void VSTargetingWobbleTracking::
setTarget(float& target_zn_rad, float& target_az_rad)
{
  double phi_rad = fRNG.Uniform()*(2.0*M_PI);
  Vec3D axis(0,0,1);
  axis.Rotate(Vec3D(fWobbleThetaRad,0,0));
  axis.Rotate(Vec3D(0,0,phi_rad));
  axis.Rotate(Vec3D(-target_zn_rad,0,0));
  axis.Rotate(Vec3D(0,0,-target_az_rad));
  target_zn_rad=atan2(sqrt(axis.x*axis.x+axis.y*axis.y), axis.z);
  target_az_rad=atan2(axis.x,axis.y);
}

void VSTargetingWobbleTracking::writeToDatabase(VSDBParameterTable* param_db,
						const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  VSDataConverter::toString(param["offset_theta_rad"],fWobbleThetaRad);
  param_db->storeParameterSet(collection, param);
}

// ----------------------------------------------------------------------------
// VSTargetingIsotropicTracking
// ----------------------------------------------------------------------------

VSTargetingIsotropicTracking::
VSTargetingIsotropicTracking(RandomNumbers& rng, const VSDBParameterSet& param)
  : fRNG(rng), fInnerRadiusRad(), fOuterRadiusRad()
{
  VSDBParameterSet::const_iterator iparam;
  
  iparam = param.find("inner_radius_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fInnerRadiusRad, iparam->second);

  iparam = param.find("outer_radius_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fOuterRadiusRad, iparam->second);
}

VSTargetingIsotropicTracking::~VSTargetingIsotropicTracking()
{
  // nothing to see here
}

void VSTargetingIsotropicTracking::
setTarget(float& target_zn_rad, float& target_az_rad)
{
  double phi_rad = fRNG.Uniform()*(2.0*M_PI);
  double theta_rad = sqrt(std::pow(fInnerRadiusRad,2) + 
			  fRNG.Uniform()*(std::pow(fOuterRadiusRad,2) - 
					  std::pow(fInnerRadiusRad,2)));
  Vec3D axis(0,0,1);
  axis.Rotate(Vec3D(theta_rad,0,0));
  axis.Rotate(Vec3D(0,0,phi_rad));
  axis.Rotate(Vec3D(-target_zn_rad,0,0));
  axis.Rotate(Vec3D(0,0,-target_az_rad));
  target_zn_rad=atan2(sqrt(axis.x*axis.x+axis.y*axis.y), axis.z);
  target_az_rad=atan2(axis.x,axis.y);
}

void VSTargetingIsotropicTracking::
writeToDatabase(VSDBParameterTable* param_db,
		const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  VSDataConverter::toString(param["inner_radius_rad"],fInnerRadiusRad);
  VSDataConverter::toString(param["outer_radius_rad"],fOuterRadiusRad);
  param_db->storeParameterSet(collection, param);
}

// ----------------------------------------------------------------------------
// VSTargetingFactory
// ----------------------------------------------------------------------------

VSTargetingFactory* VSTargetingFactory::sInstance(0);
std::string VSTargetingFactory::sCollection("Targeting");

VSTargetingFactory::~VSTargetingFactory()
{
  // nothing to see here
}

VSTargetingFactory::VSTargetingFactory(): fCollection(sCollection)
{
  // nothing to see here
}

VSTargetingFactory* VSTargetingFactory::getInstance()
{
  if(sInstance==0)sInstance=new VSTargetingFactory;
  return sInstance;
}

VSTargeting* VSTargetingFactory::getTargeting(RandomNumbers& rng,
					      const VSTokenList& tokens,
					      std::string* error_string)
{
  if(tokens.size()<2)
    {
      if(error_string)
	*error_string = 
	  std::string("directive \"")
	  + tokens[0].lower()
	  + std::string("\" requires <Type> field");
      return 0;
    }

  std::map<std::string, unsigned> types;
  types[VSTargetingFixed::name()]=4;
  types[VSTargetingTracking::name()]=2;
  types[VSTargetingOffsetTracking::name()]=4;
  types[VSTargetingWobbleTracking::name()]=3;
  types[VSTargetingIsotropicTracking::name()]=4;

  std::string type = tokens[1].lower();
  if(types.find(type) == types.end())
    {
      if(error_string)
	{
	  std::ostringstream stream;
	  stream << "directive \"" << tokens[0].lower() << "\" "
		 << "type \"" << tokens[1].string() 
		 << "\" not known. Known values are";
	  for(std::map<std::string, unsigned>::const_iterator itype = 
		types.begin(); itype != types.end(); itype++)
	    stream << ' ' << '"' << itype->first << '"';
	  *error_string = stream.str();
	}
      return 0;
    }
  
  if(tokens.size()!=types[type])
    {
      if(error_string)
	{
	  std::ostringstream stream;
	  stream << "directive \"" << tokens[0].lower() << "\" "
		 << "type \"" << tokens[1].string() 
		 << "\" requires " << types[type]-2 << " values (found "
		 << tokens.size()-2 << ")";
	  *error_string = stream.str();
	}
      return 0;
    }

  if(type == VSTargetingFixed::name())
    {
      float zn;
      float az;
      tokens[2].convertTo(zn); zn *= M_PI/180.0;
      tokens[3].convertTo(az); az *= M_PI/180.0;
      return new VSTargetingFixed(zn,az);
    }
  else if(type == VSTargetingTracking::name())
    {
      return new VSTargetingTracking;
    }
  else if(type == VSTargetingOffsetTracking::name())
    {
      float theta;
      float phi;
      tokens[2].convertTo(theta); theta *= M_PI/180.0;
      tokens[3].convertTo(phi); phi *= M_PI/180.0;
      return new VSTargetingOffsetTracking(theta,phi);
    }
  else if(type == VSTargetingWobbleTracking::name())
    {
      float theta;
      tokens[2].convertTo(theta); theta *= M_PI/180.0;
      return new VSTargetingWobbleTracking(rng,theta);
    }
  else if(type == VSTargetingIsotropicTracking::name())
    {
      float outer_radius = 0;
      float inner_radius = 0;
      tokens[2].convertTo(inner_radius); inner_radius *= M_PI/180.0;     
      tokens[3].convertTo(outer_radius); outer_radius *= M_PI/180.0;
      return new VSTargetingIsotropicTracking(rng,outer_radius,inner_radius);
    }


  // Should never get here;
  assert(0);
  return 0;
}

VSTargeting* VSTargetingFactory::getTargeting(RandomNumbers& rng,
					      VSDBParameterTable* param_db,
					      const std::string& collection)
{
  std::string the_collection = collection;
  if(the_collection.empty())the_collection=fCollection;
  
  VSDBParameterSet param;
  if(param_db->retrieveParameterSet(the_collection, param) <= 0)return 0;
  if(param.find("type") == param.end())return 0;
  
  std::string type = param["type"];

  if(type == VSTargetingFixed::name())
    return new VSTargetingFixed(param);
  else if(type == VSTargetingTracking::name())
    return new VSTargetingTracking(param);
  else if(type == VSTargetingOffsetTracking::name())
    return new VSTargetingOffsetTracking(param);
  else if(type == VSTargetingWobbleTracking::name())
    return new VSTargetingWobbleTracking(rng,param);
  else if(type == VSTargetingIsotropicTracking::name())
    return new VSTargetingIsotropicTracking(rng,param);

  return 0;
}

void VSTargetingFactory::writeTargeting(VSTargeting* pointing,
					VSDBParameterTable* param_db,
					const std::string& collection)
{
  std::string the_collection = collection;
  if(the_collection.empty())the_collection=fCollection;
  pointing->writeToDatabase(param_db,the_collection);
}

// ----------------------------------------------------------------------------
// VSArrayTracking
// ----------------------------------------------------------------------------

VSArrayTracking::~VSArrayTracking()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// VSATParallel
// ----------------------------------------------------------------------------

VSATParallel::VSATParallel()
{
  // nothing to see here
}

VSATParallel::VSATParallel(const VSDBParameterSet& param)
{
  // nothing to see here
}

VSATParallel::~VSATParallel()
{
  // nothing to see here
}

void VSATParallel::
pointArrayForEvent(VSOTelescopeArray* array,
		   float tar_zn_rad, float tar_az_rad)
{
  array->pointTelescopesAzEl(tar_az_rad, M_PI_2-tar_zn_rad);
}

void  VSATParallel::writeToDatabase(VSDBParameterTable* param_db,
				    const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  param_db->storeParameterSet(collection, param);
}

// ----------------------------------------------------------------------------
// VSATFlysEye
// ----------------------------------------------------------------------------


VSATFlysEye::VSATFlysEye(const VSDBParameterSet& param):
  fAngularSpreadRadPerM()
{
  VSDBParameterSet::const_iterator iparam;
  
  iparam = param.find("angular_spread_rad_per_m");
  if(iparam != param.end())
    VSDataConverter::fromString(fAngularSpreadRadPerM, iparam->second);
}

VSATFlysEye::~VSATFlysEye()
{
  // nothing to see here
}

void VSATFlysEye::pointArrayForEvent(VSOTelescopeArray* array,
				     float tar_zn_rad, float tar_az_rad)
{
  unsigned ntel = array->numTelescopes();
  
  float x0 = array->telescope(0)->position().x;
  float y0 = array->telescope(0)->position().y;

  for(unsigned itel = 0; itel!=ntel; itel++)
    {
      float x = array->telescope(itel)->position().x;
      float y = array->telescope(itel)->position().y;
      float dx = (x-x0)/100;
      float dy = (y-y0)/100;
      float theta = fAngularSpreadRadPerM*sqrt(dx*dx+dy*dy);
      float phi = atan2(dx,dy);
      Vec3D axis(0,0,1);
      axis.Rotate(Vec3D(-theta,0,0));
      axis.Rotate(Vec3D(0,0,-phi));
      axis.Rotate(Vec3D(-tar_zn_rad,0,0));
      axis.Rotate(Vec3D(0,0,-tar_az_rad));
      array->telescope(itel)->pointTelescope(axis);
    }
}

void VSATFlysEye::writeToDatabase(VSDBParameterTable* param_db,
				  const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  VSDataConverter::toString(param["angular_spread_rad_per_m"],
			    fAngularSpreadRadPerM);
  param_db->storeParameterSet(collection, param);
};

// ----------------------------------------------------------------------------
// VSATFactory
// ----------------------------------------------------------------------------

VSATFactory* VSATFactory::sInstance(0);
std::string VSATFactory::sCollection("ArrayPointing");


VSATFactory::~VSATFactory()
{
  // nothing to see here
}

VSATFactory::VSATFactory()
  : fCollection(sCollection)
{
  // nothing to see here
}

VSATFactory* VSATFactory::getInstance()
{
  if(sInstance==0)sInstance=new VSATFactory;
  return sInstance;
}

VSArrayTracking* VSATFactory::getArrayTracking(const VSTokenList& tokens,
					       std::string* error_string)
{
  if(tokens.size()<2)
    {
      if(error_string)
	*error_string = 
	  std::string("directive \"")
	  + tokens[0].lower()
	  + std::string("\" requires <Type> field");
      return 0;
    }

  std::map<std::string, unsigned> types;
  types[VSATParallel::name()]=2;
  types[VSATFlysEye::name()]=3;

  std::string type = tokens[1].lower();
  if(types.find(type) == types.end())
    {
      if(error_string)
	{
	  std::ostringstream stream;
	  stream << "directive \"" << tokens[0].lower() << "\" "
		 << "type \"" << tokens[1].string() 
		 << "\" not known. Known values are";
	  for(std::map<std::string, unsigned>::const_iterator itype = 
		types.begin(); itype != types.end(); itype++)
	    stream << ' ' << '"' << itype->first << '"';
	  *error_string = stream.str();
	}
      return 0;
    }
  
  if(tokens.size()!=types[type])
    {
      if(error_string)
	{
	  std::ostringstream stream;
	  stream << "directive \"" << tokens[0].lower() << "\" "
		 << "type \"" << tokens[1].string() 
		 << "\" requires " << types[type]-2 << " values (found "
		 << tokens.size()-2 << ")";
	  *error_string = stream.str();
	}
      return 0;
    }

  if(type == VSATParallel::name())
    {
      return new VSATParallel;
    }
  else if(type == VSATFlysEye::name())
    {
      float spread_rad_per_m;
      tokens[2].convertTo(spread_rad_per_m); spread_rad_per_m *= M_PI/180.0;
      return new VSATFlysEye(spread_rad_per_m);
    }

  // Should never get here;
  assert(0);
  return 0;
}

VSArrayTracking* VSATFactory::getArrayTracking(VSDBParameterTable* param_db,
					       const std::string& collection)
{
  std::string the_collection = collection;
  if(the_collection.empty())the_collection=fCollection;
  
  VSDBParameterSet param;
  if(param_db->retrieveParameterSet(the_collection, param) <= 0)return 0;
  if(param.find("type") == param.end())return 0;
  
  std::string type = param["type"];

  if(type == VSATParallel::name())
    return new VSATParallel(param);
  else if(type == VSATFlysEye::name())
    return new VSATFlysEye(param);

  return 0;
}

void VSATFactory::writeArrayTracking(VSArrayTracking* tracking,
				     VSDBParameterTable* param_db,
				     const std::string& collection)
{
  std::string the_collection = collection;
  if(the_collection.empty())the_collection=fCollection;
  tracking->writeToDatabase(param_db,the_collection);
}

// ----------------------------------------------------------------------------
// VSPADIsotropic
// ----------------------------------------------------------------------------

VSPrimaryArrivalDistribution::~VSPrimaryArrivalDistribution()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// VSPADIsotropic
// ----------------------------------------------------------------------------

VSPADIsotropic::VSPADIsotropic(float out_radius_rad, float in_radius_rad)
  : fInnerRadiusRad(in_radius_rad), fOuterRadiusRad(out_radius_rad)
{
  // nothing to see here
}

VSPADIsotropic::VSPADIsotropic(const VSDBParameterSet& param)
  : fInnerRadiusRad(0), fOuterRadiusRad(0)
{
  VSDBParameterSet::const_iterator iparam;
  
  iparam = param.find("inner_radius_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fInnerRadiusRad, iparam->second);

  iparam = param.find("outer_radius_rad");
  if(iparam != param.end())
    VSDataConverter::fromString(fOuterRadiusRad, iparam->second);
}

VSPADIsotropic::~VSPADIsotropic()
{
  // nothing to see here
}

void VSPADIsotropic::setCORSIKAParameters(float& zn_lo, float& zn_hi,
					  float& az_lo, float& az_hi,
					  float& vc_lo, float& vc_hi)
{
  vc_lo = fInnerRadiusRad/M_PI*180.0;
  vc_hi = fOuterRadiusRad/M_PI*180.0;
}

void VSPADIsotropic::setTarget(float& target_zn_rad, float& target_az_rad)
{
  // nothing to see here
}

void VSPADIsotropic::writeToDatabase(VSDBParameterTable* param_db,
				     const std::string& collection)
{
  VSDBParameterSet param;
  param["type"]=name();
  VSDataConverter::toString(param["inner_radius_rad"],fInnerRadiusRad);
  VSDataConverter::toString(param["outer_radius_rad"],fOuterRadiusRad);
  param_db->storeParameterSet(collection, param);  
}

// ----------------------------------------------------------------------------
// VSPADFactory
// ----------------------------------------------------------------------------

VSPADFactory* VSPADFactory::sInstance(0);
std::string VSPADFactory::sCollection("PrimaryArrivalDistribution");


VSPADFactory::~VSPADFactory()
{
  // nothing to see here
}

VSPADFactory::VSPADFactory()
  : fCollection(sCollection)
{
  // nothing to see here
}

VSPADFactory* VSPADFactory::getInstance()
{
  if(sInstance==0)sInstance=new VSPADFactory;
  return sInstance;
}

VSPrimaryArrivalDistribution* 
VSPADFactory::getPrimaryArrivalDistribution(const VSTokenList& tokens,
					    std::string* error_string)
{
  if(tokens.size()<2)
    {
      if(error_string)
	*error_string = 
	  std::string("directive \"")
	  + tokens[0].lower()
	  + std::string("\" requires <Type> field");
      return 0;
    }

  std::map<std::string, unsigned> types;
  types[VSPADIsotropic::name()]=4;

  std::string type = tokens[1].lower();
  if(types.find(type) == types.end())
    {
      if(error_string)
	{
	  std::ostringstream stream;
	  stream << "directive \"" << tokens[0].lower() << "\" "
		 << "type \"" << tokens[1].string() 
		 << "\" not known. Known values are";
	  for(std::map<std::string, unsigned>::const_iterator itype = 
		types.begin(); itype != types.end(); itype++)
	    stream << ' ' << '"' << itype->first << '"';
	  *error_string = stream.str();
	}
      return 0;
    }
  
  if(tokens.size()!=types[type])
    {
      if(error_string)
	{
	  std::ostringstream stream;
	  stream << "directive \"" << tokens[0].lower() << "\" "
		 << "type \"" << tokens[1].string() 
		 << "\" requires " << types[type]-2 << " values (found "
		 << tokens.size()-2 << ")";
	  *error_string = stream.str();
	}
      return 0;
    }

  if(type == VSPADIsotropic::name())
    {
      float outer_radius = 0;
      float inner_radius = 0;
      tokens[2].convertTo(inner_radius); inner_radius *= M_PI/180.0;     
      tokens[3].convertTo(outer_radius); outer_radius *= M_PI/180.0;
      return new VSPADIsotropic(outer_radius,inner_radius);
    }

  // Should never get here;
  assert(0);
  return 0;
}

VSPrimaryArrivalDistribution* 
VSPADFactory::getPrimaryArrivalDistribution(VSDBParameterTable* param_db,
					    const std::string& collection)
{
  std::string the_collection = collection;
  if(the_collection.empty())the_collection=fCollection;
  
  VSDBParameterSet param;
  if(param_db->retrieveParameterSet(the_collection, param) <= 0)return 0;
  if(param.find("type") == param.end())return 0;
  
  std::string type = param["type"];

  if(type == VSPADIsotropic::name())
    return new VSPADIsotropic(param);

  return 0;
}

void VSPADFactory::
writePrimaryArrivalDistribution(VSPrimaryArrivalDistribution* pad,
				VSDBParameterTable* param_db,
				const std::string& collection)
{
  std::string the_collection = collection;
  if(the_collection.empty())the_collection=fCollection;
  pad->writeToDatabase(param_db,the_collection);
}
