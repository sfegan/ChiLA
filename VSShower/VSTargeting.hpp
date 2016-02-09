//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTargeting.hpp

  Array and primary pointing algorithms

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/19/2005
*/

#ifndef VSTRACKING_HPP
#define VSTRACKING_HPP

#include <VSOTelescopeArray.hpp>
#include <VSDBParameterTable.hpp>
#include <VSLineTokenizer.hpp>
#include <RandomNumbers.hpp>

namespace VERITAS
{
  
  // --------------------------------------------------------------------------
  // Targeting
  // --------------------------------------------------------------------------
  
  class VSTargeting
  {
  public:
    virtual ~VSTargeting();
    virtual void setTarget(float& target_zn_rad, float& target_az_rad) = 0;    

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection) = 0;
  };
  
  class VSTargetingFixed: public VSTargeting
  {
  public:
    VSTargetingFixed(float zn, float az): 
      fFixedZnRad(zn), fFixedAzRad(az) { /* nothing to see here */ }
    VSTargetingFixed(const VSDBParameterSet& param);
    virtual ~VSTargetingFixed();
    virtual void setTarget(float& target_zn_rad, float& target_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "fixed"; }
  private:
    float fFixedZnRad;
    float fFixedAzRad;
  };

  class VSTargetingTracking: public VSTargeting
  {
  public:
    VSTargetingTracking() { /* nothing to see here */ }
    VSTargetingTracking(const VSDBParameterSet& param);
    virtual ~VSTargetingTracking();
    virtual void setTarget(float& target_zn_rad, float& target_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "tracking"; }
  };

  class VSTargetingOffsetTracking: public VSTargeting
  {
  public:
    VSTargetingOffsetTracking(const VSDBParameterSet& param);
    VSTargetingOffsetTracking(float theta, float phi):
      fOffsetThetaRad(theta), fOffsetPhiRad(phi) { /* nothing to see here */ }
    virtual ~VSTargetingOffsetTracking();
    virtual void setTarget(float& target_zn_rad, float& target_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "offsettracking"; }
  private:
    float fOffsetThetaRad;
    float fOffsetPhiRad;
  };

  class VSTargetingWobbleTracking: public VSTargeting
  {
  public:
    VSTargetingWobbleTracking(RandomNumbers& rng,
			      const VSDBParameterSet& param);
    VSTargetingWobbleTracking(RandomNumbers& rng, float theta):
      fRNG(rng), fWobbleThetaRad(theta) { /* nothing to see here */ }
    virtual ~VSTargetingWobbleTracking();
    virtual void setTarget(float& target_zn_rad, float& target_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "wobble"; }
  private:
    RandomNumbers& fRNG;
    float fWobbleThetaRad;
  };
    
  class VSTargetingIsotropicTracking: public VSTargeting
  {
  public:
    VSTargetingIsotropicTracking(RandomNumbers& rng,
				 const VSDBParameterSet& param);
    VSTargetingIsotropicTracking(RandomNumbers& rng, 
				 float out_radius_rad, 
				 float in_radius_rad=0):
      fRNG(rng), fInnerRadiusRad(in_radius_rad),
      fOuterRadiusRad(out_radius_rad) { /* nothing to see here */ }
    virtual ~VSTargetingIsotropicTracking();
    virtual void setTarget(float& target_zn_rad, float& target_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "isotropic"; }
  private:
    RandomNumbers& fRNG;
    float fInnerRadiusRad;
    float fOuterRadiusRad;
  };

  class VSTargetingFactory
  {
  public:
    virtual ~VSTargetingFactory();
    static VSTargetingFactory* getInstance();
    VSTargeting* getTargeting(RandomNumbers& rng,
			      const VSTokenList& tokens,
			      std::string* error_string = 0);
    VSTargeting* getTargeting(RandomNumbers& rng,
			      VSDBParameterTable* param_db,
			      const std::string& collection = "");
    void writeTargeting(VSTargeting* pointing,
			VSDBParameterTable* param_db,
			const std::string& collection = "");
  protected:
    VSTargetingFactory();
  private:
    std::string                fCollection;
    static std::string         sCollection;
    static VSTargetingFactory* sInstance;
  };

  // --------------------------------------------------------------------------
  // Array tracking
  // --------------------------------------------------------------------------
  
  class VSArrayTracking
  {
  public:
    virtual ~VSArrayTracking();
    virtual void pointArrayForEvent(VSOTelescopeArray* array,
				    float tar_zn_rad, float tar_az_rad) = 0;

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection) = 0;
  };

  class VSATParallel: public VSArrayTracking
  {
  public:
    VSATParallel();
    VSATParallel(const VSDBParameterSet& param);
    virtual ~VSATParallel();
    virtual void pointArrayForEvent(VSOTelescopeArray* array,
				    float tar_zn_rad, float tar_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "parallel"; }
  };

  class VSATFlysEye: public VSArrayTracking
  {
  public:
    VSATFlysEye(float spread_rad_per_m):
      fAngularSpreadRadPerM(spread_rad_per_m) { /* nothing to see here */ }
    VSATFlysEye(const VSDBParameterSet& param);
    virtual ~VSATFlysEye();
    virtual void pointArrayForEvent(VSOTelescopeArray* array,
				    float tar_zn_rad, float tar_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "flyseye"; }
  private:
    float fAngularSpreadRadPerM;
  };

  class VSATFactory
  {
  public:
    virtual ~VSATFactory();
    static VSATFactory* getInstance();
    VSArrayTracking* getArrayTracking(const VSTokenList& tokens,
				      std::string* error_string = 0);
    VSArrayTracking* getArrayTracking(VSDBParameterTable* param_db,
				      const std::string& collection = "");
    void writeArrayTracking(VSArrayTracking* tracking,
			    VSDBParameterTable* param_db,
			    const std::string& collection = "");
  protected:
    VSATFactory();
  private:
    std::string                fCollection;
    static std::string         sCollection;
    static VSATFactory*        sInstance;
  };

  // --------------------------------------------------------------------------
  // Primary arrival distribution
  // --------------------------------------------------------------------------

  class VSPrimaryArrivalDistribution
  {
  public:
    virtual ~VSPrimaryArrivalDistribution();
    virtual void setCORSIKAParameters(float& zn_lo, float& zn_hi,
				      float& az_lo, float& az_hi,
				      float& vc_lo, float& vc_hi) = 0;
    virtual void setTarget(float& target_zn_rad, float& target_az_rad) = 0;

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection) = 0;
  };

  class VSPADIsotropic: public VSPrimaryArrivalDistribution
  {
  public:
    VSPADIsotropic(float out_radius_rad, float in_radius_rad=0);
    VSPADIsotropic(const VSDBParameterSet& param);
    virtual ~VSPADIsotropic();
    virtual void setCORSIKAParameters(float& zn_lo, float& zn_hi,
				      float& az_lo, float& az_hi,
				      float& vc_lo, float& vc_hi);
    virtual void setTarget(float& target_zn_rad, float& target_az_rad);

    // Database I/O
    virtual void writeToDatabase(VSDBParameterTable* param_db,
				 const std::string& collection);
    static std::string name() { return "isotropic"; }
  private:
    float fInnerRadiusRad;
    float fOuterRadiusRad;
  };

  class VSPADFactory
  {
  public:
    virtual ~VSPADFactory();
    static VSPADFactory* getInstance();
    VSPrimaryArrivalDistribution* 
    getPrimaryArrivalDistribution(const VSTokenList& tokens,
				  std::string* error_string = 0);
    VSPrimaryArrivalDistribution* 
    getPrimaryArrivalDistribution(VSDBParameterTable* param_db,
				  const std::string& collection = "");
    void writePrimaryArrivalDistribution(VSPrimaryArrivalDistribution* pad,
					 VSDBParameterTable* param_db,
					 const std::string& collection = "");
  protected:
    VSPADFactory();
  private:
    std::string                fCollection;
    static std::string         sCollection;
    static VSPADFactory*       sInstance;
  };

} // namespace VERITAS

#endif // VSTRACKING_HPP
