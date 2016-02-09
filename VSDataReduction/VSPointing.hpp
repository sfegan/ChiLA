//-*-mode:c++; mode:font-lock;-*-

/*! \file VSPointing.hpp

  Various classes to supply pointing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/03/2006

  $Id: VSPointing.hpp,v 3.12 2008/07/21 21:35:13 matthew Exp $

*/

#ifndef VSPOINTING_HPP
#define VSPOINTING_HPP

#include<stdint.h>
#include<vector>
#include<memory>

#include <Astro.h>
#include <CorrectionParameters.h>
#include <VSCentralizedDBAccessData.hpp>
#include <VSOctaveIO.hpp>
#include <VSSimpleVBF.hpp>
#include <VSSimpleStat.hpp>

namespace VERITAS
{

  class VSPointing
  {
  public:
    VSPointing() { /* nothing ro see here */ }
    virtual ~VSPointing();
    virtual bool getMeanAzZn(SEphem::SphericalCoords& mean_az_zn,
			     double& rms_az_rad, double& rms_zn_rad,
			     double& mean_zn_rad,
			     const VSTime& lo_time, const VSTime& hi_time) = 0;
    virtual bool getMeanTarget(SEphem::SphericalCoords& mean_ra_dec,
			       double& rms_ra_rad, double& rms_dec_rad,
			       const VSTime& lo_time, const VSTime& hi_time, 
		       const SEphem::SphericalCoords& earth_position) = 0;
    virtual bool getAzZn(double& az, double& zn, 
			 unsigned iscope, const VSTime& time) = 0;

  private:
    VSPointing(const VSPointing&);
    VSPointing& operator= (const VSPointing&);
  };

  // --------------------------------------------------------------------------
  // INTERPOLATED POINTING
  // --------------------------------------------------------------------------

  template<class DataSource>
  class VSInterpolatedPointing: public VSPointing
  {
  public:
    typedef typename DataSource::Data Data;

    VSInterpolatedPointing(DataSource* source)
      : m_source(source), m_data(), m_ptr() { /* nothing to see here */ }
    VSInterpolatedPointing(const Data& data, DataSource* source=0)
      : m_source(source), m_data(data), m_ptr(data.scope.size()) 
    { /* nothing to see here */ }
    virtual ~VSInterpolatedPointing();
    virtual bool getMeanAzZn(SEphem::SphericalCoords& mean_zn_el,
			     double& rms_az_rad, double& rms_zn_rad,
			     double& mean_zn_rad,
			     const VSTime& lo_time, const VSTime& hi_time);
    virtual bool getMeanTarget(SEphem::SphericalCoords& mean_ra_dec,
			       double& rms_ra_rad, double& rms_dec_rad,
			       const VSTime& lo_time, const VSTime& hi_time, 
			       const SEphem::SphericalCoords& earth_position);
    virtual bool getAzZn(double& az, double& zn, 
			 unsigned iscope, const VSTime& time);

  private:
    VSInterpolatedPointing(const VSInterpolatedPointing&);
    VSInterpolatedPointing operator= (const VSInterpolatedPointing&);

    DataSource*           m_source;
    Data                  m_data;
    std::vector<unsigned> m_ptr;
  };

  // --------------------------------------------------------------------------
  // GET POINTING INFORMATION FROM DATABASE
  // --------------------------------------------------------------------------

  class VSDBPointingDataSource
  {
  public:

    struct Pointing
    {
      Pointing()
	: timestamp(), elevation_raw(), azimuth_raw(),
	  elevation_meas(), azimuth_meas(),
	  elevation_target(), azimuth_target() { }

      VSTime timestamp;
      float elevation_raw;
      float azimuth_raw;
      float elevation_meas;
      float azimuth_meas;
      float elevation_target;
      float azimuth_target;
    };

    typedef std::vector<Pointing> Pointings;

    struct Data
    {
      Data(): scope() { /* nothing to see here */ }

      typedef Pointing Datum;

      std::vector<Pointings> scope;

      void clear() { scope.clear(); }
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;

      void recorrect(const SEphem::CorrectionParameters& cp, unsigned iscope,
		     double* offset_x = 0, double* offset_y = 0);
    };

    VSDBPointingDataSource() { /* nothing to see here */ }
    virtual ~VSDBPointingDataSource();

    static void loadPointings(unsigned iscope, 
			      VSTime start, VSTime stop, Data& data);

  private:
    VSDBPointingDataSource(VSDBPointingDataSource&);
    VSDBPointingDataSource& operator= (VSDBPointingDataSource&);
  };

  typedef VSDBPointingDataSource::Data VSDBPointingData;

  class VSDBPointing: 
    public VSDBPointingDataSource,
    public VSInterpolatedPointing<VSDBPointingDataSource>
  {
  public:
    typedef VSDBPointingDataSource::Data Data;
    VSDBPointing()
      : VSDBPointingDataSource(),
	VSInterpolatedPointing<VSDBPointingDataSource>(this) { }
    VSDBPointing(const Data& data)
      : VSDBPointingDataSource(),
	VSInterpolatedPointing<VSDBPointingDataSource>(data) { }
    virtual ~VSDBPointing();

  private:
    VSDBPointing(const VSDBPointing&);
    VSDBPointing& operator= (const VSDBPointing&);
  };

  // --------------------------------------------------------------------------
  // GET POINTING INFORMATION FROM L3
  // --------------------------------------------------------------------------

  class VSL3PointingDataSource: public VSSimpleVBFVisitor
  {
  public:

    struct Pointing
    {
      Pointing()
	: timestamp(), elevation_meas(), azimuth_meas() { }

      VSTime timestamp;
      float elevation_meas;
      float azimuth_meas;
    };

    typedef std::vector<Pointing> Pointings;

    struct Data
    {
      Data(): scope() { /* nothing to see here */ }

      std::vector<Pointings> scope;

      void clear() { scope.clear(); }
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;
    };

    VSL3PointingDataSource();
    virtual ~VSL3PointingDataSource();

    virtual void visitSimulationHeader(void* user_data,
				       uint32_t         run_number,
				       const VSTime&    date_of_sims,
				       uint32_t         simulation_package,
				       const std::string& simulator_name,
				       const VSTime&    date_of_array,
				       uint32_t         corsika_atm_model,
				       float            obs_altitude_m,
				       const std::vector<VArrayConfiguration>&
				                        array,
				       const std::string& stringified_config);

    virtual void visitArrayEvent(bool& veto_array_event, void* user_data,
				 uint32_t               num_scope_events,
				 bool                   has_array_trigger,
				 bool                   has_event_number,
				 uint32_t               event_num,
 				 bool                   has_good_event_time,
				 const VSTime&          best_event_time,
				 EventType              event_type,
				 uint32_t               l2_trigger_mask,
				 const VArrayEvent*     array_event);

    virtual void visitArrayTriggerScope(bool& veto_array_event,
					void* user_data,
					uint32_t        telescope_num,
					bool            triggered,
					uint32_t        event_type,
					float           altitude,
					float           azimuth,
					uint32_t        tdc,
					uint32_t        shower_delay,
					uint32_t        comp_delay,
					const uint32_t* l2_counts,
					const uint32_t* cal_counts);

    void setIgnorePeds(bool ignore=true) { m_ignore_peds=ignore; }
    void getData(Data& data);
    void loadPointings(unsigned iscope, VSTime start, VSTime stop, Data& data);

  private:
    VSL3PointingDataSource(const VSL3PointingDataSource&);
    VSL3PointingDataSource& operator= (const VSL3PointingDataSource&);

    bool   m_ignore_peds;
    VSTime m_time;
    Data   m_data;
  };

  typedef VSL3PointingDataSource::Data VSL3PointingData;

  class VSL3Pointing: 
    public VSL3PointingDataSource,
    public VSInterpolatedPointing<VSL3PointingDataSource>
  {
  public:
    typedef VSL3PointingDataSource::Data Data;
    VSL3Pointing()
      : VSL3PointingDataSource(),
	VSInterpolatedPointing<VSL3PointingDataSource>(this) { }
    VSL3Pointing(const Data& data)
      : VSL3PointingDataSource(),
	VSInterpolatedPointing<VSL3PointingDataSource>(data) { }
    virtual ~VSL3Pointing();

  private:
    VSL3Pointing(const VSL3Pointing&);
    VSL3Pointing& operator= (const VSL3Pointing&);
  };

  // --------------------------------------------------------------------------
  // GET POINTING INFORMATION DIRECTLY FROM L3 RECORDS
  // --------------------------------------------------------------------------

  class VSDirectL3Pointing: public VSPointing
  {
  public:
    virtual ~VSDirectL3Pointing();
    virtual bool getMeanAzZn(SEphem::SphericalCoords& mean_az_zn,
			     double& rms_az_rad, double& rms_zn_rad,
			     double& mean_zn_rad,
			     const VSTime& lo_time, const VSTime& hi_time);
    virtual bool getMeanTarget(SEphem::SphericalCoords& mean_ra_dec,
			       double& rms_ra_rad, double& rms_dec_rad,
			       const VSTime& lo_time, const VSTime& hi_time, 
			       const SEphem::SphericalCoords& earth_position);
    virtual bool getAzZn(double& az, double& zn, 
			 unsigned iscope, const VSTime& time);
    void setAzElRad(unsigned iscope, double az, double el);
    static VSDirectL3Pointing* getInstance(); 

  protected:
    VSDirectL3Pointing();

  private:
    VSDirectL3Pointing(const VSDirectL3Pointing&);
    VSDirectL3Pointing& operator= (const VSDirectL3Pointing&);

    static std::auto_ptr<VSDirectL3Pointing> s_instance;

    std::vector<std::pair<double, double> > m_pointing;
  };

  // ==========================================================================
  // TARGET TABLE
  // ==========================================================================

  class VSTargetTable
  {
  public:

    struct Coord
    {
      Coord(): name(), radec(), epoch() 
      { /* nothing to see here */ }

      Coord(const VSCentralizedDBAccessData::TargetTableCoord& o):
	name(o.name), radec(), epoch(o.epoch)
      { 
	radec.setLatLongRad(o.dec_rad,o.ra_rad); 
	for(unsigned ichar=0;ichar<name.size();ichar++)
	  name[ichar] = tolower(name[ichar]);
      }

      std::string       name;
      SEphem::SphericalCoords radec;
      double            epoch;
    };

    class Data
    {
    public:
      Data(): coord() { /* nothing to see here */ }

      bool empty() const { return coord.empty(); }
    
      std::vector<Coord> coord;

      void clear() { coord.clear(); }
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;
    };

    struct Target
    {
      Target(): name(), coord() { /* nothing to see here */ }

      std::string name;
      SEphem::SphericalCoords coord;      
    };

    struct Observation
    {
      Observation(): 
	name("unknown"), 
	src_radec(), obs_radec(), src_radec_J2000(), obs_radec_J2000(), 
	mode(OM_UNKNOWN), mode_string("unknown"), 
	offset_min(), wobble_theta_rad(), wobble_phi_rad(),
	drift_azzn()
      { /* nothing to see here */ }
      enum ObservationMode { OM_UNKNOWN, OM_ON, OM_OFF, OM_WOBBLE,
			     OM_TRACKING, OM_DRIFT, OM_SIM_WOBBLE, 
			     OM_SIM_ON  };
      std::string name;
      SEphem::SphericalCoords src_radec;
      SEphem::SphericalCoords obs_radec;
      SEphem::SphericalCoords src_radec_J2000;
      SEphem::SphericalCoords obs_radec_J2000;
      ObservationMode mode;
      std::string mode_string;
      double offset_min;
      double wobble_theta_rad;
      double wobble_phi_rad;
      SEphem::SphericalCoords drift_azzn;

      void clear();
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;
    };

    VSTargetTable(const Data* data = 0);
    VSTargetTable(const std::vector<VSCentralizedDBAccessData::TargetTableCoord>& coords);
    
    bool empty() const { return m_coord.empty(); }
    const Data& data() { return m_coord; }
    
    SEphem::SphericalCoords
    getTarget(const std::string& name, const VSTime& approximate_start_time)
      const;

    SEphem::SphericalCoords
    getTarget(const std::string& name, const std::string& mode,
	      const VSTime& approximate_start_time) const;

    Observation getObservation(const SEphem::SphericalCoords& mean_radec,
			       double rms_ra_rad, double rms_dec_rad,
			       const SEphem::SphericalCoords& mean_azzn,
			       double rms_az_rad, double rms_zn_rad,
			       const std::string& demand_source_name,
			       const VSTime& approximate_start_time) const;

    std::vector<Target>
    getTargetForObservation(const Observation& observation, double theta_cut)
      const;

    void mergeCompiledTargets();

    static Data getCompiledTargets();

    static inline void rotate(const double theta_rad, double& x, double& y)
    {
      const double s = sin(theta_rad);
      const double c = cos(theta_rad);
      const double _x = x*c - y*s;
      const double _y = y*c + x*s;
      x = _x;
      y = _y;
    }

    static inline void wobble(const double theta_rad, const double phi_rad,
			      double& ra_rad, double& dec_rad)
    {
      static const double TWO_PI = 2.0*M_PI;  
      double x = 0;
      double y = 0;
      double z = 1;
      rotate(-theta_rad,z,x);
      rotate(phi_rad,y,x);
      rotate(M_PI_2-dec_rad,z,x);
      rotate(ra_rad,x,y);
      ra_rad = atan2(y,x);
      if(ra_rad<0)ra_rad+=TWO_PI;
      dec_rad = atan2(z,sqrt(x*x+y*y));
    }

    static inline void 
    wobble_inverse(const double src_ra_rad, double src_dec_rad,
		   const double trk_ra_rad, double trk_dec_rad,
		   double& theta_rad, double& phi_rad)
    {
      static const double TWO_PI = 2.0*M_PI;  
      double x = 0;
      double y = 0;
      double z = 1;
      rotate(M_PI_2-trk_dec_rad,z,x);
      rotate(trk_ra_rad,x,y);
      rotate(src_ra_rad,y,x);
      rotate(M_PI_2-src_dec_rad,x,z);
      phi_rad=atan2(y,-x);
      if(phi_rad<0)phi_rad+=TWO_PI;
      theta_rad=atan2(sqrt(x*x+y*y),z);
    }

  private:
    VSTargetTable(const VSTargetTable&);
    VSTargetTable& operator= (const VSTargetTable&);

    struct CoordTxt
    {
      std::string name;
      std::string ra;
      std::string dec;
      double epoch;
    };

    Data m_coord;

    static CoordTxt s_coord[];
  };

  typedef VSTargetTable::Data VSTargetTableData;

  // ==========================================================================
  // VSTargetPointing
  // ==========================================================================

  class VSTargetPointing: public VSPointing
  {
  public:
    typedef VSCentralizedDBAccessData::TrackingTargetDatum TargetDatum;
    typedef std::vector<TargetDatum> TargetData;

    typedef VSTargetTable::Observation Observation;

    VSTargetPointing(const std::vector<TargetData>& targets,
		     const SEphem::SphericalCoords& earth_position);
    virtual ~VSTargetPointing();
    virtual bool getMeanAzZn(SEphem::SphericalCoords& mean_az_zn,
			     double& rms_az_rad, double& rms_zn_rad,
			     double& mean_zn_rad,
			     const VSTime& lo_time, const VSTime& hi_time);
    virtual bool getMeanTarget(SEphem::SphericalCoords& mean_ra_dec,
			       double& rms_ra_rad, double& rms_dec_rad,
			       const VSTime& lo_time, const VSTime& hi_time, 
			       const SEphem::SphericalCoords& earth_position);

    virtual bool getAzZn(double& az, double& zn, 
			 unsigned iscope, const VSTime& time);

    bool getRaDec(double& ra, double& dec, 
		  unsigned iscope, const VSTime& time);

    bool getTargetAzZn(double& az, double& zn, 
		       unsigned iscope, const VSTime& time,
		       const TargetDatum& datum);

    bool getTargetRaDec(double& az, double& zn, 
			unsigned iscope, const VSTime& time,
			const TargetDatum& datum);

    bool getObservation(Observation& obs,
			const VSTime& run_start_time, 
			const VSTime& run_end_time) const;

  private:
    VSTargetPointing(const VSTargetPointing&);
    VSTargetPointing& operator= (const VSTargetPointing&);

    std::vector<TargetData> m_scope;
    SEphem::SphericalCoords m_earth_position;
  };
  
  // ==========================================================================
  // ==========================================================================
  // 
  // TEMPLATED CODE STARTS HERE
  //
  // ==========================================================================
  // ==========================================================================

  template<class DataSource> VSInterpolatedPointing<DataSource>::
  ~VSInterpolatedPointing()
  {
    // nothing to see here
  }

  template<class DataSource> bool VSInterpolatedPointing<DataSource>::
  getAzZn(double& az, double& zn, unsigned iscope, const VSTime& time)
  {
    if(iscope>=m_data.scope.size())
      {
	m_data.scope.resize(iscope+1);
	m_ptr.resize(iscope+1);
      }

    unsigned nscopedata = m_data.scope[iscope].size();

    if((m_source)
       &&((m_data.scope[iscope].empty())
	  ||(time<m_data.scope[iscope][0].timestamp)
	  ||(time>=m_data.scope[iscope][nscopedata-1].timestamp)))
      {
	VSTime start = time;
	if((!m_data.scope[iscope].empty())
	   &&(time>m_data.scope[iscope][nscopedata-1].timestamp))
	  start=m_data.scope[iscope][0].timestamp;
	start -= INT64_C(10000000000); // Push start time forward 10 seconds

	VSTime stop = time;
	if((!m_data.scope[iscope].empty())
	   &&(time<m_data.scope[iscope][0].timestamp))
	  stop=m_data.scope[iscope][nscopedata-1].timestamp;
	stop += INT64_C(10000000000); // Push stop time back 10 seconds

	// Get at least 35 mins and 35 min advanced past requested time
	static const int64_t MINGET = INT64_C(2100000000000);

	if(stop-start < MINGET) 
	  { 
	    stop = start; 
	    stop += MINGET;
	  } 

	if(stop-time < MINGET) 
	  { 
	    stop = time; 
	    stop += MINGET;
	  } 

	m_source->loadPointings(iscope, start, stop, m_data);

	m_ptr[iscope]=0;
	nscopedata = m_data.scope[iscope].size();
      }

    // Try to advance or retract the cached index by one, if
    // appropriate, hoping it will find the correct sample

    if((m_ptr[iscope]+2 < nscopedata)
       &&(time>m_data.scope[iscope][m_ptr[iscope]+1].timestamp))
      m_ptr[iscope]++;
    else if((m_ptr[iscope] > 0)
	    &&(time<m_data.scope[iscope][m_ptr[iscope]].timestamp))
      m_ptr[iscope]--;      

    // If the timestamp is still outside of the cached sample then do
    // a binary search

    if((m_ptr[iscope]+1 >= nscopedata)
       ||(time<m_data.scope[iscope][m_ptr[iscope]].timestamp)
       ||(time>m_data.scope[iscope][m_ptr[iscope]+1].timestamp))
      {
	unsigned lo = 0;
	unsigned hi = nscopedata;
	if(time<m_data.scope[iscope][lo].timestamp)
	  {
	    m_ptr[iscope]=0;
	    return false;
	  }
	while(hi-lo > 1)
	  {
	    unsigned test = (lo+hi)/2;
#ifdef DEBUG
	    std::cerr << lo << ' ' << test << ' ' << hi << ' '
		      << timestamp << ' ' 
		      << m_data.scope[iscope][test].timestamp << std::endl;
#endif
	    if(time>=m_data.scope[iscope][test].timestamp)lo=test;
	    else hi=test;
	  }
	m_ptr[iscope]=lo;
      }

    if(m_ptr[iscope]+1 >= nscopedata)return false;

    const VSTime t1 = m_data.scope[iscope][m_ptr[iscope]].timestamp;
    const VSTime t2 = m_data.scope[iscope][m_ptr[iscope]+1].timestamp;
    int64_t dT = t2-t1;
    int64_t dt = time-t1;
    double x = double(dt)/double(dT);

    SEphem::Angle del_az(m_data.scope[iscope][m_ptr[iscope]+1].azimuth_meas
			 - m_data.scope[iscope][m_ptr[iscope]].azimuth_meas);
    SEphem::Angle interp_az(m_data.scope[iscope][m_ptr[iscope]].azimuth_meas
			    + del_az.radPM180()*x);

    az = interp_az.rad();

    zn = m_data.scope[iscope][m_ptr[iscope]].elevation_meas*(1-x)
      + m_data.scope[iscope][m_ptr[iscope]+1].elevation_meas*x;
    zn = M_PI_2 - zn;
  
#ifdef DEBUG
    if(iscope == 0)
      std::cerr << timestamp << std::fixed << std::setprecision(7) 
		<< ' ' << x << ' ' << az/M_PI*180 << ' ' << zn/M_PI*180
		<< std::endl;
#endif

    return true;
  }

  template<class DataSource> bool VSInterpolatedPointing<DataSource>::
  getMeanAzZn(SEphem::SphericalCoords& mean_az_zn,
	      double& rms_az_rad, double& rms_zn_rad, double& mean_zn_rad,
	      const VSTime& lo_time, const VSTime& hi_time)
  {
    bool got_one = false;
    double phi_zero = 0;
    double theta_zero = 0;

    VSSimpleStat1<double> x_stat;
    VSSimpleStat1<double> y_stat;
    VSSimpleStat1<double> zn_stat;

    VSSimpleStat1<double> var_az_stat;
    VSSimpleStat1<double> var_zn_stat;

    unsigned nscope = m_data.scope.size();

    for(unsigned iscope=0;iscope<nscope;iscope++)
      {
	bool scope_got_one = false;
	double scope_phi_zero = 0;
	double scope_theta_zero = 0;
	VSSimpleStat2<double> scope_x_stat;
	VSSimpleStat2<double> scope_y_stat;

	unsigned npointing = m_data.scope[iscope].size();
	for(unsigned ipointing=0;ipointing<npointing;ipointing++)
	  {
	    if((m_data.scope[iscope][ipointing].timestamp < lo_time)
	       ||(m_data.scope[iscope][ipointing].timestamp > hi_time))
	      continue;

	    SEphem::Angle zn;
	    SEphem::Angle az;
	    zn.setCoAngleRad(m_data.scope[iscope][ipointing].elevation_meas);
	    az.setRad(m_data.scope[iscope][ipointing].azimuth_meas);
	    SEphem::SphericalCoords azzn(zn,az);

	    if(!got_one)
	      {
		phi_zero = azzn.phi();
		theta_zero = azzn.theta();
		got_one = true;
	      }

	    azzn.rotate(-phi_zero,-theta_zero,0);
	    
	    x_stat.accumulate(azzn.thetaRad()*cos(azzn.phiRad()));
	    y_stat.accumulate(azzn.thetaRad()*sin(azzn.phiRad()));
	    zn_stat.accumulate(zn.rad());

	    azzn.setThetaPhi(zn,az);

	    if(!scope_got_one)
	      {
		scope_phi_zero = azzn.phi();
		scope_theta_zero = azzn.theta();
		scope_got_one = true;
	      }

	    azzn.rotate(-scope_phi_zero,-scope_theta_zero,0);

	    scope_x_stat.accumulate(azzn.thetaRad()*cos(azzn.phiRad()));
	    scope_y_stat.accumulate(azzn.thetaRad()*sin(azzn.phiRad()));
	  }
	
	var_zn_stat.accumulate(scope_x_stat.dev(), scope_x_stat.count());
	var_az_stat.accumulate(scope_y_stat.dev(), scope_y_stat.count());
      }

    if(!got_one)return false;

    double x = x_stat.mean();
    double y = y_stat.mean();
    mean_az_zn.setRad(sqrt(x*x+y*y),atan2(y,x));
    mean_az_zn.rotate(0,theta_zero,phi_zero);

    rms_zn_rad = sqrt(var_zn_stat.mean());
    rms_az_rad = sqrt(var_az_stat.mean());

    mean_zn_rad = zn_stat.mean();

#ifdef DEBUG
    std::cout << mean_az_zn.thetaDeg() << ' ' 
	      << mean_az_zn.phiDeg() << ' ' 
	      << SEphem::Angle::toDeg(rms_zn_rad) << ' '
	      << SEphem::Angle::toDeg(rms_az_rad) << std::endl;
#endif
      
    return true;
  }

  template<class DataSource> bool VSInterpolatedPointing<DataSource>::
  getMeanTarget(SEphem::SphericalCoords& mean_ra_dec,
		double& rms_ra_rad, double& rms_dec_rad,
		const VSTime& lo_time, const VSTime& hi_time, 
		const SEphem::SphericalCoords& earth_position)
  {
    bool got_one = false;

    double phi_zero = 0;
    double theta_zero = 0;

    VSSimpleStat2<double> x_stat;
    VSSimpleStat2<double> y_stat;

    for(unsigned iscope=0;iscope<m_data.scope.size();iscope++)
      {
	unsigned npointing = m_data.scope[iscope].size();
	for(unsigned ipointing=0;ipointing<npointing;ipointing++)
	  {
	    if((m_data.scope[iscope][ipointing].timestamp < lo_time)
	       ||(m_data.scope[iscope][ipointing].timestamp > hi_time))
	      continue;

	    double mjd = m_data.scope[iscope][ipointing].timestamp.getMJDDbl();
	    SEphem::Angle lmst = 
	      SEphem::Astro::mjdToLMST(mjd, earth_position.longitudeRad());

	    SEphem::Angle zn;
	    SEphem::Angle az;
	    zn.setCoAngleRad(m_data.scope[iscope][ipointing].elevation_meas);
	    az.setRad(m_data.scope[iscope][ipointing].azimuth_meas);
	    SEphem::SphericalCoords radec(zn,az);
	    SEphem::Astro::azElToMeanRaDec(lmst, mjd, earth_position, radec);
	    
	    if(!got_one)
	      {
		phi_zero = radec.phi();
		theta_zero = radec.theta();
		got_one = true;
	      }

	    radec.rotate(-phi_zero,-theta_zero,0);
	    
	    x_stat.accumulate(radec.thetaRad()*cos(radec.phiRad()));
	    y_stat.accumulate(radec.thetaRad()*sin(radec.phiRad()));
	  }
      }

    if(!got_one)return false;

    double x = x_stat.mean();
    double y = y_stat.mean();
    mean_ra_dec.setRad(sqrt(x*x+y*y),atan2(y,x));
    mean_ra_dec.rotate(0,theta_zero,phi_zero);

    rms_dec_rad = y_stat.dev();
    rms_ra_rad = x_stat.dev();

    return true;
  }

} // namespace VERITAS

#endif // VSPOINTING_HPP
