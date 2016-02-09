//-*-mode:c++; mode:font-lock;-*-

/*! \file VSPointing.cpp

  Various classes to supply pointing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/03/2006

  $Id: VSPointing.cpp,v 3.27 2008/09/18 20:24:15 matthew Exp $

*/

#include<VSPointing.hpp>
#include<VSDataConverter.hpp>
#include<VSCentralizedDBAccess.hpp>

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// POINTING
// ============================================================================

VSPointing::~VSPointing()
{
  // nothing to see here
}

// ============================================================================
// GET POINTING INFORMATION FROM DATABASE
// ============================================================================

VSDBPointingDataSource::~VSDBPointingDataSource()
{
  // nothing to see here
}

void VSDBPointingDataSource::
loadPointings(unsigned iscope, VSTime start, VSTime stop, Data& data)
{
  if(data.scope.size() <= iscope)data.scope.resize(iscope+1);

  VSCentralizedDBAccess::getInstance()->
    getPointing(iscope, start, stop, data.scope[iscope]);
}

void VSDBPointingDataSource::Data::load(VSOctaveH5ReaderStruct* reader)
{
  scope.clear();

  VSOctaveH5ReaderCellVector* c = reader->readCellVector("scope");
  vsassert(c);
  unsigned nscope = c->dimensions();

  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(c->isStruct(iscope))
      {
	VSOctaveH5ReaderStruct* s = c->readStruct(iscope);
	
	std::vector<uint32_t> timestamp_mjd;
	std::vector<uint64_t> timestamp_dns;
	std::vector<float> elevation_raw;
	std::vector<float> azimuth_raw;
	std::vector<float> elevation_meas;
	std::vector<float> azimuth_meas;
	std::vector<float> elevation_target;
	std::vector<float> azimuth_target;
      
	s->readVector("timestamp_mjd",timestamp_mjd);
	s->readVector("timestamp_dns",timestamp_dns);
	s->readVector("elevation_raw",elevation_raw);
	s->readVector("azimuth_raw",azimuth_raw);
	s->readVector("elevation_meas",elevation_meas);
	s->readVector("azimuth_meas",azimuth_meas);
	s->readVector("elevation_target",elevation_target);
	s->readVector("azimuth_target",azimuth_target);

	unsigned np = timestamp_mjd.size();

	if(np!=0)
	  {
	    scope[iscope].resize(np);
	    Pointings& p(scope[iscope]);
	    for(unsigned ip=0;ip<np;ip++)
	      {
		p[ip].timestamp.setFromMJDIntAndNS(timestamp_mjd[ip],
						   timestamp_dns[ip]);
		p[ip].elevation_raw     = elevation_raw[ip];
		p[ip].azimuth_raw       = azimuth_raw[ip];
		p[ip].elevation_meas    = elevation_meas[ip];
		p[ip].azimuth_meas      = azimuth_meas[ip];
		p[ip].elevation_target  = elevation_target[ip];
		p[ip].azimuth_target    = azimuth_target[ip];
	      }
	  }

	delete s;
      }

  delete c;
}

void VSDBPointingDataSource::Data::save(VSOctaveH5WriterStruct* writer) const
{
  unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* c = writer->writeCellVector("scope",nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const Pointings& p(scope[iscope]);
      unsigned np = p.size();
      if(np)
	{
	  VSOctaveH5WriterStruct* s = c->writeStruct(iscope);
	  
	  std::vector<uint32_t> timestamp_mjd(np);
	  std::vector<uint64_t> timestamp_dns(np);
	  std::vector<float> elevation_raw(np);
	  std::vector<float> azimuth_raw(np);
	  std::vector<float> elevation_meas(np);
	  std::vector<float> azimuth_meas(np);
	  std::vector<float> elevation_target(np);
	  std::vector<float> azimuth_target(np);
	  
	  for(unsigned ip=0;ip<np;ip++)
	    {
	      timestamp_mjd[ip]     = p[ip].timestamp.getMJDInt();
	      timestamp_dns[ip]     = p[ip].timestamp.getDayNS();
 	      elevation_raw[ip]     = p[ip].elevation_raw;
	      azimuth_raw[ip]       = p[ip].azimuth_raw;
	      elevation_meas[ip]    = p[ip].elevation_meas;
	      azimuth_meas[ip]      = p[ip].azimuth_meas;
	      elevation_target[ip]  = p[ip].elevation_target;
	      azimuth_target[ip]    = p[ip].azimuth_target;
	    }

	  s->writeVector("timestamp_mjd",timestamp_mjd);
	  s->writeVector("timestamp_dns",timestamp_dns);
	  s->writeVector("elevation_raw",elevation_raw);
	  s->writeVector("azimuth_raw",azimuth_raw);
	  s->writeVector("elevation_meas",elevation_meas);
	  s->writeVector("azimuth_meas",azimuth_meas);
	  s->writeVector("elevation_target",elevation_target);
	  s->writeVector("azimuth_target",azimuth_target);

	  delete s;
	}
    }
  delete c;
}

void VSDBPointing::Data::
recorrect(const SEphem::CorrectionParameters& cp, unsigned iscope,
	  double* offset_x, double* offset_y)
{
  unsigned nscope = scope.size();
  if(iscope>=nscope)return;

  if(offset_x && offset_y)
    {
      *offset_x = 0;
      *offset_y = 0;
    }

  Pointings& p(scope[iscope]);
  unsigned np = p.size();
  for(unsigned ip=0;ip<np;ip++)
    {
      double el = p[ip].elevation_raw;
      double az = p[ip].azimuth_raw;
      cp.undoAzElCorrections(az, el, true);

      if(offset_x && offset_y)
	{
	  SphericalCoords c(M_PI_2-el,az);
	  c.rotate(0,p[ip].elevation_meas-M_PI_2,-p[ip].azimuth_meas);
	  *offset_y += c.theta()*cos(c.phi());
	  *offset_x -= c.theta()*sin(c.phi());
	}

      p[ip].elevation_meas = el;
      p[ip].azimuth_meas   = az;
    }

  if(offset_x && offset_y)
    {
      *offset_x /= double(np);
      *offset_y /= double(np);
    }
}

VSDBPointing::~VSDBPointing()
{
  // nothing to see here
}

// ============================================================================
// GET POINTING INFORMATION FROM L3
// ============================================================================

VSL3PointingDataSource::VSL3PointingDataSource()
  : VSSimpleVBFVisitor(), m_ignore_peds(), m_time(), m_data()
{
  // nothing to see here
}

VSL3PointingDataSource::~VSL3PointingDataSource()
{
  // nothing to see here
}

void VSL3PointingDataSource::
visitSimulationHeader(void* user_data,
		      uint32_t         run_number,
		      const VSTime&    date_of_sims,
		      uint32_t         simulation_package,
		      const std::string& simulator_name,
		      const VSTime&    date_of_array,
		      uint32_t         corsika_atm_model,
		      float            obs_altitude_m,
		      const std::vector<VArrayConfiguration>& array,
		      const std::string& stringified_config)
{
  m_ignore_peds = true;
}

void VSL3PointingDataSource::
visitArrayEvent(bool& veto_array_event, void* user_data,
		uint32_t               num_scope_events,
		bool                   has_array_trigger,
		bool                   has_event_number,
		uint32_t               event_num,
		bool                   has_good_event_time,
		const VSTime&          best_event_time,
		EventType              event_type,
		uint32_t               l2_trigger_mask,
		const VArrayEvent*     array_event)
{
  if((!has_good_event_time) 
     ||(m_ignore_peds && event_type == VSSimpleVBFVisitor::ET_PED))
    {
      veto_array_event=true;
      return;
    }
  m_time = best_event_time;
}

void VSL3PointingDataSource::
visitArrayTriggerScope(bool& veto_array_event,
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
		       const uint32_t* cal_counts)
{
  // This algorithm eliminates duplicate entries, since the L3
  // pointing is updated only periodically and keeps an extra
  // termination entry on the end of the list so that the
  // interpolation works correctly at the end point

  if(telescope_num >= m_data.scope.size())
    m_data.scope.resize(telescope_num+1);

  float az = Angle::frDeg(azimuth);
  float el = Angle::frDeg(altitude);

  Pointing p;
  p.timestamp      = m_time;
  p.elevation_meas = el;
  p.azimuth_meas   = az;
  
  if(m_data.scope[telescope_num].empty())
    {
      m_data.scope[telescope_num].push_back(p);
      m_data.scope[telescope_num].push_back(p);
    }    
  else if((m_data.scope[telescope_num].back().elevation_meas == el)
	  &&(m_data.scope[telescope_num].back().azimuth_meas == az))
    {
      m_data.scope[telescope_num].back() = p;
    }
  else
    {
      m_data.scope[telescope_num].back() = p;
      m_data.scope[telescope_num].push_back(p);
    }
}

void VSL3PointingDataSource::
getData(Data& data)
{
  data=m_data;
}

void VSL3PointingDataSource::
loadPointings(unsigned iscope, VSTime start, VSTime stop, Data& data)
{
  data=m_data;
}

void VSL3PointingDataSource::Data::load(VSOctaveH5ReaderStruct* reader)
{
  scope.clear();

  VSOctaveH5ReaderCellVector* c = reader->readCellVector("scope");
  vsassert(c);
  unsigned nscope = c->dimensions();

  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(c->isStruct(iscope))
      {
	VSOctaveH5ReaderStruct* s = c->readStruct(iscope);
	
	std::vector<uint32_t> timestamp_mjd;
	std::vector<uint64_t> timestamp_dns;
	std::vector<float> elevation_meas;
	std::vector<float> azimuth_meas;
      
	s->readVector("timestamp_mjd",timestamp_mjd);
	s->readVector("timestamp_dns",timestamp_dns);
	s->readVector("elevation_meas",elevation_meas);
	s->readVector("azimuth_meas",azimuth_meas);

	unsigned np = timestamp_mjd.size();
	if(np!=0)
	  {
	    scope[iscope].resize(np);
	    Pointings& p(scope[iscope]);
	    for(unsigned ip=0;ip<np;ip++)
	      {
		p[ip].timestamp.setFromMJDIntAndNS(timestamp_mjd[ip],
						   timestamp_dns[ip]);
		p[ip].elevation_meas    = elevation_meas[ip];
		p[ip].azimuth_meas      = azimuth_meas[ip];
	      }
	  }

	delete s;
      }

  delete c;
}

void VSL3PointingDataSource::Data::save(VSOctaveH5WriterStruct* writer) const
{
  unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* c = writer->writeCellVector("scope",nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const Pointings& p(scope[iscope]);
      unsigned np = p.size();
      if(np)
	{
	  VSOctaveH5WriterStruct* s = c->writeStruct(iscope);
	  
	  std::vector<uint32_t> timestamp_mjd(np);
	  std::vector<uint64_t> timestamp_dns(np);
	  std::vector<float> elevation_meas(np);
	  std::vector<float> azimuth_meas(np);
	  
	  for(unsigned ip=0;ip<np;ip++)
	    {
	      timestamp_mjd[ip]     = p[ip].timestamp.getMJDInt();
	      timestamp_dns[ip]     = p[ip].timestamp.getDayNS();
	      elevation_meas[ip]    = p[ip].elevation_meas;
	      azimuth_meas[ip]      = p[ip].azimuth_meas;
	    }
	  s->writeVector("timestamp_mjd",timestamp_mjd);
	  s->writeVector("timestamp_dns",timestamp_dns);
	  s->writeVector("elevation_meas",elevation_meas);
	  s->writeVector("azimuth_meas",azimuth_meas);

	  delete s;
	}
    }
  delete c;
}

VSL3Pointing::~VSL3Pointing()
{
  // nothing to see here
}

// ============================================================================
// GET POINTING INFORMATION DIRECTLY FROM L3 RECORDS
// ============================================================================

std::auto_ptr<VSDirectL3Pointing> VSDirectL3Pointing::s_instance;

VSDirectL3Pointing* VSDirectL3Pointing::getInstance()
{
  // WARNING: THIS CODE IS NOT THREAD SAFE
  if(!s_instance.get())s_instance.reset(new VSDirectL3Pointing);
  return s_instance.get();
}

VSDirectL3Pointing::VSDirectL3Pointing(): m_pointing()
{
  // nothing to see here
}

VSDirectL3Pointing::~VSDirectL3Pointing()
{
  // nothing to see here
}

bool VSDirectL3Pointing::
getMeanAzZn(SphericalCoords& mean_az_zn,
	    double& rms_az_rad, double& rms_zn_rad, double& mean_zn_rad,
	    const VSTime& lo_time, const VSTime& hi_time)
{
  return false;
  
#if 0  
  bool got_one = false;
  double phi_zero = 0;
  double theta_zero = 0;

  VSSimpleStat1<double> x_stat;
  VSSimpleStat1<double> y_stat;

  VSSimpleStat1<double> var_az_stat;
  VSSimpleStat1<double> var_zn_stat;

  unsigned nscope = m_pointing.size();

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      bool scope_got_one = false;
      double scope_phi_zero = 0;
      double scope_theta_zero = 0;
      VSSimpleStat2<double> scope_x_stat;
      VSSimpleStat2<double> scope_y_stat;

      unsigned npointing = 100;
      for(unsigned ipointing=0;ipointing<npointing;ipointing++)
	{
	  VSTime t(lo_time);
	  t += (hi_time-lo_time)/INT64_C(100)*int64_t(ipointing);
	  
	  double zn;
	  double az;
	  
	  if(getAzZn(az, zn, iscope, t))
	    {
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
	}

      if(scope_got_one)
	{
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

#ifdef DEBUG
      std::cout << mean_az_zn.thetaDeg() << ' ' 
		<< mean_az_zn.phiDeg() << ' ' 
		<< SEphem::Angle::toDeg(rms_zn_rad) << ' '
		<< SEphem::Angle::toDeg(rms_az_rad) << std::endl;
#endif
    }
      
  return true;
#endif
}

bool VSDirectL3Pointing::
getMeanTarget(SphericalCoords& mean_ra_dec,
	      double& rms_ra_rad, double& rms_dec_rad,
	      const VSTime& lo_time, const VSTime& hi_time, 
	      const SphericalCoords& earth_position)
{
  return false;
}

bool VSDirectL3Pointing::getAzZn(double& az, double& zn, 
			   unsigned iscope, const VSTime& time)
{
  if(iscope >= m_pointing.size())return false;
  az = m_pointing[iscope].first;
  zn = M_PI_2 - m_pointing[iscope].second;
  return true;
}

void VSDirectL3Pointing::setAzElRad(unsigned iscope, double az, double el)
{
  if(iscope >= m_pointing.size())m_pointing.resize(iscope+1);
  m_pointing[iscope].first = az;
  m_pointing[iscope].second = el;  
}

// ============================================================================
// TARGET TABLE
// ============================================================================

// perl -ne 'chomp; @bins=split /\s+/; $bins[1]=lc $bins[1]; $ls=15-length($bins[1]); $space=" "x($ls>0?$ls:0); print "    { \"$bins[1]\",$space \"$bins[2]\", \"$bins[3]\", $bins[4] },\n"' ~/tracking/VTracking/sources.trk

VSTargetTable::CoordTxt VSTargetTable::s_coord[] =
  {
    { "tycho",           "00:25:08.1", "+64:09:56", 2000 },
    { "1es0033+595",     "00:35:52.6", "+59:50:05", 2000 },
    { "1es0120+340",     "01:23:08.5", "+34:20:48", 2000 },
    { "m33",             "01:33:51.0", "+30:39:37", 2000 },
    { "3c58",            "02:05:38.0", "+64:49:42", 2000 },
    { "rgbj0214+517",    "02:14:17.9", "+51:44:52", 2000 },
    { "3c66a",           "02:22:39.6", "+43:02:08", 2000 },
    { "bl0224+014",      "02:27:16.6", "+02:01:58", 2000 },
    { "1es0229+200",     "02:32:48.4", "+20:17:16", 2000 },
    { "lsi+61303",       "02:40:31.7", "+61:13:46", 2000 },
    { "rxj0314.0+2445",  "03:14:02.7", "+24:44:33", 2000 },
    { "b20321+33b",      "03:24:41.2", "+34:10:46", 2000 },
    { "1es0323+022",     "03:26:14.0", "+02:25:15", 2000 },
    { "ic342",           "03:46:48.5", "+68:05:46", 2000 },
//  { "crab-0.5sw",      "05:32:52.4", "+21:41:43", 2000 },
    { "crab",            "05:34:31.9", "+22:00:52", 2000 },
//  { "crab-0.5n",       "05:34:31.9", "+22:30:52", 2000 },
//  { "crab-0.8n",       "05:34:31.9", "+22:48:52", 2000 },
//  { "crab-0.5ne",      "05:36:11.9", "+22:19:58", 2000 },
    { "theta1oric",      "05:35:16.5", "-05:23:23", 2000 },
    { "pks0607+17",      "06:09:24.0", "+17:19:00", 2000 },
    { "ic443",           "06:17:05.3", "+22:21:27", 2000 },
    { "psrj0631+1036",   "06:31:27.7", "+10:36:58", 2000 },
    { "geminga",         "06:33:54.2", "+17:46:13", 2000 },
    { "monoceros",       "06:34:51.6", "+05:56:05", 2000 },
    { "pks0646+06",      "06:48:41.0", "+06:26:42", 2000 },
    { "1es0647+250",     "06:50:46.5", "+25:03:00", 2000 },
    { "psrb0656+14",     "06:59:48.1", "+14:14:22", 2000 },
    { "1es0806+524",     "08:09:49.2", "+52:18:58", 2000 },
    { "mrk1218",         "08:38:10.9", "+24:53:43", 2000 },
    { "m82",             "09:55:52.2", "+69:40:49", 2000 },
    { "1es1101-232",     "11:03:37.6", "-23:29:30", 2000 },
    { "1es1011+496",     "10:15:04.1", "+49:26:01", 2000 },
    { "mrk421",          "11:04:27.3", "+38:12:32", 2000 },
    { "rxj1117.1+2014",  "11:17:06.2", "+20:14:07", 2000 },
    { "mrk180",          "11:36:26.4", "+70:09:27", 2000 },
    { "rxj1136.5+6737",  "11:36:30.1", "+67:37:04", 2000 },
    { "1es1218+304",     "12:21:21.9", "+30:10:37", 2000 },
    { "wcomae",          "12:21:31.7", "+28:13:58", 2000 },
    { "3c273",           "12:29:06.7", "+02:03:09", 2000 },
    { "m87",             "12:30:49.4", "+12:23:28", 2000 },
    { "ngc4589",         "12:37:25.2", "+74:11:31", 2000 },
    { "3c279",           "12:56:11.1", "-05:47:22", 2000 },
    { "coma cluster",    "12:59:48.7", "+27:58:50", 2000 },
    { "3egj1337+5029",   "13:35:00.0", "+51:00:00", 2000 },
    { "ngc5322",         "13:49:15.3", "+60:11:26", 2000 },
    { "rgbj1413+436",    "14:13:43.8", "+43:39:45", 2000 },
    { "h1426+428",       "14:28:32.6", "+42:40:29", 2000 },
    { "ursa minor",      "15:09:11.3", "+67:12:52", 2000 },
    { "1es1553+113",     "15:55:43.0", "+11:11:24", 2000 },
    { "1es1627+402",     "16:29:01.3", "+40:08:00", 2000 },
    { "m13",             "16:41:41.4", "+36:27:37", 2000 },
    { "mrk501",          "16:53:52.2", "+39:45:36", 2000 },
    { "draco",           "17:20:12.4", "+57:54:55", 2000 },
    { "rgbj1725+118",    "17:25:04.3", "+11:52:15", 2000 },
    { "1es1727+502",     "17:28:18.6", "+50:13:10", 2000 },
    { "psrj1740+1000",   "17:40:26.0", "+10:00:06", 2000 },
    { "1es1741+196",     "17:43:57.8", "+19:35:09", 2000 },
    { "psrj1837-0604",   "18:37:43.5", "-06:04:49", 2000 },
    { "psrj1841-0345",   "18:41:38.7", "-03:48:43", 2000 },
    { "psrj1841-0524",   "18:41:49.3", "-05:24:29", 2000 },
    { "psrj1846-0258",   "18:46:24.5", "-02:58:28", 2000 },
    { "3c391",           "18:49:22.4", "-00:55:21", 2000 },
    { "w44",             "18:56:08.3", "+01:17:58", 2000 },
    { "mgroj1853+01",    "18:53:34.0", "+01:02:23", 2000 },
    { "psrb1853+01",     "18:56:10.7", "+01:13:21", 2000 },
    { "psrj1857+0143",   "18:57:33.0", "+01:43:47", 2000 },
    { "g39.2-0.3",       "19:04:05.4", "+05:25:32", 2000 },
    { "psrj1907+0918",   "19:07:22.4", "+09:18:31", 2000 },
    { "mgroj1909+06",    "19:08:52.0", "+06:16:17", 2000 },
    { "psrj1909+0912",   "19:09:19.9", "+09:12:54", 2000 },
    { "psrj1913+1011",   "19:13:20.3", "+10:11:23", 2000 },
    { "w51",             "19:23:42.0", "+14:30:33", 2000 },
    { "psrj1928+1746",   "19:28:42.8", "+17:46:27", 2000 },
    { "psrj1930+1852",   "19:30:30.1", "+18:52:14", 2000 },
    { "psrb1930+22",     "19:32:22.7", "+22:20:52", 2000 },
    { "g54.4-0.3",       "19:33:38.2", "+18:55:14", 2000 },
    { "psrb1937+21",     "19:39:38.6", "+21:34:59", 2000 },
    { "psrb1951+32",     "19:52:58.2", "+32:52:41", 2000 },
    { "cygnus x-1",      "19:58:21.7", "+35:12:06", 2000 },
    { "psrb1957+20",     "19:59:36.8", "+20:48:15", 2000 },
    { "1es1959+650",     "19:59:59.9", "+65:08:55", 2000 },
    { "ctb87",           "20:16:04.8", "+37:13:33", 2000 },
    { "w63",             "20:19:00.0", "+45:42:00", 2000 },
    { "psrj2021+3651",   "20:21:04.5", "+36:51:27", 2000 },
    { "gamma cygni",     "20:22:13.7", "+40:15:24", 2000 },
    { "psrj2021+3651",   "20:21:04.5", "+36:51:27", 2000 },
    { "cygnus x-3",      "20:32:25.8", "+40:57:28", 2000 },
    { "tevj2032+4130",   "20:32:07.0", "+41:30:30", 2000 },
    { "wr147",           "20:36:43.7", "+40:21:07", 2000 },
    { "g89.0+4.7",       "20:45:43.1", "+50:40:26", 2000 },
    { "m15",             "21:29:58.4", "+12:10:01", 2000 },
    { "bllac",           "22:02:43.3", "+42:16:40", 2000 },
    { "pg2209+184",      "22:11:53.7", "+18:41:51", 2000 },
    { "psrj2229+6114",   "22:29:05.0", "+61:14:13", 2000 },
    { "mgroj2232+60",    "22:32:35.0", "+60:02:13", 2000 },
    { "1es2344+514",     "23:47:04.8", "+51:42:18", 2000 },
    { "casa",            "23:23:24.0", "+58:48:54", 2000 },
#if 1
    { "dark_001_8.58",   "16:26:52.8", "-05:42:00", 2000 },
    { "dark_002_8.55",   "12:11:49.9", "+36:48:00", 2000 },
    { "dark_003_8.55",   "00:33:36.6", "+01:24:00", 2000 },
    { "dark_004_8.43",   "03:15:43.9", "-53:48:00", 2000 },
    { "dark_005_8.43",   "15:26:13.7", "-03:36:00", 2000 },
    { "dark_006_8.39",   "19:29:36.0", "-38:06:00", 2000 },
    { "dark_007_8.37",   "00:54:15.7", "-40:36:00", 2000 },
    { "dark_008_8.36",   "00:13:48.9", "-29:42:00", 2000 },
    { "dark_009_8.31",   "22:21:17.5", "+03:12:00", 2000 },
    { "dark_010_8.30",   "23:54:55.0", "-07:12:00", 2000 },
    { "dark_011_8.27",   "00:35:43.8", "-13:06:00", 2000 },
    { "dark_012_8.27",   "15:11:26.9", "+14:54:00", 2000 },
    { "dark_013_8.26",   "07:10:12.2", "+42:06:00", 2000 },
    { "dark_014_8.26",   "23:31:22.8", "-18:18:00", 2000 },
    { "dark_015_8.26",   "10:41:57.8", "-04:00:00", 2000 },
    { "dark_016_8.26",   "12:03:51.1", "-16:24:00", 2000 },
    { "dark_017_8.26",   "11:38:31.0", "+15:06:00", 2000 },
    { "dark_018_8.25",   "02:15:28.1", "+10:30:00", 2000 },
    { "dark_019_8.25",   "12:20:59.8", "+11:42:00", 2000 },
    { "dark_020_8.24",   "17:07:18.2", "-08:12:00", 2000 },
    { "dark_021_8.24",   "13:33:17.0", "+28:54:00", 2000 },
    { "dark_022_8.23",   "05:07:21.6", "-28:42:00", 2000 },
    { "dark_023_8.21",   "10:50:49.0", "+20:48:00", 2000 },
    { "dark_024_8.21",   "20:41:11.3", "-53:42:00", 2000 },
    { "dark_025_8.20",   "10:05:15.4", "+26:00:00", 2000 },
    { "dark_026_8.20",   "14:28:19.0", "-18:18:00", 2000 },
    { "dark_027_8.20",   "08:46:27.8", "+23:30:00", 2000 },
    { "dark_028_8.19",   "13:48:12.7", "+02:12:00", 2000 },
    { "dark_029_8.19",   "00:54:10.7", "+10:54:00", 2000 },
    { "dark_030_8.16",   "22:19:47.0", "-37:00:00", 2000 },
    { "dark_031_8.16",   "16:57:48.0", "+04:24:00", 2000 },
    { "dark_032_8.16",   "12:56:57.1", "+41:48:00", 2000 },
    { "dark_033_8.12",   "01:17:13.0", "-27:00:00", 2000 },
    { "dark_034_8.11",   "14:22:45.6", "+23:18:00", 2000 },
    { "dark_035_8.10",   "10:54:15.6", "-24:42:00", 2000 },
    { "dark_036_8.10",   "03:41:38.6", "+12:00:00", 2000 },
    { "dark_037_8.10",   "00:44:27.8", "-67:48:00", 2000 },
    { "dark_038_8.10",   "13:41:07.0", "-20:48:00", 2000 },
    { "dark_039_8.09",   "16:00:29.0", "+61:54:00", 2000 },
    { "dark_040_8.09",   "18:30:21.8", "-12:42:00", 2000 },
    { "dark_041_8.08",   "01:50:02.2", "-08:36:00", 2000 },
    { "dark_042_8.06",   "15:39:03.8", "+28:12:00", 2000 },
    { "dark_043_8.06",   "23:30:11.3", "-48:00:00", 2000 },
    { "dark_044_8.05",   "16:04:37.7", "+19:48:00", 2000 },
    { "dark_045_8.05",   "14:25:48.7", "+73:18:00", 2000 },
    { "dark_046_8.04",   "13:21:23.8", "+10:48:00", 2000 },
    { "dark_047_8.03",   "08:04:18.5", "+30:48:00", 2000 },
    { "dark_048_8.02",   "23:12:22.8", "-37:18:00", 2000 },
    { "dark_049_8.02",   "16:40:06.7", "-15:42:00", 2000 },
    { "dark_050_8.02",   "13:25:51.6", "-08:36:00", 2000 },
    { "dark_051_8.01",   "02:51:40.7", "-17:12:00", 2000 },
    { "dark_052_8.01",   "14:28:22.1", "+33:48:00", 2000 },
    { "dark_053_7.97",   "23:12:41.0", "+22:24:00", 2000 },
    { "dark_054_7.97",   "16:34:37.0", "-31:54:00", 2000 },
    { "dark_055_7.96",   "13:15:16.1", "+60:36:00", 2000 },
    { "dark_056_7.93",   "03:12:16.5", "+04:48:00", 2000 },
    { "dark_057_7.92",   "07:55:57.1", "+52:24:00", 2000 },
    { "dark_058_7.92",   "02:22:26.5", "-32:18:00", 2000 },
    { "dark_059_7.92",   "21:17:12.0", "-34:06:00", 2000 },
    { "dark_060_7.91",   "09:37:11.5", "+18:18:00", 2000 },
    { "dark_061_7.91",   "04:20:42.2", "-27:54:00", 2000 },
    { "dark_062_7.89",   "10:47:05.0", "+35:30:00", 2000 },
    { "dark_063_7.89",   "04:31:42.6", "-65:06:00", 2000 },
    { "dark_064_7.88",   "12:45:43.0", "+01:00:00", 2000 },
    { "dark_065_7.88",   "15:00:14.4", "+42:42:00", 2000 },
    { "dark_066_7.87",   "17:07:44.9", "+17:54:00", 2000 },
    { "dark_067_7.85",   "01:27:20.9", "-55:24:00", 2000 },
    { "dark_068_7.85",   "14:40:08.9", "+05:00:00", 2000 },
    { "dark_069_7.85",   "00:44:43.6", "+36:24:00", 2000 },
    { "dark_070_7.84",   "15:51:00.7", "+51:06:00", 2000 },
    { "dark_071_7.84",   "04:28:33.4", "+48:36:00", 2000 },
    { "dark_072_7.84",   "06:58:56.9", "+65:24:00", 2000 },
    { "dark_073_7.83",   "04:33:12.0", "+03:06:00", 2000 },
    { "dark_074_7.82",   "02:05:44.5", "-86:54:00", 2000 },
    { "dark_075_7.82",   "04:03:16.2", "-14:48:00", 2000 },
    { "dark_076_7.82",   "03:53:29.1", "+28:42:00", 2000 },
    { "dark_077_7.81",   "11:14:31.7", "+04:06:00", 2000 },
    { "dark_078_7.80",   "20:47:11.3", "-07:24:00", 2000 },
    { "dark_079_7.79",   "10:00:36.2", "+01:30:00", 2000 },
    { "dark_080_7.79",   "02:24:03.5", "+21:06:00", 2000 },
    { "dark_081_7.79",   "22:21:15.4", "-18:36:00", 2000 },
    { "dark_082_7.79",   "21:39:21.4", "-02:12:00", 2000 },
    { "dark_083_7.79",   "09:51:31.2", "+43:12:00", 2000 },
    { "dark_084_7.78",   "05:22:45.6", "+26:12:00", 2000 },
    { "dark_085_7.77",   "05:06:48.1", "-46:36:00", 2000 },
    { "dark_086_7.77",   "19:38:40.8", "+65:12:00", 2000 },
    { "dark_087_7.76",   "10:50:00.0", "+48:18:00", 2000 },
    { "dark_088_7.76",   "17:45:49.7", "-11:24:00", 2000 },
    { "dark_089_7.75",   "20:16:32.9", "-23:36:00", 2000 },
    { "dark_090_7.75",   "19:25:23.8", "+06:48:00", 2000 },
    { "dark_091_7.74",   "16:57:55.0", "+37:12:00", 2000 },
    { "dark_092_7.74",   "10:43:19.0", "+61:48:00", 2000 },
    { "dark_093_7.74",   "12:02:38.6", "+24:48:00", 2000 },
    { "dark_094_7.74",   "08:35:03.8", "+44:24:00", 2000 },
    { "dark_095_7.72",   "09:21:30.2", "+30:18:00", 2000 },
    { "dark_096_7.72",   "19:59:36.2", "+54:42:00", 2000 },
    { "dark_097_7.71",   "14:11:51.8", "-38:24:00", 2000 },
    { "dark_098_7.71",   "03:34:22.2", "-38:24:00", 2000 },
    { "dark_099_7.71",   "10:15:18.2", "+11:00:00", 2000 },
    { "dark_100_7.71",   "18:18:37.7", "+09:06:00", 2000 },
    { "dark_101_7.71",   "18:41:57.6", "-68:30:00", 2000 },
    { "dark_102_7.71",   "09:40:34.8", "-07:12:00", 2000 },
    { "dark_103_7.71",   "12:32:16.3", "-37:06:00", 2000 },
    { "dark_104_7.70",   "17:42:17.0", "-25:24:00", 2000 },
    { "dark_105_7.70",   "21:38:04.1", "-60:30:00", 2000 },
    { "dark_106_7.70",   "23:48:38.2", "+15:18:00", 2000 },
    { "dark_107_7.69",   "03:25:45.4", "+37:42:00", 2000 },
    { "dark_108_7.69",   "05:32:16.5", "+46:06:00", 2000 },
    { "dark_109_7.67",   "17:48:40.8", "+58:42:00", 2000 },
    { "dark_110_7.67",   "08:54:48.0", "-00:06:00", 2000 },
    { "dark_111_7.66",   "12:40:54.0", "-10:24:00", 2000 },
    { "dark_112_7.66",   "15:20:07.0", "-34:24:00", 2000 },
    { "dark_113_7.66",   "08:07:46.3", "+07:30:00", 2000 },
    { "dark_114_7.66",   "11:34:00.2", "+71:06:00", 2000 },
    { "dark_115_7.65",   "13:52:18.2", "+44:42:00", 2000 },
    { "dark_116_7.65",   "22:29:13.0", "+22:12:00", 2000 },
    { "dark_117_7.64",   "07:55:09.4", "-07:24:00", 2000 },
    { "dark_118_7.64",   "06:51:06.0", "-75:00:00", 2000 },
    { "dark_119_7.64",   "04:41:55.2", "-10:30:00", 2000 },
    { "dark_120_7.64",   "21:36:59.8", "-48:18:00", 2000 },
    { "dark_121_7.63",   "02:15:01.5", "+36:36:00", 2000 },
    { "dark_122_7.62",   "23:06:43.0", "-02:18:00", 2000 },
    { "dark_123_7.62",   "16:19:52.3", "-60:24:00", 2000 },
    { "dark_124_7.62",   "09:58:37.7", "-74:06:00", 2000 },
    { "dark_125_7.62",   "11:26:15.6", "-33:36:00", 2000 },
    { "dark_126_7.61",   "19:34:34.8", "-79:48:00", 2000 },
    { "dark_127_7.61",   "00:09:38.0", "+70:36:00", 2000 },
    { "dark_128_7.61",   "21:27:29.8", "+73:48:00", 2000 },
    { "dark_129_7.61",   "01:46:35.1", "-43:54:00", 2000 },
    { "dark_130_7.60",   "11:52:23.5", "-01:54:00", 2000 },
    { "dark_131_7.60",   "05:35:50.1", "-13:54:00", 2000 },
    { "dark_132_7.60",   "23:08:16.1", "+12:18:00", 2000 },
    { "dark_133_7.59",   "08:36:17.8", "+62:36:00", 2000 },
    { "dark_134_7.59",   "03:05:00.9", "+71:36:00", 2000 },
    { "dark_135_7.59",   "07:15:05.5", "-02:06:00", 2000 },
    { "dark_136_7.58",   "14:22:43.0", "-04:12:00", 2000 },
    { "dark_137_7.57",   "08:59:43.4", "-30:48:00", 2000 },
    { "dark_138_7.57",   "21:14:57.4", "+14:12:00", 2000 },
    { "dark_139_7.56",   "19:46:36.2", "+15:24:00", 2000 },
    { "dark_140_7.56",   "06:35:07.7", "+50:36:00", 2000 },
    { "dark_141_7.54",   "19:14:40.8", "-03:36:00", 2000 },
    { "dark_142_7.54",   "05:55:46.1", "+81:24:00", 2000 },
    { "dark_143_7.54",   "21:15:06.0", "+28:06:00", 2000 },
    { "dark_144_7.51",   "01:52:10.5", "-19:06:00", 2000 },
    { "dark_145_7.51",   "04:42:43.8", "+35:42:00", 2000 },
    { "dark_146_7.50",   "12:34:58.8", "-25:24:00", 2000 },
    { "dark_147_7.49",   "09:57:51.4", "-37:48:00", 2000 },
    { "dark_148_7.48",   "00:25:00.3", "+21:54:00", 2000 },
    { "dark_149_7.48",   "18:37:00.5", "+50:00:00", 2000 },
    { "dark_150_7.48",   "06:55:57.6", "-63:30:00", 2000 },
    { "dark_151_7.48",   "13:03:37.7", "-45:18:00", 2000 },
    { "dark_152_7.47",   "20:59:51.6", "-22:24:00", 2000 },
    { "dark_153_7.47",   "13:24:50.2", "-30:12:00", 2000 },
    { "dark_154_7.45",   "05:19:06.2", "+14:36:00", 2000 },
    { "dark_155_7.45",   "01:41:02.4", "+23:18:00", 2000 },
    { "dark_156_7.45",   "06:15:16.1", "-32:18:00", 2000 },
    { "dark_157_7.45",   "06:58:25.2", "-37:24:00", 2000 },
    { "dark_158_7.45",   "22:36:27.8", "+33:48:00", 2000 },
    { "dark_159_7.44",   "14:23:55.2", "-71:54:00", 2000 },
    { "dark_160_7.44",   "23:46:15.4", "+33:06:00", 2000 },
    { "dark_161_7.43",   "10:44:52.6", "-45:06:00", 2000 },
    { "dark_162_7.43",   "18:36:36.5", "-26:00:00", 2000 },
    { "dark_163_7.43",   "22:06:41.3", "-73:30:00", 2000 },
    { "dark_164_7.43",   "17:38:05.8", "+82:12:00", 2000 },
    { "dark_165_7.42",   "12:05:17.3", "+49:24:00", 2000 },
    { "dark_166_7.42",   "03:38:34.5", "-23:00:00", 2000 },
    { "dark_167_7.42",   "11:52:53.3", "-54:18:00", 2000 },
    { "dark_168_7.41",   "15:19:10.1", "-50:06:00", 2000 },
    { "dark_169_7.41",   "16:03:14.6", "+02:24:00", 2000 },
    { "dark_170_7.41",   "01:35:52.9", "+04:24:00", 2000 },
    { "dark_171_7.41",   "18:04:11.3", "-53:54:00", 2000 },
    { "dark_172_7.40",   "07:20:26.2", "+13:12:00", 2000 },
    { "dark_173_7.40",   "23:15:27.4", "+55:24:00", 2000 },
    { "dark_174_7.38",   "00:28:31.6", "+46:36:00", 2000 },
    { "dark_175_7.36",   "15:06:57.6", "-13:54:00", 2000 },
    { "dark_176_7.36",   "09:00:13.2", "-13:54:00", 2000 },
    { "dark_177_7.35",   "17:49:48.0", "+42:24:00", 2000 },
    { "dark_178_7.35",   "18:13:20.2", "-00:54:00", 2000 },
    { "dark_179_7.32",   "21:49:46.6", "-25:36:00", 2000 },
    { "dark_180_7.32",   "20:20:57.4", "-44:06:00", 2000 },
    { "dark_181_7.30",   "23:39:21.6", "+03:42:00", 2000 },
    { "dark_182_7.29",   "15:23:19.2", "-23:48:00", 2000 },
    { "dark_183_7.29",   "06:19:38.5", "+20:30:00", 2000 },
    { "dark_184_7.29",   "21:01:13.4", "+61:12:00", 2000 },
    { "dark_185_7.29",   "23:53:46.6", "+81:24:00", 2000 },
    { "dark_186_7.28",   "21:58:21.4", "+50:36:00", 2000 },
    { "dark_187_7.27",   "06:51:35.0", "+29:54:00", 2000 },
    { "dark_188_7.27",   "02:59:23.6", "+14:36:00", 2000 },
    { "dark_189_7.27",   "08:25:26.4", "-61:36:00", 2000 },
    { "dark_190_7.26",   "18:28:07.2", "+32:54:00", 2000 },
    { "dark_191_7.26",   "18:09:37.4", "-39:06:00", 2000 },
    { "dark_192_7.25",   "20:18:25.4", "+06:54:00", 2000 },
    { "dark_193_7.25",   "12:47:22.3", "+29:18:00", 2000 },
    { "dark_194_7.24",   "03:54:45.1", "-76:24:00", 2000 },
    { "dark_195_7.24",   "16:32:11.3", "+28:18:00", 2000 },
    { "dark_196_7.24",   "02:43:08.2", "+51:18:00", 2000 },
    { "dark_197_7.24",   "19:57:46.6", "-12:00:00", 2000 },
    { "dark_198_7.23",   "19:18:45.8", "-58:42:00", 2000 },
    { "dark_199_7.22",   "23:07:05.0", "-64:54:00", 2000 },
    { "dark_200_7.22",   "11:42:40.8", "+83:12:00", 2000 },
    { "dark_201_7.21",   "01:20:47.3", "+53:54:00", 2000 },
    { "dark_202_7.20",   "09:31:37.0", "+53:48:00", 2000 },
    { "dark_203_7.20",   "12:42:49.2", "-81:48:00", 2000 },
    { "dark_204_7.20",   "18:16:41.0", "+19:54:00", 2000 },
    { "dark_205_7.19",   "08:14:43.0", "-25:06:00", 2000 },
    { "dark_206_7.18",   "21:32:28.1", "-16:00:00", 2000 },
    { "dark_207_7.17",   "10:18:39.8", "-16:06:00", 2000 },
    { "dark_208_7.17",   "07:27:11.0", "-45:54:00", 2000 },
    { "dark_209_7.15",   "00:15:38.0", "-57:30:00", 2000 },
    { "dark_210_7.15",   "01:11:25.8", "-04:36:00", 2000 },
    { "dark_211_7.15",   "02:54:22.4", "-05:24:00", 2000 },
    { "dark_212_7.13",   "04:30:55.5", "-38:06:00", 2000 },
    { "dark_213_7.10",   "14:34:04.3", "-28:18:00", 2000 },
    { "dark_214_7.09",   "11:10:03.6", "-15:24:00", 2000 },
    { "dark_215_7.09",   "05:15:34.1", "-04:00:00", 2000 },
    { "dark_216_7.08",   "09:19:36.0", "+09:12:00", 2000 },
    { "dark_217_7.07",   "19:51:06.5", "+43:30:00", 2000 },
    { "dark_218_7.07",   "04:02:32.0", "-03:48:00", 2000 },
    { "dark_219_7.07",   "13:29:21.8", "-54:30:00", 2000 },
    { "dark_220_7.06",   "17:35:46.1", "+07:48:00", 2000 },
    { "dark_221_7.05",   "20:26:42.2", "+28:42:00", 2000 },
    { "dark_222_7.04",   "02:12:32.4", "+62:30:00", 2000 },
    { "dark_223_7.04",   "06:21:42.5", "+37:54:00", 2000 },
    { "dark_224_7.03",   "06:08:00.0", "+00:06:00", 2000 },
    { "dark_225_7.03",   "04:45:49.7", "+20:48:00", 2000 },
    { "dark_226_7.03",   "14:49:02.2", "+54:54:00", 2000 },
    { "dark_227_7.02",   "19:12:49.4", "-14:12:00", 2000 },
    { "dark_228_7.01",   "15:55:24.5", "-12:36:00", 2000 },
    { "dark_229_7.01",   "14:03:46.3", "+13:24:00", 2000 },
    { "dark_230_7.00",   "20:51:26.2", "+38:06:00", 2000 },
    { "dark_231_6.99",   "16:47:05.8", "+70:48:00", 2000 },
    { "dark_232_6.99",   "21:58:46.1", "+40:30:00", 2000 },
    { "dark_233_6.98",   "19:15:50.6", "+24:54:00", 2000 },
    { "dark_234_6.98",   "09:01:21.6", "+75:48:00", 2000 },
    { "dark_235_6.97",   "07:09:04.8", "-14:24:00", 2000 },
    { "dark_236_6.97",   "10:04:41.0", "-27:54:00", 2000 },
    { "dark_237_6.96",   "11:20:05.0", "+28:18:00", 2000 },
    { "dark_238_6.92",   "05:07:32.7", "+56:48:00", 2000 },
    { "dark_239_6.89",   "17:06:07.4", "-40:54:00", 2000 },
    { "dark_240_6.89",   "16:55:49.4", "+51:54:00", 2000 },
    { "dark_241_6.86",   "16:39:43.4", "-75:24:00", 2000 },
    { "dark_242_6.86",   "12:49:29.3", "-65:00:00", 2000 },
    { "dark_243_6.85",   "22:28:24.7", "-08:30:00", 2000 },
    { "dark_244_6.84",   "19:18:18.2", "-48:36:00", 2000 },
    { "dark_245_6.84",   "06:39:17.5", "-06:54:00", 2000 },
    { "dark_246_6.84",   "06:16:32.4", "-47:00:00", 2000 },
    { "dark_247_6.83",   "22:52:24.0", "-27:48:00", 2000 },
    { "dark_248_6.79",   "02:24:17.1", "-62:36:00", 2000 },
    { "dark_249_6.78",   "02:47:41.0", "-44:30:00", 2000 },
    { "dark_250_6.77",   "06:48:12.7", "+05:24:00", 2000 },
    { "dark_251_6.76",   "05:37:29.7", "-56:48:00", 2000 },
    { "dark_252_6.76",   "03:11:38.5", "+24:12:00", 2000 },
    { "dark_253_6.76",   "02:14:00.1", "+00:18:00", 2000 },
    { "dark_254_6.74",   "20:27:56.4", "-33:48:00", 2000 },
    { "dark_255_6.72",   "16:05:12.0", "-43:12:00", 2000 },
    { "dark_256_6.72",   "19:36:30.7", "-28:00:00", 2000 },
    { "dark_257_6.71",   "21:56:32.4", "+15:48:00", 2000 },
    { "dark_258_6.71",   "19:58:35.5", "-01:48:00", 2000 },
    { "dark_259_6.70",   "20:56:23.8", "+02:30:00", 2000 },
    { "dark_260_6.70",   "08:54:46.8", "-50:18:00", 2000 },
    { "dark_261_6.69",   "16:06:19.4", "+40:54:00", 2000 },
    { "dark_262_6.67",   "22:33:49.2", "-52:42:00", 2000 },
    { "dark_263_6.64",   "17:41:58.1", "+29:24:00", 2000 },
    { "dark_264_6.62",   "04:08:05.0", "-48:00:00", 2000 },
    { "dark_265_6.62",   "11:51:17.0", "-43:24:00", 2000 },
    { "dark_266_6.60",   "07:39:15.4", "-31:24:00", 2000 },
    { "dark_267_6.59",   "10:04:03.4", "-64:00:00", 2000 },
    { "dark_268_6.58",   "06:27:24.6", "-22:00:00", 2000 },
    { "dark_269_6.57",   "05:25:58.1", "+70:48:00", 2000 },
    { "dark_270_6.54",   "07:30:07.0", "+24:00:00", 2000 },
    { "dark_271_6.43",   "06:09:43.8", "+10:06:00", 2000 },
    { "dark_272_6.42",   "12:50:17.8", "+18:48:00", 2000 },
    { "dark_273_6.31",   "08:31:15.4", "-40:42:00", 2000 },
    { "dark_274_6.30",   "19:16:02.4", "+35:42:00", 2000 },
    { "dark_275_6.29",   "03:40:30.1", "+57:24:00", 2000 },
    { "dark_276_6.25",   "23:19:45.1", "+44:06:00", 2000 },
    { "dark_277_6.23",   "20:31:17.8", "+16:24:00", 2000 },
    { "dark_278_6.21",   "08:07:38.6", "+19:06:00", 2000 },
    { "dark_279_6.17",   "05:32:10.9", "+36:06:00", 2000 },
    { "dark_280_6.16",   "10:10:54.0", "-53:48:00", 2000 },
    { "dark_281_6.11",   "16:24:33.8", "+10:54:00", 2000 },
    { "dark_282_6.10",   "14:35:29.3", "-59:06:00", 2000 },
    { "dark_283_6.03",   "00:02:04.8", "-39:42:00", 2000 },
    { "dark_284_5.98",   "01:28:41.6", "+42:18:00", 2000 },
    { "dark_285_5.86",   "16:10:54.7", "-23:12:00", 2000 },
    { "dark_286_5.74",   "00:22:28.1", "+60:06:00", 2000 },
    { "dark_287_5.74",   "09:37:51.6", "-18:30:00", 2000 },
    { "dark_288_5.62",   "08:41:46.3", "+13:00:00", 2000 },
    { "dark_289_5.62",   "16:46:25.7", "-51:06:00", 2000 },
    { "dark_290_5.60",   "20:48:51.1", "+48:12:00", 2000 },
    { "dark_291_5.52",   "22:24:21.6", "+64:48:00", 2000 },
    { "dark_292_5.45",   "05:31:38.3", "-37:42:00", 2000 },
    { "dark_293_5.27",   "18:55:26.2", "+13:36:00", 2000 },
    { "dark_294_5.14",   "20:30:21.1", "-67:18:00", 2000 },
    { "dark_295_4.20",   "05:33:22.3", "+05:12:00", 2000 },
#endif
  };

VSTargetTable::VSTargetTable(const Data* data): m_coord()
{
  if(data)m_coord = *data;
}

VSTargetTable::VSTargetTable(const std::vector<VSCentralizedDBAccessData::TargetTableCoord>& coords)
{
  unsigned ncoord = coords.size();
  for(unsigned icoord=0;icoord<ncoord;icoord++)
    if(coords[icoord].epoch)m_coord.coord.push_back(coords[icoord]);
}

SphericalCoords VSTargetTable::
getTarget(const std::string& name, const VSTime& approximate_start_time) const
{
  static const char sep = '@';

  bool matched = false;
  SphericalCoords radec;
  double epoch;

  for(unsigned icoord=0; icoord<m_coord.coord.size(); icoord++)
    {
      if(name == m_coord.coord[icoord].name)
	{
	  radec = m_coord.coord[icoord].radec;
	  epoch = m_coord.coord[icoord].epoch;
	  matched = true;
	  break;
	}
    }
  
  if(!matched)
    {
      // This mouthfull extracts user-specified coordinates
      std::string::const_iterator find1 = name.begin();
      while((find1!=name.end())&&(*find1!=sep))find1++;
      if(find1==name.end())goto user_specified_target_exit;
      std::string ra_string(name.begin(),find1);
      find1++;
      std::string::const_iterator find2 = find1;
      while((find2!=name.end())&&(*find2!=sep))find2++;
      std::string dec_string(find1,find2);
      epoch = 2000.0; // default
      matched = true;
      if(find2==name.end())goto user_specified_target_exit;
      find2++;
      std::string epoch_string(find2,name.end());
      VSDataConverter::fromString(epoch, epoch_string);

      Angle ra;
      Angle dec;
      if((!ra.setFromHMSString(ra_string))
	 ||(!dec.setCoAngleFromDMSString(dec_string)))
	{
	  std::cerr << "Could not extract RA and DEC from: " 
		    << ra_string << " and " << dec_string
		    << std::endl;
	  exit(EXIT_FAILURE);	      
	}

      radec.setLatLong(dec, ra);
    }
 user_specified_target_exit:
      
  if(!matched)
    {
      std::cerr << "Unknown target: " << name
		<< std::endl;
      exit(EXIT_FAILURE);
    }
      
  Astro::precess(approximate_start_time.getMJDDbl(), radec, 
		 Astro::julianEpochToMJD(epoch));
  
  return radec;
}

SphericalCoords VSTargetTable::
getTarget(const std::string& name, const std::string& mode,
	  const VSTime& approximate_start_time) const
{
  SphericalCoords c = getTarget(name, approximate_start_time);

  if(mode == "on")
    {
      // nothing to see here
    }
  else if(mode.substr(0,6) == "wobble")
    {
      // This junk takes the offset theta and phi out of the offset
      static const char sep = '@';
      const std::string offset_data(mode.substr(6));
      std::string::const_iterator find1 = offset_data.begin();
      while((find1!=offset_data.end())&&(*find1!=sep))find1++;
      if(find1==offset_data.end())
	{
	  std::cerr << "Could not understand offset: " << offset_data
		    << std::endl;
	  exit(EXIT_FAILURE);	      
	}
      std::string theta_string(offset_data.begin(),find1);
      find1++;
      std::string phi_string(find1,offset_data.end());
      double theta = 0;
      double phi = 0;
      VSDataConverter::fromString(theta, theta_string);
      VSDataConverter::fromString(phi, phi_string);
      theta *= M_PI/180.0;
      phi *= M_PI/180.0;
      double ra_rad = c.longitudeRad();
      double dec_rad = c.latitudeRad();
      wobble(theta,phi,ra_rad,dec_rad);
      c.setLatLongRad(dec_rad,ra_rad);
    }
  else if(mode.substr(0,6) == "woboff")
    {
      // This junk takes the wobble theta and phi and psi out of the offset
      static const char sep = '@';
      const std::string wobble_data(mode.substr(6));
      std::string::const_iterator find1 = wobble_data.begin();
      while((find1!=wobble_data.end())&&(*find1!=sep))find1++;
      if(find1==wobble_data.end())
	{
	  std::cerr << "Could not understand offset: " << wobble_data
		    << std::endl;
	  exit(EXIT_FAILURE);	      
	}
      std::string theta_string(wobble_data.begin(),find1);
      find1++;
      std::string::const_iterator find2 = find1;
      while((find2!=wobble_data.end())&&(*find2!=sep))find2++;
      std::string phi_string(find1,find2);
      double theta = 0;
      double phi = 0;
      VSDataConverter::fromString(theta, theta_string);
      VSDataConverter::fromString(phi, phi_string);
      double off = 0;      
      if(find2!=wobble_data.end())
	{
	  find2++;
	  std::string off_string(find2,wobble_data.end());
	  VSDataConverter::fromString(off, off_string);
	}
      theta *= M_PI/180.0;
      phi *= M_PI/180.0;
      off *= M_PI/180.0;
      double ra_rad = c.longitudeRad();
      double dec_rad = c.latitudeRad();
      wobble(theta,phi,ra_rad,dec_rad);
      wobble(theta,phi-off,ra_rad,dec_rad);
      c.setLatLongRad(dec_rad,ra_rad);
    }
  else if(mode.substr(0,3) == "off")
    {
      double offset = 30;
      if(mode.length() > 3)
	VSDataConverter::fromString(offset,
				    mode.substr(3));
      c.setPhi(c.phi() +
	       Angle::makeHrs(Astro::siderealRate()*double(offset)/60.0));
    }
  else
    {
      std::cerr << "Could not understand offset: " << mode
		<< std::endl;
      exit(EXIT_FAILURE);	      
    }

  return c;
}

VSTargetTable::Observation VSTargetTable::
getObservation(const SEphem::SphericalCoords& mean_radec,
	       double rms_ra_rad, double rms_dec_rad,
	       const SEphem::SphericalCoords& mean_azzn,
	       double rms_az_rad, double rms_zn_rad,
	       const std::string& demand_source_name,
	       const VSTime& approximate_start_time) const
{
  VSTargetTable::Observation obs;
  
  double on_sep = -1;

#ifdef DEBUG
  std::cout << Angle::toDeg(rms_az_rad) << ' '
	    << Angle::toDeg(rms_zn_rad) << ' '
	    << Angle::toDeg(rms_ra_rad) << ' '
	    << Angle::toDeg(rms_dec_rad) << std::endl;
#endif

  if(sqrt(rms_az_rad*rms_az_rad+rms_zn_rad*rms_zn_rad)<Angle::frDeg(0.1))
    {
      obs.mode        = Observation::OM_DRIFT;
      obs.mode_string = "drift";
      obs.drift_azzn  = mean_azzn;
      if(mean_azzn.theta() < Angle::frDeg(5.0))
	obs.name      = "zenith";
      else
	obs.name      = "drift";
      return obs;
    }
  else if(sqrt(rms_ra_rad*rms_ra_rad+rms_dec_rad*rms_dec_rad)
	  >= Angle::frDeg(0.1))
    {
      obs.clear();
      return obs;
    }

#if 0
  std::cout << demand_source_name << ' ' << approximate_start_time.toString()
	    << mean_radec.phi().hmsString() << ' ' 
	    << mean_radec.latitude().dmsString() << std::endl;
#endif

  if(!demand_source_name.empty())
    {
      for(unsigned icoord=0; icoord<m_coord.coord.size(); icoord++)
	if(m_coord.coord[icoord].name == demand_source_name)
	  {
	    double epoch        = m_coord.coord[icoord].epoch;
	    obs.src_radec       = m_coord.coord[icoord].radec;
	    obs.src_radec_J2000 = m_coord.coord[icoord].radec;

	    Astro::precess(approximate_start_time.getMJDDbl(), 
			   obs.src_radec, 
			   Astro::julianEpochToMJD(epoch));

	    Astro::precess(Astro::julianEpochToMJD(2000.0),
			   obs.src_radec_J2000, 
			   Astro::julianEpochToMJD(epoch));

	    on_sep          = mean_radec.separation(obs.src_radec).deg();
	    obs.name        = m_coord.coord[icoord].name;
	  }
    }
  else
    {
      double min_sep = -1;
      for(unsigned icoord=0; icoord<m_coord.coord.size(); icoord++)
	{
	  SphericalCoords radec = m_coord.coord[icoord].radec;
	  SphericalCoords radec_J2000 = m_coord.coord[icoord].radec;
	  Astro::precess(approximate_start_time.getMJDDbl(), 
			 radec, 
			 Astro::julianEpochToMJD(m_coord.coord[icoord].epoch));
	  Astro::precess(Astro::julianEpochToMJD(2000.0),
			 radec_J2000, 
			 Astro::julianEpochToMJD(m_coord.coord[icoord].epoch));
	  
	  Angle _on_sep     = mean_radec.separation(radec);
	  Angle sep         = _on_sep;

	  double offsepdeg = fabs(radec.thetaDeg()-mean_radec.thetaDeg());

	  if(offsepdeg < 0.05)
	    {
	      Angle radiff = radec.phi() - mean_radec.phi();
	      if((fabs(radiff.hrsPM())>=0.25)&&(fabs(radiff.hrsPM())<=0.75))
		sep = Angle::frDeg(offsepdeg);
	    }

#if 0
	  std::cerr << m_coord.coord[icoord].name << ' '
		    << _on_sep.deg() << ' '
		    << offsepdeg << ' ' << sep.deg() << '\n';
#endif

	  if((min_sep < 0)||(sep.deg() < min_sep))
	    {
	      min_sep             = sep.deg();
	      on_sep              = _on_sep.deg();
	      obs.src_radec       = radec;
	      obs.src_radec_J2000 = radec_J2000;
	      obs.name            = m_coord.coord[icoord].name;
	    }
	}
    }

  Angle radiff = obs.src_radec.phi() - mean_radec.phi();

#if 0
  std::cerr 
    << on_sep << ' ' << ' ' << obs.name << '\n'
    << obs.src_radec.phiDeg() << ' ' << mean_radec.phiDeg() << '\n'
    << obs.src_radec.thetaDeg() << ' ' << mean_radec.thetaDeg() << '\n'
    << fabs(obs.src_radec.thetaDeg()-mean_radec.thetaDeg()) << '\n'
    << fabs(radiff.hrsPM()) << ' ' << fabs(radiff.hrsPM()) << '\n';
#endif

  if(on_sep < 0.05)
    {
      // ON -------------------------------------------------------------------

      obs.mode              = Observation::OM_ON;
      obs.mode_string       = "on";
      obs.obs_radec         = obs.src_radec;
    }
  else if((fabs(obs.src_radec.thetaDeg()-mean_radec.thetaDeg())<0.05)
	  &&(fabs(radiff.hrsPM())>0.25) 
	  &&(fabs(radiff.hrsPM())<0.75))
    {
      // OFF ------------------------------------------------------------------

      // Calculate offset rounded to nearest minute
      obs.offset_min        = 
	round(-radiff.hrsPM()/Astro::siderealRate()*60);

      obs.mode              = Observation::OM_OFF;
      obs.mode_string       = "off";
      if(obs.offset_min>=0)
	obs.mode_string    += std::string("+");
      else 
	obs.mode_string    += std::string("-");
      obs.mode_string      +=
	VSDataConverter::toString(unsigned(fabs(obs.offset_min)),true);

      obs.obs_radec.setTheta(obs.src_radec.theta());
      obs.obs_radec.setPhi(obs.src_radec.phi() 
	  + Angle::makeHrs(Astro::siderealRate()*double(obs.offset_min)/60.0));
    }
  else if((on_sep > 0.10)&&(on_sep < 0.85)
	  &&(fabs(on_sep-round(on_sep/0.1)*0.1)<0.05))
    {
      // WOBBLE ---------------------------------------------------------------

      const double offset_round = Angle::frDeg(0.1);
      const double phi_round = Angle::frDeg(15.0);

      double source_ra = obs.src_radec.longitudeRad();
      double source_dec = obs.src_radec.latitudeRad();

      wobble_inverse(source_ra, source_dec,
		     mean_radec.longitudeRad(), mean_radec.latitudeRad(),
		     obs.wobble_theta_rad, obs.wobble_phi_rad);
      
      obs.wobble_theta_rad  = 
	round(obs.wobble_theta_rad/offset_round)*offset_round;
      obs.wobble_phi_rad    = 
	round(obs.wobble_phi_rad/phi_round)*phi_round;

      obs.mode              = Observation::OM_WOBBLE;
      obs.mode_string       = "wobble";
      obs.mode_string      += 
	VSDataConverter::toString(Angle::toDeg(obs.wobble_theta_rad),true);
      obs.mode_string      += std::string("@");
      obs.mode_string      += 
	VSDataConverter::toString(Angle::toDeg(obs.wobble_phi_rad),true);

      wobble(obs.wobble_theta_rad, obs.wobble_phi_rad, source_ra, source_dec);

      obs.obs_radec.setLatLongRad(source_dec,source_ra);
    }
  else
    {
      obs.clear();
      obs.mode        = Observation::OM_TRACKING;
      obs.mode_string = "tracking";
    }

  obs.obs_radec_J2000 = obs.obs_radec;
  Astro::precess(Astro::julianEpochToMJD(2000.0), obs.obs_radec_J2000,
		 approximate_start_time.getMJDDbl());

  obs.src_radec_J2000 = obs.src_radec;
  Astro::precess(Astro::julianEpochToMJD(2000.0), obs.src_radec_J2000,
		 approximate_start_time.getMJDDbl());     

  return obs;
}

std::vector<VSTargetTable::Target> VSTargetTable::
getTargetForObservation(const Observation& observation, double theta_cut) const
{
  static const double TWO_PI = 2.0*M_PI;  
  std::vector<VSTargetTable::Target> target;

  switch(observation.mode)
    {
    case Observation::OM_ON:
      {
	VSTargetTable::Target x;
	x.name     = observation.name + std::string("/on");
	x.coord    = observation.obs_radec;
	target.push_back(x);
      }
      break;

    case Observation::OM_SIM_ON:
      {
	VSTargetTable::Target x;
	x.name     = observation.name + std::string("/on");
	x.coord    = observation.src_radec;
	target.push_back(x);
      }
      break;

    case Observation::OM_OFF:
      {
	VSTargetTable::Target x;
	x.name     = observation.name+std::string("/")+observation.mode_string;
	x.coord    = observation.obs_radec;
	target.push_back(x);
      }
      break;

    case Observation::OM_WOBBLE:
    case Observation::OM_SIM_WOBBLE:
      {
	double region_sep = 
	  asin(Angle::frDeg(theta_cut)/observation.wobble_theta_rad)*2.0;
	double nregion = TWO_PI/region_sep;
	if(nregion>1)
	  {
	    unsigned offregion = unsigned(floor(nregion))-1;

	    VSTargetTable::Target xon;
	    xon.name       = observation.name + std::string("/on");
	    xon.coord      = observation.src_radec;
	    target.push_back(xon);

	    for(unsigned iregion=0;iregion<offregion;iregion++)
	      {
		double region_phi = 
		  (offregion%2==0)
		  ?region_sep*(double(iregion/2)+0.5)
		  :region_sep*(double((iregion+1)/2));
		if(iregion%2==1)region_phi = -region_phi;

		double off_ra = observation.obs_radec.longitudeRad();
		double off_dec = observation.obs_radec.latitudeRad();
		wobble(observation.wobble_theta_rad,
		       observation.wobble_phi_rad-region_phi,off_ra,off_dec);

		SphericalCoords woboff;
		woboff.setLatLongRad(off_dec,off_ra);

		VSTargetTable::Target xoff;
		xoff.name  = observation.name+std::string("/woboff");
		xoff.name += 
   VSDataConverter::toString(Angle::toDeg(observation.wobble_theta_rad),true);
		xoff.name += std::string("@");
		xoff.name += 
   VSDataConverter::toString(Angle::toDeg(observation.wobble_phi_rad),true);
		xoff.name += std::string("@");
		xoff.name += 
   VSDataConverter::toString(int(round(Angle::toDeg(region_phi))));

		xoff.coord = woboff;
		target.push_back(xoff);
	      }
	  }
      }
      break;

    case Observation::OM_UNKNOWN:
    default:
      break;
    };

  return target;
}

void VSTargetTable::mergeCompiledTargets()
{
  unsigned ntarget = sizeof(s_coord)/sizeof(*s_coord);
  for(unsigned icoord=0;icoord<ntarget;icoord++)
    {
      std::string ra_string   = s_coord[icoord].ra;
      std::string dec_string  = s_coord[icoord].dec;

      Angle ra;
      Angle dec;
      if((!ra.setFromHMSString(ra_string))
	 ||(!dec.setCoAngleFromDMSString(dec_string)))
	{
	  std::cerr << "Could not extract RA and DEC from: " 
		    << ra_string << " and " << dec_string
		    << std::endl;
	  exit(EXIT_FAILURE);	      
	}

      Coord coord;      
      coord.name        = s_coord[icoord].name;
      coord.radec.setLatLong(dec, ra);
      coord.epoch       = s_coord[icoord].epoch;

      m_coord.coord.push_back(coord);
    }
}

VSTargetTable::Data VSTargetTable::getCompiledTargets()
{
  VSTargetTable table;
  table.mergeCompiledTargets();
  return table.data();
}

void VSTargetTable::Data::load(VSOctaveH5ReaderStruct* reader)
{
  std::vector<double> ra_rad;
  std::vector<double> dec_rad;
  std::vector<double> epoch;
  reader->readVector("ra_rad",ra_rad);
  reader->readVector("dec_rad",dec_rad);

  if(reader->isVector("epoch"))
    reader->readVector("epoch",epoch);
  else
    reader->readVector("epoch_rad",epoch);

  VSOctaveH5ReaderCellVector* c = reader->readCellVector("name");
  unsigned ntarget = ra_rad.size();
  vsassert(ntarget == c->dimensions());
  coord.resize(ntarget);
  for(unsigned itarget=0;itarget<ntarget;itarget++)
    {
      coord[itarget].radec.setLatLongRad(dec_rad[itarget], ra_rad[itarget]);
      coord[itarget].epoch = epoch[itarget];
      c->readString(itarget,coord[itarget].name);
    }
  delete c;
}

void VSTargetTable::Data::save(VSOctaveH5WriterStruct* writer) const
{
  unsigned ntarget = coord.size();
  std::vector<double> ra_rad(ntarget);
  std::vector<double> dec_rad(ntarget);
  std::vector<double> epoch(ntarget);
  VSOctaveH5WriterCellVector* c = writer->writeCellVector("name",ntarget);
  for(unsigned icoord=0;icoord<ntarget;icoord++)
    {
      ra_rad[icoord] = coord[icoord].radec.phiRad();
      dec_rad[icoord] = coord[icoord].radec.latitudeRad();
      epoch[icoord] = coord[icoord].epoch;
      c->writeString(icoord,coord[icoord].name);
    }
  delete c;
  writer->writeVector("ra_rad",ra_rad);
  writer->writeVector("dec_rad",dec_rad);
  writer->writeVector("epoch",epoch);
}

void VSTargetTable::Observation::clear() 
{
  name              = "unknown";
  src_radec         = SphericalCoords();
  src_radec_J2000   = SphericalCoords();
  obs_radec         = SphericalCoords();
  obs_radec_J2000   = SphericalCoords();
  mode              = OM_UNKNOWN;
  mode_string       = "unknown";
  offset_min        = 0;
  wobble_theta_rad  = 0;
  wobble_phi_rad    = 0;
}

void VSTargetTable::Observation::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  reader->readString("name",name);
  int imode = mode;
  reader->readScalar("mode",imode);
  mode = ObservationMode(imode);
  reader->readString("mode_string",mode_string);
  if((mode != OM_DRIFT)&&(mode != OM_TRACKING))
    {
      double src_ra_rad;
      double src_dec_rad;
      reader->readScalar("src_ra_rad",src_ra_rad);
      reader->readScalar("src_dec_rad",src_dec_rad);
      src_radec.setLatLongRad(src_dec_rad, src_ra_rad);
      double src_ra_J2000_rad;
      double src_dec_J2000_rad;
      reader->readScalar("src_ra_J2000_rad",src_ra_J2000_rad);
      reader->readScalar("src_dec_J2000_rad",src_dec_J2000_rad);
      src_radec_J2000.setLatLongRad(src_dec_J2000_rad, src_ra_J2000_rad);
      double obs_ra_rad;
      double obs_dec_rad;
      reader->readScalar("obs_ra_rad",obs_ra_rad);
      reader->readScalar("obs_dec_rad",obs_dec_rad);
      obs_radec.setLatLongRad(obs_dec_rad, obs_ra_rad);
      double obs_ra_J2000_rad;
      double obs_dec_J2000_rad;
      reader->readScalar("obs_ra_J2000_rad",obs_ra_J2000_rad);
      reader->readScalar("obs_dec_J2000_rad",obs_dec_J2000_rad);
      obs_radec_J2000.setLatLongRad(obs_dec_J2000_rad, obs_ra_J2000_rad);
    }
  if(mode==OM_OFF)
    reader->readScalar("offset_min",offset_min);
  else if(mode==OM_WOBBLE)
    reader->readScalar("wobble_theta_rad",wobble_theta_rad),
      reader->readScalar("wobble_phi_rad",wobble_phi_rad);
  else if(mode==OM_DRIFT)
    {
      double drift_az_rad;
      double drift_zn_rad;
      reader->readScalar("drift_az_rad",drift_az_rad);
      reader->readScalar("drift_zn_rad",drift_zn_rad);
      drift_azzn.setThetaPhi(drift_zn_rad,drift_az_rad);
    }
}

void VSTargetTable::Observation::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeString("name",name);
  writer->writeScalar("mode",int(mode));
  writer->writeString("mode_string",mode_string);
  if((mode != OM_DRIFT)&&(mode != OM_TRACKING))
    {
      writer->writeScalar("src_ra_rad",src_radec.longitudeRad());
      writer->writeScalar("src_dec_rad",src_radec.latitudeRad());
      writer->writeScalar("src_ra_J2000_rad",src_radec_J2000.longitudeRad());
      writer->writeScalar("src_dec_J2000_rad",src_radec_J2000.latitudeRad());
      writer->writeString("src_ra_string",src_radec.longitude().hmsString(1));
      writer->writeString("src_dec_string",src_radec.latitude().dmsString(1));
      writer->writeString("src_ra_J2000_string",
			  src_radec_J2000.longitude().hmsString(1));
      writer->writeString("src_dec_J2000_string",
			  src_radec_J2000.latitude().dmsString(1));
      writer->writeScalar("obs_ra_rad",obs_radec.longitudeRad());
      writer->writeScalar("obs_dec_rad",obs_radec.latitudeRad());
      writer->writeScalar("obs_ra_J2000_rad",obs_radec_J2000.longitudeRad());
      writer->writeScalar("obs_dec_J2000_rad",obs_radec_J2000.latitudeRad());
      writer->writeString("obs_ra_string",obs_radec.longitude().hmsString(1));
      writer->writeString("obs_dec_string",obs_radec.latitude().dmsString(1));
      writer->writeString("obs_ra_J2000_string",
			  obs_radec_J2000.longitude().hmsString(1));
      writer->writeString("obs_dec_J2000_string",
			  obs_radec_J2000.latitude().dmsString(1));
    }
  if(mode==OM_OFF)
    writer->writeScalar("offset_min",offset_min);
  else if(mode==OM_WOBBLE)
    writer->writeScalar("wobble_theta_rad",wobble_theta_rad),
      writer->writeScalar("wobble_phi_rad",wobble_phi_rad);
  else if(mode==OM_DRIFT)
    writer->writeScalar("drift_az_rad",drift_azzn.phiRad()),
      writer->writeScalar("drift_zn_rad",drift_azzn.thetaRad());    
}

// ============================================================================
// VSTargetPointing
// ============================================================================

VSTargetPointing::
VSTargetPointing(const std::vector<TargetData>& targets,
		 const SEphem::SphericalCoords& earth_position):
  VSPointing(),
  m_scope(targets), m_earth_position(earth_position)
{
  // nothing to see here
}

VSTargetPointing::~VSTargetPointing()
{
  // nothing to see here
}

bool VSTargetPointing::
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

  unsigned nscope = m_scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!m_scope[iscope].empty())
      {
	bool scope_got_one = false;
	double scope_phi_zero = 0;
	double scope_theta_zero = 0;
	VSSimpleStat2<double> scope_x_stat;
	VSSimpleStat2<double> scope_y_stat;

	unsigned npointing = 100;
	for(unsigned ipointing=0;ipointing<npointing;ipointing++)
	  {
	    VSTime t(lo_time);
	    t += ((hi_time-lo_time)/npointing)*ipointing;
	  
	    double az;
	    double zn;

	    bool found = getAzZn(az, zn, iscope, t);
	    if(!found)continue;

	    SEphem::SphericalCoords azzn(zn,az);
	    if(!got_one)
	      {
		phi_zero = az;
		theta_zero = zn;
		got_one = true;
	      }

	    azzn.rotate(-phi_zero,-theta_zero,0);

	    x_stat.accumulate(azzn.thetaRad()*cos(azzn.phiRad()));
	    y_stat.accumulate(azzn.thetaRad()*sin(azzn.phiRad()));
	    zn_stat.accumulate(zn);

	    azzn.setThetaPhi(zn,az);

	    if(!scope_got_one)
	      {
		scope_phi_zero = az;
		scope_theta_zero = zn;
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

bool VSTargetPointing::
getMeanTarget(SEphem::SphericalCoords& mean_ra_dec,
	      double& rms_ra_rad, double& rms_dec_rad,
	      const VSTime& lo_time, const VSTime& hi_time, 
	      const SEphem::SphericalCoords& earth_position)
{
  bool              has_one_target = false; 

  std::string       target_type;
  double            target_angle1_rad = 0;
  double            target_angle2_rad = 0;
  double            target_epoch = 0;

  unsigned nscope = m_scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!m_scope[iscope].empty())
      {
	unsigned ntarget = m_scope[iscope].size();
	if(ntarget>1)return false;

	for(unsigned itarget=0; itarget<ntarget; itarget++)
	  {
	    const TargetDatum& t(m_scope[iscope][itarget]);
	    if(!has_one_target)
	      {
		has_one_target      = true;
		if(t.target_type != "tracking")return false;
		target_type         = t.target_type;
		target_angle1_rad   = t.target_angle1_rad;
		target_angle2_rad   = t.target_angle2_rad;
		target_epoch        = t.target_epoch;
	      }

	    if((target_type           != t.target_type)
	       ||(target_angle1_rad   != t.target_angle1_rad)
	       ||(target_angle2_rad   != t.target_angle2_rad)
	       ||(target_epoch        != t.target_epoch))
	      return false;
	  }
      }

  if(!has_one_target)return false;

  mean_ra_dec.setLatLongRad(target_angle2_rad,target_angle1_rad);
  rms_ra_rad = 0;
  rms_dec_rad = 0;

  return true;
}

bool VSTargetPointing::
getAzZn(double& az, double& zn, unsigned iscope, const VSTime& time)
{
  if((iscope>=m_scope.size())
     ||(m_scope[iscope].empty())
     ||(time<m_scope[iscope].front().db_start_time))return false;
  
  for(TargetData::const_iterator itarget=m_scope[iscope].begin();
      itarget != m_scope[iscope].end(); itarget++)
    if((time>=itarget->db_start_time)
       &&((itarget->db_end_time_is_null)||(time<=itarget->db_end_time)))
      return getTargetAzZn(az, zn, iscope, time, *itarget);
  
  return false;
}

bool VSTargetPointing::
getRaDec(double& ra, double& dec, unsigned iscope, const VSTime& time)
{
  if((iscope>=m_scope.size())
     ||(m_scope[iscope].empty())
     ||(time<m_scope[iscope].front().db_start_time))return false;
  
  for(TargetData::const_iterator itarget=m_scope[iscope].begin();
      itarget != m_scope[iscope].end(); itarget++)
    if((time>=itarget->db_start_time)
       &&((itarget->db_end_time_is_null)||(time<=itarget->db_end_time)))
      return getTargetRaDec(ra, dec, iscope, time, *itarget);

  return false;
}

bool VSTargetPointing::
getTargetAzZn(double& az, double& zn, unsigned iscope, const VSTime& time,
	      const TargetDatum& datum)
{
  if(datum.target_type == "tracking")
    {
      // Approximation to positioner target -- offset is not applied
      // to "real time" coodinates

      SphericalCoords coords;
      coords.setLatLongRad(datum.target_angle2_rad, datum.target_angle1_rad);

      double mjd = time.getMJDDbl();
      Angle lmst = Astro::mjdToLMST(mjd, m_earth_position.longitudeRad());

      Astro::precess(mjd, coords, Astro::julianEpochToMJD(datum.target_epoch));
      Astro::raDecMeanToApparent(mjd, coords);
      Astro::apparentRaDecToAzEl(mjd, lmst, m_earth_position, coords);

      az = coords.longitudeRad();
      zn = coords.thetaRad();
    }
  else if(datum.target_type == "fixed")
    {
      az = datum.target_angle1_rad;
      zn = M_PI_2 - datum.target_angle2_rad;
    }
  else
    {
      return false;
    }

  return true;
}

bool VSTargetPointing::
getTargetRaDec(double& ra, double& dec, unsigned iscope, const VSTime& time,
	       const TargetDatum& datum)
{
  if(datum.target_type == "tracking")
    {
      SphericalCoords coords;
      coords.setLatLongRad(datum.target_angle2_rad, datum.target_angle1_rad);

      double mjd = time.getMJDDbl();
      Astro::precess(mjd, coords, Astro::julianEpochToMJD(datum.target_epoch));
 
      ra = coords.longitudeRad();
      dec = coords.latitudeRad();
    }
  else if(datum.target_type == "fixed")
    {
      SphericalCoords coords;
      coords.setLatLongRad(datum.target_angle2_rad, datum.target_angle1_rad);

      double mjd = time.getMJDDbl();
      Angle lmst = Astro::mjdToLMST(mjd, m_earth_position.longitudeRad());

      Astro::azElToMeanRaDec(lmst, mjd, m_earth_position, coords);

      ra = coords.longitudeRad();
      dec = coords.latitudeRad();
      return true;
    }
  else
    {
      return false;
    }

  return true;
}

class TDCmp
{
public:
  bool operator() (const VSTargetPointing::TargetDatum& a,
		   const VSTargetPointing::TargetDatum& b)
  {
    if((a.target_type           < b.target_type)
       ||(a.target_name         < b.target_name)
       ||(a.tracking_mode       < b.tracking_mode)
       ||(a.tracking_param1     < b.tracking_param1)
       ||(a.tracking_param2     < b.tracking_param2)
       ||(a.pointing_mode       < b.pointing_mode)
       ||(a.pointing_array_mask < b.pointing_array_mask)
       ||(a.pointing_param      < b.pointing_param))return true;
    return false;
  }
};

bool VSTargetPointing::
getObservation(Observation& obs, 
	       const VSTime& run_start_time, const VSTime& run_end_time) const
{
  bool              has_one_target = false; 
  TargetDatum       the_one_target;

  unsigned nscope = m_scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!m_scope[iscope].empty())
      {
	std::map<TargetDatum, int64_t, TDCmp> targets;

	unsigned ntarget = m_scope[iscope].size();
	for(unsigned itarget=0; itarget<ntarget; itarget++)
	  if(m_scope[iscope][itarget].has_extended_target_info)
	    {
	      const TargetDatum& t(m_scope[iscope][itarget]);

	      VSTime start = t.db_start_time;
	      if(start < run_start_time)start = run_start_time;
	      VSTime end = t.db_end_time;
	      if(end > run_end_time)end = run_end_time;

	      if(start<end)
		targets[t] += (end-start);
	    }

	if(targets.empty())continue;

	std::map<TargetDatum, int64_t, TDCmp>::iterator itmax = 
	  targets.begin();

	for(std::map<TargetDatum, int64_t, TDCmp>::iterator 
	      itarget = targets.begin(); itarget != targets.end(); itarget++)
	  if(itarget->second > itmax->second)itmax = itarget;

	if(double(itmax->second)/double(run_end_time-run_start_time) < 0.9)
	  return false;

//std::cout << 'T' << iscope+1 << ": " << itmax->first.target_name << ' ' << double(itmax->second)/double(run_end_time-run_start_time) << '\n';
	if(!has_one_target)
	  {
	    the_one_target = itmax->first;
	    has_one_target = true;
	    if((the_one_target.target_type != "tracking")
	       &&(the_one_target.target_type != "fixed"))
	      return false;
	  }
	else
	  {
	    TDCmp cmp;
	    if(cmp(the_one_target, itmax->first)||
	       cmp(itmax->first, the_one_target))
//{ std::cout << 'T' << iscope+1 << ": UNEQUAL\n"; 
	      return false; //}
	  }
      }

  if(!has_one_target)return false;

  for(unsigned ichar=0;ichar<the_one_target.target_name.size();ichar++)
    the_one_target.target_name[ichar] = 
      tolower(the_one_target.target_name[ichar]);

  obs.clear();

  if(the_one_target.target_type == "fixed")
    {
      obs.mode           = Observation::OM_DRIFT;
      obs.mode_string    = "drift";
      obs.drift_azzn     .setLatLongRad(the_one_target.target_angle2_rad,
					the_one_target.target_angle1_rad);
      if(obs.drift_azzn.theta() < Angle::frDeg(5.0))
	obs.name          = "zenith";
      else
	obs.name          = "drift";
    }
  else
    {
      obs.name            = the_one_target.target_name;
      if(obs.name.empty())
	obs.name          = "unknown";
      obs.obs_radec       .setLatLongRad(the_one_target.target_angle2_rad,
					 the_one_target.target_angle1_rad);
      obs.src_radec       = obs.obs_radec;
      if(the_one_target.tracking_mode == "on")
	{
	  obs.mode        = Observation::OM_ON;
	  obs.mode_string = "on";
	}
      else if(the_one_target.tracking_mode == "off")
	{
	  obs.mode        = Observation::OM_OFF;
	  obs.mode_string = "off";
	  obs.offset_min  = 
	    the_one_target.tracking_param1/M_PI*12/Astro::siderealRate()*60;
	  if(obs.offset_min>=0)
	    obs.mode_string += std::string("+");
	  else 
	    obs.mode_string += std::string("-");
	  obs.mode_string   +=
	    VSDataConverter::toString(unsigned(fabs(round(obs.offset_min))),true);
	  obs.src_radec   .rotateRad(-the_one_target.tracking_param1, 0, 0);
	}
      else
	{
	  obs.mode        = Observation::OM_WOBBLE;
	  obs.mode_string = "wobble";
	  obs.mode_string += 
	    VSDataConverter::toString(round(Angle::toDeg(the_one_target.tracking_param1)*100)/100,true);
	  obs.mode_string += std::string("@");
	  obs.mode_string += 
	    VSDataConverter::toString(Angle::toDeg(the_one_target.tracking_param2),true);

	  obs.wobble_theta_rad = 
	    Angle::frDeg(round(Angle::toDeg(the_one_target.tracking_param1)*100)/100);
	  obs.wobble_phi_rad   = the_one_target.tracking_param2;
	  obs.src_radec   .setThetaPhi(the_one_target.tracking_param1,
				       -the_one_target.tracking_param2);
	  obs.src_radec   .rotateRad(obs.obs_radec.phi().rad(),
				     obs.obs_radec.theta().rad(), 0);
	}

      obs.obs_radec_J2000 = obs.obs_radec;
      Astro::precess(Astro::julianEpochToMJD(2000.0), obs.obs_radec_J2000,
		     run_start_time.getMJDDbl());

      obs.src_radec_J2000 = obs.src_radec;
      Astro::precess(Astro::julianEpochToMJD(2000.0), obs.src_radec_J2000,
		     run_start_time.getMJDDbl());     
    }

  return true;
}
