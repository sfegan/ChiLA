//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisData.cpp
  
  Misc analysis running information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/23/2006

  $Id: VSAnalysisData.cpp,v 3.3 2007/12/04 18:05:03 sfegan Exp $

*/

#include <unistd.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/types.h>
#include <pwd.h>

#include <VSDataConverter.hpp>
#include <VSAnalysisData.hpp>

using namespace VERITAS;

VSAnalysisData::VSAnalysisData():
  time_start(), time_end(), 
  revision_string(), revision_int(), version_string(),
  utime_msec(), stime_msec(), rtime_msec(), user_name(), host_name(),
  m_running(false), m_tms_start()
{
  // nothing to see here
}

VSAnalysisData::
VSAnalysisData(const char*const revision, const char*const version):
  time_start(VSTime::now()), time_end(), 
  revision_string(revision), revision_int(), version_string(version),
  utime_msec(), stime_msec(), rtime_msec(), user_name(), host_name(),
  m_running(true), m_tms_start()
{
  // Execution times at start -------------------------------------------------

  times(&m_tms_start);

  // Revision integer ---------------------------------------------------------

  revision_string = revision_string.substr(11,revision_string.size()-13);
  std::string::size_type rev_str_dot = revision_string.find('.');
  unsigned rev_int_minor = 0;
  unsigned rev_int_major = 0;
  if(rev_str_dot != std::string::npos)
    {
      VSDataConverter::fromString(rev_int_major,
				  revision_string.substr(0,rev_str_dot));
      VSDataConverter::fromString(rev_int_minor,
				  revision_string.substr(rev_str_dot+1));
    }
  revision_int = rev_int_major*1000+rev_int_minor;

  // Login name ---------------------------------------------------------------

  char* login = getlogin();
  if(login)user_name = login;
  else 
    {
      struct passwd* pwd = getpwuid(getuid());
      if(pwd)user_name = pwd->pw_name;
      else user_name = "unknown";
    }
  
  // Host name ----------------------------------------------------------------

  char hostname_buffer[81];
  gethostname(hostname_buffer,80);
  host_name = hostname_buffer;
}

void VSAnalysisData::stop()
{
  vsassert(m_running);

  time_end = VSTime::now();
  
  const long ticks_per_sec = sysconf(_SC_CLK_TCK);

  struct tms tms_end;
  times(&tms_end);

  utime_msec = (tms_end.tms_utime-m_tms_start.tms_utime)*1000/ticks_per_sec;
  stime_msec = (tms_end.tms_stime-m_tms_start.tms_stime)*1000/ticks_per_sec;
  rtime_msec = unsigned((time_end-time_start)/1000000ULL);

  m_running = false;
}

void VSAnalysisData::clear()
{
  time_start                 = VSTime();
  time_end                   = VSTime();
  revision_string.clear();
  revision_int               = 0;
  version_string.clear();
  utime_msec                 = 0;
  stime_msec                 = 0;
  rtime_msec                 = 0;
  user_name.clear();
  host_name.clear();
  m_running                  = false;
  memset(&m_tms_start,0,sizeof(m_tms_start));
}

void VSAnalysisData::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  reader->readCompositeHere(*this);
//   H5READ_VATIME(reader,"start_time",time_start);
//   H5READ_VATIME(reader,"end_time",time_end);
//   reader->readString("revision_string",revision_string);
//   reader->readScalar("revision",revision_int);
//   reader->readString("version",version_string);
//   reader->readScalar("run_utime_ms",utime_msec);
//   reader->readScalar("run_stime_ms",stime_msec);
//   reader->readScalar("run_rtime_ms",rtime_msec);
//   reader->readString("user",user_name);
//   reader->readString("host",host_name);
}

void VSAnalysisData::save(VSOctaveH5WriterStruct* writer) const
{
  vsassert(!m_running);
  writer->writeCompositeHere(*this);
//   writer->writeString("start_time", time_start.getString());
//   H5WRITE_VATIME(writer,"start_time",time_start);
//   H5WRITE_VATIME(writer,"end_time",time_end);
//   writer->writeString("revision_string",revision_string);
//   writer->writeScalar("revision",revision_int);
//   writer->writeString("version",version_string);
//   writer->writeScalar("run_utime_ms",utime_msec);
//   writer->writeScalar("run_stime_ms",stime_msec);
//   writer->writeScalar("run_rtime_ms",rtime_msec);
//   writer->writeString("user",user_name);
//   writer->writeString("host",host_name);
}

void VSAnalysisData::
printTiming(std::ostream& stream, const std::string& filename,
	    unsigned nevents) const
{
  stream << "Start time: " << time_start << std::endl
	 << "Stop time:  " << time_end << std::endl
	 << "Exec time:  " << std::fixed << std::setprecision(2)
	 << double(utime_msec)/1000 << " (user) " 
	 << double(stime_msec)/1000 << " (sys) " 
	 << double(rtime_msec)/1000 << " (real) "
	 << double(utime_msec+stime_msec)/double(rtime_msec)*100 << "% CPU"
	 << std::endl;
  if(!filename.empty())
    {
      static const double kB = 1024.0;
      static const double MB = kB*kB;
      static const double GB = kB*MB;
      struct stat filestat;
      if(stat(filename.c_str(),&filestat)==0)
	{
	  double size = double(filestat.st_size);
	  double rate_u = size/(double(utime_msec+stime_msec)/1000);
	  double rate_r = size/(double(rtime_msec)/1000);
	  stream << "Data rate:  ";
	  if(size>GB)
	    stream << std::setprecision(3) << size/GB << " GB";
	  else if(size>MB)
	    stream << std::setprecision(3) << size/MB << " MB";
	  else if(size>kB)
	    stream << std::setprecision(3) << size/kB << " kB";
	  else
	    stream << std::setprecision(3) << size/MB << " B ";
	  stream << "  - ";
	  if(rate_u>GB)
	    stream << std::setprecision(3) << rate_u/GB << " GB/s";
	  else if(rate_u>MB)
	    stream << std::setprecision(3) << rate_u/MB << " MB/s";
	  else if(rate_u>kB)
	    stream << std::setprecision(3) << rate_u/kB << " kB/s";
	  else
	    stream << std::setprecision(3) << rate_u    << " B/s ";
	  stream << " (user+sys) ";
	  if(rate_r>GB)
	    stream << std::setprecision(3) << rate_r/GB << " GB/s";
	  else if(rate_r>MB)
	    stream << std::setprecision(3) << rate_r/MB << " MB/s";
	  else if(rate_r>kB)
	    stream << std::setprecision(3) << rate_r/kB << " kB/s";
	  else
	    stream << std::setprecision(3) << rate_r    << " B/s ";
	  stream << " (real)" << std::endl;
	}
    }
  if(nevents)
    {
      static const double kEv = 1000.0;
      static const double MEv = kEv*kEv;
      static const double GEv = kEv*MEv;

      double events = double(nevents);
      double rate_u = events/(double(utime_msec+stime_msec)/1000);
      double rate_r = events/(double(rtime_msec)/1000);
      stream << "Event rate: ";
      if(events>GEv)
	stream << std::setprecision(3) << events/GEv << " GEv";
      else if(events>MEv)
	stream << std::setprecision(3) << events/MEv << " MEv";
      else if(events>kEv)
	stream << std::setprecision(3) << events/kEv << " kEv";
      else
	stream << std::setprecision(3) << events     << " Ev ";
      stream << " - ";
      if(rate_u>GEv)
	stream << std::setprecision(3) << rate_u/GEv << " GHz";
      else if(rate_u>MEv)
	stream << std::setprecision(3) << rate_u/MEv << " MHz";
      else if(rate_u>kEv)
	stream << std::setprecision(3) << rate_u/kEv << " kHz";
      else
	stream << std::setprecision(3) << rate_u     << " Hz ";
      stream << " (user+sys) ";
      if(rate_r>GEv)
	stream << std::setprecision(3) << rate_r/GEv << " GHz";
      else if(rate_r>MEv)
	stream << std::setprecision(3) << rate_r/MEv << " MHz";
      else if(rate_r>kEv)
	stream << std::setprecision(3) << rate_r/kEv << " kHz";
      else
	stream << std::setprecision(3) << rate_r     << " Hz ";
      stream << " (real)" << std::endl;
    }
}
