//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisData.hpp
  
  Misc analysis running information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/23/2006

  $Id: VSAnalysisData.hpp,v 3.2 2007/04/25 15:06:32 sfegan Exp $

*/

#ifndef VSANALYSISDATA_HPP
#define VSANALYSISDATA_HPP

#include <iostream>

#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include <pwd.h>

#include <VSTime.hpp>
#include <VSOctaveIO.hpp>

namespace VERITAS
{

  class VSAnalysisData
  {
  public:
    VSAnalysisData();
    VSAnalysisData(const char*const revision, const char*const version);

    void stop();
    void printTiming(std::ostream& stream, 
		     const std::string& filename="", unsigned nevents=0) const;

    void clear();
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    VSTime        time_start;
    VSTime        time_end;
    std::string   revision_string;
    unsigned      revision_int;
    std::string   version_string;
    unsigned      utime_msec;
    unsigned      stime_msec;
    unsigned      rtime_msec;
    std::string   user_name;
    std::string   host_name;

  static void _compose(VSOctaveH5CompositeDefinition& c)
  {
    // Do this by hand since we don't want to use the same variable
    // names in the HDF file
    H5_ADDNAMEDSIMPLECOMPOSITE(c,VSAnalysisData,time_start,"start_time");
    H5_ADDNAMEDSIMPLECOMPOSITE(c,VSAnalysisData,time_end,"end_time");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,revision_string,"revision_string");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,revision_int,"revision");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,version_string,"version");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,utime_msec,"run_utime_ms");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,stime_msec,"run_stime_ms");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,rtime_msec,"run_rtime_ms");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,user_name,"user");
    H5_ADDNAMEDMEMBER(c,VSAnalysisData,host_name,"host");
  }

  private:
    bool          m_running;
    struct tms    m_tms_start;
  };

}

#endif // VSANALYSISDATA_HPP
