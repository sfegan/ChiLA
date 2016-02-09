// vatime.cpp - Stephen Fegan - 2006-05-07
// $Id: vatime.cpp,v 1.1 2008/10/02 11:21:57 sfegan Exp $

// compile with: g++ -I. -DNOROOT -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -o vatime vatime.cpp VSOption.cpp VSTime.cpp -lm

// g++ -I. -DNOROOT -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -o vatime vatime.cpp VSTime.cpp -I ~/ChiLA/VSUtility -L ~/ChiLA/VSUtility/ -lVSUtility -lm

#include<string>
#include<iostream>
#include<iomanip>
#include<sstream>
#include<fstream>

#include"vsassert"
#include"VSAssert.hpp"
#include"VSDataConverter.hpp"
#include"VSOptions.hpp"
#include"VSTime.hpp"

using namespace std;
using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname << " [options] date_time" << std::endl
         << std::endl
	 << "Convert time between a number of different formats using VSTime."
	 << std::endl
	 << "The time can be specified as an MJD, UNIX or in a number of"
	 << std::endl
	 << "UTC formats. The time will be converted to a UTC string, UTC "
	 << std::endl
	 << "integer time-stamp (as used by the tracking program), UNIX time "
	 << std::endl
	 << "and broken down calendar time (with day of week and day of year)."
	 << std::endl
	 << "Additionally a reference time can be set and the offset between "
	 << std::endl
	 << "the input time and the reference is displayed in seconds, hours "
	 << std::endl
	 << "days and Julian years (overflow will occur at +/- 292 years)."
	 << std::endl << std::endl;

  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char *argv[])
{
  std::string progname(*argv);
  VSOptions options(argc, argv, true);
  bool print_usage = false;

  enum DateType { DT_NOW, DT_MJD, DT_UTC_STR, DT_UTC_INT, DT_UNIX,
		  DT_REF_SECOND, DT_REF_HOUR, DT_REF_DAY, DT_REF_YEAR };
  
  DateType type = DT_NOW;

  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  std::string ref_string("");
  options.findWithValue("ref", ref_string,
			"Set the reference time. Must be given as a valid UTC "
			"string (see \"utc\" option below). If an empty "
			"string is given the reference point is set to the "
			"current time.");
  
  if(options.find("now","Convert from time now [default].")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_NOW;
  
  if(options.find("mjd","Convert from MJD.")!=VSOptions::FS_NOT_FOUND)
    type = DT_MJD;

  if(options.find("utc",
		  "Convert from UTC string. The date should be given in the "
		  "form \"YYYY-MM-DD HH:MM:SS.NNNNNNNNN\". The time component "
		  "can be left out, as can the nanoseconds component.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_UTC_STR;

  if(options.find("int",
		  "Convert from UTC integer. The date should be given in the "
		  "form YYYYMMDDhhmmssmmm. The format is quite flexible, for "
		  "example the time portion can be eliminated completely, or "
		  "components can. A 2-digit year can be given. "
		  "For example 0605071121 is a valid date/time.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_UTC_INT;

  if(options.find("unix",
		  "Convert from UNIX time. Two parameters can be given, the "
		  "first (mandatory) is the UNIX time, the second (optional) "
		  "is the number of nanoseconds.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_UNIX;
  
  if(options.find("ref_sec",
		  "Convert from time in seconds with respect to reference. "
		  "The date should be given as the number of seconds after "
		  "(positive) or before (negative) the reference time.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_REF_SECOND;

  if(options.find("ref_hour",
		  "Convert from time in hours with respect to reference. "
		  "The date should be given as the number of hours after "
		  "(positive) or before (negative) the reference time.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_REF_HOUR;

  if(options.find("ref_day",
		  "Convert from time in days with respect to reference. "
		  "The date should be given as the number of days after "
		  "(positive) or before (negative) the reference time.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_REF_DAY;

  if(options.find("ref_year",
		  "Convert from time in years with respect to reference. "
		  "The date should be given as the number of years after "
		  "(positive) or before (negative) the reference time.")
     !=VSOptions::FS_NOT_FOUND)
    type = DT_REF_YEAR;

  bool print_short = false;
  if(options.find("s","Print the converted date in short form."))
    print_short = true;

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  if((argc==0)&&(type != DT_NOW))
    {
      std::cerr << progname << ": need at least one argumenent" 
		<< std::endl << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }
  
  VSTime now = VSTime::now();

  VSTime ref = now;
  if(!ref_string.empty())ref.setFromString(ref_string);

  VSTime t;
  bool good = false;
  if(type == DT_NOW)
    {
      t = now;
      good=true;
    }
  else if(type == DT_MJD)
    {
      double mjd;
      good = VSDataConverter::fromString(mjd,*argv);
      if(good)good = t.setFromMJDDbl(mjd);
    }
  else if(type == DT_UTC_STR)
    {
      good = t.setFromString(*argv);
    }
  else if(type == DT_UTC_INT)
    {
      uint64_t utc;
      good = VSDataConverter::fromString(utc,*argv);
      if(good)
	{
	  if(utc < UINT64_C(500000))utc += UINT64_C(20000000);
	  else if(utc < UINT64_C(999999))utc += UINT64_C(19000000);
//	  else if(utc < UINT64_C(50000000))utc += UINT64_C(2000000000);
//	  else if(utc < UINT64_C(99999999))utc += UINT64_C(1900000000);
//	  else if(utc < UINT64_C(99999999))utc += 0;
//	  else if(utc < UINT64_C(5000000000))utc += UINT64_C(200000000000);
//	  else if(utc < UINT64_C(9999999999))utc += UINT64_C(190000000000);
//	  else if(utc < UINT64_C(500000000000))utc += UINT64_C(20000000000000);
//	  else if(utc < UINT64_C(999999999999))utc += UINT64_C(19000000000000);
//	  else if(utc < UINT64_C(500000000000000))utc += UINT64_C(20000000000000000);
//	  else if(utc < UINT64_C(999999999999999))utc += UINT64_C(19000000000000000);
	  
	  if(utc < UINT64_C(100000000))utc = utc*UINT64_C(1000000000);
	  else if(utc < UINT64_C(10000000000))utc = utc*UINT64_C(10000000);
	  else if(utc < UINT64_C(1000000000000))utc = utc*UINT64_C(100000);
	  else if(utc < UINT64_C(100000000000000))utc = utc*UINT64_C(1000);
	  good = t.setFromMSTimeStamp(utc);
	}
    }
  else if(type == DT_UNIX)
    {
      time_t s = 0;
      uint32_t ns = 0;
      good = VSDataConverter::fromString(s,*argv);
      if((good)&&(argc>1))
	{
	  good = VSDataConverter::fromString(ns,*(argv+1));
	  std::istringstream stream1(*argv);
	}
      if(good)good = t.setFromPOSIXTimeT(&s,ns);
    }
  else if(type == DT_REF_SECOND)
    {
      t = ref;
      int32_t s = 0;
      good = VSDataConverter::fromString(s,*argv);
      if(good)
	{
	  int64_t ns = s;
	  ns *= VATime::SECOND_NS_I;
	  t += ns;
	}
    }
  else if(type == DT_REF_HOUR)
    {
      t = ref;
      double h = 0;
      good = VSDataConverter::fromString(h,*argv);
      if(good)
	{
	  int64_t ns = llround(h * double(VATime::HOUR_NS_I));
	  t += ns;
	}
    }
  else if(type == DT_REF_DAY)
    {
      t = ref;
      double d = 0;
      good = VSDataConverter::fromString(d,*argv);
      if(good)
	{
	  int64_t ns = llround(d * double(VATime::JULIAN_NS_I));
	  t += ns;
	}
    }
  else if(type == DT_REF_YEAR)
    {
      t = ref;
      double y = 0;
      good = VSDataConverter::fromString(y,*argv);
      if(good)
	{
	  int64_t ns = llround(y * 365.25 * double(VATime::JULIAN_NS_I));
	  t += ns;
	}
    }


  if(!good)
    {
      std::cerr << progname << ": Could not convert date/time: " << *argv
		<< std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);      
    }

  if(!print_short)std::cout << "UTC:          ";
  std::cout << t.getString();
  if(!print_short)std::cout << std::endl;
  else std::cout << ' ';

  if(!print_short)std::cout << "Timestamp:    ";
  std::cout << t.getMSTimeStamp();
  if(!print_short)std::cout << std::endl;
  else std::cout << ' ';

  if(!print_short)std::cout << "MJD:          ";
  std::cout << std::fixed << std::setprecision(9) << t.getMJDDbl();
  if(!print_short)std::cout << std::endl;
  else std::cout << ' ';

  if(!print_short)std::cout << "UNIX:         ";
  std::cout << t.getPOSIXTimeT();
  if(!print_short)std::cout << std::endl;
  else std::cout << ' ';

  if(!print_short)std::cout << "Components:   ";
  std::cout << t.getYear() << ' ' << t.getMonth() << ' ' << t.getDay() << ' '
	    << t.getHour() << ' ' << t.getMin() << ' ' << t.getSec() << ' '
	    << t.getNanoSec() << ' ' << t.getDayOfYear() << ' ' 
	    << t.getDayOfWeek();
  if(!print_short)std::cout << std::endl;
  else std::cout << ' ';

  if(!print_short)std::cout << "Ref UTC:      ";
  std::cout << ref.getString();
  if(!print_short)std::cout << std::endl;
  else std::cout << ' ';

  if(!print_short)std::cout << "From ref:     ";
  std::cout << (t-ref)/INT64_C(1000000000) << ' '
	    << setprecision(2) << double((t-ref)/INT64_C(1000000000))/3600.0 << ' '
	    << setprecision(2) << double((t-ref)/INT64_C(1000000000))/86400.0 << ' '
	    << setprecision(2) << double((t-ref)/INT64_C(1000000000))/31557600.00;
  if(!print_short)std::cout << std::endl;
  else ' ';

  if(!print_short)std::cout << "Julian epoch: ";
  std::cout << setprecision(3) << t.getJulianEpoch();
  if(!print_short)std::cout << std::endl;
  else std::cout << std::endl;
}
