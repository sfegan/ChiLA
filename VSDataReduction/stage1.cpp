/*! \file stage1.cpp

  Driver for stage 1 of analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/13/2006

  $Id: stage1.cpp,v 3.10 2009/10/14 22:03:30 matthew Exp $

*/

#include <VSOptions.hpp>
#include <VSTime.hpp>
#include <VSDBFactory.hpp>
#include <VSSimpleVBF.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSFileUtility.hpp>

using namespace VERITAS;
using namespace SEphem;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] vbf_file" 
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  try
    {
      std::string progname(*argv);
      std::vector<std::string> command_line(argc);
      for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
	command_line[iarg]=argv[iarg];
 
      VSOptions options(argc, argv, true);

      bool print_usage = false;
      if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
	print_usage=true;
      if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
	print_usage=true;

      std::string profile = "veritas";
      options.findWithValue("profile", profile,
			    "Set the default values for a given profile. The "
			    "profiles presently defined are: \"veritas\", "
			    "\"simulations\" and \"whipple\"");

      bool no_db = false;
      bool no_l3 = false;
      bool no_threads = false;
      bool do_use_overflow = false;
      std::string laser_file("");
      std::string hilo_file("");
      std::string output_file("");

#ifdef NOVDB
      VSDBFactory::configure(&options);
#endif

      VSAnalysisStage1::configure(options,profile,"");

      triple<std::string, std::string, double> 
	array_pos("31d40m29.04s","-110d57m10.08s",1270);

      options.findWithValue("array_pos", array_pos,
			    "Set the reference position for the array. The "
			    "coordinate triple should be specified as "
			    "latitude/longitude/elevation.");

      Angle latitude;
      Angle longitude;
      latitude.setFromDMSString(array_pos.first);
      longitude.setFromDMSString(array_pos.second);
      SphericalCoords earth_position;
      earth_position.setLatLong(latitude,longitude);
      //  double earth_elevation = array_pos.third;

      options.findBoolValue("no_db", no_db, true,
			    "Do not collect information from the database.");
      options.findBoolValue("no_l3", no_l3, true, 
			    "Do not use L3 time.");

      VSSimCoordTransform::configure(options,profile);

#ifndef NOTHREADS
      options.findBoolValue("no_threads", no_threads, true,
			"Disable multi-threaded operation. This may slow the "
			"analysis on machines with multiple processors, "
			"multiple cores, hyperthreading and when reading data "
			"from a slow source (such as NFS).");
#endif

      options.findBoolValue("do_merge_overflow", do_use_overflow, true,
			"Do attempt to merge events from the overflow bank.");

      options.findWithValue("o", output_file, 
			"Set the output file name. If the output file is "
			"empty a name defined by the run number will be used, "
			"for example: x12345_s1.h5.");

      options.
	findWithValue("laser", laser_file, 
		      "Set a comma-separated list of LASER calibration "
		      "files to be used in the analysis. "
		      "Calibration files can be preceded by a telescope "
		      "index (e.g. 3/x30001.h5) specifying the telescope "
		      "to which that laser calibration should be applied. "
		      "A calibration file without an index will be used for "
		      "all telescopes for which no specific calibration "
		      "was chosen (e.g. x30000.h5,3/x30001.h5 will "
		      "use x30000.h5 for T1-T3 and x30001.h5 for T4).");

      options.findWithValue("hilo", hilo_file, 
			"Set the name of the Hi-Lo gain calibration file to "
			"be used in the analysis.");

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

      int arg_req = 1;
      if(argc != arg_req)
	{
	  std::cerr << progname << ": need " << arg_req
		    << " arguments, got " << argc << std::endl;
	  usage(progname, options, std::cout);
	  exit(EXIT_SUCCESS);
	}

      std::string vbf_file = *argv;
      argv++,argc--;

      if(!hilo_file.empty())VSFileUtility::expandFilename(hilo_file);
      if(!output_file.empty())VSFileUtility::expandFilename(output_file);

      // ----------------------------------------------------------------------
      // Calculate stage 1 information
      // ----------------------------------------------------------------------

      VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(SIGINT);
      VSAnalysisStage1 stage1_calc;
      VSAnalysisStage1Data stage1_data;
      stage1_calc.runStage1(vbf_file, stage1_data, earth_position,
			    no_db, no_l3, no_threads, do_use_overflow);

      // ----------------------------------------------------------------------
      // Laser
      // ----------------------------------------------------------------------

      if(!laser_file.empty())
	{
	  if(stage1_data.laser)
	    {
	      delete stage1_data.laser;
	      stage1_data.laser = 0;
	    }

	  VSLaserData* laser = new VSLaserData;
	  if(laser->load(laser_file)) stage1_data.laser = laser;
	  else
 	    {
 	      delete laser;
 	      std::cerr << "Could not open LASER file: " << laser_file 	
			<< std::endl;
	      exit(EXIT_FAILURE);
 	    }
	}

      // ----------------------------------------------------------------------
      // Hi/Lo gain
      // ----------------------------------------------------------------------

      if(!hilo_file.empty())
	{
          if(!VSFileUtility::isFile(hilo_file))
	    {
	      std::cerr << "Hi-Lo gain file does not exist: " 
			<< hilo_file << std::endl;
	      exit(EXIT_FAILURE);
	    }

	  if(stage1_data.hilo)
	    {
	      delete stage1_data.hilo;
	      stage1_data.hilo = 0;
	    }

	  if(VSOctaveH5ReaderStruct::isHDF5(hilo_file))
	    {
	      VSOctaveH5Reader reader(hilo_file);
	      VSOctaveH5ReaderStruct* s = reader.readStruct("hilo");
	      VSHiLoData* hilo = new VSHiLoData;
	      hilo->load(s);
	      stage1_data.hilo = hilo;
	      delete s;
	    }
	  else 
	    {
	      std::cerr << "Could not open Hi-Lo file: " << hilo_file
			<< std::endl;
	      exit(EXIT_FAILURE);
	    }
	}

      // ----------------------------------------------------------------------
      // Write stage 1 file
      // ----------------------------------------------------------------------

      std::ostringstream s1_filename;
      s1_filename << 'x' << stage1_data.run_info.run_number << "_s1.h5";
      if(output_file.empty())output_file=s1_filename.str();
  
      time_t thetime = time(NULL);
      std::string thetimestring(ctime(&thetime));
      std::string comment =
	std::string("# VERITAS analysis stage 1, written by ") + progname +
	std::string(" - ") + thetimestring.substr(0,thetimestring.length()-1);

      VSOctaveH5Writer writer(output_file, true, comment);

      VSOctaveH5WriterStruct* stage1_struct = writer.writeStruct("stage1");
      stage1_data.save(stage1_struct);
      delete stage1_struct;

      // ----------------------------------------------------------------------
      // Write configuration etc..
      // ----------------------------------------------------------------------

      VSOctaveH5WriterCellVector* wc = 
	writer.writeCellVector("command_line", command_line.size());
      for(unsigned iarg=0;iarg<command_line.size();iarg++)
	wc->writeString(iarg,command_line[iarg]);
      delete wc;

      unsigned noptionrecords = options.getOptionRecords().size();
      wc = writer.writeCellVector("configuration", noptionrecords);
      for(unsigned irecord=0;irecord<noptionrecords;irecord++)
	{
	  VSOctaveH5WriterStruct* ws = wc->writeStruct(irecord);
	  ws->writeString("key", options.getOptionRecords()[irecord].key);
	  ws->writeScalar("status", 
			  int(options.getOptionRecords()[irecord].status));
	  ws->writeScalar("val_requested", 
			  options.getOptionRecords()[irecord].val_requested);
	  ws->writeString("val", options.getOptionRecords()[irecord].val);      
	  delete ws;
	}
      delete wc;  
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cout << "Caught instance of VSOctaveH5Exception" << std::endl
		<< e.message() << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(const VSAssert& a)
    {
      std::cout << a.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(const std::exception& x)
    {
      std::cout << "Caught instance of " << x.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(...)
    {
      std::cout << "Caught some exception" << std::endl;
      exit(EXIT_FAILURE);
    }

  H5close();
}
