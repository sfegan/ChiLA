/*! \file stage2.cpp

  Driver for stage 2 of analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/25/2006

  $Id: stage2.cpp,v 3.18 2009/10/14 22:03:30 matthew Exp $

*/

#include <stdlib.h>
#include <signal.h>

#include <VSOptions.hpp>
#include <VSTime.hpp>
#include <VSDBFactory.hpp>
#include <VSSimpleVBF.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSAnalysisStage2.hpp>
#include <VSFileUtility.hpp>
#include <VBFSimpleRunInfoExtractor.hpp>
#include <VSCentralizedDBAccess.hpp>

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// MAIN
// ============================================================================

void sigabort(int signal)
{
  const char* text;
  text = "Caught signal: ";
  write(2,text,strlen(text));
  text = strsignal(signal);
  write(2,text,strlen(text));
  text = " (aborting)\n";
  write(2,text,strlen(text));
  abort();
}

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
  //signal(SIGSEGV,sigabort);

  try
    {
      std::string progname(*argv);
      std::vector<std::string> command_line(argc);
      for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
	command_line[iarg]=argv[iarg];

      std::cout << "Command line:";
      for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
	std::cout << ' ' << argv[iarg];
      std::cout << std::endl;
 
      VSOptions options(argc, argv, true);

      bool print_usage = false;
      if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
	print_usage=true;
      if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
	print_usage=true;

      std::string profile = "analysis";
      options.findWithValue("profile", profile,
			    "Set the default values for a given profile. The "
			    "profiles presently defined are: \"analysis\", "
			    "\"diagnostics\" and \"simulations\"", "common");

      bool no_db = false;
      bool no_l3 = false;
      bool no_threads = false;
      bool no_verbose = false;
      bool do_use_overflow = false;

      unsigned runno = 0;

      bool compress = true;
      std::string laser_file("");
      std::string hilo_file("");
      std::string output_file("");
      std::string stage1_file("");
      std::string pad_stage1_file("");

#ifdef NOVDB
      options.addCatagory("db","Configure connection to database");
      VSDBFactory::configure(&options,"db");
#endif
      VSAnalysisStage1::configure(options,profile);
      VSAnalysisStage2::configure(options,profile);

      triple<std::string, std::string, double> 
	array_pos("31d40m29.04s","-110d57m10.08s",1270);

      options.addCatagory("common", 
			  "Options common to multiple stages.");

      options.addCatagory("file", 
			  "Input and output files");

      options.findWithValue("array_pos", array_pos,
			    "Set the reference position for the array. The "
			    "coordinate triple should be specified as "
			    "latitude/longitude/elevation.",
			    "common");

      Angle latitude;
      Angle longitude;
      latitude.setFromDMSString(array_pos.first);
      longitude.setFromDMSString(array_pos.second);
      SphericalCoords earth_position;
      earth_position.setLatLong(latitude,longitude);
      double earth_elevation = array_pos.third;

      options.findBoolValue("no_db", no_db, true,
			    "Do not collect information from the database.",
			    "common");
      options.findBoolValue("no_l3", no_l3, true, 
			    "Do not use L3 time.",
			    "common");

      VSSimCoordTransform::configure(options,profile);

#ifndef NOTHREADS
      options.findBoolValue("no_threads", no_threads, true,
			"Disable multi-threaded operation. This may slow the "
			"analysis on machines with multiple processors, "
			"multiple cores, hyperthreading and when reading data "
			"from a slow source (such as NFS).");
#endif

      options.findBoolValue("no_verbose", no_verbose, true, 
			"Do not write anything to the terminal.");

      options.findBoolValue("do_merge_overflow", do_use_overflow, true,
			"Do attempt to merge events from the overflow bank.",
			"common");

      options.findBoolValue("compress", compress, true,
			"Compress data in output files where possible.",
			"file"); 
  
      options.findWithValue("o", output_file, 
			"Set the output file name. If the output file is "
			"empty a name defined by the run number will be used, "
			"for example: \"x12345_s2.h5\". "
			"The first instance of a question mark in the name "
			"will be replaced with the run number (there is no "
			"way to override this), for example: \"x?_s2.h5\" "
			"might be expanded to \"x12345_s2.h5\".",
			"file"); 

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
		      "use x30000.h5 for T1-T3 and x30001.h5 for T4). "
		      "The LASER file(s) specified will "
		      "override any LASER data in the Stage 1 input file. "
		      "If this option is not given then the LASER data "
		      "embedded in the Stage 1 file will be used, or no "
		      "LASER calibration will be done if that data does "
		      "not exist. Whether the LASER file is actually used "
		      "in the analysis is determined by the stage 2 option "
		      "\"s2_no_laser\".",
		      "file");

      options.findWithValue("hilo", hilo_file, 
			    "Set the name of the Hi-Lo gain calibration file "
			    "to be used in the analysis.",
			    "file");

      options.findWithValue("stage1", stage1_file, 
			"Set the name of the stage 1 file corresponding to "
			"the VBF file. If no file is speicified then the "
			"VBF file is used to determine the run number and a "
			"filename based on the run number is used. If the "
			"file is not available the stage 1 analysis is run "
			"directly here, and its output written before the "
			"stage 2 analysis starts.",
			"file");

      options.findWithValue("pad_stage1", pad_stage1_file, 
			"Set the name of the stage 1 file to use for "
			"padding.",
			"file");

      options.findWithValue("runno", runno, 
			"Set the run number of the file. If not explictly "
			"set, the number will be determined from the VBF file "
			"or from its name. If no number can be determined "
			"the analysis will not proceed.",
			"common");

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
      VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(SIGINT);

      if(!hilo_file.empty())VSFileUtility::expandFilename(hilo_file);
      if(!output_file.empty())VSFileUtility::expandFilename(output_file);
      if(!stage1_file.empty())VSFileUtility::expandFilename(stage1_file);
      if(!pad_stage1_file.empty())
	VSFileUtility::expandFilename(pad_stage1_file);

      if(!no_verbose)
	std::cout 
	  << "Input file:   " << vbf_file << std::endl;

      VSOctaveH5Writer::setDefaultCompress(compress);

      // ----------------------------------------------------------------------
      // Pass zero through the (start of the) file.. extract the run number etc
      // ----------------------------------------------------------------------

      VBFSimpleRunInfoExtractor simple_run_info;
      VSSimpleVBFDispatcher dispatcher(&simple_run_info);
      dispatcher.processFile(vbf_file.c_str());

      if(!no_verbose)
	{
	  if(simple_run_info.hasRunNumber())
	    std::cout 
	      << "VBF run num:  " << simple_run_info.runNumber() << std::endl;
	  else
	    std::cout 
	      << "VBF run num:  unavailable from VBF file" << std::endl;
	  
	  if(simple_run_info.hasApproximateStartTime())
	    std::cout 
	      << "Approx time:  " 
	      << simple_run_info.approximateStartTime().toString()
	      << std::endl;
	  else
	    std::cout 
	      << "Approx time:  unavailable from VBF file" << std::endl;

	  if(simple_run_info.hasSubarrayTelescopeList())
	    {
	      std::cout << "Telescopes:  ";
	      if(simple_run_info.subarrayTelescopeList().empty())
		std::cout << " none";
	      else
		for(unsigned iscope=0;
		    iscope<simple_run_info.subarrayTelescopeList().size();
		    iscope++)
		  std::cout << ' ' 
			    << simple_run_info.subarrayTelescopeList()[iscope];
	      std::cout << std::endl;
	    }
	  else
	    std::cout 
	      << "Telescopes:   unavailable from VBF file" << std::endl;
	}

      if(runno==0)
	{
	  if(simple_run_info.hasRunNumber())
	    runno = simple_run_info.runNumber();
	  if(runno == 0)runno = 
	    VSFileUtility::extractNumberFromFilename(vbf_file);

	  if(runno==0)
	    {
	      std::cerr << "FATAL: Run number not available from VBF file "
		"and not set explicitly."
			<< std::endl;
	      exit(EXIT_FAILURE);
	    }
	}
      
      std::cout << "Run number:   " << runno << std::endl;


      
      // ----------------------------------------------------------------------
      // Load or calculate stage 1 files
      // ----------------------------------------------------------------------

      VSAnalysisStage1Data stage1_data;

      if(stage1_file.empty())stage1_file="x?_s1.h5";

      if(stage1_file.find('?') != stage1_file.npos)
	{
	  std::string::size_type iquestion = stage1_file.find('?');
	  stage1_file.replace(iquestion,1, VSDataConverter::toString(runno));
	}
  
      if(!no_verbose)
	std::cout 
	  << "Stage1 file:  " << stage1_file << std::endl;

      bool save_stage1 = false;
      if(!VSFileUtility::isFile(stage1_file))
	{
	  if((!no_db)
	     &&(simple_run_info.hasApproximateStartTime())
	     &&(simple_run_info.hasSubarrayTelescopeList()))
	    VSCentralizedDBAccess::getInstance()->
	      prefetchFromDB(runno,
			     simple_run_info.approximateStartTime(),
			     simple_run_info.subarrayTelescopeList(),
			     no_threads);

	  VSAnalysisStage1 stage1;
	  stage1.runStage1(vbf_file, stage1_data, earth_position,
			   no_db, no_l3, no_threads, no_verbose, 
			   do_use_overflow);
	  save_stage1 = true;
	}
      else
	{
	  VSOctaveH5Reader reader(stage1_file);
	  VSOctaveH5ReaderStruct* s = reader.readStruct("stage1");
	  stage1_data.load(s);
	  delete s;
	}

      // ----------------------------------------------------------------------
      // Laser
      // ----------------------------------------------------------------------

      if(!laser_file.empty())
	{
	  if(!no_verbose)
	    std::cout << "Laser file:   " << laser_file << std::endl;

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
	  if(!no_verbose)
	    std::cout 
	      << "Hi/Lo file:   " << hilo_file << std::endl;

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
      // Write stage 1 file if necessary
      // ----------------------------------------------------------------------

      if(save_stage1 == true)
	{
	  time_t thetime = time(NULL);
	  std::string thetimestring(ctime(&thetime));
	  std::string comment =
	    std::string("# VERITAS analysis stage 1, written by ") + progname 
	    + std::string(" - ") 
	    + thetimestring.substr(0,thetimestring.length()-1);
      
	  VSOctaveH5Writer writer(stage1_file, true, comment);
	  VSOctaveH5WriterStruct* stage1_struct = writer.writeStruct("stage1");
	  stage1_data.save(stage1_struct);
	  delete stage1_struct;
	}

      VSAnalysisStage1Data* pad_stage1_data = 0;
      if(!pad_stage1_file.empty())
	{
	  pad_stage1_data = new VSAnalysisStage1Data;

	  VSOctaveH5Reader reader(pad_stage1_file);
	  VSOctaveH5ReaderStruct* s = reader.readStruct("stage1");
	  pad_stage1_data->load(s);
	  delete s;
	}

      // ----------------------------------------------------------------------
      // Stage 2
      // ----------------------------------------------------------------------

      if(output_file.empty())output_file="x?_s2.h5";

      if(output_file.find('?') != output_file.npos)
	{
	  std::string::size_type iquestion = output_file.find('?');
	  output_file.replace(iquestion,1,VSDataConverter::toString(runno));
	}
  
      if(!no_verbose)
	std::cout 
	  << "Output file:  " << output_file << std::endl;

      time_t thetime = time(NULL);
      std::string thetimestring(ctime(&thetime));
      std::string comment =
	std::string("# VERITAS analysis stage 2, written by ") + progname +
	std::string(" - ") + thetimestring.substr(0,thetimestring.length()-1);

      VSOctaveH5Writer writer(output_file, true, comment);
      VSAnalysisStage2 stage2;

      stage2.runStage2(vbf_file, &stage1_data, &writer, 
		       earth_position, earth_elevation, pad_stage1_data,
		       no_db, no_l3, no_threads, no_verbose, do_use_overflow);

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
  catch(const std::string& x)
    {
      std::cout << "Caught error message: " << x << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(...)
    {
      std::cout << "Caught some exception" << std::endl;
      exit(EXIT_FAILURE);
    }

  H5close();
};
