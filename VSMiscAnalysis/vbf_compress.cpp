#include <iostream>
#include <vector>

#include <VSOptions.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSSimpleVBF.hpp>
#include <VSPixelCleaningVisitor.hpp>
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
  std::string progname(*argv);
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  uint32_t dispatch_flags = 0;
  bool one_packet = false;
  unsigned packet_num = 0;
  unsigned npackets = 0;

  bool no_db = true;
  bool no_l3 = false;
  bool no_threads = false;
  bool do_use_overflow = false;
  std::string output_file("");
  std::string stage1_file("");

  options.findWithValue("stage1", stage1_file,
			"Set the name of the stage1 output file.  If this "
			"options is not defined no stage1 "
			"file will be generated.");

  options.findWithValue("npackets", npackets,
			"Set the number of packets to process.");

  options.findBoolValue("no_db", no_db, true,
			"Do not collect information from the database.");

  VSAnalysisStage1::configure(options,"","s1_");
  VSPixelCleaningVisitor::configure(options,"","");

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

//   if(no_l3)dispatcher.setDontUseL3Time(true);
//   if(no_channels)dispatcher.setVisitChannels(false);
  
  if(!stage1_file.empty())VSFileUtility::expandFilename(stage1_file);

  // --------------------------------------------------------------------------
  // Calculate stage 1 information
  // --------------------------------------------------------------------------
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
  

  VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(SIGINT);
  VSAnalysisStage1 stage1_calc;
  VSAnalysisStage1Data stage1_data;

  bool save_stage1 = false;

  if(!stage1_file.empty() && VSFileUtility::isFile(stage1_file))
    {
      VSOctaveH5Reader reader(stage1_file);
      VSOctaveH5ReaderStruct* s = reader.readStruct("stage1");
      vsassert(s);
      stage1_data.load(s);
      delete s;
    }
  else
    {
      stage1_calc.runStage1(vbf_file, stage1_data, earth_position,
			    no_db, no_l3, no_threads, do_use_overflow);    
      save_stage1 = true;
    }

  if(!stage1_file.empty() && save_stage1)
    {
      VSFileUtility::expandFilename(stage1_file);
      VSOctaveH5Writer writer(stage1_file, true);
      VSOctaveH5WriterStruct* stage1_struct = writer.writeStruct("stage1");
      stage1_data.save(stage1_struct);
      delete stage1_struct;
    }
  
  std::cout << "Finished stage1 analysis." << std::endl;

  VSPixelCleaningVisitor pixel_visitor(stage1_data);
  VSSimpleVBFDispatcher dispatcher(&pixel_visitor);
  
  try 
    {
      dispatcher.openFile(vbf_file.c_str());
      dispatcher.dispatchAllPackets(npackets, dispatch_flags);	 
      dispatcher.closeFile();
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
