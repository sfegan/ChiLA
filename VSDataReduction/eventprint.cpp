//-*-mode:c++; mode:font-lock;-*-

/*! \file eventprint.cpp
  Dump some events

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005

  $Id: eventprint.cpp,v 3.4 2009/11/05 23:02:20 matthew Exp $

*/

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <exception>

#include <VSOptions.hpp>
#include <VBFDumper.hpp>

using namespace VERITAS;

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

  bool no_l3 = false;
  bool no_channels = false;
  bool dump_raw_clock = false;
  uint32_t dispatch_flags = 0;
  bool one_packet = false;
  unsigned packet_num = 0;
  bool integrate_overflow = false;
  bool lo_gain = false;
  bool one_scope = false;
  unsigned scope_id = 0;

  options.findBoolValue("no_channels", no_channels, true,
			"Do not dump individual channel data");

  options.findBoolValue("dump_raw_clock", dump_raw_clock, true,
			"Dump the raw BCD (hex) from the clock structure");

  options.findBoolValue("integrate_overflow", integrate_overflow, false,
			"Integrate the overflow events into the array event");

  if(integrate_overflow)
    dispatch_flags |= VSSimpleVBFDispatcher::DISPATCH_INTEGRATE_OVERFLOW;

  options.findWithValue("dispatch_flags", dispatch_flags,
			"Set the flags for the VBF dispatcher, see "
			"\"VSSimpleVBF.hpp\" for details.");

  options.findBoolValue("no_l3", no_l3, true, "Do not use L3 time.");

  if(options.findWithValue("packet", packet_num,
			   "Dispach only single, given packet.")
     !=VSOptions::FS_NOT_FOUND)one_packet=true;

  if(options.findWithValue("scope", scope_id,
			   "Dispach only data from a given telescope.")
     !=VSOptions::FS_NOT_FOUND)one_scope=true;

  options.findBoolValue("lo_gain", lo_gain, true,
			"Dispach only data from a given telescope.");
  

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

  if(argc != 1)
    {
      std::cerr << progname << ": need exctly one file name" << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  VBFDumper dumper(std::cout, lo_gain, one_scope, scope_id, 
		   no_channels, dump_raw_clock);
  VSSimpleVBFDispatcher dispatcher(&dumper);
  if(no_l3)dispatcher.setDontUseL3Time(true);
  if(no_channels)dispatcher.setVisitChannels(false);

  try 
    {
      dispatcher.openFile(*argv);
      argc--,argv++;

      if(one_packet)
	dispatcher.dispatchPacket(packet_num, 0, dispatch_flags);
      else
	dispatcher.dispatchAllPackets(0, dispatch_flags);
	 
      dispatcher.closeFile();
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
