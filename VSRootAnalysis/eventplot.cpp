#include <iostream>
#include <vector>

#include <VSOptions.hpp>
#include <VSAnalysisStage1.hpp>


// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TH1F.h>
#include <TApplication.h>
#include <TPDF.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TMultiGraph.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSREventVisitor.hpp"
#include "VSRHistogram.hpp"
#include "VSRHistogram1D.hpp"
#include "VSRHistogram2D.hpp"
#include "VSRCanvas.hpp"

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

  VSRCanvas::configure(options);

  uint32_t dispatch_flags = 
    VSSimpleVBFDispatcher::DISPATCH_IN_ORDER |
    VSSimpleVBFDispatcher::DISPATCH_VERBOSE;
  bool one_packet = false;
  unsigned packet_num = 1;
  unsigned npackets = 0;
  std::string stage1_file = "";

  options.findWithValue("stage1", stage1_file,
			"");

  options.findWithValue("packet", packet_num,
			"");

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

  gROOT->SetStyle("Plain");
  TApplication *app = new TApplication("app",NULL,NULL);


  VSRCanvas* canvas = new VSRCanvas("trace");
  

  VSOctaveH5Reader reader(stage1_file);
  VSOctaveH5ReaderStruct* s = reader.readStruct("stage1");
  VSAnalysisStage1Data stage1_data;
  stage1_data.load(s);
  delete s;

  VSREventVisitor event_visitor(stage1_data);
  VSSimpleVBFDispatcher dispatcher(&event_visitor);
//   if(no_l3)dispatcher.setDontUseL3Time(true);
//   if(no_channels)dispatcher.setVisitChannels(false);

  try 
    {
      dispatcher.openFile(*argv);
      argc--,argv++;

      dispatcher.dispatchPacket(packet_num, 0, dispatch_flags);
      //dispatcher.dispatchAllPackets(0, dispatch_flags);	 
      dispatcher.closeFile();
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      return EXIT_FAILURE;
    }

  event_visitor.draw();
//   for(std::vector< std::pair< unsigned, unsigned > >::iterator itr =
// 	channel_list.begin(); itr != channel_list.end(); ++itr)
//     {
//       VSRHistogram1D* hist = new VSRHistogram1D(24,0,24);
//       std::vector< double > trace = 
// 	trace_visitor.getAverageTrace(itr->first,itr->second);
//       const unsigned nsample = trace.size();
//       for(unsigned isample = 0; isample < nsample; isample++)
// 	hist->set(isample,-trace[isample]);

//       canvas->add(0,hist);
//     }

//  canvas->draw();

  app->Run(true);

  return EXIT_SUCCESS;
}
