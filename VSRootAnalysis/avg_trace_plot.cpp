#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSOptions.hpp>
#include <VSFileUtility.hpp>
#include <VSOctaveIO.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSLineTokenizer.hpp>
#include <Vec3D.hpp>
#include <VSTimeDepPedData.hpp>
#include <VSRunInfoData.hpp>
#include <VSLaserData.hpp>
#include <VSEventData.hpp>
#include <VSNSpace.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSPointing.hpp>
#include <VSResultsData.hpp>
#include <VSChannelMap.hpp>
#include <VSAnalysisStage1.hpp>
#include <VBFLaserCalc.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TH1F.h>
#include <TRint.h>
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
#include "VSRHistogram.hpp"
#include "VSRHistogram1D.hpp"
#include "VSRHistogram2D.hpp"
#include "VSRCanvas.hpp"
#include "VSRGraph1D.hpp"

using namespace VERITAS;
using namespace Physics;

template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname << " [options] <vbf file>" << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}


int main(int argc, char **argv) 
{
  std::string progname(*argv);
  VSOptions options(argc, argv, true);


  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;


  VSRCanvas::configure(options);

  unsigned ichan = 0;

  std::vector<std::string> labels;

  options.findWithValue("labels",labels,"");

  options.findWithValue("chan", ichan, 
			"");

  if(print_usage)
    {
      usage(progname, options, std::cerr);
      exit(EXIT_SUCCESS);
    }

  if(options.argC() == 0)
    {
      std::cerr << "Error: No input files specified."
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  std::list< std::string > file_list;
  for(unsigned iarg = 1; iarg < argc; iarg++)
    file_list.push_back(argv[iarg]);

  

  gROOT->SetStyle("Plain");
  TApplication *app = new TApplication("app",NULL,NULL);


  std::vector<VSRGraph1D*> mean_trace;

  VSRCanvas mean_trace_canvas("mean_trace");

  
  std::vector<std::string>::iterator label_itr = labels.begin();
  // --------------------------------------------------------------------------
  // Loop on files
  // --------------------------------------------------------------------------
  for(std::list<std::string>::iterator fileItr = file_list.begin();
      fileItr != file_list.end(); ++fileItr) 
    {
      if(!VSFileUtility::isFile(*fileItr))
	{
	  std::cerr << "Error: " << *fileItr
		    << " is not a valid file." << std::endl;
	  exit(EXIT_FAILURE);
	}

      VSOctaveH5Reader *reader = new VSOctaveH5Reader(*fileItr);

      VSArrayDiagnosticsData diagnostics;

      diagnostics.load(reader->readStruct("diagnostics"));

      VSScopeDiagnosticsData* scope = diagnostics.scope[0];
      
      VSRGraph1D* mean_trace_graph = new VSRGraph1D;

      for(unsigned isamp = 0; isamp < 24; isamp++)
	mean_trace_graph->set(isamp,scope->channeMeanHitrace(ichan,isamp));

      //	std::cout << scope->channeMeanHitrace(0,isamp) << std::endl;

      mean_trace_canvas.add(0,mean_trace_graph,*label_itr);

      label_itr++;
    }

  mean_trace_canvas.draw();

  app->Run(true);
}

