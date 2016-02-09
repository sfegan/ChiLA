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
#include "VSRHistogramCamera.hpp"
#include "VSRCanvas.hpp"

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

struct PixelStatus
{
  unsigned m_scopeID;
  unsigned m_pixelID;
  float    m_absgain;
  float    m_relgain;
  float    m_releff;
  float    m_l1rate;
  float    m_pedrms;
  float    m_current;
};

int main(int argc, char **argv) 
{
  const unsigned max_scopes = 4;
  const unsigned max_pixels = 500;

  std::string progname(*argv);

  VSOptions options(argc, argv, true);

  std::string output_file;
  std::string sort_mode;
  unsigned window = 24;

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  options.findWithValue("o", output_file, 
			"Dump channel data to a text file.");

  options.findWithValue("sort", sort_mode, 
			"Sort channels (absgain, l1rate, ).");

  options.findWithValue("window", window, 
			"Sort channels (absgain, l1rate, ).");

  if(print_usage)
    {
      usage(progname, options, std::cerr);
      exit(EXIT_SUCCESS);
    }

  if(options.argC() < 2)
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


  VSRCanvas ped_canvas;
  VSRCanvas current_canvas;
  VSRCanvas current_camera_canvas;

  // --------------------------------------------------------------------------
  // Loop on files
  // --------------------------------------------------------------------------
  for(std::list<std::string>::iterator fileItr = file_list.begin();
      fileItr != file_list.end(); ++fileItr) 
    {
      VSAnalysisStage1Data stage1_data;
 

      if(!VSFileUtility::isFile(*fileItr))
	{
	  std::cerr << "Error: " << *fileItr
		    << " is not a valid file." << std::endl;
	  exit(EXIT_FAILURE);
	}
      else if(VSOctaveH5ReaderStruct::isHDF5(*fileItr))
	{
	  VSOctaveH5Reader *reader = new VSOctaveH5Reader(*fileItr);
	  stage1_data.load(reader->readStruct("stage1"));
	  delete reader;
	}

      VSMiscellaneousDBData* db_data      = stage1_data.misc_db;
      VSTimeDepPedData*      ped_data     = &stage1_data.pedestals;

      if(db_data == NULL)
	db_data = new VSMiscellaneousDBData;

      std::vector< VSRHistogram1D* > ped_hist(max_scopes);

      for(unsigned iscope = 0; iscope < max_scopes; iscope++)
	ped_hist[iscope] = new VSRHistogram1D(200,0,20);

      //      VSRHistogram1D* ped_hist = new VSRHistogram1D(200,0,20);
      //     VSRHistogramCamera* current_camera_hist = new VSRHistogramCamera;
      VSRHistogram1D* current_hist = new VSRHistogram1D(25,0,10);

      std::vector< double >       mean_ped(max_scopes);

      // ----------------------------------------------------------------------
      // Pedestals
      // ----------------------------------------------------------------------
      unsigned nscope = ped_data->nscope();
      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  std::vector< double >       peds;
	  double mean = 0;
	  
	  const unsigned nchan = ped_data->nchan(iscope);
	  for(unsigned ichan = 0; ichan < nchan; ichan++)
	    { 
	      float pedrms = ped_data->dev(iscope,ichan,window);
	      ped_hist[iscope]->fill(pedrms);
	      mean += pedrms;
	      
	      peds.push_back(pedrms);
	    }
	  
	  mean /= (double)nchan;
	  mean_ped[iscope]=mean;
	  
	  std::sort(peds.begin(),peds.end());
	  
	  std::cout << iscope << " " 
		    << *(peds.begin()+peds.size()/2) << std::endl;
	  ped_canvas.add(iscope,ped_hist[iscope]);
	}

      // ----------------------------------------------------------------------
      // Currents 
      // ----------------------------------------------------------------------
      std::vector< float > avg_current;

      
      unsigned nrecord = 0;

 //      if(db_data->scope.size() != 0)
// 	nrecord = db_data->scope[iscope].hv_status.size();

//       for(unsigned irecord = 0; irecord < nrecord; irecord++)
// 	{
// 	  const unsigned nchan = 
// 	    db_data->scope[iscope].hv_status[irecord].chan.size();
// 	  for(unsigned ichan = 0; ichan < nchan; ichan++)
// 	    {
// 	      if(ichan >= avg_current.size())
// 		avg_current.resize(ichan+1);
	      
// 	      avg_current[ichan] += 
// 		db_data->scope[iscope].
// 		hv_status[irecord].chan[ichan].current;		  
// 	    }
// 	}
      
//       for(unsigned ichan = 0; ichan < avg_current.size(); ichan++)
// 	{
// 	  avg_current[ichan] /= (float)nrecord;
// 	  current_camera_hist->setPixel(ichan,avg_current[ichan]); 
// 	  current_hist->fill(avg_current[ichan]);
// 	}

 
//       current_canvas.addElement(0,current_hist);
//       current_camera_canvas.addElement(0,current_camera_hist);
    }

  //  ped_canvas.setErrorOpt(0,false);
  //  current_canvas.setErrorOpt(0,false);


  ped_canvas.draw();
//   current_canvas.createCanvas("currents")->Draw();
//   current_camera_canvas.createCanvas("currents_camera")->Draw();

  app->Run(true);
}

