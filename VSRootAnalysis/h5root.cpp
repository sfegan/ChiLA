#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSOptions.hpp>
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
#include <VSFileUtility.hpp>
#include <VSScaledParameterCalc.hpp>
#include <VSAnalysisStage1.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TROOT.h>
#include <TH1F.h>
#include <TRint.h>
#include <TPDF.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram.hpp"
#include "VSRHistogram1D.hpp"
#include "VSRHistogram2D.hpp"
#include "VSRH5Loader.hpp"
#include "VSRCanvas.hpp"
#include "VSREventVisitor.hpp"

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
  stream << "Usage: " << progname << " [options] <h5 file> <path>" << std::endl
         << std::endl;
  stream << "Description: " 
	 << "Plot a histogram, graph, or histogram in an H5 file." 
	 << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char **argv) 
{
  // --------------------------------------------------------------------------
  // Parse options
  // --------------------------------------------------------------------------
  std::string progname(*argv);
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  VSScaledParameterCalc::configure(options,"","");

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  argv++,argc--;
  std::list<std::string> file_list;
  for(int i=0;i<argc;i++) file_list.push_back(argv[i]);

  try
    {
      VSREventVisitor* visitor = new VSREventVisitor;

      VSEventDataDispatcher* dispatcher = 
	new VSEventDataDispatcher(visitor);
      
      dispatcher->processFiles(file_list);  

      delete visitor;
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

#if 0
  std::list< std::string > files;

  std::string path = "events";

  TFile f("test.root","recreate");

  TTree *t = new TTree("t1","tree");

  // --------------------------------------------------------------------------
  // Insert the list of files to process into files
  // --------------------------------------------------------------------------
  for(int iarg = 0; iarg < argc; iarg++)
    files.push_back(argv[iarg]);

  for(std::list<std::string>::iterator itr = files.begin();
      itr != files.end(); )
    {
      VSFileUtility::expandFilename(*itr);
      if(!VSFileUtility::isFile(*itr)) 
	{
	  std::cerr << *itr << " is not a valid file." << std::endl;
	  exit(EXIT_FAILURE);
	}
      else if(!VSOctaveH5ReaderStruct::isHDF5(*itr))
	{
	  std::ifstream fileStream(itr->c_str());
	  itr = files.erase(itr);
	  VSLineTokenizer tokenizer;

	  while( !fileStream.eof() )
	    {
	      VSTokenList tokens;
	      std::string tmp;
	      getline(fileStream,tmp);
	      tokenizer.tokenize(tmp,tokens);
	      if(tokens.size() == 0)continue;

	      if(tokens[0].string().find('#') == std::string::npos)
		{
		  std::string h5_file = tokens[0].string();		  
		  VSFileUtility::expandFilename(h5_file);
		  if(!VSFileUtility::isFile(h5_file))
		    {
		      std::cerr << h5_file 
				<< " is not a valid file." << std::endl;
		      exit(EXIT_FAILURE);
		    } 
		  else if(!VSOctaveH5ReaderStruct::isHDF5(h5_file))
		    {
		      std::cerr << h5_file 
				<< " is not a valid HDF5 file." << std::endl;
		      exit(EXIT_FAILURE);
		    }
		  else
		    files.push_back(h5_file);	      
		}
	    }
	}
      else
	++itr;
    }

  VSEventArrayDatum                evt;

  unsigned run_id;

  t->Branch("run",&run_id,"run_id/I");
  t->Branch("rec_mask",&evt.used_in_reconstruction_mask,"rec_mask/I");
  t->Branch("msc_width",&evt.msc_width,"msc_width/D");
  t->Branch("msc_length",&evt.msc_length,"msc_length/D");
  t->Branch("msc_disp",&evt.msc_disp,"msc_disp/D");
  t->Branch("N2",&evt.N2,"N2/D");
  t->Branch("R",&evt.R,"R/D");
  t->Branch("Rx",&evt.Rx,"Rx/D");
  t->Branch("Ry",&evt.Ry,"Ry/D");
  t->Branch("zn",&evt.zn,"zn/D");
  t->Branch("az",&evt.az,"az/D");
  t->Branch("fovx_derotated",&evt.mean_derotated_fov_x,"fovx_derotated/D");
  t->Branch("fovy_derotated",&evt.mean_derotated_fov_y,"fovy_derotated/D");
  t->Branch("fovx",&evt.mean_fov_x,"fovx/D");
  t->Branch("fovy",&evt.mean_fov_y,"fovy/D");
  t->Branch("array_zn",&evt.mean_array_zn,"array_zn/D");
  t->Branch("array_az",&evt.mean_array_az,"array_az/D");
  t->Branch("theta0",&evt.theta0,"theta0/D");

  for(std::list<std::string>::iterator itr = files.begin();
      itr != files.end(); ++itr)
    {
      std::string file = *itr;      
      VSOctaveH5Reader *reader = new VSOctaveH5Reader(file);
      vsassert(reader);

      VSOctaveH5ReaderStruct* s = reader->readStruct(path);
      vsassert(s);

      VSAnalysisStage1Data         stage1;
      VSTargetTable::Observation   obs;
      VSArrayMergedCalibrationData cal;
      VSArrayDiagnosticsData       diagnostics;
      
      stage1.load(reader->readStruct("stage1"));
      obs.load(reader->readStruct("observation"));
      cal.load(reader);
      diagnostics.load(reader->readStruct("diagnostics"));

      run_id = stage1.run_info.run_number;

      VSScaledParameterCalc* sp_calc = new VSScaledParameterCalc;
      sp_calc->load(stage1.run_info.nchan,stage1.run_info.zn_mean_deg,
		    stage1.run_info.az_mean_deg,cal.mean_scaled_dev);


      VSEventDataReader* ed = new VSEventDataReader(s);
      const unsigned nrow = ed->rows();

      std::cout << file << ": " << nrow << std::endl;

      for(unsigned irow=0;irow<nrow;irow++)
	{
	  ed->element(evt,irow);
	  sp_calc->calcSP(evt);

	  if(!std::isfinite(evt.msc_width) || evt.msc_width > 20) continue;
	  else if(evt.theta0 > 1.7) continue;
	  
	  t->Fill();
	}

      delete ed;

      delete sp_calc;
      delete s;
      delete reader;
    }

  t->Write();
#endif
}


