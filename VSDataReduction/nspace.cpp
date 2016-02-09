/*! \file nspace.cpp
  NSpace manipulation application. When will it ever stop?

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/17/2006

  $Id: nspace.cpp,v 3.10 2008/04/07 03:41:52 matthew Exp $

*/

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <exception>

#include <VSLineTokenizer.hpp>
#include <VSOctaveIO.hpp>
#include <VSOptions.hpp>
#include <VSNSpace.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSNSpaceCutsCalc.hpp>
#include <VSNSpaceUtility.hpp>
#include <VSDatumElementExtractor.hpp>
#include <VSNSpaceDataSource.hpp>
#include <VSNSpaceEventData.hpp>
#include <VSScaledParameterCalc.hpp>
#include <VSRunInfoData.hpp>

using namespace VERITAS;

void main_load_data(const std::vector< std::string >& files,
	       VSNSpace::Space& space,
	       const std::pair< double, double >& theta0_cut,
	       const std::pair< double, double >& theta1_cut,
	       bool has_sp_lookup,
	       const std::string& output_file);

void main_create_cuts(const std::string& on_file,
		      const std::string& off_file,
		      const std::string& output);

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  std::string pad(7+progname.length(),' ');
  stream << "Usage: " << progname 
	 << " -load_data [options] data_file [data_file...]" << std::endl
	 << pad << " -index_data [options] data_file" << std::endl
	 << pad << " -cut_data [options] data_file" << std::endl
	 << pad << " -define_space [options]" << std::endl
	 << pad << " -produce_ordering [options] sig_hist_file bkg_hist_file"
	 << std::endl
    	 << pad << " -create_cuts [options] sig_hist_file bkg_hist_file"
	 << std::endl
	 << pad << " -project_1d [options] hist_file/filter_file dim" 
	 << std::endl
	 << pad << " -project_2d [options] hist_file/filter_file dim1 dim2" 
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  std::vector<std::string> command_line(argc);
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)command_line[iarg]=argv[iarg];
  
  VSOptions options(argc, argv, true);

  bool print_usage              = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  // Variables to hold options ------------------------------------------------

  bool mode_load_data           = false;
  bool mode_create_cuts         = false;

  std::vector<VSNSpaceUtility::AxisDefinition> space_axes;
  std::string space_file        = "";
  std::pair< double, double > theta0_cut;
  std::pair< double, double > theta1_cut;
  unsigned theta_set            = 0;

  std::string ordering_file     = "";
  std::string filter_file       = "";
  std::string output_file       = "";

  std::vector< std::string > scope_filters;
  std::string array_filter = "";
  bool has_sp_lookup = false;

  VSScaledParameterCalc::configure(options,"","");
  VSNSpaceUtility::configure(options);

  if(options.wasOptionFound("sc_parameter_lookup")) has_sp_lookup = true;

  options.findBoolValue("load_data", mode_load_data, true,
                        "Load data from input file(s) into histogram and "
			"write histogram to output file.");

  options.findBoolValue("create_cuts", mode_create_cuts, true,
			"Create a set of cuts from a set of telescope and "
			"array filters.");

  options.findWithValue("space", space_axes,
			"Specify definition of space over which NSpace will "
			"operate, when creating a new histogram. See below "
			"for details.");

  options.findWithValue("theta0", theta0_cut, "Theta cut.");
  options.findWithValue("theta1", theta1_cut, "Theta cut.");

  options.findWithValue("o", output_file, "Set the output file name.");


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

  std::vector< std::string > args;
  for(unsigned iarg = 0; iarg < (unsigned)argc; iarg++)
    args.push_back(argv[iarg]);

  // Define the space ---------------------------------------------------------
  VSNSpace::Space space;
  VSNSpaceUtility::SpaceDefinition space_definition = space_axes;
  if(!space_axes.empty())
    VSNSpaceUtility::createSpace(space_definition, space);

  if(mode_load_data)
    main_load_data(args,space,theta0_cut,theta1_cut,has_sp_lookup,output_file);
  else if(mode_create_cuts)
    main_create_cuts(args[0],args[1],output_file);
}

void main_load_data(const std::vector< std::string >& files,
	       VSNSpace::Space& space,
	       const std::pair< double, double >& theta0_cut,
	       const std::pair< double, double >& theta1_cut,
	       bool has_sp_lookup,
	       const std::string& output_file)
{
  VSNSpaceIO* io = new VSNSpaceOctaveH5IO;

  VSNSpace hist(space);

  // Loop over files ----------------------------------------------------------
  for(std::vector< std::string >::const_iterator itr = files.begin();
      itr != files.end(); ++itr)
    {
      std::cerr << "opening file \"" << *itr << '"' << std::endl;
      
      VSScaledParameterCalc* sp_calc = NULL;

      if(has_sp_lookup)
	{
	  VSOctaveH5Reader *reader = new VSOctaveH5Reader(*itr);
	  vsassert(reader);
	  double zn_mean_deg;
	  double az_mean_deg;

	  reader->readScalar("stage1.run_info.zn_mean_deg",zn_mean_deg);
	  reader->readScalar("stage1.run_info.az_mean_deg",az_mean_deg);
	  
	  VSRunInfoData run_info;
	  VSOctaveH5ReaderStruct* stg1_struct = reader->readStruct("stage1");
	  vsassert(stg1_struct);
	  run_info.load(stg1_struct->readStruct("run_info"));
	  
	  delete reader;
	  
	  sp_calc = new VSScaledParameterCalc(run_info.nchan,zn_mean_deg,
					      az_mean_deg,0);
	}

      VSNSpaceDataSource* data_source = 
	new VSNSpaceEventData(*itr, space, sp_calc);
      
      // Read data, transform and enter it into histogram ---------------------
      unsigned nevents = 0;
      unsigned nloaded = 0;
      std::vector< VSNSpace::Point > points;
      double theta0,theta1;
      while(data_source->getData(points,theta0,theta1))
	{	  
	  nevents++;
	  for(std::vector< VSNSpace::Point >::iterator itr = points.begin();
	      itr != points.end(); ++itr)
	    {
	      if(theta0<theta0_cut.first || 
		 (theta0_cut.second&&(theta0>theta0_cut.second))) continue;
	      else if(theta1<theta1_cut.first || 
		      (theta1_cut.second&&(theta1>theta1_cut.second))) 
		continue;

	      if(hist.accumulate(*itr))nloaded++;
	    }
	  points.clear();
	}
      // Close the file -------------------------------------------------------
      std::cerr << *itr << ": " << nevents << " events found, "
		<< nloaded << " loaded" << std::endl;
      
      delete data_source;
      delete sp_calc;
    }

  // Write the histogram ------------------------------------------------------
  std::string filename(output_file);
  if(filename.empty())filename = "output.h5";
  std::cerr << "main_load_data: writing histogram to \"" << filename << '"'
	    << std::endl;
  io->writeHistogram(filename,hist);
  delete io;
}

void main_create_cuts(const std::string& on_file,
		      const std::string& off_file,
		      const std::string& output)
{
  const char* fn = "main_create_cuts: ";
  VSNSpaceIO* io = new VSNSpaceOctaveH5IO;

  VSNSpaceUtility* nspace_calc = new VSNSpaceUtility;

  // Get the ON and OFF histograms --------------------------------------------
  VSNSpace on;
  VSNSpace off;

  VSNSpaceUtility::loadHist(on_file,on);
  VSNSpaceUtility::loadHist(off_file,off);

  
  VSNSpace::Ordering ordering;
  nspace_calc->createOrdering(on,off,ordering);

  VSNSpace::Volume filter;
  nspace_calc->createFilter(ordering,filter);

  // Set up various elements --------------------------------------------------
  bool is_array_filter = false;
  bool is_scope_filter_default = false;

  VSNSpace::Volume* array_filter = NULL;
  VSNSpace::Volume* scope_filter_default = NULL;
  std::vector<const VSNSpace::Volume*> scope_filters;

  for(unsigned idim = 0; idim < filter.space().ndim; idim++)
    {
      std::string index = 
	VSH5DatumElementParser::getIndex(filter.space().axes[idim].name);

      if(index == "*")
	is_scope_filter_default = true;
      else if(index.empty())
	is_array_filter = true;      
    }

  if(is_scope_filter_default)
    scope_filter_default = &filter;
  else if(is_array_filter)
    array_filter = &filter;

  std::string filename = output;
  if(filename.empty()) filename = "nspace_cuts.h5";

  std::cerr << fn << " Loading cuts to " << filename << std::endl;

  VSOctaveH5Writer* file = new VSOctaveH5Writer(filename,true);
  VSNSpaceCutsCalc* nspace_cuts_calc =
    new VSNSpaceCutsCalc(array_filter, scope_filter_default, scope_filters);

  nspace_cuts_calc->save(file);

  delete nspace_cuts_calc;
  delete file;
  delete io;
}
