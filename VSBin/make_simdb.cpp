//-*-mode:c++; mode:font-lock;-*-

/*! \file make_simdb.cpp

  Load configuration file and initialize the database

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       08/15/2005
  \note
*/

#include <cerrno>
#include <cmath>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <VSOptions.hpp>
#include <VSLineTokenizer.hpp>
#include <VSDatabase.hpp>
#include <VSDBFactory.hpp>
#include <VSDBParameterTable.hpp>
#include <VSSimDB.hpp>
#include <VSSimDBTables.hpp>
#include <VSSimDBCORSIKADatasets.hpp>
#include <VSTargeting.hpp>
#include <VSOArrayParameters.hpp>
#include <VSOTelescopeArray.hpp>
#include <RandomNumbers.hpp>

using namespace VERITAS;
using namespace Physics;

#if 0
struct DataStorage
{
  enum Storage { S_DATABASE, S_HDF5 };
  Storage storage;
  std::string directory;
  DataStorage(): storage(S_DATABASE), directory() { }
};
#endif 

struct PowerLawSpectrum
{ 
  double index;
  double e_lo;
  double e_hi;
  double flux_e_lo;
};

struct MonoChromaticSpectrum
{
  double e;
  double flux_e;
};

struct SamplingRadiusBase
{
  double zn_lo;
  double zn_hi;
  double radius_m;
};

struct SamplingRadiusLogEnergy
{ 
  double zn_lo;
  double zn_hi;
  double e_lo;
  double e_hi;
  double radius_factor;
};

struct EventCountPerPatch
{
  double zn_lo;
  double zn_hi;
  unsigned event_count;
};

struct WorkunitSize
{
  double e_lo;
  double e_hi;
  unsigned n;
};

struct WorkunitPriority
{
  enum Mode { M_ENERGY, M_ZENITH };
  Mode mode;
  double lo_value;
  double hi_value;
  int adjustment;
};

struct Optics
{
  double zn_lo;
  double zn_hi;
  unsigned optics_id;
};

struct PointingGrid
{
  std::string mode;
  std::string target_mode;
  double zn_lo;
  double zn_hi;
  unsigned zn_n;
  double az_lo;
  double az_hi;
  unsigned az_n;
};

struct GridPoint
{
  std::string target_mode;
  double cos_zn_lo;
  double cos_zn_hi;
  double az_lo;
  double az_hi;
  double zn_mid;
  double az_mid;
  unsigned optics_id;
};

struct DataFile
{
  const char* filename;
  const char* data;
  unsigned data_len;
};

void writeDataFileToStream(const DataFile& df, std::ostream& stream)
{
  stream << df.data;
}

void print_sample_config(std::ostream& stream)
{
#include "db_config_template.h"
  writeDataFileToStream(DB_CONFIG_TEMPLATE,stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  
  VSOptions options(argc,argv);

  // --------------------------------------------------------------------------
  // PROCESS COMMAND LINES
  // --------------------------------------------------------------------------

  int exit_value = EXIT_SUCCESS;
  bool print_usage = false;
  if(options.find("h","Print this help message")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this help message")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  VSDBFactory::configure(&options);

  std::string rngstatefile(RandomNumbers::defaultFilename());
  options.findWithValue("rng_state_file",rngstatefile,
                        "Set random number generator state file");

  bool no_db = false;
  if(options.find("no_db","Do not write configuration to the database, just "
		  "parse the configuration file and write information to "
		  "the console if the \"verbose\" option is given")
     !=VSOptions::FS_NOT_FOUND)
    no_db=true;  

  bool verbose = false;
  if(options.find("verbose","Write additional information to the console "
		  "during processing")!=VSOptions::FS_NOT_FOUND)
    verbose=true;  
  if(options.find("v","Write additional information to the console "
		  "during processing")!=VSOptions::FS_NOT_FOUND)
    verbose=true;  

  bool print_sample = false;
  if(options.find("sample","Print sample configuration file to the console "
		  "and exit")!=VSOptions::FS_NOT_FOUND)
    print_sample=true;  

  bool drop_database = false;
  if(options.find("drop_database","Drop the database (if it exists) before "
		  "loading new data. This option is potentially dangerous "
		  "since ALL EXISTING DATA WILL BE DELETED")
     !=VSOptions::FS_NOT_FOUND)
    drop_database=true;  

  // --------------------------------------------------------------------------
  // ALL KNOWN OPTIONS HAVE BEEN PROCESSED
  // --------------------------------------------------------------------------

  argc--,argv++;

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": Unknown command line options:";
      for(int iopt=0;iopt<argc;iopt++)
	if(*argv[iopt]=='-')std::cerr << ' ' << argv[iopt];
      std::cerr << std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }
  else if(print_sample)
    {
      if(argc==0)
	{
	  print_sample_config(std::cout);
	  return EXIT_SUCCESS;
	}
      else 
	{
	  std::cerr << progname 
		    << ": \"-sample\" option has illegal extra arguments:";
	  for(int iopt=0;iopt<argc;iopt++)std::cerr << ' ' << *argv;
	  std::cerr << std::endl << std::endl;
	  print_usage=true;
	}
    }
  else if(argc!=2)
    {
      std::cerr << progname 
		<< ": two arguments required (found " << argc << ")"
		<< std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }


  if(print_usage)
    {
      std::cerr << "Usage: " << progname 
		<< " [options] database_name config_file" << std::endl
		<< "   or: " << progname << " -sample" << std::endl
		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_value;
    }

  std::string database(*argv);
  argc--,argv++;

  std::string filename(*argv);
  argc--,argv++;

  // --------------------------------------------------------------------------
  // PREPARE PARAMETER SET
  // --------------------------------------------------------------------------

  VSDBParameterSet config_param;
  VSDBParameterSet input_file_param;

  // --------------------------------------------------------------------------
  // RANDOM NUMBER GENERATOR USED WHEN GENERATING OPTICS
  // --------------------------------------------------------------------------

  RandomNumbers rng(rngstatefile.c_str());

  // --------------------------------------------------------------------------
  // OPEN CONFIGURATION FILE AND SET UP PARSER
  // --------------------------------------------------------------------------

  std::ifstream config_file(filename.c_str());
  if(!config_file.good())
    {
      std::cerr << progname << ": could not open " << filename << std::endl
		<< strerror(errno) << std::endl;
      return EXIT_FAILURE;
    }

  VSLineTokenizer tokenizer;

  // --------------------------------------------------------------------------
  // LOCAL VARIABLES TO HOLD PARAMETERS /AND/ LIST OF KNOWN DIRECTIVES
  // --------------------------------------------------------------------------

  VSDBParameterSet                        data_storage;
  unsigned                                primary_id                  = 1;
  VSPrimaryArrivalDistribution*           pad                         = 0;
  double                                  spectrum_bins_per_decade    = 16;
  double                                  spectrum_e0_gev             = 10;
  std::vector<PowerLawSpectrum>           spectrum_power_law;
  std::vector<MonoChromaticSpectrum>      spectrum_monochromatic;
  unsigned                                spectrum_minimum            = 0;
  std::vector<SamplingRadiusBase>         sampling_radius_base;
  std::vector<SamplingRadiusLogEnergy>    sampling_radius_log_energy;
  std::vector<PointingGrid>               pointing_grid;
  VSTargeting*                            targeting                   = 0;
  VSArrayTracking*                        tracking                    = 0;
  double                                  event_count_multiplier      = 1;
  std::vector<EventCountPerPatch>         event_count_per_patch;
  double                                  work_unit_size_fraction     = 0.01;
  std::vector<WorkunitSize>               work_unit_size;
  std::vector<WorkunitPriority>           work_unit_priority;
  std::vector<Optics>                     optics;
  std::map<unsigned, VSOTelescopeArray*>  optics_array;
  std::map<unsigned, VSOArrayParameters*> optics_param;
  VSDBParameterSet                        make_steering_param;
  std::string                             load_modtran_profile;
  std::string                             load_atmospheric_absorption;
  std::string                             load_quantum_efficiency;
  double                                  scale_quantum_efficiency    = 1.0;
  std::string                             load_mirror_reflectivity;

  std::map<std::string, unsigned> directives;
  directives["datastorage"]=0;
  directives["primary"]=2;
  directives["primaryarrivaldistribution"]=0;
  directives["spectrumresolution"]=3;
  directives["spectrumpowerlaw"]=5;
  directives["spectrummonochromatic"]=3;
  directives["spectrumminimumperenergybin"]=2;
  directives["samplingradiusbase"]=4;
  directives["samplingradiuslogenergy"]=6;
  directives["pointinggrid"]=9;
  directives["targeting"]=0;
  directives["arraypointing"]=0;
  directives["eventcountperpatch"]=4;
  directives["eventcountmultiplier"]=2;
  directives["workunitsizefraction"]=2;
  directives["workunitsize"]=4;
  directives["workunitpriority"]=5;
  directives["optics"]=4;
  directives["loadoptics"]=4;
  directives["makesteeringcommandline"]=0;
  directives["loadmodtranprofile"]=2;
  directives["loadatmosphericabsorption"]=2;
  directives["loadquantumefficiency"]=2;
  directives["scalequantumefficiency"]=2;
  directives["loadmirrorreflectivity"]=2;

  // --------------------------------------------------------------------------
  // PARSE CONFIGURATION FILE
  // --------------------------------------------------------------------------

  if(verbose)std::cout << "Configuration options:" << std::endl;

  unsigned line_no = 0;
  while(!config_file.eof())
    {
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // TOKENIZE LINE
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      VSTokenList tokens;
      tokenizer.tokenize(config_file,tokens);
      if(tokens.size() == 0)continue;

      line_no++;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // ADD LINE TO THE PARAMETER SET (FOR LATER WRITING TO THE DATABASE)
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      std::string key = VSDataConverter::toString(line_no);
      std::string val = tokens[0].escaped();
      for(unsigned itoken=1; itoken<tokens.size(); itoken++)
	val += std::string(" ")+tokens[itoken].escaped();
      input_file_param[key]=val;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // WRITE TOKENS TO STDOUT IF VERBOSE
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      if(verbose)
	{
	  std::cout << '"' << tokens[0].string() << '"';
	  for(unsigned i=1;i<tokens.size();i++)
	    std::cout << ' ' << '"' <<  tokens[i].string() << '"';
	  std::cout << std::endl;
	}

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // TEST IF DIRECTIVE IS KNOWN AND CONTAINS REQUIRED NUMBER OF ARGUMENTS
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      key = tokens[0].lower();
      if(directives.find(key) == directives.end())
	{
	  std::cerr << filename << ": directive \"" << tokens[0].string() 
		    << "\" not known" << std::endl;
	  return EXIT_FAILURE;
	}
	  
      if((directives[key]>0)&&(tokens.size()!=directives[key]))
	{
	  std::cerr << filename << ": directive \"" << tokens[0].string() 
		    << "\" requires " << directives[key]-1 << " values (found "
		    << tokens.size()-1 << ")" << std::endl;
	  return EXIT_FAILURE;
	}

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // USE TOKENS TO SET PARAMATERS
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      if(key == "datastorage")
	{
	  if(tokens.size() == 1)
	    {
	      std::cerr << progname 
			<< ": directive \"datastorage\" requires at least "
			<< "1 value" << std::endl;
	      return EXIT_FAILURE;
	    }
	  else if(tokens[1].lower() == "database")
	    {
	      if(tokens.size() != 2)
		{
		  std::cerr << progname 
			    << ": data storage mode \"database\" must not "
			    << "have any parameters" << std::endl;
		  return EXIT_FAILURE;
		}
	      data_storage["mode"]="database";
	    }
	  else if(tokens[1].lower() == "hdf5")
	    {
	      if(tokens.size() != 3)
		{
		  std::cerr << progname 
			    << ": data storage mode \"hdf5\" must have one"
			    << "parameter, a directory to store the files" 
			    << std::endl;
		  return EXIT_FAILURE;
		}
	      data_storage["mode"]="hdf5";
	      data_storage["directory"]=tokens[2].string();
	    }
	  else
	    {
	      std::cerr << progname 
			<< ": unknown data storage mode \"" 
			<< tokens[1].string() << '"' << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      else if(key == "primary")
	{
#define NUCLEUS(A,Z) ((A)*100+(Z))
	  std::string pri_key = tokens[1].lower();
	  if(pri_key == "gamma")primary_id=1;
	  else if(pri_key == "proton")primary_id=14;
	  else if(pri_key == "electron")primary_id=3;
	  else if(pri_key == "muon")primary_id=6;
	  else if(pri_key == "helium")primary_id=NUCLEUS(4,2);
	  else if(pri_key == "iron")primary_id=NUCLEUS(56,26);
	  else tokens[1].convertTo(primary_id);
	}
      else if(key == "primaryarrivaldistribution")
	{
	  std::string error_string;
	  if(pad)delete pad;
	  pad = 
	    VSPADFactory::getInstance()->
	    getPrimaryArrivalDistribution(tokens,&error_string);
	  if(pad==0)
	    {
	      std::cerr << filename << ": " << error_string << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      else if(key == "spectrumresolution")
	{
	  tokens[1].convertTo(spectrum_bins_per_decade);
	  tokens[2].convertTo(spectrum_e0_gev);
	}
      else if(key == "spectrumpowerlaw")
	{
	  PowerLawSpectrum pls;
	  tokens[1].convertTo(pls.index);
	  tokens[2].convertTo(pls.e_lo);
	  tokens[3].convertTo(pls.e_hi);
	  tokens[4].convertTo(pls.flux_e_lo);
	  spectrum_power_law.push_back(pls);
	}
      else if(key == "spectrummonochromatic")
	{
	  MonoChromaticSpectrum mcs;
	  tokens[1].convertTo(mcs.e);
	  tokens[2].convertTo(mcs.flux_e);
	  spectrum_monochromatic.push_back(mcs);
	}
      else if(key == "spectrumminimumperenergybin")
	{
	  tokens[1].convertTo(spectrum_minimum);
	}
      else if(key == "samplingradiusbase")
	{
	  SamplingRadiusBase srb;	
	  tokens[1].convertTo(srb.zn_lo);
	  tokens[2].convertTo(srb.zn_hi);
	  tokens[3].convertTo(srb.radius_m);
	  sampling_radius_base.push_back(srb);
	}
      else if(key == "samplingradiuslogenergy")
	{
	  SamplingRadiusLogEnergy srle;
	  tokens[1].convertTo(srle.zn_lo);
	  tokens[2].convertTo(srle.zn_hi);
	  tokens[3].convertTo(srle.e_lo);
	  tokens[4].convertTo(srle.e_hi);
	  tokens[5].convertTo(srle.radius_factor);
	  sampling_radius_log_energy.push_back(srle);
	}
      else if(key == "pointinggrid")
	{
	  PointingGrid ptgd;

	  ptgd.mode = tokens[1].lower();
	  ptgd.target_mode = tokens[2].lower();
	  tokens[3].convertTo(ptgd.zn_lo);
	  tokens[4].convertTo(ptgd.zn_hi);
	  tokens[5].convertTo(ptgd.zn_n);
	  tokens[6].convertTo(ptgd.az_lo);
	  tokens[7].convertTo(ptgd.az_hi);
	  tokens[8].convertTo(ptgd.az_n);

	  pointing_grid.push_back(ptgd);
	}
      else if(key == "targeting")
	{
	  std::string error_string;
	  if(targeting)delete targeting;
	  targeting =
	    VSTargetingFactory::getInstance()->
	    getTargeting(rng,tokens,&error_string);
	  if(targeting==0)
	    {
	      std::cerr << filename << ": " << error_string << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      else if(key == "arraypointing")
	{
	  std::string error_string;
	  if(tracking)delete tracking;
	  tracking =
	    VSATFactory::getInstance()->getArrayTracking(tokens,&error_string);
	  if(tracking==0)
	    {
	      std::cerr << filename << ": " << error_string << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      else if(key == "eventcountperpatch")
	{
	  EventCountPerPatch ecpp;
	  tokens[1].convertTo(ecpp.zn_lo);
	  tokens[2].convertTo(ecpp.zn_hi);
	  tokens[3].convertTo(ecpp.event_count);
	  event_count_per_patch.push_back(ecpp);
	}
      else if(key == "eventcountmultiplier")
	{
	  tokens[1].convertTo(event_count_multiplier);
	}
      else if(key == "workunitsizefraction")
	{
	  tokens[1].convertTo(work_unit_size_fraction);
	}
      else if(key == "workunitsize")
	{
	  WorkunitSize wus;
	  tokens[1].convertTo(wus.e_lo);
	  tokens[2].convertTo(wus.e_hi);
	  tokens[3].convertTo(wus.n);
	  work_unit_size.push_back(wus);
	}
      else if(key == "workunitpriority")
	{
	  WorkunitPriority wup;
	  if(tokens[1].lower() == "energy")
	    wup.mode = WorkunitPriority::M_ENERGY;
	  else if(tokens[1].lower() == "zenith")
	    wup.mode = WorkunitPriority::M_ZENITH;
	  else
	    {
	      std::cerr << progname 
			<< ": unknown work unit priority mode \"" 
			<< tokens[1].string() << '"' << std::endl;
	      return EXIT_FAILURE;
	    }
	  tokens[2].convertTo(wup.lo_value);
	  tokens[3].convertTo(wup.hi_value);
	  tokens[4].convertTo(wup.adjustment);
	  work_unit_priority.push_back(wup);
	}
      else if(key == "optics")
	{
	  Optics opt;
	  tokens[1].convertTo(opt.zn_lo);
	  tokens[2].convertTo(opt.zn_hi);
	  tokens[3].convertTo(opt.optics_id);
	  optics.push_back(opt);
	}
      else if(key == "loadoptics")
	{
	  unsigned optics_id = 0;
	  tokens[1].convertTo(optics_id);

	  if(tokens[2].lower() == "frominifile")
	    {
	      VSOArrayParameters* param = new VSOArrayParameters;
	      if(param->readFromArrayINIFile(tokens[3].string()) == false)
		{
		  std::cerr << progname 
			    << ": could not load optics ArrayINI file: " 
			<< '"' << tokens[3].string() << '"' << std::endl;
		  delete param;
		  return EXIT_FAILURE;
		}
	      
	      VSOTelescopeArray* array = new VSOTelescopeArray;
	      array->generateFromArrayParameters(*param, rng);

	      optics_param[optics_id] = param;
	      optics_array[optics_id] = array;
	    }
	  else if(tokens[2].lower() == "fromdumpfile")
	    {
	      VSOTelescopeArray* array = new VSOTelescopeArray;
	      if(array->readFromShortDump(tokens[3].string()) == false)
		{
		  std::cerr << progname 
			    << ": could not load optics dump file: " 
			<< '"' << tokens[3].string() << '"' << std::endl;
		  delete array;		  
		  return EXIT_FAILURE;
		} 	      

	      optics_array[optics_id] = array;
	    }
	  else
	    {
	      std::cerr << filename 
			<< ": directive LoadOptics: unrecognised source type: "
			<< '"' << tokens[2].string() << '"' << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      else if(key == "makesteeringcommandline")
	{
	  if(tokens.size()==2)
	    {
	      make_steering_param[tokens[1].string()]="";
	    }
	  else if(tokens.size()==3)
	    {
	      make_steering_param[tokens[1].string()]=tokens[2].string();
	    }
	  else
	    {
	      std::cerr << filename 
			<< ": directive MakeSteeringCommandLine requires "
			<< "one or two arguments (found " << tokens.size()-1
			<< ")" << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      else if(key == "loadmodtranprofile")
	{
	  load_modtran_profile = tokens[1].string();
	}
      else if(key == "loadatmosphericabsorption")
	{
	  load_atmospheric_absorption = tokens[1].string();
	}
      else if(key == "loadquantumefficiency")
	{
	  load_quantum_efficiency = tokens[1].string();
	}
      else if(key == "scalequantumefficiency")
	{
	  tokens[1].convertTo(scale_quantum_efficiency);
	}
      else if(key == "loadmirrorreflectivity")
	{
	  load_mirror_reflectivity = tokens[1].string();
	}
      else
	{
	  std::cerr << progname << ": stupid programmer error: unhandled key: "
		    << key << std::endl;
	  assert(0);
	}
    }

  // --------------------------------------------------------------------------
  // TEST THE DATA STORAGE HAS BEEN SPECIFIED
  // ------------------------------------------------------------------------

  if(data_storage.size() == 0)
    {
      std::cerr << filename
		<< ": must specify \"datastorage\" directive"
		<< std::endl;
      return EXIT_FAILURE;
    }

  // --------------------------------------------------------------------------
  // LOAD THE CORSIKA DATA SETS
  // --------------------------------------------------------------------------

  VSSimDBWavelengthDataset*          ds_quantum_efficiency(0);
  VSSimDBWavelengthDataset*          ds_mirror_reflectivity(0);
  VSSimDBWavelengthAltitudeDataset*  ds_atmospheric_absorption(0);
  VSSimDBModtranProfileDataset*      ds_modtran_profile(0);

  if(!load_quantum_efficiency.empty())
    {
      ds_quantum_efficiency = 
	VSSimDBWavelengthDataset::createFromCORSIKA(load_quantum_efficiency);

      if(!ds_quantum_efficiency)
	{
	  std::cerr << progname 
		    << ": could not load quantum efficiency file: "
		    << load_quantum_efficiency
		    << std::endl;
	  return EXIT_FAILURE;
	}

      if((scale_quantum_efficiency > 0)&& (scale_quantum_efficiency != 1.0))
	{
	  std::string scale_str = 
	    VSDataConverter::toString(scale_quantum_efficiency);

	  config_param["ScaleQuantumEfficiency"] = scale_str;
	  ds_quantum_efficiency->comment += 
	    std::string(" [scaled by make_simdb: ") 
	    + scale_str + std::string("]");

	  for(VSSimDBWavelengthData::iterator iwl = 
		ds_quantum_efficiency->data.begin(); 
	      iwl != ds_quantum_efficiency->data.end(); iwl++)
	    {
	      iwl->second *= scale_quantum_efficiency;
	      if(iwl->second > 1.0)
		{
		  std::cerr << progname 
			    << ": scaled quantum efficiency exceeds 1.0: "
			    << iwl->first << '/' << iwl->second
			    << std::endl;
		  return EXIT_FAILURE;
		}
	    }
	  
	}
    }
  
  if(!load_mirror_reflectivity.empty())
    {
      ds_mirror_reflectivity = 
	VSSimDBWavelengthDataset::createFromCORSIKA(load_mirror_reflectivity);

      if(!ds_mirror_reflectivity)
	{
	  std::cerr << progname 
		    << ": could not load mirror reflectivity file: "
		    << load_mirror_reflectivity
		    << std::endl;
	  return EXIT_FAILURE;
	}
    }

  if(!load_atmospheric_absorption.empty())
    {
      ds_atmospheric_absorption = 
	VSSimDBWavelengthAltitudeDataset::
	createFromCORSIKA(load_atmospheric_absorption);

      if(!ds_atmospheric_absorption)
	{
	  std::cerr << progname 
		    << ": could not load atmospheric absorption file: "
		    << load_atmospheric_absorption
		    << std::endl;
	  return EXIT_FAILURE;
	}
    }

  if(!load_modtran_profile.empty())
    {
      ds_modtran_profile = 
	VSSimDBModtranProfileDataset::
	createFromCORSIKA(load_modtran_profile);

      if(!ds_modtran_profile)
	{
	  std::cerr << progname 
		    << ": could not load modtran atmospheric profile: "
		    << load_modtran_profile
		    << std::endl;
	  return EXIT_FAILURE;
	}
    }

  typedef std::map<unsigned, float>              VSSimDBAltitudeData;
  typedef std::map<unsigned,VSSimDBAltitudeData> VSSimDBWavelengthAltitudeData;
  
  // --------------------------------------------------------------------------
  // CHECK ARRAY POINTING AND PRIMARY ARRIVAL DISTIBUTION
  // --------------------------------------------------------------------------

  if(targeting==0)
    {
      std::cerr << filename << ": must contain \"Targeting\" directive"
		<< std::endl;
      return EXIT_FAILURE;
    }

  if(tracking==0)
    {
      std::cerr << filename << ": must contain \"ArrayPointing\" directive"
		<< std::endl;
      return EXIT_FAILURE;
    }

  if(pad==0)
    {
      std::cerr << filename 
		<< ": must contain \"PrimaryArrivalDistribution\" directive"
		<< std::endl;
      return EXIT_FAILURE;
    }

  // --------------------------------------------------------------------------
  // BUILD SPECTRUM
  // --------------------------------------------------------------------------

  std::map<int, double> spectrum;
  double log_e0 = log10(spectrum_e0_gev);
  
  for(std::vector<PowerLawSpectrum>::const_iterator ipl = 
	spectrum_power_law.begin(); ipl!=spectrum_power_law.end(); ipl++)
    {
      int elo = int(round((log10(ipl->e_lo)-log_e0)*spectrum_bins_per_decade));
      int ehi = int(round((log10(ipl->e_hi)-log_e0)*spectrum_bins_per_decade));
      
      for(int ebin=elo;ebin<=ehi;ebin++)
	{
	  double df = 
	    double(ipl->flux_e_lo)
	    *pow10(ipl->index*double(ebin-elo)/spectrum_bins_per_decade);
	  if(df>spectrum[ebin])spectrum[ebin]=df;
	}
    }
  
  for(std::vector<MonoChromaticSpectrum>::const_iterator imono = 
	spectrum_monochromatic.begin(); 
      imono!=spectrum_monochromatic.end(); imono++)
    {
      int ebin = int(round((log10(imono->e)-log_e0)*spectrum_bins_per_decade));
      if(imono->flux_e>spectrum[ebin])spectrum[ebin]=imono->flux_e;
    }

  // --------------------------------------------------------------------------
  // BUILD TARGET GRID
  // --------------------------------------------------------------------------

  std::vector<GridPoint> grid;

  for(std::vector<PointingGrid>::iterator itr = pointing_grid.begin();
      itr != pointing_grid.end(); ++itr)
    {
      if((itr->target_mode != "center") &&(itr->target_mode != "patch")
	 &&(itr->target_mode != "fixedzenith")
	 &&(itr->target_mode != "fixedzazimuth"))
	{
	  std::cerr << "Unknown target distribution mode \"" 
		    << itr->target_mode << '"' << std::endl;
	  return EXIT_FAILURE;
	}

      bool fix_az = 
	itr->target_mode=="center" || itr->target_mode=="fixedazimuth";

      if(itr->mode == "equalareaequalazimuth")
	{
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	  // EQUAL AREA BINS EQUAL AZIMUTH STEPS VARIABLE ZENITH STEPS
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	  double dcoszn = 
	    cos(itr->zn_lo/180.0*M_PI) 
	    - cos(itr->zn_hi/180.0*M_PI);
	  double daz = 
	    itr->az_hi - itr->az_lo;
      
	  for(unsigned znbin=0;znbin<itr->zn_n;znbin++)
	    {
	      double y = double(znbin)/double(itr->zn_n);
	      double Y = double(znbin+1)/double(itr->zn_n);

	      double cos_zn_lo = cos(itr->zn_lo/180.0*M_PI) - dcoszn*y;
	      double cos_zn_hi = cos(itr->zn_lo/180.0*M_PI) - dcoszn*Y;

	      for(unsigned azbin=0;azbin<itr->az_n;azbin++)
		{
		  double x = double(azbin)/double(itr->az_n);
		  double X = double(azbin+1)/double(itr->az_n);
	      
		  double az_lo = itr->az_lo + daz*x;
		  double az_hi = itr->az_lo + daz*X;
	      
		  GridPoint patch;
		  patch.target_mode = itr->target_mode;
		  patch.cos_zn_lo = cos_zn_lo;
		  patch.cos_zn_hi = cos_zn_hi;
		  patch.az_lo = az_lo;
		  patch.az_hi = az_hi;

		  if((cos_zn_lo == 1.0)&&(az_lo==0)&&(az_hi==360)&&(fix_az))
		    {
		      patch.zn_mid=0;
		      patch.az_mid=0;
		    }
		  else
		    {
		      patch.zn_mid=
			acos(0.5*(patch.cos_zn_hi+patch.cos_zn_lo))/M_PI*180.0;
		      patch.az_mid=0.5*(az_lo+az_hi);
		    }

		  patch.optics_id = 0;

		  grid.push_back(patch);
		}
	    }
	}
      else if(itr->mode == "equalareaequalzenith")
	{
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	  // EQUAL AREA BINS (APPROX) EQUAL ZENITH STEPS VARIABLE AZIMUTH STEPS
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	  double mean_step_zn = (itr->zn_hi - itr->zn_lo)/double(itr->zn_n);
      
	  double mean_area_max = 
	    cos((itr->zn_hi-mean_step_zn)/180.0*M_PI)
	    - cos(itr->zn_hi/180.0*M_PI);

	  double daz = itr->az_hi - itr->az_lo;

	  if(itr->az_n == 0)
	    {
	      double n_equator = 
		double(itr->zn_n)* (itr->az_hi-itr->az_lo)
		/ (itr->zn_hi-itr->zn_lo);
	      double zn_mean = (itr->zn_hi-mean_step_zn/2)/180.0*M_PI;
	      itr->az_n = unsigned(round(n_equator * (1.0-cos(zn_mean))));
	    }

	  double zn_lo = itr->zn_lo;

	  for(unsigned znbin=0;znbin<itr->zn_n;znbin++)
	    {
	      double zn_hi = zn_lo + mean_step_zn;

	      double mean_area = cos(zn_lo/180.0*M_PI) - cos(zn_hi/180.0*M_PI);
	      double mean_n_az = double(itr->az_n) * mean_area/mean_area_max;

	      unsigned n_az = unsigned(round(mean_n_az));
	      if(n_az==0)n_az=1;

	      double area = double(n_az)*mean_area_max/double(itr->az_n);
	  
	      zn_hi = acos((cos(zn_lo/180.0*M_PI) - area))/M_PI*180.0;

	      for(unsigned azbin=0;azbin<n_az;azbin++)
		{
		  double x = double(azbin)/double(n_az);
		  double X = double(azbin+1)/double(n_az);
	      
		  double az_lo = itr->az_lo + daz*x;
		  double az_hi = itr->az_lo + daz*X;
	      
		  GridPoint patch;
		  patch.target_mode = itr->target_mode;
		  patch.cos_zn_lo = cos(zn_lo/180.0*M_PI);
		  patch.cos_zn_hi = cos(zn_hi/180.0*M_PI);
		  patch.az_lo = az_lo;
		  patch.az_hi = az_hi;
		  if((zn_lo == 0)&&(az_hi-az_lo==360)&&(fix_az))
		    {
		      patch.zn_mid=0;
		      patch.az_mid=0;
		    }
		  else
		    {
		      patch.zn_mid=
			acos(0.5*(patch.cos_zn_hi+patch.cos_zn_lo))/M_PI*180.0;
		      patch.az_mid=0.5*(az_lo+az_hi);
		    }
		  patch.optics_id = 0;

		  grid.push_back(patch);
		}
	  
	      zn_lo = zn_hi;
	    }
	}
      else if(itr->mode == "equalzenithequalazimuth")
	{
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	  // EQUAL ZENITH STEPS EQUAL AZIMUTH STEPS
	  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

	  double mean_step_zn = (itr->zn_hi - itr->zn_lo)/double(itr->zn_n);
	  double daz = itr->az_hi - itr->az_lo;
      
	  double zn_lo = itr->zn_lo;

	  for(unsigned znbin=0;znbin<itr->zn_n;znbin++)
	    {
	      double zn_hi = zn_lo + mean_step_zn;
	      unsigned n_az;

	      if(zn_lo == 0) n_az = 1;
	      else n_az = itr->az_n;
	    
	      for(unsigned azbin=0;azbin<n_az;azbin++)
		{
		  double x = double(azbin)/double(n_az);
		  double X = double(azbin+1)/double(n_az);

		  double az_lo = itr->az_lo + daz*x;
		  double az_hi = itr->az_lo + daz*X;

		  GridPoint patch;
		  patch.target_mode = itr->target_mode;
		  patch.cos_zn_lo = cos(zn_lo/180.0*M_PI);
		  patch.cos_zn_hi = cos(zn_hi/180.0*M_PI);
		  patch.az_lo = az_lo;
		  patch.az_hi = az_hi;

		  if((zn_lo == 0)&&(az_hi-az_lo==360)&&(fix_az))
		    {
		      patch.zn_mid=0;
		      patch.az_mid=0;
		    }
		  else
		    {
		      patch.zn_mid=0.5*(zn_hi+zn_lo);
		      patch.az_mid=0.5*(az_lo+az_hi);
		    }

		  patch.optics_id = 0;

		  grid.push_back(patch);
		}

	      zn_lo = zn_hi;
	    }
	}
      else
	{
	  std::cerr << "Unknown grid distribution \"" << itr->mode 
		    << std::endl;
	  return EXIT_FAILURE;
	}
    }

  // --------------------------------------------------------------------------
  // WRITE GRID TO STDOUT IF VERBOSE
  // --------------------------------------------------------------------------

  if(verbose)
    {
      std::cout << std::endl
		<< "Grid:" << std::endl;
      for(std::vector<GridPoint>::const_iterator igrid = grid.begin();
	  igrid!=grid.end(); igrid++)
	{
	  std::cout 
	    << acos(igrid->cos_zn_lo)/M_PI*180.0 << ' '
	    << acos(igrid->cos_zn_hi)/M_PI*180.0 << ' '
	    << igrid->zn_mid << ' '
	    << igrid->az_lo << ' '
	    << igrid->az_hi << ' '
	    << igrid->az_mid << ' '
	    << igrid->optics_id
	    << std::endl;
	}
    }

  // --------------------------------------------------------------------------
  // SET OPTICS
  // --------------------------------------------------------------------------

  for(std::vector<Optics>::const_iterator ioptics = optics.begin();
      ioptics!=optics.end(); ioptics++)
    {
      for(std::vector<GridPoint>::iterator igrid = grid.begin();
	  igrid!=grid.end(); igrid++)
	{
	  if((igrid->zn_mid>=ioptics->zn_lo)&&(igrid->zn_mid<ioptics->zn_hi))
	    igrid->optics_id=ioptics->optics_id;
	}
    }

  // --------------------------------------------------------------------------
  // BUILD WORKUNIT SIZES
  // --------------------------------------------------------------------------

  std::map<int, unsigned> spectrum_workunit;
  for(std::vector<WorkunitSize>::const_iterator iwus = work_unit_size.begin();
      iwus!=work_unit_size.end(); iwus++)
    {
      int elo = 
	int(ceil((log10(iwus->e_lo)-log_e0)*spectrum_bins_per_decade)-0.1);
      int ehi = 
	int(floor((log10(iwus->e_hi)-log_e0)*spectrum_bins_per_decade)+0.1);

      for(int ebin=elo;ebin<=ehi;ebin++)
	spectrum_workunit[ebin]=iwus->n;
    }

  // --------------------------------------------------------------------------
  // BUILD ENERGY DEPENDENT PRIORITIES
  // --------------------------------------------------------------------------

  std::map<int, int> spectrum_priority;
  for(std::vector<WorkunitPriority>::const_iterator iwup = 
	work_unit_priority.begin(); iwup!=work_unit_priority.end(); iwup++)
    if(iwup->mode == WorkunitPriority::M_ENERGY)
      {
	int elo = int(ceil((log10(iwup->lo_value)-log_e0)
			   *spectrum_bins_per_decade)-0.1);
	int ehi = int(floor((log10(iwup->hi_value)-log_e0)
			    *spectrum_bins_per_decade)+0.1);
	
	for(int ebin=elo;ebin<=ehi;ebin++)
	  spectrum_priority[ebin] += iwup->adjustment;
      }

  // --------------------------------------------------------------------------
  // BUILD THE TABLES AND WRITE SPECTRUM FOR THIS ELEVATION
  // --------------------------------------------------------------------------

  std::vector<VSSimDBTableParam> tables;
  std::set<double> spectrum_written;

  for(std::vector<GridPoint>::const_iterator igrid = grid.begin();
      igrid!=grid.end(); igrid++)
    {
      bool write_spectrum = false;
      if((verbose)
	 &&(spectrum_written.find(igrid->zn_mid) == spectrum_written.end()))
	{
	  spectrum_written.insert(igrid->zn_mid);
	  write_spectrum=true;
	  std::cout << std::endl
		    << "Spectrum for patch at zn=" << igrid->zn_mid 
		    << std::endl;
	}

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // BUILD SAMPLING RADIUS
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#define RADIUSFN(E) log10(log10(E))
      double radius_base = 0;
      for(std::vector<SamplingRadiusBase>::const_iterator isrb =
	    sampling_radius_base.begin(); isrb!=sampling_radius_base.end();
	  isrb++)
	if((igrid->zn_mid>=isrb->zn_lo)&&(igrid->zn_mid<=isrb->zn_hi))
	  radius_base = isrb->radius_m;
      
      std::map<int, double> spectrum_sampling_radius;
      for(std::map<int, double>::const_iterator ispec = spectrum.begin();
	  ispec!=spectrum.end(); ispec++)
	spectrum_sampling_radius[ispec->first] = radius_base;

      for(std::vector<SamplingRadiusLogEnergy>::const_iterator isrle =
	    sampling_radius_log_energy.begin(); 
	  isrle!=sampling_radius_log_energy.end(); isrle++)
	if((igrid->zn_mid>=isrle->zn_lo)&&(igrid->zn_mid<=isrle->zn_hi))
	  for(std::map<int, double>::const_iterator ispec = spectrum.begin();
	      ispec!=spectrum.end(); ispec++)
	    {
	      double logelo = log10(isrle->e_lo);
	      double logehi = log10(isrle->e_hi);
	      
	      double x0 = RADIUSFN(isrle->e_lo);
	      double x1 = RADIUSFN(isrle->e_hi);
	      
	      int elo = 
		int(ceil((logelo-log_e0)*spectrum_bins_per_decade)-0.1);
	      int ehi = 
		int(floor((logehi-log_e0)*spectrum_bins_per_decade)+0.1);
	      
	      for(int ebin=elo;ebin<=ehi;ebin++)
		{
		  double x = 
		    RADIUSFN(pow10(double(ebin)
				   /spectrum_bins_per_decade+log_e0));
		  double f = (x-x0)/(x1-x0);
		  double r = radius_base*(1.0+(isrle->radius_factor-1.0)*f);
		  if(r>spectrum_sampling_radius[ebin])
		    spectrum_sampling_radius[ebin]=r;
		}
	    }
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // FIND TOTAL NUMBER OF COUNTS AT THIS ZENITH ANGLE
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      unsigned event_count = 0;
      for(std::vector<EventCountPerPatch>::const_iterator iecpp = 
	    event_count_per_patch.begin();
	  iecpp!=event_count_per_patch.end(); iecpp++)
	if((igrid->zn_mid>=iecpp->zn_lo)&&(igrid->zn_mid<=iecpp->zn_hi))
	  event_count = iecpp->event_count;
      event_count = 
	unsigned(round(double(event_count)*event_count_multiplier));
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // BUILD SPECTRUM COUNTS
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      double integral_flux = 0;
      for(std::map<int, double>::const_iterator ispec = spectrum.begin();
	  ispec!=spectrum.end(); ispec++)
	integral_flux += 
	  ispec->second
	  *spectrum_sampling_radius[ispec->first]
	  *spectrum_sampling_radius[ispec->first];

      std::map<int, unsigned> spectrum_counts;
      for(std::map<int, double>::const_iterator ispec = spectrum.begin();
	  ispec!=spectrum.end(); ispec++)
	{
	  unsigned n = 
	    lrint(double(event_count)
		  *ispec->second
		  *spectrum_sampling_radius[ispec->first]
		  *spectrum_sampling_radius[ispec->first]
		  /integral_flux);
	  if(n<spectrum_minimum)n=spectrum_minimum;
	  spectrum_counts[ispec->first] = n;
#if 0
	  std::cerr << event_count << ' '
	 	    << ispec->second << ' '
		    << spectrum_sampling_radius[ispec->first] << ' '
		    << integral_flux << ' ' << n << ' '
		    << lrint(double(event_count)
			   *ispec->second
			   *spectrum_sampling_radius[ispec->first]
			   *spectrum_sampling_radius[ispec->first]
			   /integral_flux) << '\n';
#endif
	}

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // BUILD WORKUNIT COUNTS
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      std::map<int, unsigned> workunit_counts;
      for(std::map<int, unsigned>::const_iterator iscnt = 
	    spectrum_counts.begin(); iscnt!=spectrum_counts.end(); iscnt++)
	{
	  unsigned n = 
	    unsigned(ceil(double(iscnt->second)
			  *work_unit_size_fraction));
	  if((spectrum_workunit.find(iscnt->first) != spectrum_workunit.end())
	     &&(spectrum_workunit[iscnt->first]<n))
	    n=spectrum_workunit[iscnt->first];
	  if(n==0)n=1;
	  workunit_counts[iscnt->first]=n;
	}	

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // WORKUNIT PRIORITY
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      int base_priority = 0;
      for(std::vector<WorkunitPriority>::const_iterator iwup = 
	    work_unit_priority.begin(); iwup!=work_unit_priority.end(); iwup++)
	if((iwup->mode == WorkunitPriority::M_ZENITH)
	   &&(igrid->zn_mid>=iwup->lo_value)&&(igrid->zn_mid<=iwup->hi_value))
	  base_priority += iwup->adjustment;

      std::map<int, int> priority;
      for(std::map<int, double>::const_iterator ispec = spectrum.begin();
	  ispec!=spectrum.end(); ispec++)
	priority[ispec->first] = base_priority+spectrum_priority[ispec->first];

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // BUILD SPECTRUM COUNTS
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      for(std::map<int, double>::const_iterator ispec = spectrum.begin();
	  ispec!=spectrum.end(); ispec++)
	{
	  int ebin = ispec->first;
	  double e = pow10(double(ebin)/spectrum_bins_per_decade+log_e0);

	  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
	  // WRITE SPECTRUM
	  // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
	  if(write_spectrum)
	    std::cout << ispec->first << ' ' 
		      << e << ' ' 
		      << spectrum_counts[ebin] << ' '
		      << workunit_counts[ebin] << ' '
		      << spectrum_sampling_radius[ebin] << ' '
		      << priority[ebin] << std::endl;

	  bool fix_az = 
	    igrid->target_mode=="center" || 
	    igrid->target_mode=="fixedazimuth";

	  bool fix_zn = 
	    igrid->target_mode=="center" || 
	    igrid->target_mode=="fixedzenith";

	  VSSimDBTableParam table;
	  table.fPrimaryID            = primary_id;
	  table.fEnergyGeV            = e;
	  table.fZenithMinRad         = 
	    fix_zn?igrid->zn_mid:acos(igrid->cos_zn_lo);
	  table.fZenithMaxRad         =
	    fix_zn?igrid->zn_mid:acos(igrid->cos_zn_hi);
	  table.fAzimuthMinRad        =
	    fix_az?igrid->az_mid:igrid->az_lo;
	  table.fAzimuthMaxRad        =
	    fix_az?igrid->az_mid:igrid->az_hi;
	  table.fOpticsID = igrid->optics_id;
	  table.fSamplingRadiusM      = spectrum_sampling_radius[ebin];
	  table.fTargetEventCount     = spectrum_counts[ebin];
	  table.fWorkunitEventCount   = workunit_counts[ebin];
	  table.fWorkunitPriority     = priority[ebin];

	  if(fix_zn)table.fZenithMinRad *= M_PI/180.0;
	  if(fix_zn)table.fZenithMaxRad *= M_PI/180.0;
	  table.fAzimuthMinRad *= M_PI/180.0;
	  table.fAzimuthMaxRad *= M_PI/180.0;

	  tables.push_back(table);
	}
    }

  // --------------------------------------------------------------------------
  // WRITE TABLE COUNT / EVENT COUNT INFORMATION TO STDOUT IF VERBOSE
  // --------------------------------------------------------------------------

  if(verbose)
    {
      std::cout << std::endl
		<< "Total number of tables: " << tables.size()
		<< std::endl;
      
      unsigned nevents=0;
      for(std::vector<VSSimDBTableParam>::const_iterator itab =
	    tables.begin(); itab!=tables.end(); itab++)
	nevents += itab->fTargetEventCount;

      std::cout << "Total number of events: " << nevents
		<< std::endl;
    }

  // --------------------------------------------------------------------------
  // EXIT NOW IF WE ARE NOT ALLOWED USE THE DATABASE
  // --------------------------------------------------------------------------

  if(no_db)
    {
      delete tracking;
      delete pad;
      return EXIT_SUCCESS;
    }

  // --------------------------------------------------------------------------
  // CONNECT TO THE DATABASE
  // --------------------------------------------------------------------------

  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  if(db==0)
    {
      std::cerr << progname << ": could not connect to database server" 
		<< std::endl;
      return EXIT_FAILURE;
    }

  if(drop_database)
    db->dropDatabase(database,VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);

  db->createDatabase(database,VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(db->useDatabase(database) < 0)
    {
      std::cerr << progname << ": could not connect to database " 
		<< database << std::endl;
      return EXIT_FAILURE;
    }

  VSDBParameterTable db_param(db);
  VSSimDB db_sim(db);

  db_param.createParameterTable();
  
  // --------------------------------------------------------------------------
  // LOAD OPTICS
  // --------------------------------------------------------------------------

  VSOTelescopeArray::createTelescopeArrayTable(db);
  VSOTelescope::createTelescopeTable(db);
  VSOMirror::createMirrorTable(db);
  VSOPixel::createPixelTable(db);

  for(std::map<unsigned, VSOTelescopeArray*>::const_iterator iarray =
	optics_array.begin(); iarray!=optics_array.end(); iarray++)
    {
      unsigned optics_id = iarray->first;
      std::string cond(std::string("OpticsID=")
		       + VSDataConverter::toString(optics_id));
      db->deleteFromTable(VSIMDB_TABLE_NAME_ARRAY, cond);
      db->deleteFromTable(VSIMDB_TABLE_NAME_TELESCOPE, cond);
      db->deleteFromTable(VSIMDB_TABLE_NAME_MIRROR, cond);
      db->deleteFromTable(VSIMDB_TABLE_NAME_PIXEL, cond);
      iarray->second->writeToDatabase(db, optics_id);
    }

  // --------------------------------------------------------------------------
  // LOAD CORSIKA DATAFILES
  // --------------------------------------------------------------------------

  if((ds_quantum_efficiency)||(ds_mirror_reflectivity)
     ||(ds_atmospheric_absorption)||(ds_modtran_profile))
    {
      VSSimDBCORSIKADatasets db_datasets(db);
      
      if(ds_quantum_efficiency)
	db_datasets.storeWavelengthDataset(VSIMDB_TABLE_NAME_QUANEFF,
					   *ds_quantum_efficiency);
      
      if(ds_mirror_reflectivity)
	db_datasets.storeWavelengthDataset(VSIMDB_TABLE_NAME_MIRRREF,
					   *ds_mirror_reflectivity);
      
      if(ds_atmospheric_absorption)
	db_datasets.storeWavelengthAltitudeDataset(VSIMDB_TABLE_NAME_ATMOABS,
						   *ds_atmospheric_absorption);
      
      if(ds_modtran_profile)
	db_datasets.storeModtranProfileDataset(VSIMDB_TABLE_NAME_MODTRAN,
					       *ds_modtran_profile);
    }

  // --------------------------------------------------------------------------
  // WRITE PARAMETERS
  // --------------------------------------------------------------------------

  db_param.deleteParameterSet("MakeSimDB");
  db_param.storeParameterSet("MakeSimDB",config_param);
  db_param.deleteParameterSet("MakeSimDBConfigFile");
  db_param.storeParameterSet("MakeSimDBConfigFile",input_file_param);

  db_param.deleteParameterSet("DataStorage");
  db_param.storeParameterSet("DataStorage",data_storage);

  if(make_steering_param.size() != 0)
    {
      db_param.deleteParameterSet("MakeSteering");
      db_param.storeParameterSet("MakeSteering",make_steering_param);
    }
  
  VSTargetingFactory::getInstance()->writeTargeting(targeting, &db_param);
  VSATFactory::getInstance()->writeArrayTracking(tracking, &db_param);
  VSPADFactory::getInstance()->
    writePrimaryArrivalDistribution(pad, &db_param);

  for(std::map<unsigned, VSOArrayParameters*>::const_iterator iparam =
	optics_param.begin(); iparam!=optics_param.end(); iparam++)
    {
      unsigned optics_id = iparam->first;
      std::string name = std::string(VSOArrayParameters::scCollection)
	+ std::string("[")
	+ VSDataConverter::toString(optics_id)
	+ std::string("]");
      db_param.deleteParameterSet(name);
      iparam->second->writeToDatabase(db, optics_id);
    }

  // --------------------------------------------------------------------------
  // GENERATE TABLES
  // --------------------------------------------------------------------------

  db_sim.createInfrastructureTables();

  for(std::vector<VSSimDBTableParam>::iterator itab =
	tables.begin(); itab!=tables.end(); itab++)
    {
      itab->fTableName            = db_sim.tableName(*itab);
      db_sim.createDataTable(*itab, data_storage["mode"]!="database");
    }
  
  // --------------------------------------------------------------------------
  // CLEAN UP AND EXIT
  // --------------------------------------------------------------------------

  delete ds_quantum_efficiency;
  delete ds_mirror_reflectivity;
  delete ds_atmospheric_absorption;
  delete ds_modtran_profile;

  for(std::map<unsigned, VSOArrayParameters*>::iterator iparam =
	optics_param.begin(); iparam!=optics_param.end(); iparam++)
    delete iparam->second;

  for(std::map<unsigned, VSOTelescopeArray*>::iterator iarray =
	optics_array.begin(); iarray!=optics_array.end(); iarray++)
    delete iarray->second;
  
  delete db;
  delete tracking;
  delete pad;
  return EXIT_SUCCESS;
}
