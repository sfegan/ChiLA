//-*-mode:c++; mode:font-lock;-*-

/*! \file make_optics.cpp

  Program to create random array in DB from array.ini file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    12/05/2004
  \version 0.2
  \note
*/

#include <iostream>
#include <string>
#include <cstdlib>

#include <VSOptions.hpp>
#include <VSDBFactory.hpp>

#include <VSDBParameterTable.hpp>

#include "VSOTelescopeArray.hpp"

using namespace VERITAS;

enum DataSource { DS_NONE, DS_INIFILE, DS_DUMP, DS_DB };

int main(int argc, char ** argv)
{
  std::string progname(*argv);
  VSOptions options(argc,argv);

  // --------------------------------------------------------------------------
  // Options
  // --------------------------------------------------------------------------

  int exit_value = EXIT_SUCCESS;
  bool print_usage = false;
  if(options.find("h","Print this help message")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this help message")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  DataSource ds_from = DS_NONE;
  std::string ds_from_name = "array.ini";
  if(options.findWithValue("from_ini", ds_from_name,
			   "Create random optics instance from specified "
			   "parameter (INI) file. If the parameter file "
			   "name is given as an empty string, the "
			   "canonical values are used.") != 
     VSOptions::FS_NOT_FOUND)ds_from = DS_INIFILE;

  if(options.findWithValue("from_dump", ds_from_name,
			   "Load optics instance from dump file") != 
     VSOptions::FS_NOT_FOUND)ds_from = DS_DUMP;

  if(options.findWithValue("from_db", ds_from_name,
			   "Load optics instance from database. Argument "
			   "should be in the form of database[/optics_id]") != 
     VSOptions::FS_NOT_FOUND)ds_from = DS_DB;


  DataSource ds_to = DS_NONE;
  std::string ds_to_name = "array.dump";
  if(options.findWithValue("to_db", ds_to_name,
			   "Write optics instance and array parameters to "
			   "database. Argument should be in the form of "
			   "database[/optics_id]") != 
     VSOptions::FS_NOT_FOUND)ds_to = DS_DB;

  if(options.findWithValue("to_dump", ds_to_name,
			   "Write optics instance to dump file") != 
     VSOptions::FS_NOT_FOUND)ds_to = DS_DUMP;

  if(options.findWithValue("to_ini", ds_to_name,
			   "Write optics parameters to ini file") != 
     VSOptions::FS_NOT_FOUND)ds_to = DS_INIFILE;

  VSOArrayParameters::setCanonicalValuesFromOptions(options);
  VSDBFactory::configure(&options);

  std::string rngstatefile(RandomNumbers::defaultFilename());
  options.findWithValue("rng_state_file",rngstatefile,
			"Set random number generator state file");
  
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
  else if(argc!=0)
    {
      std::cerr << progname 
		<< ": unwanted arguments required (found " << argc << ")"
		<< std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }  

  if(ds_from == DS_NONE)
    {
      std::cerr << progname 
		<< ": must specify one of \"-from_ini\", \"-from_dump\" and "
		<< "\"-from_db\"" << std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;      
    }

  if(ds_to == DS_NONE)
    {
      std::cerr << progname 
		<< ": must specify one of \"-to_ini\", \"-to_dump\" and "
		<< "\"-to_db\"" << std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;      
    }

  if(print_usage)
    {
      std::cerr << "Usage: " << progname 
		<< " [options]" << std::endl
		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_value;
    }

  // --------------------------------------------------------------------------
  // Connect to the database if required
  // --------------------------------------------------------------------------
  
  VSDatabase* db = 0;
  if((ds_from == DS_DB)||(ds_to == DS_DB))
    db = VSDBFactory::getInstance()->createVSDB();

  // --------------------------------------------------------------------------
  // Create random number generator
  // --------------------------------------------------------------------------

  RandomNumbers rng(rngstatefile.c_str());

  // --------------------------------------------------------------------------
  // Read/Create the array and INI parameters
  // --------------------------------------------------------------------------

  VSOArrayParameters param;
  VSOTelescopeArray array;
  bool have_param = false;

  if(ds_from == DS_INIFILE)
    {
      if(ds_from_name.empty()) 
	{
	  param.reset(true); 
	  have_param=true; 
	}
      else 
	{
	  have_param=param.readFromArrayINIFile(ds_from_name);
	  if(!have_param)
	    {
	      std::cerr << "Could not load array parameters (INI) file " 
			<< ds_from_name << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      array.generateFromArrayParameters(param, rng);
    }
  else if(ds_from == DS_DUMP)
    {
      array.readFromShortDump(ds_from_name);
    }
  else if(ds_from == DS_DB)
    {
      unsigned optics_id = 0;
      std::string database = ds_from_name;
      std::string::size_type ichar = ds_from_name.find('/');
      if((ichar != std::string::npos)
	 &&(ichar < ds_from_name.size()-1))
	{
	  database = ds_from_name.substr(0,ichar);
	  VSDataConverter::fromString(optics_id,ds_from_name.substr(ichar+1));
	}
									
      db->useDatabase(database);
      have_param = param.readFromDatabase(db,optics_id);
      array.readFromDatabase(db,optics_id);
    }

  // --------------------------------------------------------------------------
  // Write/Create the array and INI parameters
  // --------------------------------------------------------------------------

  if(ds_to == DS_INIFILE)
    {
      if(!have_param)
	{
	  std::cerr << "Could not load array parameters (INI) therefore they "
		    << "cannot be saved" << std::endl;
	  return EXIT_FAILURE;
	}
      param.writeToArrayINIFile(ds_to_name);
    }
  else if(ds_to == DS_DUMP)
    {
      array.dumpShort(ds_to_name);
    }
  else if(ds_to == DS_DB)
    {
      unsigned optics_id = 0;
      std::string database = ds_to_name;
      std::string::size_type ichar = ds_to_name.find('/');
      if((ichar != std::string::npos)
	 &&(ichar < ds_to_name.size()-1))
	{
	  database = ds_to_name.substr(0,ichar);
	  VSDataConverter::fromString(optics_id,ds_to_name.substr(ichar+1));
	}
									
      db->createDatabase(database,
			 VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
      db->useDatabase(database);

      // The tables are (possibly) shared by others components/optics
      // in the simulations chain. Do not just delete the tables here.
      // 1) Try to create the table (which does nothing if it already exists)
      // 2) Delete all "ArrayINI" entries in the table
      
      VSOTelescopeArray::createTelescopeArrayTable(db);
      VSOTelescope::createTelescopeTable(db);
      VSOMirror::createMirrorTable(db);
      VSOPixel::createPixelTable(db);
      VSOArrayParameters::createSimulationParametersTable(db);
      
      std::string cond(std::string("OpticsID=")
		       + VSDataConverter::toString(optics_id));
      
      db->deleteFromTable(VSIMDB_TABLE_NAME_ARRAY, cond);
      db->deleteFromTable(VSIMDB_TABLE_NAME_TELESCOPE, cond);
      db->deleteFromTable(VSIMDB_TABLE_NAME_MIRROR, cond);
      db->deleteFromTable(VSIMDB_TABLE_NAME_PIXEL, cond);
      
      cond = std::string(VSOArrayParameters::scCollection)
	+ std::string("[") 
	+ VSDataConverter::toString(optics_id) 
	+ std::string("]");
      
      VSDBParameterTable parameter_table(db);
      parameter_table.deleteParameterSet(cond);

      if(have_param)param.writeToDatabase(db,optics_id);
      array.writeToDatabase(db,optics_id);
    }

  // --------------------------------------------------------------------------
  // Done
  // --------------------------------------------------------------------------
  return EXIT_SUCCESS;
}
