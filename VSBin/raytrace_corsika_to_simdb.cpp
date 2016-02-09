//-*-mode:c++; mode:font-lock;-*-

/*! \file raytrace_corsika_to_simdb.cpp
  Process CORSIKA files through the optics and into the database

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       07/26/2005
*/

#include <VSSimDBVisitor.hpp>

#include <VSOptions.hpp>
#include <VSDBFactory.hpp>
#include <VSTargeting.hpp>
#include <VSHDFEventStore.hpp>
#include <VSHDFEventStoreManyTables.hpp>
#include <VSSimDBCORSIKADatasets.hpp>

#include "RandomNumbers.hpp"

using namespace VERITAS;

int main(int argc, char** argv)
{
  std::string progname(*argv);

  VSOptions options(argc,argv);

  // --------------------------------------------------------------------------
  // Options
  // --------------------------------------------------------------------------

  VSDBFactory::configure(&options);

  std::string rngstatefile(RandomNumbers::defaultFilename());
  options.findWithValue("rng_state_file",rngstatefile,
			"Set random number generator state file");

  bool compress = true;
  bool fake_camera = false;
  float fake_camera_radius_deg = 0;
  std::string tablename;
  uint32_t workunit_runid=0;

  options.findBoolValue("compress", compress, true,
			"Compress data in output files where possible.");
  if(options.findWithValue("fake_camera",fake_camera_radius_deg,
			   "Use a fake camera, filled with pixels out to "
			   "a specified angular radius") == 
     VSOptions::FS_FOUND)fake_camera=true;
  
  options.findWithValue("table",tablename,
			"Database table to which the contents of the "
			"CORSIKA file will be written.");
  options.findWithValue("workunit_runid",workunit_runid,
			"Unique ID identifying this workunit.  If "
			"the workunit ID is specified then the table "
			"option must be defined as well.");
  
  // --------------------------------------------------------------------------
  // Get the database and table names from the command line arguments
  // --------------------------------------------------------------------------

  bool print_usage = false;
  int exit_code = EXIT_SUCCESS;

  if(options.find("h") != VSOptions::FS_NOT_FOUND)print_usage=true;
  if(options.find("help") != VSOptions::FS_NOT_FOUND)print_usage=true;
  
  if(!options.assertNoOptions())
    {
      std::cerr << "Unknown options:";
      for(int iopt=1;iopt<argc;iopt++)
	if(argv[iopt][0]=='-')std::cerr << ' ' << argv[iopt];
      std::cerr << std::endl;
      print_usage=true;
      exit_code=EXIT_FAILURE;
    }
  
  argc--,argv++;

  VSOctaveH5Writer::setDefaultCompress(compress);

  if(options.wasOptionFound("workunit_runid") && 
     !options.wasOptionFound("table"))
    {
      print_usage=true;
      exit_code=EXIT_FAILURE;
    }
  else if(argc < 2)
    {
      print_usage=true;
      exit_code=EXIT_FAILURE;
    }
  
  if(print_usage)
    {
      std::cerr << "Usage: " << std::endl << std::endl
		<< progname << " [options] database corsika_file"
		<< std::endl << std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_code;
    }

  std::string database;
  std::list<std::string> filenames;

  database = *argv;
  argc--; argv++;

  while(argc)
    {
      filenames.push_back(*argv);
      argc--; argv++;      
    }

  // --------------------------------------------------------------------------
  // Create the RNG
  // --------------------------------------------------------------------------

  RandomNumbers rng(rngstatefile.c_str());

  // --------------------------------------------------------------------------
  // Create the database
  // --------------------------------------------------------------------------

  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  if(db->useDatabase(database) < 0)
    {
      std::cerr << "Could not connect to database " << database << std::endl;
      return EXIT_FAILURE;
    };
 
  // --------------------------------------------------------------------------
  // Simulation database
  // --------------------------------------------------------------------------
  
  VSSimDB* db_sim = new VSSimDB(db);
  VSDBParameterTable* db_param = new VSDBParameterTable(db);

  // --------------------------------------------------------------------------
  // Load simulation datasets
  // --------------------------------------------------------------------------

  VSSimDBCORSIKADatasets db_datasets(db);
  VSSimDBWavelengthData* quaneff_data =
    db_datasets.retrieveWavelengthData(VSIMDB_TABLE_NAME_QUANEFF);
  VSSimDBWavelengthData* mirrref_data =
    db_datasets.retrieveWavelengthData(VSIMDB_TABLE_NAME_MIRRREF);
  VSSimDBWavelengthAltitudeData* atmoabs_data =
    db_datasets.retrieveWavelengthAltitudeData(VSIMDB_TABLE_NAME_ATMOABS);

  // --------------------------------------------------------------------------
  // Set up the event store
  // --------------------------------------------------------------------------

  VSSimDBEventStore* event_store = db_sim;
  VSSimDBEventStore* hdf5_event_store = 0;
  VSSimDBTableParam* param = 0;
  unsigned opticsID = 0;

  if(tablename.empty()) 
    {
      hdf5_event_store = new VSHDFEventStoreManyTables(db_sim);
      event_store = hdf5_event_store;
    }
  else
    {
      // ----------------------------------------------------------------------
      // Get the table parameters from the DB and prepare for INSERTs
      // ----------------------------------------------------------------------
      param = db_sim->getDataTableByName(tablename);
      if(!param)
	{
	  std::cerr << "Could not get information for table " << tablename 
		    << std::endl;
	  return EXIT_FAILURE;      
	}

      VSDBParameterSet data_storage_param;
      db_param->retrieveParameterSet("DataStorage", data_storage_param);
      if(data_storage_param["mode"]=="hdf5")
	{
	  hdf5_event_store = 
	    new VSHDFEventStore(db_sim, workunit_runid, param);
	  event_store = hdf5_event_store;
	}

      db_sim->setInsertTableName(tablename,event_store!=db_sim);
      opticsID = param->fOpticsID;
    }

  // --------------------------------------------------------------------------
  // Get the definition of the array from the database
  // --------------------------------------------------------------------------

  VSOTelescopeArray array;
  array.readFromDatabase(db, opticsID);

  // --------------------------------------------------------------------------
  // Set up the CORSIKA visitor
  // --------------------------------------------------------------------------

  VSTargeting* targeting;
  targeting = VSTargetingFactory::getInstance()->getTargeting(rng,db_param);
  if(!targeting)
    {
      std::cerr << "Could not get array pointing for database"
		<< std::endl;
      return EXIT_FAILURE;      
    }

  VSArrayTracking* tracking;
  tracking = VSATFactory::getInstance()->getArrayTracking(db_param);
  if(!tracking)
    {
      std::cerr << "Could not get array pointing for database"
		<< std::endl;
      return EXIT_FAILURE;      
    }

  VSPrimaryArrivalDistribution* pad;
  pad = VSPADFactory::getInstance()->getPrimaryArrivalDistribution(db_param);
  if(!pad)
    {
      std::cerr << "Could not get primary arrival distribution for database"
		<< std::endl;
      return EXIT_FAILURE;      
    }
  
  VSSimDBVisitor* visitor = 
    new VSSimDBVisitor(workunit_runid, event_store, false,
		       &array, &rng, targeting, tracking, pad,
		       quaneff_data, mirrref_data, atmoabs_data);
  if(fake_camera)visitor->setFakeCamera(fake_camera_radius_deg/180.0*M_PI);

  // --------------------------------------------------------------------------
  // INSERT the data
  // --------------------------------------------------------------------------

  VSCORSIKAEventDispatcher dispatcher(visitor);

  for(std::list<std::string>::iterator itr = filenames.begin();
      itr != filenames.end(); ++itr)
    {
      std::cout << itr->c_str() << std::endl;
      dispatcher.processFile(itr->c_str());
    }

  // --------------------------------------------------------------------------
  // Clean up
  // --------------------------------------------------------------------------

  delete visitor;
  delete pad;
  delete tracking;
  delete targeting;
  delete quaneff_data;
  delete mirrref_data;
  delete atmoabs_data;
  delete hdf5_event_store;
  delete param;
  delete db_param;
  delete db_sim;
  delete db;
}
