//-*-mode:c++; mode:font-lock;-*-

/*! \file manage_simdb.cpp

  Dump the simulations database

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
#include <VSTargeting.hpp>
#include <VSHDFEventStore.hpp>

using namespace VERITAS;
using namespace Physics;

class TablePriority
{
public:
  TablePriority(unsigned id=0, int p=0): fTableID(id), fPriority(p) { }  
  unsigned fTableID;
  int fPriority;
  bool operator < (const TablePriority& o) const
  { return (fPriority>o.fPriority)
      ||((fPriority==o.fPriority)&&(fTableID<o.fTableID)); }
};

int main(int argc, char** argv)
{
  std::string progname(*argv);
  
  VSOptions options(argc,argv,true);

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
  int mode_selected=0;
  bool create_db = false;

  bool get_workunit = false;
  if(options.find("get_workunit",
		  "Find a table which is not complete, print its name to the "
		  "console and exit")!=VSOptions::FS_NOT_FOUND)
    get_workunit=true, mode_selected++;

  bool start_workunit_run = false;
  if(options.find("start_workunit_run",
		  "Register start of workunit on node, print unique "
		  "workunit run ID to console, then exit")
     !=VSOptions::FS_NOT_FOUND)
    start_workunit_run=true, mode_selected++;

  bool finish_workunit_run = false;
  if(options.find("finish_workunit_run",
		  "Register completion of workunit then "
		  "exit")!=VSOptions::FS_NOT_FOUND)
    finish_workunit_run=true, mode_selected++;

  bool set_disable_run = false;
  if(options.find("set_disable_run",
		  "Set database run terminate flag, commanding all jobs "
		  "inserting events into this database running to terminate, "
		  "then exit")!=VSOptions::FS_NOT_FOUND)
    set_disable_run=true, mode_selected++;

  bool clear_disable_run = false;
  if(options.find("clear_disable_run",
		  "Clear database run terminate flag, allowing jobs to insert "
		  "events into this database, then exit")
     !=VSOptions::FS_NOT_FOUND)
    clear_disable_run=true, mode_selected++;

  bool set_min_priority = false;
  int min_priority = 0;
  if(options.findWithValue("set_min_priority",min_priority,
			   "Set minimum priority that table must have to be "
			   "processed, then exit")!=VSOptions::FS_NOT_FOUND)
    set_min_priority=true, mode_selected++;

  bool clear_min_priority = false;
  if(options.find("clear_min_priority",
		  "Clear minimum priority that table must have to be "
		  "processed, allowing jobs to insert events into any table, "
		  "then exit")
     !=VSOptions::FS_NOT_FOUND)
    clear_min_priority=true, mode_selected++;

  bool get_parameter = false;
  if(options.find("get_parameter",
		  "Get the value of a parameter from the database")
     !=VSOptions::FS_NOT_FOUND)
    get_parameter=true, mode_selected++;

  bool mark_all_workunit_events_as_complete = false;
  if(options.find("mark_all_workunit_events_as_complete",
		  "Mark all events in workunit as being complete")
     !=VSOptions::FS_NOT_FOUND)
    mark_all_workunit_events_as_complete=true, mode_selected++;

  bool get_hdf_name = false;
  if(options.find("get_hdf_name", "Get the name of the output HDF file")
     !=VSOptions::FS_NOT_FOUND)
    get_hdf_name=true, mode_selected++;

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

  if(mode_selected>1)
    {
      std::string spaces(' ',progname.size()+2);
      std::cerr << progname << ": Only one option from:" << std::endl
		<< spaces << "-get_workunit"
		<< std::endl
		<< spaces << "-start_workunit_run"
		<< std::endl
		<< spaces << "-finish_workunit_run"
		<< std::endl
		<< spaces << "-set_disable_run"
		<< std::endl
		<< spaces << "-clear_disable_run"
		<< std::endl
		<< spaces << "-set_min_priority"
		<< std::endl
		<< spaces << "-clear_min_priority"
		<< std::endl
		<< spaces << "-get_parameter"
		<< std::endl
		<< spaces << "-mark_all_workunit_events_as_complete"
		<< std::endl
		<< spaces << "-get_hdf_name"
		<< std::endl
		<< "may be selected at any time." << std::endl
		<< std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }
  
  if(((mode_selected==0)&&(argc!=1))
     ||((mode_selected==1)&&(get_workunit)&&(argc!=1))
     ||((mode_selected==1)&&(start_workunit_run)&&(argc!=4))
     ||((mode_selected==1)&&(finish_workunit_run)&&(argc!=2))
     ||((mode_selected==1)&&(set_disable_run)&&(argc!=1))
     ||((mode_selected==1)&&(clear_disable_run)&&(argc!=1))
     ||((mode_selected==1)&&(set_min_priority)&&(argc!=1))
     ||((mode_selected==1)&&(clear_min_priority)&&(argc!=1))
     ||((mode_selected==1)&&(get_parameter)&&(argc!=3))
     ||((mode_selected==1)&&(mark_all_workunit_events_as_complete)&&(argc!=3))
     ||((mode_selected==1)&&(get_hdf_name)&&(argc!=1)))
    {
      std::cerr << progname 
		<< ": incorrect number of arguments found (" << argc << ')'
		<< std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }

  if(print_usage)
    {
      std::cerr << "Usage: " << progname 
		<< " [options] database_name" << std::endl
		<< "   or: " << progname
		<< " -start_workunit_run [options] database_name table_name host_name job_id" << std::endl
		<< "   or: " << progname
		<< " -finish_workunit_run [options] database_name workunit_run_id" << std::endl

 		<< "   or: " << progname
		<< " -set_disable_run [options] database_name" << std::endl
		<< "   or: " << progname
		<< " -clear_disable_run [options] database_name" << std::endl

 		<< "   or: " << progname
		<< " -set_min_priority=priority [options] database_name" << std::endl
		<< "   or: " << progname
		<< " -clear_min_priority [options] database_name" << std::endl

		<< "   or: " << progname
		<< " -get_parameter [options] database_name collection_name parameter_name" << std::endl

		<< "   or: " << progname
		<< " -mark_all_workunit_events_as_complete [options] database_name table_name workunit_run_id" << std::endl

		<< "   or: " << progname
		<< " -get_hdf_name [options] workunit_run_id" << std::endl

		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_value;
    }


  // --------------------------------------------------------------------------
  // GET THE DEFAULT NAME OF THE OUTPUT HDF NAME
  // --------------------------------------------------------------------------

  if(get_hdf_name)
    {
      unsigned workunit_run_id = 0;
      VSDataConverter::fromString(workunit_run_id,*argv);
      argc--,argv++;

      std::cout 
	<< VSHDFEventStore::nameForWorkunitID(workunit_run_id)
	<< std::endl;

      exit(EXIT_SUCCESS);
    }

  // --------------------------------------------------------------------------
  // GET THE DATABASE NAME
  // --------------------------------------------------------------------------

  std::string database(*argv);
  argc--,argv++;

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

  if(create_db)
    db->createDatabase(database,
		       VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(db->useDatabase(database) < 0)
    {
      std::cerr << progname << ": could not connect to database " 
		<< database << std::endl;
      return EXIT_FAILURE;
    }
     
  VSDBParameterTable db_param(db);
  VSSimDB db_sim(db);

  // --------------------------------------------------------------------------
  // REGISTER START / FINISH OF WORKUNIT IF REQUESTED
  // --------------------------------------------------------------------------

  if(start_workunit_run)
    {
      std::string tablename(*argv);
      argc--,argv++;

      VSSimDBTableParam* table = db_sim.getDataTableByName(tablename);
      if(!table)
	{
	  std::cerr << progname << ": could not find table " << tablename
		    << std::endl;
	  return EXIT_FAILURE;
	}

      std::string hostname(*argv);
      argc--,argv++;
      if(hostname.empty())
	{
	  char hostbuffer[1000];
	  gethostname(hostbuffer,1000);
	  hostname=hostbuffer;
	}

      unsigned jobid=getppid();
      VSDataConverter::fromString(jobid,*argv);
      argc--,argv++;

      uint64_t workunit_run_id = 
	db_sim.registerWorkunitRunStart(table->fTableID,hostname,jobid);
      
      std::cout << workunit_run_id << std::endl;
      
      delete table;
      
      return EXIT_SUCCESS;
    }

  if(finish_workunit_run)
    {
      uint64_t workunit_run_id;
      VSDataConverter::fromString(workunit_run_id,*argv);
      argc--,argv++;

      db_sim.registerWorkunitRunFinish(workunit_run_id);
      
      return EXIT_SUCCESS;
    }

  // --------------------------------------------------------------------------
  // ENABLE/DISABLE DATABASE AND SET/DELETE MINIMUM PRIORITY
  // --------------------------------------------------------------------------

  if(set_disable_run)
    {
      db_param.createParameterTable();
      db_param.setOrUpdateParameterValue("ManageSimDB","DisableRun","1");
      return EXIT_SUCCESS;
    }

  if(clear_disable_run)
    {
      db_param.deleteParameter("ManageSimDB","DisableRun");
      return EXIT_SUCCESS;
    }

  if(set_min_priority)
    {
      std::string min_priority_string;
      VSDataConverter::toString(min_priority_string,min_priority);
      db_param.createParameterTable();
      db_param.setOrUpdateParameterValue("ManageSimDB","MinimumPriority",
					 min_priority_string);
      return EXIT_SUCCESS;
    }

  if(clear_min_priority)
    {
      db_param.deleteParameter("ManageSimDB","MinimumPriority");
      return EXIT_SUCCESS;
    }

  // --------------------------------------------------------------------------
  // GET PARAMETER
  // --------------------------------------------------------------------------

  if(get_parameter)
    {
      VSDBParameterSet parameters;
      db_param.retrieveParameterSet(*argv++, parameters);
      if(parameters.empty())return EXIT_FAILURE;
      
      std::string key(*argv++);
      if(parameters.find(key) == parameters.end())return EXIT_FAILURE;
      
      std::cout << parameters[key] << std::endl;
      return EXIT_SUCCESS;
    }

  // --------------------------------------------------------------------------
  // MARK ALL WORKUNIT EVENTS AS COMPLETE
  // --------------------------------------------------------------------------

  if(mark_all_workunit_events_as_complete)
    {
      std::string table_prefix(*argv);
      argc--,argv++;

      std::string table_name = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ table_prefix
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);
      
      unsigned workunit_run_id;
      VSDataConverter::fromString(workunit_run_id,*argv);
      argc--,argv++;

      std::ostringstream query;
      query << "UPDATE " << table_name
	    << " SET EventComplete=1 WHERE WorkunitRunID=?";

      VSDBStatement* stmt = 
	db->createQuery(query.str(),VSDatabase::FLAG_NO_SERVER_PS);
      stmt->bindToParam(workunit_run_id);
      stmt->execute();

      return EXIT_SUCCESS;
    }

  // --------------------------------------------------------------------------
  // GET ALL TABLE INFORMATION
  // --------------------------------------------------------------------------

  std::vector<VSSimDBTableParam> tables = db_sim.getAllDataTables();
  if(tables.size()==0)
    {
      std::cerr << progname << ": no tables found in database" << std::endl;
      return EXIT_FAILURE;
    }

  std::vector<uint32_t> counts = db_sim.getEventCountOfAllTables(tables,true);
  if(counts.size()!=tables.size())
    {
      std::cerr << progname << ": number of table counts (" 
		<< counts.size() << ") not equal number of table entries ("
		<< tables.size() << ")" << std::endl;
      return EXIT_FAILURE;
    }

  VSDBParameterSet parameters;
  db_param.retrieveParameterSet("ManageSimDB", parameters);

  // --------------------------------------------------------------------------
  // IF REQUESTED: FIND A WORKUNIT AT RANDOM AND EXIT
  // --------------------------------------------------------------------------

  if(get_workunit)
    {
      if(tables.size()==0)return EXIT_SUCCESS;
      if(parameters.find("DisableRun") != parameters.end())
	{
	  bool disable_run = true;
	  VSDataConverter::fromString(disable_run,
				      parameters["DisableRun"]);
	  if(disable_run)return EXIT_SUCCESS;
	}

      RandomNumbers rng(rngstatefile.c_str());

      std::vector<unsigned> cumulative_units;
      cumulative_units.reserve(tables.size());
      unsigned total_units = 0;

      bool found_one_unfilled_table = false;
      int highest_priority = 0;
      for(unsigned itab=0;itab<tables.size();itab++)
	{
	  int events = tables[itab].fTargetEventCount-counts[itab];
	  int priority = tables[itab].fWorkunitPriority;

	  if((events>0)
	     &&((!found_one_unfilled_table)||(priority>highest_priority)))
	    {
	      highest_priority=tables[itab].fWorkunitPriority;
	      found_one_unfilled_table = true;
	    }
	}
      
      if(!found_one_unfilled_table)return EXIT_SUCCESS;

      if(parameters.find("MinimumPriority") != parameters.end())
	{
	  VSDataConverter::fromString(min_priority,
				      parameters["MinimumPriority"]);
	  if(highest_priority < min_priority)return EXIT_SUCCESS;
	}

      for(unsigned itab=0;itab<tables.size();itab++)
	{
	  int events = tables[itab].fTargetEventCount-counts[itab];
	  int priority = tables[itab].fWorkunitPriority;
	  if((events<=0)||(priority < highest_priority))events=0;
	  unsigned workunits = 
	    unsigned(ceil(double(events)
			  /double(tables[itab].fWorkunitEventCount)));
	  total_units += (unsigned)std::pow(workunits,2.0);    
	  cumulative_units.push_back(total_units);
	}

      assert(total_units!=0);
      
      unsigned x = unsigned(floor(rng.Uniform()*double(total_units)));

      for(unsigned itab=0;itab<tables.size();itab++)
	if(x < cumulative_units[itab])
	  {
	    std::cout << tables[itab].fTableName << std::endl;
	    return EXIT_SUCCESS;	    
	  }
      
      // should never get here
      assert(0);
    }
  
  // --------------------------------------------------------------------------
  // PRINT LIST OF TABLES
  // --------------------------------------------------------------------------

  for(unsigned itab=0;itab<tables.size();itab++)
    {
      int events = tables[itab].fTargetEventCount-counts[itab];
      if(events<=0)events=0;
      unsigned workunits = 
	unsigned(ceil(double(events)
		      /double(tables[itab].fWorkunitEventCount)));
      
      std::cout << tables[itab].fTableName << ' ' 
		<< std::setw(7) << counts[itab] << ' '
		<< std::setw(7) << tables[itab].fTargetEventCount << ' '
		<< std::setw(4) << workunits << std::endl;
    }
}
