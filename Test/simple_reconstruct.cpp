//-*-mode:c++; mode:font-lock;-*-

/*! \file simple_reconstruct.cpp
  Simple reconstruction program to illustrate use of various DB components

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       09/16/2005
*/

#include <list>
#include <iostream>
#include <iterator>

#include <VSOptions.hpp>
#include <VSDatabase.hpp>
#include <VSDBFactory.hpp>
#include <VSDBParameterTable.hpp>
#include <VSSimDB.hpp>
#include <VSSimDBTables.hpp>
#include <VSOTelescopeArray.hpp>
#include <VSReconstruction.hpp>

#include <VSTestTimer.hpp>

using namespace VERITAS;
using namespace Physics;

struct ChannelInfo
{
  ChannelInfo(): valid(false), vec() { } 
  bool valid;
  Vec3D vec;
};

typedef std::list<VSReconstruction::Chi2MapElement> Chi2Map;

int main(int argc, char** argv)
{
  std::string progname(*argv);

  VSOptions options(argc,argv,true);

  // --------------------------------------------------------------------------
  // Options
  // --------------------------------------------------------------------------

  VSDBFactory::configure(&options);

#if 0
  std::ostringstream stream;
  stream << "/tmp/rng_state_uid" << getuid() << ".dat";
  std::string rngstatefile(stream.str());
  options.findWithValue("rng_state_file",rngstatefile,
			"Set random number generator state file");
#endif

  bool print_usage = false;
  int exit_code = EXIT_SUCCESS;

  if(options.find("h","Print this message") != VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help", "Print this message") != VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool single_event_mode = false;
  unsigned event_no = 0;
  if(options.findWithValue("ev",event_no,
			   "Process only one event, writing details of chi^2 "
			   "grid search. Number of event to process should be "
			   "given as an argument to this option") 
     != VSOptions::FS_NOT_FOUND)single_event_mode = true;

  unsigned max_events = 0;
  options.findWithValue("max_events",max_events,
			"Maximum number of events to process, or zero to "
			"process all events");
  
  double plate_scale = 1.0;
  options.findWithValue("platescale",plate_scale,
			"Set the plate scale correction "
			"(1.0 gives no correction)");

  std::vector<unsigned> scope_suppress_id;
  options.findWithValue("scope_suppress",scope_suppress_id,
			"Comma separated list of telescope numbers to "
			"suppress in analysis");

  double fov = 3.5;
  options.findWithValue("fov",fov,
			"Set the camera field of view");

  unsigned l1_pe = 6;
  options.findWithValue("l1_pe",l1_pe,
			"Set the L1 trigger level in PE");

  unsigned l2_mult = 3;
  options.findWithValue("l2_mult",l2_mult,
			"Set the L2 trigger multiplicity in channels");

  bool l2_force = false;
  if(options.find("l2_force",
		  "Force L2 trigger on all telescopes when L3 "
		  "trigger occurs") != VSOptions::FS_NOT_FOUND)
    l2_force = true;

  unsigned l3_mult = 2;
  options.findWithValue("l3_mult",l3_mult,
			"Set the L3 trigger multiplicity in telescopes");

  unsigned cleaning_pe = unsigned(floor(0.1878*10*3.0));
  options.findWithValue("cleaning_pe",cleaning_pe,
			"Set the cleaning level in PE");

  // --------------------------------------------------------------------------
  // Get the database and table names from the command line arguments
  // --------------------------------------------------------------------------

  argc--,argv++;

  if (argc!=2)
    {
      print_usage=true;
      exit_code=EXIT_FAILURE;
    }

  if(print_usage)
    {
      std::cerr << "Usage: " << progname 
		<< " [options] database table_name" 
		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_code;
    }

  std::string database(*argv);
  argc--; argv++;

  std::string tablename(*argv);
  argc--; argv++;

  // --------------------------------------------------------------------------
  // Create the database
  // --------------------------------------------------------------------------

  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  if(db->useDatabase(database) < 0)
    {
      std::cerr << "Could not connect to database " << database << std::endl;
      return EXIT_FAILURE;
    };

  VSSimDB db_sim(db);
  VSDBParameterTable db_param(db);

  // --------------------------------------------------------------------------
  // Get the table parameters from the DB
  // --------------------------------------------------------------------------

  VSSimDBTableParam* param;
  param = db_sim.getDataTableByName(tablename);
  if(!param)
    {
      std::cerr << "Could not get information for table " << tablename 
		<< std::endl;
      return EXIT_FAILURE;      
    }

  // --------------------------------------------------------------------------
  // Get the definition of the array from the database
  // --------------------------------------------------------------------------

  VSOTelescopeArray array;
  array.readFromDatabase(db, param->fOpticsID);

  std::vector<bool> scope_suppress;
  std::vector<std::vector<ChannelInfo> > channel_info;
  std::vector<unsigned> max_channel;

  scope_suppress.resize(array.numTelescopes());
  channel_info.resize(array.numTelescopes());
  max_channel.resize(array.numTelescopes());
  
  for(std::vector<unsigned>::const_iterator iid=scope_suppress_id.begin();
      iid!=scope_suppress_id.end(); iid++)
    if(*iid<scope_suppress.size())scope_suppress[*iid]=true;

  for(unsigned iscope=0;iscope<array.numTelescopes();iscope++)
    {
      unsigned nvalid=0;
      unsigned npix=array.telescope(iscope)->numPixels();

      channel_info[iscope].resize(npix);
      for(unsigned ipix=0;ipix<npix; ipix++)
	{
	  channel_info[iscope][ipix].vec = 
	    array.telescope(iscope)->pixel(ipix)->incomingSkyVectorAtZenith(plate_scale);
	  
	  double theta = acos(channel_info[iscope][ipix].vec *
			      Vec3D(0,0,-1));
	  
	  if(theta <= fov/2/180*M_PI)
	    {
	      max_channel[iscope]=ipix;
	      channel_info[iscope][ipix].valid=true;
	      nvalid++;
	    }
	}
      
      if(npix==0)
	{
	  max_channel[iscope] = 0xFFFFFFFFU;
	}

      std::cerr << progname << ": scope " << iscope << " has " << nvalid 
		<< " pixels (max id:" << max_channel[iscope] << ')'
		<< std::endl;
    }

  // --------------------------------------------------------------------------
  // Set up the reconstruction
  // --------------------------------------------------------------------------

  VSReconstruction::ScopeSpecs scopes;
  for(unsigned iscope=0;iscope<array.numTelescopes();iscope++)
    {
      VSCORSIKATelescopeSpec scope;
      scope.scope_num = iscope;
      scope.scope_x   = array.telescope(iscope)->pos().x/100;
      scope.scope_y   = array.telescope(iscope)->pos().y/100;
      scope.scope_z   = array.telescope(iscope)->pos().z/100;
      scope.scope_r   = array.telescope(iscope)->aperture()/100;
      scopes.push_back(scope);
    }  
  VSReconstruction reconstruction(scopes);

  // --------------------------------------------------------------------------
  // Get a list of all (or single, if requested) complete events
  // --------------------------------------------------------------------------
 
  std::string tablename_ev = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    +tablename
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

  std::string tablename_tel = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    +tablename
    + std::string(VSIMDB_TABLE_POSTFIX_SCOPES);

  std::string tablename_pe = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    +tablename
    + std::string(VSIMDB_TABLE_POSTFIX_PES);

  std::string ev_condition = "EventComplete=1";
  if(single_event_mode)
    ev_condition += " AND EventID=" + VSDataConverter::toString(event_no);
  else if(max_events!=0)
    ev_condition += " AND EventID<" + VSDataConverter::toString(max_events);    
  VSDBStatement* stmt = 
    db->createSelectQuery(tablename_ev, ev_condition, "*",
			  VSDatabase::FLAG_NO_BUFFER);

  if(!stmt)
    {
      std::cerr << progname
		<< ": could no query table: " << tablename_ev << std::endl;
      return EXIT_FAILURE;
    }
  
  VSSimDBEventData event;
  std::list<VSSimDBEventData> events;

  stmt->bindToResult(event.fEventID);
  stmt->bindToResult(event.fTargetZenithRad);
  stmt->bindToResult(event.fTargetAzimuthRad);
  stmt->bindToResult(event.fPrimaryZenithRad);
  stmt->bindToResult(event.fPrimaryAzimuthRad);
  stmt->bindToResult(event.fPrimaryCoreEastM);
  stmt->bindToResult(event.fPrimaryCoreNorthM);
  stmt->bindToResult(event.fPrimaryCoreUpASLM);
  stmt->bindToResult(event.fNumHitScopes);
  stmt->bindToResult(event.fEventComplete);

  if(stmt->execute() < 0)
    {
      std::cerr << progname
		<< ": error executing query on: " << tablename_ev << std::endl;
      return EXIT_FAILURE;
    }

  VSTestTimer<VSTimerCoreIA32> timer_query;
  VSTestTimer<VSTimerCoreIA32> timer_reconstruction;

  timer_query.start();
  while(stmt->retrieveNextRow())events.push_back(event);
  timer_query.stop();
  delete stmt;

  std::cerr << progname << ": loaded " << events.size() 
	    << " event headers in "  << timer_query << " seconds" << std::endl;

  // --------------------------------------------------------------------------
  // Create the database queries and bind variables
  // --------------------------------------------------------------------------

  VSDBStatement* stmt_tel = 
    db->createSelectQuery(tablename_tel, "EventID=?", 
			  "ScopeID,ScopeZenithRad,ScopeAzimuthRad");

  if(!stmt_tel)
    {
      std::cerr << progname
		<< ": could no query table: " << tablename_tel << std::endl;
      return EXIT_FAILURE;
    }

  VSDBStatement* stmt_pix = 
    db->createQuery(std::string("SELECT PixelID,COUNT(PixelID) FROM ")
		    + tablename_pe
		    + std::string(" WHERE EventID=? AND ScopeID=? "
				  "AND PixelID<=? GROUP BY PixelID"));

  if(!stmt_pix)
    {
      std::cerr << progname
		<< ": could no query table: " << tablename_pe << std::endl;
      return EXIT_FAILURE;
    }

  uint32_t event_id;
  uint16_t scope_id;

  float scope_zn;
  float scope_az;

  uint32_t max_pixel_id;
  uint32_t pixel_id;
  uint32_t pe_count;

  stmt_tel->bindToParam(event_id);
  stmt_tel->bindToResult(scope_id);
  stmt_tel->bindToResult(scope_zn);
  stmt_tel->bindToResult(scope_az);

  stmt_pix->bindToParam(event_id);
  stmt_pix->bindToParam(scope_id);
  stmt_pix->bindToParam(max_pixel_id);
  stmt_pix->bindToResult(pixel_id);
  stmt_pix->bindToResult(pe_count);
  
  // --------------------------------------------------------------------------
  // Query database and reconstruct events
  // --------------------------------------------------------------------------

  std::vector<std::vector<std::pair<uint32_t, uint32_t> > > data;
  data.resize(array.numTelescopes());
  for(unsigned iscope=0;iscope<array.numTelescopes();iscope++)
    data[iscope].reserve(max_channel[scope_id]+1);

  std::vector<double> scopes_zn;
  std::vector<double> scopes_az;
  scopes_az.resize(array.numTelescopes());
  scopes_zn.resize(array.numTelescopes());

  std::vector<VSReconstruction::GridSpecs> grid;
  grid.push_back(VSReconstruction::GridSpecs(fov/2.0,0.100));
  grid.push_back(VSReconstruction::GridSpecs(fov/20.0,0.010));
  grid.push_back(VSReconstruction::GridSpecs(fov/200.0,0.001));

  timer_query.reset();
  unsigned ntrig_ev = 0;
  for(std::list<VSSimDBEventData>::const_iterator ievent = events.begin();
      ievent!=events.end(); ievent++)
    {
      unsigned ntrig_l3 = 0;

      std::vector<bool> scopes_l2;
      scopes_l2.resize(array.numTelescopes());

      event_id = ievent->fEventID;

      for(unsigned iscope=0;iscope<array.numTelescopes();iscope++)
	data[iscope].clear();

      timer_query.start();
      if(stmt_tel->execute()<0)
	{
	  std::cerr << progname
		    << ": error executing query on: " << tablename_tel 
		    << std::endl;
	  return EXIT_FAILURE;
	}
      
      while(stmt_tel->retrieveNextRow())
	{
	  if(scope_suppress[scope_id])continue;

	  uint32_t ntrig_l2 = 0;
	  max_pixel_id = max_channel[scope_id];

	  scopes_zn[scope_id] = scope_zn;
	  scopes_az[scope_id] = scope_az;
	  
	  if(stmt_pix->execute()<0)
	    {
	      std::cerr << progname
			<< ": error executing query on: " << tablename_pe 
			<< std::endl;
	      return EXIT_FAILURE;
	    }

	  while(stmt_pix->retrieveNextRow())
	    {
	      if((pixel_id<channel_info[scope_id].size())&&
		 (!channel_info[scope_id][pixel_id].valid))continue;
	      if(pe_count>=l1_pe)ntrig_l2++;
	      if(pe_count>=cleaning_pe)
		data[scope_id].
		  push_back(std::make_pair<uint32_t,uint32_t>(pixel_id,
							      pe_count));
	    }

	  if(ntrig_l2>=l2_mult)
	    {
	      ntrig_l3++;
	      scopes_l2[scope_id]=true;
	    }
	}

      timer_query.stop();

#ifndef NO_OUTPUT
      if(!single_event_mode)
	std::cout 
	  << ievent->fEventID
	  << ' ' << ievent->fTargetZenithRad
	  << ' ' << ievent->fTargetAzimuthRad
	  << ' ' << ievent->fPrimaryZenithRad
	  << ' ' << ievent->fPrimaryAzimuthRad
	  << ' ' << ievent->fPrimaryCoreEastM
	  << ' ' << ievent->fPrimaryCoreNorthM
	  << ' ' << ievent->fPrimaryCoreUpASLM
	  << ' ' << ntrig_l3;
#endif
	    
      if(ntrig_l3>=l3_mult)
	{
#ifndef NO_OUTPUT
	  if(!single_event_mode)
	    std::cout << ' ' << 1;
#endif

	  // ******************************************************************
	  // R E C O N S T R U C T I O N
	  // ******************************************************************

	  timer_reconstruction.start();
	 
	  Vec3D axis(0,0,-1);
	  axis.Rotate(Vec3D(-ievent->fTargetZenithRad,0,0));
	  axis.Rotate(Vec3D(0,0,-ievent->fTargetAzimuthRad));
	  VSReconstruction::Ray optical_axis(axis.x, axis.y);

	  unsigned nray = 0;

	  VSReconstruction::AllScopeRays array_rays;
	  array_rays.resize(array.numTelescopes());

	  for(unsigned iscope=0;iscope<array.numTelescopes();iscope++)
	    {
	      unsigned npix = data[iscope].size();
	      if((npix==0)||((l2_force==false)&&(scopes_l2[iscope]==false)))
		 continue;

	      array_rays[iscope].reserve(npix);

	      Vec3D scope_rotation(-scopes_zn[iscope],0,0);
	      scope_rotation &= Vec3D(0,0,-scopes_az[iscope]);

	      for(unsigned ipix=0;ipix<npix;ipix++)
		{
		  unsigned ipixel = data[iscope][ipix].first;
		  unsigned npe = data[iscope][ipix].second;
		  Vec3D r;

		  if(ipixel<channel_info[iscope].size())
		    r=channel_info[iscope][ipixel].vec;
		  else
		    {
		      const VSOTelescope* scope = array.telescope(iscope);
		      int hexid = ipixel;
		      Vec3D pos;
		      nh_to_xy(&hexid, &pos.x, &pos.z); 

		      if(scope->pixelParity())pos.x=-pos.x;
		      pos.x *= scope->pixelSpacing();
		      pos.z *= scope->pixelSpacing();
		      pos.y = 0;
		      pos.Rotate(scope->focalPlaneRotion());
		      
		      VSOPixel* pixel = 
			new VSOPixel(scope, hexid, hexid, false, pos);
		      
		      r = pixel->incomingSkyVectorAtZenith(plate_scale);
		    }

		  r.Rotate(scope_rotation);
		  array_rays[iscope]
		    .push_back(VSReconstruction::Ray(r.x, r.y, npe));
		  nray++;
		}
	    }
	  
	  VSReconstruction::MinimizationResults res;

	  if(single_event_mode)
	    {
	      Chi2Map chi2map;
	      reconstruction.gridSearch(array_rays, optical_axis, grid, res,
				std::back_insert_iterator<Chi2Map>(chi2map));

	      for(Chi2Map::iterator i = chi2map.begin(); 
		  i!=chi2map.end(); i++)
		{
		  Vec3D p(i->ec_x, i->ec_y, i->ec_z);
		  p.Rotate(Vec3D(0,0,ievent->fTargetAzimuthRad));
		  p.Rotate(Vec3D(ievent->fTargetZenithRad,0,0));
		  double x_fp = p.y/p.z/M_PI*180;
		  double y_fp = p.x/p.z/M_PI*180;
		  std::cout << x_fp << ' ' << y_fp << ' ' << i->chi2 << ' '
		    << i->ec_x << ' ' <<  i->ec_y << ' ' <<  i->ec_z << ' '
		    << i->rc_x << ' ' <<  i->rc_y << ' ' <<  i->rc_z 
			    << std::endl;
		}
	    }
	  else
	    {
	      reconstruction.gridSearch(array_rays, optical_axis, grid, res);
	    }

	  timer_reconstruction.stop();

#ifndef NO_OUTPUT
	  if(!single_event_mode)
	    std::cout
	      << ' ' << nray
	      << ' ' << res.chi2
	      << ' ' << res.weight
	      << ' ' << atan2(sqrt(res.ec_x*res.ec_x+res.ec_y*res.ec_y),-res.ec_z)*180/M_PI
	      << ' ' << atan2(res.ec_x,res.ec_y)*180/M_PI
	      << ' ' << res.rc_x
	      << ' ' << res.rc_y
	      << ' ' << res.rc_z
	      << ' ' << res.chi2_curv_l
	      << ' ' << res.chi2_curv_w;
#endif

	  ntrig_ev++;
	}
      else 
	{
#ifndef NO_OUTPUT
	  if(!single_event_mode)
	    std::cout << ' ' << 0;
#endif
	}
      
#ifndef NO_OUTPUT
      if(!single_event_mode)
	std::cout << std::endl;
#endif      
    }
  
  std::cerr << progname << ": " << events.size() << " events read in "
	    << timer_query << " seconds" << std::endl;

  std::cerr << progname << ": " << ntrig_ev 
	    << " triggering events reconstructed in " 
	    << timer_reconstruction << " seconds" << std::endl;

}
