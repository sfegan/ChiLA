//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierFactory.hpp

  Factory for creating event suppliers

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       19/07/2006
  \note
*/

#include <VSDBParameterTable.hpp>

#include "VSEventSupplierLaser.hpp"
#include "VSEventSupplierPedestal.hpp"
#include "VSEventSupplierOneTableDB.hpp"
#include "VSEventSupplierOneTableHDF.hpp"
#include "VSEventSupplierManyTables.hpp"

#include "VSEventSupplierFactory.hpp"

using namespace VERITAS;

VSEventSupplierFactory::
VSEventSupplierFactory(VSSimDB* sim_db, RandomNumbers* rng,
                       const std::string& dir):
  fSimDB(sim_db), fRNG(rng), fOverrideDataDirectory(dir), 
  fMaxBufferSize(16*1024), fCheckHDFEventCount(),
  fMinEnergyGeV(), fMaxEnergyGeV(),
  fTargetZenithRad(), fTargetAzimuthRad(),
  fHasTargetZenith(), fHasTargetAzimuth()
{
  // nothing to see here
}

VSEventSupplierFactory::~VSEventSupplierFactory()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// EMPTY (PEDESTAL) SUPPLIER
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newEmptySupplier(unsigned nevents, unsigned nscope, unsigned nchan)
{
  return new VSEventSupplierPedestal(nevents, nscope, nchan);
}

// ----------------------------------------------------------------------------
// LASER SUPPLIER
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newLaserSupplier(unsigned nevents,
		 double mean_pe_per_pixel, double mean_pe_per_pixel_dev, 
		 double mean_pulse_fwhm_ns,
		 unsigned nscope, unsigned nchan)
{
  assert(fRNG);
  return new VSEventSupplierLaser(fRNG, nevents, 
				  mean_pe_per_pixel, mean_pe_per_pixel_dev,
				  mean_pulse_fwhm_ns,
				  nscope, nchan);
}

// ----------------------------------------------------------------------------
// ONE TABLE SUPPLIER
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newOneTableSupplier(const std::string& table_name)
{
  assert(fSimDB);

  VSDBParameterTable db_param(fSimDB->db());

  VSDBParameterSet parameters;
  db_param.retrieveParameterSet("DataStorage", parameters);

  if((!parameters.empty())
     &&(parameters.find("mode") != parameters.end())
     &&(parameters["mode"] == "hdf5"))
    {
      std::string data_directory = parameters["directory"];
      if(!fOverrideDataDirectory.empty())
	data_directory = fOverrideDataDirectory;
      return 
	new VSEventSupplierOneTableHDF(table_name, fSimDB, data_directory,
				       fCheckHDFEventCount);
    }

  return new VSEventSupplierOneTableDB(table_name, fSimDB, fMaxBufferSize);
}      

VSEventSupplier* VSEventSupplierFactory::
newOneTableSupplier(const std::string& table_name,
		    uint32_t nevent_to_supply, double ievent_to_start_fraction)
{
  assert(fSimDB);

  VSDBParameterTable db_param(fSimDB->db());

  VSDBParameterSet parameters;
  db_param.retrieveParameterSet("DataStorage", parameters);

  if((!parameters.empty())
     &&(parameters.find("mode") != parameters.end())
     &&(parameters["mode"] == "hdf5"))
    {
      std::string data_directory = parameters["directory"];
      if(!fOverrideDataDirectory.empty())
	data_directory = fOverrideDataDirectory;
      return 
	new VSEventSupplierOneTableHDF(table_name, fSimDB, 
				       nevent_to_supply,
				       ievent_to_start_fraction,
				       data_directory,
				       fCheckHDFEventCount);
    }

  return new VSEventSupplierOneTableDB(table_name, fSimDB, 
				       nevent_to_supply,
				       ievent_to_start_fraction,
				       fMaxBufferSize);
}

// ----------------------------------------------------------------------------
// ALL EVENTS SUPPLIER
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newAllEventsSupplier(bool randomize_if_available)
{
  std::vector<VSSimDBTableParam> tables = fSimDB->getAllDataTables();
  if(tables.empty())return 0;

  std::list<VSEventSupplierManyTables::TableName> tables_list;
  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	tables.begin(); itr != tables.end(); ++itr)
    tables_list.push_back(itr->fTableName);

  return new VSEventSupplierManyTables(this, tables_list,
				       randomize_if_available?fRNG:0);
}

// ----------------------------------------------------------------------------
// ALL EVENTS WEIGHTED TO HIGH ENERGIES
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newHEWeightedSupplier(uint32_t nevent_to_supply,
		      bool randomize_if_available)
{
  std::vector<VSSimDBTableParam> tables = fSimDB->getAllDataTables();
  if(tables.empty())return 0;

  selectByEnergy(tables);
  selectByZenith(tables);
  selectByAzimuth(tables);

  std::vector< std::pair< std::string, uint32_t > > tables_count;
  findHEWeighting(nevent_to_supply,
		  randomize_if_available,
		  tables,
		  tables_count);

  std::list<VSEventSupplierManyTables::TableNameAndCount> tables_list;
  for(std::vector< std::pair< std::string, uint32_t > >::iterator itr =
	tables_count.begin(); itr != tables_count.end(); ++itr)
    tables_list.push_back(*itr);

  return new VSEventSupplierManyTables(this, tables_list,
				       randomize_if_available?fRNG:0);
}

// ----------------------------------------------------------------------------
// MONOENERGETIC SPECTRUM SUPPLIER
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newMonoenergeticSupplier(double target_energy_gev,
			 uint32_t nevent_to_supply,
			 bool randomize_if_available)
{
  std::vector<VSSimDBTableParam> tables = fSimDB->getAllDataTables();
  if(tables.empty())return 0;

  selectByZenith(tables);
  selectByAzimuth(tables);

  double closest_energy_gev = findClosestEnergy(target_energy_gev,tables);

  for(std::vector<VSSimDBTableParam>::iterator itr =
	tables.begin(); itr != tables.end(); )
    {
      if(fabs(itr->fEnergyGeV-closest_energy_gev)/itr->fEnergyGeV < 1E-4)
	itr++;
      else
	tables.erase(itr);
    }

  std::list<VSEventSupplierManyTables::TableName> tables_list;
  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	tables.begin(); itr != tables.end(); ++itr)
    tables_list.push_back(itr->fTableName);

  return new VSEventSupplierManyTables(this, tables_list,
				       randomize_if_available?fRNG:0);
}

// ----------------------------------------------------------------------------
// POWER-LAW SPECTRUM SUPPLIER
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newPLSpectrumSupplier(double spectral_index,
		      uint32_t nevent_to_supply,
		      bool randomize_if_available)
{
  std::vector<VSSimDBTableParam> tables = fSimDB->getAllDataTables();
  if(tables.empty())return 0;

  selectByEnergy(tables);
  selectByZenith(tables);
  selectByAzimuth(tables);

  std::vector< std::pair< std::string, uint32_t > > tables_count;
  findPLWeighting(spectral_index,
		  nevent_to_supply,
		  randomize_if_available,
		  tables,
		  tables_count);

  std::list<VSEventSupplierManyTables::TableNameAndCount> tables_list;
  for(std::vector< std::pair< std::string, uint32_t > >::iterator itr =
	tables_count.begin(); itr != tables_count.end(); ++itr)
    tables_list.push_back(*itr);
  
  return new VSEventSupplierManyTables(this, tables_list,
				       randomize_if_available?fRNG:0);
}

// ----------------------------------------------------------------------------
// LEGACY METHODS
// ----------------------------------------------------------------------------

VSEventSupplier* VSEventSupplierFactory::
newPLSpectrumClosestZenithSupplier(double spectral_index,
				   double min_energy_gev,
				   double max_energy_gev,
				   double zn_rad,
				   uint32_t nevent_to_supply,
				   bool randomize_if_available)
{
  setEnergyRange(min_energy_gev,max_energy_gev);
  setTargetZenith(zn_rad);
  return newPLSpectrumSupplier(spectral_index,nevent_to_supply,
			       randomize_if_available);
}

VSEventSupplier* VSEventSupplierFactory::
newPLSpectrumClosestPointingSupplier(double spectral_index,
				     double min_energy_gev,
				     double max_energy_gev,
				     double zn_rad,
				     double az_rad,
				     uint32_t nevent_to_supply,
				     bool randomize_if_available)
{
  setEnergyRange(min_energy_gev,max_energy_gev);
  setTargetZenith(zn_rad);
  setTargetAzimuth(az_rad);
  return newPLSpectrumSupplier(spectral_index,nevent_to_supply,
			       randomize_if_available);
}

VSEventSupplier* VSEventSupplierFactory::
newAllTablesClosestZenithSupplier(double min_energy_gev,
				  double max_energy_gev,
				  double zn_rad, 
				  uint32_t nevent_to_supply,
				  bool randomize_if_available)
{
  setEnergyRange(min_energy_gev,max_energy_gev);
  setTargetZenith(zn_rad);
  return newHEWeightedSupplier(nevent_to_supply,randomize_if_available);
}

VSEventSupplier* VSEventSupplierFactory::
newAllTablesClosestPointingSupplier(double min_energy_gev,
				    double max_energy_gev,
				    double zn_rad, 
				    double az_rad,
				    uint32_t nevent_to_supply,
				    bool randomize_if_available)
{
  setEnergyRange(min_energy_gev,max_energy_gev);
  setTargetZenith(zn_rad);
  setTargetAzimuth(az_rad);
  return newHEWeightedSupplier(nevent_to_supply,randomize_if_available);
}

// ----------------------------------------------------------------------------
// UTILITY FUNCTIONS
// ----------------------------------------------------------------------------

double VSEventSupplierFactory::
angleDistance(double target_angle, 
	      double angle_range_min, double angle_range_max)
{
  if(((angle_range_min <= angle_range_max)
      &&((target_angle >= angle_range_min)
	 &&(target_angle <= angle_range_max)))||
     ((angle_range_min > angle_range_max)
      &&((target_angle >= angle_range_min)
	 ||(target_angle <= angle_range_max))))return 0;

  double distance_min = fabs(target_angle-angle_range_min);
  double distance_max = fabs(target_angle-angle_range_max);
  if(distance_min > M_PI)distance_min = 2*M_PI - distance_min;
  if(distance_max > M_PI)distance_max = 2*M_PI - distance_max;
  return distance_min<distance_max?distance_min:distance_max;
}

double VSEventSupplierFactory::
findClosestEnergy(double target_energy_gev,
		  const std::vector<VSSimDBTableParam>& tables)
{
  double closest_energy_gev = 0;
  double de_min = 0;
  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	tables.begin(); itr != tables.end(); itr++)
    {
      double de = fabs(itr->fEnergyGeV-target_energy_gev);
      
      if(itr == tables.begin() || de < de_min)
	{
	  closest_energy_gev = itr->fEnergyGeV;
	  de_min = de;
	}
    }

  return closest_energy_gev;
}

double VSEventSupplierFactory::
findMinZenith(double target_angle,const std::vector<VSSimDBTableParam>& tables)
{
  double closest_zn_rad = 0;
  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	tables.begin(); itr != tables.end(); itr++)
    {
      double d = 
	angleDistance(target_angle, itr->fZenithMinRad,itr->fZenithMaxRad);
      
      if(itr == tables.begin() || d < closest_zn_rad)closest_zn_rad = d;
    }

  return closest_zn_rad;
}

double VSEventSupplierFactory::
findMinAzimuth(double target_angle, double zn_angle,
	       const std::vector<VSSimDBTableParam>& tables)
{
  double closest_zn_rad = findMinZenith(zn_angle,tables);
  double closest_az_rad = 2*M_PI;

  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	tables.begin(); itr != tables.end();itr++)
    {
      double d_az = 
	angleDistance(target_angle,itr->fAzimuthMinRad,itr->fAzimuthMaxRad);
      
      double d_zn = 
	angleDistance(zn_angle,itr->fZenithMinRad,itr->fZenithMaxRad);

      if(d_zn != closest_zn_rad) continue;

      if(itr == tables.begin() || d_az < closest_az_rad)
	closest_az_rad = d_az;
    }

  return closest_az_rad;
}

void VSEventSupplierFactory::
findPLWeighting(double spectral_index,
		uint32_t nevent_to_supply,
		bool randomize_if_available,
		std::vector<VSSimDBTableParam> &table_params,
		std::vector< std::pair< std::string, uint32_t > > &tables)
{ 
  std::map<double,uint32_t> max_event_count;
  std::map<double,double> energy_effarea;

  double sum = 0;
  for(std::vector<VSSimDBTableParam>::iterator itr =
	table_params.begin(); itr != table_params.end(); itr++)
    {     
      uint32_t table_count = fSimDB->getMaxEventNumByName(itr->fTableName);
      double effarea = M_PI*std::pow(itr->fSamplingRadiusM,2);
      
      if(table_count == 0)
	std::cerr << "In VSEventSupplierFactory::findPLWeighting(): "
		  << "Table " << itr->fTableName 
		  << " has 0 events." << std::endl;

      if(max_event_count.find(itr->fEnergyGeV) == max_event_count.end())
	sum+=std::pow( itr->fEnergyGeV, (float)(spectral_index+1.0) )*effarea;
      
      if(max_event_count.find(itr->fEnergyGeV) == max_event_count.end())
	max_event_count[itr->fEnergyGeV] = table_count;
      else if(table_count < max_event_count[itr->fEnergyGeV])
	max_event_count[itr->fEnergyGeV] = table_count;

      energy_effarea[itr->fEnergyGeV] = effarea;
    }

  // --------------------------------------------------------------------------
  // Find the largest possible value of nevent_to_supply given the
  // number of events which can be provided from each energy bin.
  // --------------------------------------------------------------------------
  for(std::map<double,uint32_t>::iterator itr = max_event_count.begin();
      itr != max_event_count.end(); ++itr) 
    {
      double wegy = std::pow(itr->first,spectral_index+1 );
      double area = energy_effarea[itr->first];
      uint32_t nevents = lround(sum*itr->second/(wegy*area));

      if(nevents < nevent_to_supply) nevent_to_supply = nevents;
    }

  std::map<double,uint32_t> event_count; 
  double nevent = nevent_to_supply;
  for(std::map<double,uint32_t>::reverse_iterator itr = 
	max_event_count.rbegin(); itr != max_event_count.rend(); ++itr) 
    {
      double wegy = std::pow(itr->first,spectral_index+1 );
      double area = energy_effarea[itr->first];

      event_count[itr->first] = 
	std::min(max_event_count[itr->first],
		 (uint32_t)lround(nevent/sum*wegy*area));
      nevent -= event_count[itr->first];
      sum -= wegy*area;
    }

  // --------------------------------------------------------------------------
  // Distribute events evenly among all tables in each energy bin.
  // --------------------------------------------------------------------------
  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	table_params.begin(); itr != table_params.end(); itr++)
    {
      uint32_t nevents = event_count[itr->fEnergyGeV];
      tables.push_back(std::make_pair(itr->fTableName,nevents));
    }
}

void VSEventSupplierFactory::
findHEWeighting(uint32_t nevent_to_supply,
		bool randomize_if_available,
		std::vector<VSSimDBTableParam> &tables,
		std::vector< std::pair< std::string, 
		uint32_t > > &tables_count)
{
  std::map<double,unsigned> energy_event_count;
  std::map<double,unsigned> energy_table_count;
  
  // --------------------------------------------------------------------------
  // Find the number of events and tables at each energy
  // --------------------------------------------------------------------------
  for(std::vector<VSSimDBTableParam>::iterator itr =
	tables.begin(); itr != tables.end(); ++itr)
    {
      uint32_t event_count = fSimDB->getMaxEventNumByName(itr->fTableName);
      energy_event_count[itr->fEnergyGeV] += event_count;
      energy_table_count[itr->fEnergyGeV]++;
    }      

  unsigned nenergy = energy_event_count.size();
  unsigned nevent_per_energy = 
    (unsigned)ceil((double)nevent_to_supply/(double)nenergy);

  for(std::map<double,unsigned>::reverse_iterator itr =
	energy_event_count.rbegin(); itr != energy_event_count.rend(); ++itr)
    {
      if(nevent_per_energy <= itr->second)
	{
	  itr->second = nevent_per_energy;
	  nevent_to_supply -= nevent_per_energy;
	  nenergy--;
	}
      else
	{
	  nevent_to_supply -= itr->second;
	  nenergy--;
	  nevent_per_energy = 
	    (unsigned)ceil((double)nevent_to_supply/(double)nenergy);
	}
    }

  for(std::vector<VSSimDBTableParam>::const_iterator itr =
	tables.begin(); itr != tables.end(); itr++)
    {
      unsigned nevent = energy_event_count[itr->fEnergyGeV]/
	energy_table_count[itr->fEnergyGeV];
      uint32_t nevent_table = fSimDB->getMaxEventNumByName(itr->fTableName);

      if(nevent > nevent_table) nevent = nevent_table;

      energy_event_count[itr->fEnergyGeV] -= nevent;
      energy_table_count[itr->fEnergyGeV]--;
      tables_count.push_back(std::make_pair(itr->fTableName,nevent));
    }
}

void VSEventSupplierFactory::
selectByEnergy(std::vector<VSSimDBTableParam>& tables)
{
  if(fMinEnergyGeV == 0 && fMaxEnergyGeV == 0) return;

  for(std::vector<VSSimDBTableParam>::iterator itr =
	tables.begin(); itr != tables.end(); )
    {
      if(itr->fEnergyGeV < fMinEnergyGeV ||
	 (fMaxEnergyGeV > 0 && itr->fEnergyGeV > fMaxEnergyGeV))
	tables.erase(itr);
      else
	++itr;
    }
}

void VSEventSupplierFactory::
selectByZenith(std::vector<VSSimDBTableParam>& tables)
{
  if(!fHasTargetZenith) return;

  double closest_zn_rad = findMinZenith(fTargetZenithRad, tables);

  for(std::vector<VSSimDBTableParam>::iterator itr =
	tables.begin(); itr != tables.end(); )
    {
      if(angleDistance(fTargetZenithRad, itr->fZenithMinRad, 
		       itr->fZenithMaxRad) - closest_zn_rad < 1E-4)
	++itr;
      else
	tables.erase(itr);
    }
}

void VSEventSupplierFactory::
selectByAzimuth(std::vector<VSSimDBTableParam>& tables)
{
  if(!fHasTargetAzimuth) return;

  double closest_az_rad = 
    findMinAzimuth(fTargetAzimuthRad, fTargetZenithRad, tables);

  for(std::vector<VSSimDBTableParam>::iterator itr =
	tables.begin(); itr != tables.end(); )
    {
      if(angleDistance(fTargetAzimuthRad, itr->fAzimuthMinRad, 
		       itr->fAzimuthMaxRad) - closest_az_rad < 1E-4)
	++itr;
      else
	tables.erase(itr);
    }
}
