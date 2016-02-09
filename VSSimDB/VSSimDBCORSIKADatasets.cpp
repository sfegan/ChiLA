//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimDBCORSIKADatasets.cpp
  Database access to wavelength/value tables

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       08/19/2005
*/

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cassert>

#include <VSSimDBCORSIKADatasets.hpp>
#include <VSDBParameterTable.hpp>

using namespace VERITAS;

// ----------------------------------------------------------------------------
// VSSimDBWavelengthDataset
// ----------------------------------------------------------------------------

bool VSSimDBWavelengthDataset::
writeToCORSIKA(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  if(!stream)return false;

  stream << comment << std::endl;

  unsigned wavelength=180;
  unsigned count=0;
  while(wavelength<=700)
    {
      VSSimDBWavelengthData::const_iterator idatum = data.find(wavelength);
      assert(idatum != data.end());

      if(count==0){ if(wavelength!=180)stream << std::endl; count=1; }
      else { stream << ' '; count++; if(count==8)count=0; }

      stream << std::fixed << std::setprecision(3) << idatum->second;

      wavelength += 5;
    }
  stream << std::endl;
  return true;
}

VSSimDBWavelengthDataset* 
VSSimDBWavelengthDataset::createFromCORSIKA(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(!stream)return 0;
  
  VSSimDBWavelengthDataset* dataset = new VSSimDBWavelengthDataset;
  if(dataset==0)return 0;

  std::getline(stream,dataset->comment);
  
  unsigned wavelength=180;
  while(wavelength<=700)
    {
      float value;
      stream >> value;
      dataset->data[wavelength]=value;
      wavelength += 5;
    }    
  return dataset;
}

// ----------------------------------------------------------------------------
// VSSimDBWavelengthAltitudeDataset
// ----------------------------------------------------------------------------

bool VSSimDBWavelengthAltitudeDataset::
writeToCORSIKA(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  if(!stream)return false;
  
  stream << comment << std::endl;
  unsigned wavelength=180;
  while(wavelength<=700)
    {
      VSSimDBWavelengthAltitudeData::const_iterator iwlset = 
	data.find(wavelength);
      assert(iwlset != data.end());
      
      stream << std::fixed << std::setw(4) << wavelength << std::endl;
      unsigned altitude=0;
      unsigned count=0;
      while(altitude<=50)
	{
	  VSSimDBAltitudeData::const_iterator idatum = 
	    iwlset->second.find(altitude);
	  assert(idatum != iwlset->second.end());

	  if(count==0){ if(altitude!=0)stream << std::endl; count=1; }
	  else { stream << ' '; count++; if(count==10)count=0; }

	  stream << std::fixed 
		 << std::setw(9) << std::setprecision(3) << idatum->second;

	  altitude++;
	}
      stream << std::endl;
      wavelength += 5;
    }
  return true;
}

VSSimDBWavelengthAltitudeDataset* VSSimDBWavelengthAltitudeDataset::
createFromCORSIKA(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(!stream)return 0;

  std::map<unsigned,VSSimDBAltitudeData> temp_data;
  std::string comment;
  
  std::getline(stream,comment);
  while(stream)
    {
      unsigned wavelength = 0;
      stream >> wavelength;
      if(!stream)continue;

      unsigned altitude=0;
      while(altitude<=50)
	{
	  float value = 0;
	  stream >> value;
	  assert(stream);
	  temp_data[wavelength][altitude]=value;
	  altitude++;
	}    
    }

  if(temp_data.size() == 0)
    {
      std::cerr << filename << ": did not read any data" << std::endl;
      return 0;
    }

  VSSimDBWavelengthAltitudeDataset* dataset = 
    new VSSimDBWavelengthAltitudeDataset;
  if(dataset==0)return 0;

  dataset->comment = comment;

  unsigned wavelength = 180;
  while(wavelength < temp_data.begin()->first)
    {
      std::cerr << filename << ": replicating " << temp_data.begin()->first
		<< "nm data at " << wavelength << "nm" << std::endl;
      dataset->data[wavelength] = temp_data.begin()->second;
      wavelength += 5;
    }

  std::map<unsigned,VSSimDBAltitudeData>::const_iterator iwl_lo = 
    temp_data.begin();
  std::map<unsigned,VSSimDBAltitudeData>::const_iterator iwl_hi = 
    iwl_lo;
  iwl_hi++;

  while((wavelength <= 700)&&(iwl_hi != temp_data.end()))
    {
      if(iwl_lo->first == wavelength)
	dataset->data[wavelength] = iwl_lo->second;
      else
	{
	  std::cerr << filename << ": interpolating " 
		    << iwl_lo->first << ',' << iwl_hi->first
		    << "nm data to " << wavelength << "nm" << std::endl;
	  double x = 
	    (double(wavelength) - double(iwl_lo->first))/
	    (double(iwl_hi->first) - double(iwl_lo->first));
	  for(unsigned altitude=0;altitude<=50;altitude++)
	    {
	      std::map<unsigned, float>::const_iterator 
		ialt_lo = iwl_lo->second.find(altitude);
	      std::map<unsigned, float>::const_iterator 
		ialt_hi = iwl_hi->second.find(altitude);
	      dataset->data[wavelength][altitude] = 
		ialt_lo->second*(1-x) + ialt_hi->second*x; 
	    }
	}
      
      wavelength += 5;

      if(iwl_hi->first <= wavelength)
	{
	  iwl_lo++;
	  iwl_hi++;
	}
    }

  while(wavelength <= 700)
    {
      if(wavelength != iwl_lo->first)
	std::cerr << filename << ": replicating " << iwl_lo->first
		  << "nm data at " << wavelength << "nm" << std::endl;
      dataset->data[wavelength] = iwl_lo->second;
      wavelength += 5;
    }

  dataset->writeToCORSIKA("test.dat");
  return dataset;
}

// ----------------------------------------------------------------------------
// VSSimDBModtranProfileDataset
// ----------------------------------------------------------------------------

bool VSSimDBModtranProfileDataset::
writeToCORSIKA(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  if(!stream)return false;

  stream << comment << std::endl;
  for(VSSimDBModtranProfileData::const_iterator idata = 
	data.begin(); idata!=data.end(); idata++)
    {
      stream << std::fixed 
	     << std::setw(9) << std::setprecision(3) << idata->altitude 
	     << "     "
	     << std::scientific
	     << std::setw(11) << std::setprecision(5) << idata->rho << "  "
	     << std::setw(11) << std::setprecision(5) << idata->thick << "  "
	     << std::setw(11) << std::setprecision(5) << idata->n_minus_one 
	     << std::endl;
    }
  return true;
}

VSSimDBModtranProfileDataset* VSSimDBModtranProfileDataset::
createFromCORSIKA(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(!stream)return 0;
  
  VSSimDBModtranProfileDataset* dataset = 
    new VSSimDBModtranProfileDataset;
  if(dataset==0)return 0;

  int test = stream.peek();
  while(test=='#')
    {
      std::string line;
      std::getline(stream,line);
      if(!dataset->comment.empty())dataset->comment += '\n';
      dataset->comment += line;
      test = stream.peek();
    }

  while(stream)
    {
      std::string line;
      std::getline(stream,line);
      if(stream)
	{
	  std::istringstream linestream(line);
	  VSSimDBModtranProfileDatum datum;
	  linestream >> datum.altitude
		     >> datum.rho
		     >> datum.thick
		     >> datum.n_minus_one;
	  dataset->data.push_back(datum);
	}
    }
  return dataset;
}

// ----------------------------------------------------------------------------
// VSSimDBCORSIKADatasets
// ----------------------------------------------------------------------------

VSSimDBCORSIKADatasets::~VSSimDBCORSIKADatasets()
{
  // nothing to see here
}

//
// Wavelength Data
//

void VSSimDBCORSIKADatasets::
storeWavelengthData(const std::string& table_name,
		    const VSSimDBWavelengthData& data)
{
  unsigned wavelength;
  float value;

  fDB->createTable(table_name,
		   fDB->sqlSpecOf("Wavelength",wavelength,true,"NOT NULL")+
		   fDB->sqlSpecOf("Value",value,false,"NOT NULL")+
		   ", PRIMARY KEY (Wavelength)",
		   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  fDB->deleteFromTable(table_name,"");
  
  VSDBStatement* stmt = fDB->createInsertQuery(table_name, 2);
  assert(stmt);

  stmt->bindToParam(wavelength);
  stmt->bindToParam(value);

  for(VSSimDBWavelengthData::const_iterator idatum = data.begin();
      idatum!=data.end(); idatum++)
    {
      wavelength = idatum->first;
      value = idatum->second;
      assert(stmt->execute());
    }
  
  delete stmt;
}

VSSimDBWavelengthData* VSSimDBCORSIKADatasets::
retrieveWavelengthData(const std::string& table_name)
{
  VSSimDBWavelengthData* data(0);

  VSDBStatement* stmt = 
    fDB->createSelectQuery(table_name,"","",
			   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(!stmt)return 0;

  unsigned wavelength;
  float value;
  stmt->bindToResult(wavelength);
  stmt->bindToResult(value);
  
  if(stmt->execute() > 0)
    {
      data = new VSSimDBWavelengthData;
      while(stmt->retrieveNextRow())
	(*data)[wavelength]=value;
    }
  
  delete stmt;
  return data;
}

void VSSimDBCORSIKADatasets::
storeWavelengthDataset(const std::string& table_name,
		       const VSSimDBWavelengthDataset& dataset)
{
  VSDBParameterTable db_param(fDB);
  db_param.createParameterTable();
  db_param.deleteParameterSet(table_name);

  VSDBParameterSet param;
  param["Comment"]=dataset.comment;
  db_param.storeParameterSet(table_name, param);

  storeWavelengthData(table_name,dataset.data);
}

VSSimDBWavelengthDataset* VSSimDBCORSIKADatasets::
retrieveWavelengthDataset(const std::string& table_name)
{
  VSSimDBWavelengthDataset* dataset(0);

  VSDBParameterTable db_param(fDB);

  VSDBParameterSet param;
  db_param.retrieveParameterSet(table_name, param);
  
  VSSimDBWavelengthData* data = retrieveWavelengthData(table_name);
  if(data)
    {
      dataset = new VSSimDBWavelengthDataset;
      dataset->comment=param["Comment"];
      dataset->data=*data;
    }

  return dataset;
}

//
// Wavelength-Altitude Data
//

void VSSimDBCORSIKADatasets::
storeWavelengthAltitudeData(const std::string& table_name,
			    const VSSimDBWavelengthAltitudeData& data)
{
  unsigned wavelength;
  unsigned altitude;
  float value;

  fDB->createTable(table_name,
		   fDB->sqlSpecOf("Wavelength",wavelength,true,"NOT NULL")+
		   fDB->sqlSpecOf("Altitude",altitude,false,"NOT NULL")+
		   fDB->sqlSpecOf("Value",value,false,"NOT NULL")+
		   ", PRIMARY KEY (Wavelength, Altitude)",
		   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  fDB->deleteFromTable(table_name,"");
  
  VSDBStatement* stmt = fDB->createInsertQuery(table_name, 3);
  assert(stmt);

  stmt->bindToParam(wavelength);
  stmt->bindToParam(altitude);
  stmt->bindToParam(value);

  for(VSSimDBWavelengthAltitudeData::const_iterator iwldata = data.begin();
      iwldata!=data.end(); iwldata++)
    {
      wavelength = iwldata->first;
      for(VSSimDBAltitudeData::const_iterator idatum = iwldata->second.begin();
	  idatum!=iwldata->second.end(); idatum++)
	{
	  altitude = idatum->first;
	  value = idatum->second;
	  assert(stmt->execute());
	}
    }
  
  delete stmt;
}

VSSimDBWavelengthAltitudeData* VSSimDBCORSIKADatasets::
retrieveWavelengthAltitudeData(const std::string& table_name)
{
  VSSimDBWavelengthAltitudeData* data(0);

  VSDBStatement* stmt = 
    fDB->createSelectQuery(table_name,"","*",
			   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(!stmt)return 0;

  unsigned wavelength;
  unsigned altitude;
  float value;
  stmt->bindToResult(wavelength);
  stmt->bindToResult(altitude);
  stmt->bindToResult(value);
  
  if(stmt->execute() > 0)
    {
      data = new VSSimDBWavelengthAltitudeData;
      while(stmt->retrieveNextRow())
	(*data)[wavelength][altitude]=value;
    }
  
  delete stmt;
  return data;
}

void VSSimDBCORSIKADatasets::
storeWavelengthAltitudeDataset(const std::string& table_name,
			       const VSSimDBWavelengthAltitudeDataset& dataset)
{
  VSDBParameterTable db_param(fDB);
  db_param.createParameterTable();
  db_param.deleteParameterSet(table_name);

  VSDBParameterSet param;
  param["Comment"]=dataset.comment;
  db_param.storeParameterSet(table_name, param);

  storeWavelengthAltitudeData(table_name,dataset.data);
}

VSSimDBWavelengthAltitudeDataset* VSSimDBCORSIKADatasets::
retrieveWavelengthAltitudeDataset(const std::string& table_name)
{
  VSSimDBWavelengthAltitudeDataset* dataset(0);

  VSDBParameterTable db_param(fDB);

  VSDBParameterSet param;
  db_param.retrieveParameterSet(table_name, param);

  VSSimDBWavelengthAltitudeData* data =
    retrieveWavelengthAltitudeData(table_name);
  if(data)
    {
      dataset=new VSSimDBWavelengthAltitudeDataset;
      dataset->comment=param["Comment"];
      dataset->data = *data;
    }

  return dataset;
}

//
// Modtran Profile Data
//

void VSSimDBCORSIKADatasets::
storeModtranProfileData(const std::string& table_name,
			const VSSimDBModtranProfileData& data)
{
  float altitude;
  float rho;
  float thick;
  float n_minus_one;

  fDB->createTable(table_name,
		   fDB->sqlSpecOf("Altitude",altitude,true,"NOT NULL")+
		   fDB->sqlSpecOf("Rho",rho,false,"NOT NULL")+
		   fDB->sqlSpecOf("Thickness",thick,false,"NOT NULL")+
		   fDB->sqlSpecOf("NMinusOne",n_minus_one,false,"NOT NULL")+
		   ", PRIMARY KEY (Altitude)",
		   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  fDB->deleteFromTable(table_name,"");
  
  VSDBStatement* stmt = fDB->createInsertQuery(table_name, 4);
  assert(stmt);

  stmt->bindToParam(altitude);
  stmt->bindToParam(rho);
  stmt->bindToParam(thick);
  stmt->bindToParam(n_minus_one);

  for(VSSimDBModtranProfileData::const_iterator iprof = data.begin();
      iprof!=data.end(); iprof++)
    {
      altitude = iprof->altitude;
      rho = iprof->rho;
      thick = iprof->thick;
      n_minus_one = iprof->n_minus_one;
      stmt->execute();
    }
  
  delete stmt;
}

VSSimDBModtranProfileData* VSSimDBCORSIKADatasets::
retrieveModtranProfileData(const std::string& table_name)
{
  VSSimDBModtranProfileData* data(0);

  VSDBStatement* stmt = 
    fDB->createSelectQuery(table_name,"","*",
			   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(!stmt)return 0;

  float altitude;
  float rho;
  float thick;
  float n_minus_one;

  stmt->bindToResult(altitude);
  stmt->bindToResult(rho);
  stmt->bindToResult(thick);
  stmt->bindToResult(n_minus_one);
  
  if(stmt->execute() > 0)
    {
      data = new VSSimDBModtranProfileData;
      while(stmt->retrieveNextRow())
	{
	  VSSimDBModtranProfileDatum datum;
	  datum.altitude = altitude;
	  datum.rho = rho;
	  datum.thick = thick;
	  datum.n_minus_one = n_minus_one;
	  data->push_back(datum);
	}
    }
  std::sort(data->begin(), data->end());
  
  delete stmt;
  return data;
}

void VSSimDBCORSIKADatasets::
storeModtranProfileDataset(const std::string& table_name,
			   const VSSimDBModtranProfileDataset& dataset)
{
  VSDBParameterTable db_param(fDB);
  db_param.createParameterTable();
  db_param.deleteParameterSet(table_name);

  VSDBParameterSet param;
  param["Comment"]=dataset.comment;
  db_param.storeParameterSet(table_name, param);

  storeModtranProfileData(table_name,dataset.data);
}

VSSimDBModtranProfileDataset* VSSimDBCORSIKADatasets::
retrieveModtranProfileDataset(const std::string& table_name)
{
  VSSimDBModtranProfileDataset* dataset = 0;

  VSDBParameterTable db_param(fDB);

  VSDBParameterSet param;
  db_param.retrieveParameterSet(table_name, param);

  VSSimDBModtranProfileData* data = retrieveModtranProfileData(table_name);
  if(data)
    {
      dataset = new VSSimDBModtranProfileDataset;
      dataset->comment=param["Comment"];
      dataset->data = *data;
    }

  return dataset;
}

#ifdef TESTMAIN

#include<VSOptions.hpp>

void main()
{

}

#endif
