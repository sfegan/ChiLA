//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBWavelengthTable.cpp
  Database access to wavelength/value tables

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       08/19/2005
*/

#include <fstream>
#include <iomanip>

#include <VSDBWavelengthTable.hpp>
#include <VSDBParameterTable.hpp>

using namespace VERITAS;

// ----------------------------------------------------------------------------
// VSDBWavelengthDataset
// ----------------------------------------------------------------------------

bool VSDBWavelengthDataset::writeToCORSIKA(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  if(!stream)return false;

  stream << comment << std::endl;
  unsigned wavelength=180;
  unsigned count=0;
  while(wavelength<=700)
    {
      float value = 0;
      VSDBWavelengthData::const_iterator idatum = data.find(wavelength);
      if(idatum != data.end())value=idatum->second;
      stream << std::setprecision(3) << value;
      if(count==7){ stream << std::endl; count=0; }
      else { stream << ' '; count++; }
      wavelength += 5;
    }
  return true;
}

VSDBWavelengthDataset* 
VSDBWavelengthDataset::createFromCORSIKA(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(!stream)return 0;
  
  VSDBWavelengthDataset* dataset = new VSDBWavelengthDataset;
  if(dataset==0)return 0;

  std::getline(stream,dataset->comment);

  unsigned wavelength=180;
  while(wavelength<=700)
    {
      float value;
      stream >> value;
      dataset->data[wavelength]=value;
    }    
  return dataset;
}

// ----------------------------------------------------------------------------
// VSDBWavelengthAltitudeDataset
// ----------------------------------------------------------------------------

bool VSDBWavelengthAltitudeDataset::
writeToCORSIKA(const std::string& filename) const
{
  std::ofstream stream(filename.c_str());
  if(!stream)return false;

  stream << comment << std::endl;
  unsigned wavelength=180;
  while(wavelength<=700)
    {
      VSDBWavelengthAltitudeData::const_iterator iwlset = 
	data.find(wavelength);
      stream << std::setw(4) << wavelength << std::endl;
      unsigned altitude=0;
      unsigned count=0;
      while(altitude<=50)
	{
	  float value = 0;
	  if(iwlset != data.end())
	    {
	      VSDBAltitudeData::const_iterator idatum = 
		iwlset->second.find(altitude);
	      if(idatum != iwlset->second.end())value=idatum->second;
	    }
	  stream << std::setw(9) << std::setprecision(3) << value;
	  if(count==9) { stream << std::endl; count=0; }
	  else { stream << ' '; count++; }
	  altitude++;
	}
      wavelength += 5;
    }
  return true;
}

VSDBWavelengthAltitudeDataset* 
VSDBWavelengthAltitudeDataset::createFromCORSIKA(const std::string& filename)
{
  std::ifstream stream(filename.c_str());
  if(!stream)return 0;
  
  VSDBWavelengthAltitudeDataset* dataset = new VSDBWavelengthAltitudeDataset;
  if(dataset==0)return 0;

  std::getline(stream,dataset->comment);

  unsigned wavelength=180;
  while(wavelength<=700)
    {
      unsigned unused;
      stream >> unused;
      unsigned altitude=0;
      while(altitude<=50)
	{
	  float value;
	  stream >> value;
	  dataset->data[wavelength][altitude]=value;
	}    
    }
  return dataset;
}

// ----------------------------------------------------------------------------
// VSDBCORSIKADatasets
// ----------------------------------------------------------------------------

VSDBCORSIKADatasets::~VSDBCORSIKADatasets()
{
  // nothing to see here
}

void VSDBCORSIKADatasets::
storeWavelengthData(const std::string& table_name,
		    const VSDBWavelengthData& data)
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

  for(VSDBWavelengthData::const_iterator idatum = data.begin();
      idatum!=data.end(); idatum++)
    {
      wavelength = idatum->first;
      value = idatum->second;
      assert(stmt->execute());
    }
  
  delete stmt;
}

VSDBWavelengthData VSDBCORSIKADatasets::
retrieveWavelengthData(const std::string& table_name)
{
  VSDBWavelengthData data;

  VSDBStatement* stmt = fDB->createSelectQuery(table_name);
  if(!stmt)return data;

  unsigned wavelength;
  float value;
  stmt->bindToResult(wavelength);
  stmt->bindToResult(value);
  
  if(stmt->execute() > 0)
    {
      while(stmt->retrieveNextRow())
	data[wavelength]=value;
    }
  
  delete stmt;
  return data;
}

void VSDBCORSIKADatasets::
storeWavelengthDataset(const std::string& table_name,
		       const VSDBWavelengthDataset& dataset)
{
  VSDBParameterTable db_param(fDB);
  db_param.createParameterTable();
  db_param.deleteParameterSet(table_name);

  VSDBParameterSet param;
  param["Comment"]=dataset.comment;
  db_param.storeParameterSet(table_name, param);

  storeWavelengthData(table_name,dataset.data);
}

VSDBWavelengthDataset VSDBCORSIKADatasets::
retrieveWavelengthDataset(const std::string& table_name)
{
  VSDBWavelengthDataset dataset;

  VSDBParameterTable db_param(fDB);

  VSDBParameterSet param;
  db_param.retrieveParameterSet(table_name, param);
  dataset.comment=param["Comment"];

  dataset.data = retrieveWavelengthData(table_name);
  return dataset;
}

void VSDBCORSIKADatasets::
storeWavelengthAltitudeData(const std::string& table_name,
			    const VSDBWavelengthAltitudeData& data)
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
  
  VSDBStatement* stmt = fDB->createInsertQuery(table_name, 2);
  assert(stmt);

  stmt->bindToParam(wavelength);
  stmt->bindToParam(altitude);
  stmt->bindToParam(value);

  for(VSDBWavelengthAltitudeData::const_iterator iwldata = data.begin();
      iwldata!=data.end(); iwldata++)
    {
      wavelength = iwldata->first;
      for(VSDBAltitudeData::const_iterator idatum = iwldata->second.begin();
	  idatum!=iwldata->second.end(); idatum++)
	{
	  altitude = idatum->first;
	  value = idatum->second;
	  assert(stmt->execute());
	}
    }
  
  delete stmt;
}

// ----------------------------------------------------------------------------
// VSDBWavelengthAltitudeTable
// ----------------------------------------------------------------------------

VSDBWavelengthAltitudeData VSDBCORSIKADatasets::
retrieveWavelengthAltitudeData(const std::string& table_name)
{
  VSDBWavelengthAltitudeData data;

  VSDBStatement* stmt = fDB->createSelectQuery(table_name);
  if(!stmt)return data;

  unsigned wavelength;
  unsigned altitude;
  float value;
  stmt->bindToResult(wavelength);
  stmt->bindToResult(altitude);
  stmt->bindToResult(value);
  
  if(stmt->execute() > 0)
    {
      while(stmt->retrieveNextRow())
	data[wavelength][altitude]=value;
    }
  
  delete stmt;
  return data;
}

void VSDBCORSIKADatasets::
storeWavelengthAltitudeDataset(const std::string& table_name,
			       const VSDBWavelengthAltitudeDataset& dataset)
{
  VSDBParameterTable db_param(fDB);
  db_param.createParameterTable();
  db_param.deleteParameterSet(table_name);

  VSDBParameterSet param;
  param["Comment"]=dataset.comment;
  db_param.storeParameterSet(table_name, param);

  storeWavelengthAltitudeData(table_name,dataset.data);
}

VSDBWavelengthAltitudeDataset VSDBCORSIKADatasets::
retrieveWavelengthAltitudeDataset(const std::string& table_name)
{
  VSDBWavelengthAltitudeDataset dataset;

  VSDBParameterTable db_param(fDB);

  VSDBParameterSet param;
  db_param.retrieveParameterSet(table_name, param);
  dataset.comment=param["Comment"];

  dataset.data = retrieveWavelengthAltitudeData(table_name);
  return dataset;
}
