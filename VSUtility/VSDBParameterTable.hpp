//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBParameterTable.hpp
  Database access to parameter table

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       05/18/2005
*/

#ifndef VSPARAMETERTABLE_HPP
#define VSPARAMETERTABLE_HPP

#include <map>
#include <set>
#include <string>

#include "VSDatabase.hpp"

//! VERITAS namespace
namespace VERITAS
{
  
  typedef std::string                                 VSDBCollection;
  typedef std::string                                 VSDBParameter;
  typedef std::string                                 VSDBValue;

  typedef std::map<VSDBParameter, VSDBValue>          VSDBParameterSet;
  typedef std::map<VSDBCollection, VSDBParameterSet>  VSDBCollectionSet;
  typedef std::set<VSDBCollection>                    VSDBCollectionNameSet;
  
  class VSDBParameterTable
  {
  public:
    VSDBParameterTable(VSDatabase* db): fDB(db) { }
    virtual ~VSDBParameterTable();

    void createParameterTable();

    int deleteParameterSet(const VSDBCollection& collection);
    int storeParameterSet(const VSDBCollection& collection,
			  const VSDBParameterSet& parameters);
    int retrieveParameterSet(const VSDBCollection& collection,
			     VSDBParameterSet& parameters);

    int retrieveAllParameterSets(VSDBCollectionSet& collections);

    int deleteParameter(const VSDBCollection& collection,
			const VSDBParameter& parameter);

    int setOrUpdateParameterValue(const VSDBCollection& collection,
				  const VSDBParameter& parameter,
				  const VSDBValue& value);

  private:
    VSDatabase* fDB;
  };

}

#endif // VSPARAMETERTABLE_HPP
