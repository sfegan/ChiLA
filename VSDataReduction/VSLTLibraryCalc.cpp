//-*-mode:c++; mode:font-lock;-*-

#include <fstream>

#include "VSLTLibraryCalc.hpp"
#include <VSLineTokenizer.hpp>
#include <VSFileUtility.hpp>
#include <VSOctaveH5Reader.hpp>

using namespace VERITAS;

bool VSLTLibraryCalc::load(const std::string& file,
			   const VSTime& date,
			   const std::vector<unsigned>& nchan,
			   double zn_deg, double az_deg, double ped_rms)
{

  //
  // Here's where it decides whether the sp_parameter_lookup option is a list
  //   of files or just a single h5 file.

  //std::cout<<"\n  Loading file "<<file<<" in VSLTLibraryCalc::load()...\n\n";
  
  if(file.empty()) return false;
  else if(!VSFileUtility::isFile(file)) return false;
  else if(VSOctaveH5ReaderStruct::isHDF5(file)) 
    return load(file,nchan,zn_deg,az_deg,ped_rms);

  std::ifstream datastream(file.c_str());
  if(!datastream)return false;

  VSLineTokenizer tokenizer;
  VSTokenList tokens;
  std::string line;
  while(getline(datastream,line))
    {
      std::string line_copy = line;
      tokenizer.tokenize(line_copy, tokens);

      //std::cout<<"tokens.size(): "<<tokens.size()<<std::endl;

      if(line_copy.substr(0,1) == "#" || tokens.size() == 0) continue;
      else if(tokens.size() == 3)
	{

	  std::string lo_date_string = tokens[0].string();
	  std::string hi_date_string = tokens[1].string();
	  std::string table_file = tokens[2].string();

	  VSTime lo_date;
	  VSTime hi_date;
	  
	  uint64_t lo_date_uint64 = 0;
	  VSDataConverter::fromString(lo_date_uint64,lo_date_string);
	  if(!lo_date.setFromDBTimeStamp(lo_date_uint64))
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": "
			<< "Error parsing timestamp "
			<< lo_date_uint64 << std::endl;
	      continue;
	    }

	  uint64_t hi_date_uint64 = 0;
	  VSDataConverter::fromString(hi_date_uint64,hi_date_string);
	  if(!hi_date.setFromDBTimeStamp(hi_date_uint64))
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": "
			<< "Error parsing timestamp "
			<< hi_date_uint64 << std::endl;
	      continue;
	    }

	  if(date >= lo_date && date <= hi_date)
	    return load(table_file,nchan,zn_deg,az_deg,ped_rms);

	}
      else if(tokens.size() == 5)
	{
	  std::string lo_date_string = tokens[0].string();
	  std::string lo_time_string = tokens[1].string();
	  std::string hi_date_string = tokens[2].string();
	  std::string hi_time_string = tokens[3].string();
	  std::string table_file = tokens[4].string();

	  VSTime lo_date;
	  VSTime hi_date;

	  if(!lo_date.setFromString(lo_date_string + " " + lo_time_string))
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": "
			<< "Error parsing timestamp "
			<< lo_date_string + " " + lo_time_string << std::endl;
	      continue;
	    }

	  if(!hi_date.setFromString(hi_date_string + " " + hi_time_string))
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": "
			<< "Error parsing timestamp "
			<< hi_date_string + " " + hi_time_string << std::endl;
	      continue;
	    }

	  if(date >= lo_date && date <= hi_date)
	    return load(table_file,nchan,zn_deg,az_deg,ped_rms);

	}
      
    }

  // Table did not load correctly for the given run...Quit with error...
  std::cerr<< "\nERROR! "<< std::string(__PRETTY_FUNCTION__) + ": "
	   << "Run is out of range of "<<file<<" dates."
	   << std::endl;
  std::cerr<<"Exit failure"<<std::endl;
  exit(EXIT_FAILURE);

  return false;
}
