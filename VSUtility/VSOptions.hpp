//-*-mode:c++; mode:font-lock;-*-

// This file is shared between the Simulations and OAWG. To avoid a
// fork in development, please coordinate changes with Stephen Fegan.
// email: sfegan@astro.ucla.edu

// Old header for consistency with Simulations package

/*! \file VSOptions.hpp
  Class for processing of command line options

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

// New header for consistency with OAWG package

/**
 *
 * Original Author: Stephen Fegan
 * $Author: sfegan $
 * $Date: 2008/03/14 01:41:01 $
 * $Revision: 1.13 $
 * $Tag$
 *
 **/

#ifndef VSOPTIONS_HPP
#define VSOPTIONS_HPP

#include<ostream>
#include<sstream>
#include<vector>
#include<list>
#include<map>
#include<string>
#include<memory>
#include<vsassert>

#include<VSDataConverter.hpp>

#ifndef _OAWG

namespace VERITAS
{
#endif

  //! Class encapsulating extraction of data type from option string
  template<typename T> class VSOptionValueExtractor
  {
  public:
    static inline bool extract(const std::string& buffer, T& value);
    static inline std::string typeName();
  };

  //! Command line options processor
  class VSOptions
  {
  public:
    enum FindStatus { FS_NOT_FOUND, FS_FOUND, 
		      FS_FOUND_WITH_UNDESIRED_VALUE,
		      FS_FOUND_BUT_WITHOUT_VALUE, 
		      FS_FOUND_BUT_COULD_NOT_CONVERT_VALUE };

    struct OptionRecord
    {
      OptionRecord(): key(), status(), val_requested(), val() 
      { /* nothing to see here */ }
      std::string key;
      FindStatus  status;
      bool        val_requested;
      std::string val;
    };

    VSOptions(int& argc, char**argv, bool print_default = false);
    virtual ~VSOptions();
    
    void addOption(const std::string& opt, bool override = false);
    void addOptionWithValue(const std::string& opt, const std::string& val,
			    bool override = false);

    const std::string& arg0() const { return fArg0; }
    int argC() const { return fArgC; }
    char** argV() const { return fArgV; };
    
    inline FindStatus find(const std::string& key, const std::string& help="",
			   bool option_is_simple = true,
			   const std::string& catagory = "");

    inline FindStatus find(const std::string& key, const std::string& help,
			   const std::string& catagory,
			   bool option_is_simple = true)
    { return find(key,help,option_is_simple,catagory); }

    inline FindStatus find(const std::string& key, const std::string& help,
			   const char* catagory,
			   bool option_is_simple = true)
    { return find(key,help,option_is_simple,std::string(catagory)); }
    
    template<typename T> inline FindStatus 
    findWithValue(const std::string& key, T& value, 
		  const std::string& help="", 
		  bool option_is_simple = true,
		  bool attached_val_only = false,
		  const std::string& catagory = "");

    template<typename T> inline FindStatus 
    findWithValue(const std::string& key, T& value, 
		  const std::string& help, 
		  const std::string& catagory,
		  bool option_is_simple = true,
		  bool attached_val_only = false)
    { return findWithValue(key, value, help, option_is_simple,
			   attached_val_only, catagory); }

    template<typename T> inline FindStatus 
    findWithValue(const std::string& key, T& value, 
		  const std::string& help, 
		  const char* catagory,
		  bool option_is_simple = true,
		  bool attached_val_only = false)
    { return findWithValue(key, value, help, option_is_simple,
			   attached_val_only, std::string(catagory)); }

    inline FindStatus findBoolValue(const std::string& key, bool& value, 
				    bool value_when_not_explicitly_set,
				    const std::string& help="",
				    bool option_is_simple = true,
				    const std::string& catagory = "");

    inline FindStatus findBoolValue(const std::string& key, bool& value, 
				    bool value_when_not_explicitly_set,
				    const std::string& help,
				    const std::string& catagory,
				    bool option_is_simple = true)
    { return findBoolValue(key, value, value_when_not_explicitly_set, help, 
			   option_is_simple, catagory); }
    

    inline FindStatus findBoolValue(const std::string& key, bool& value, 
				    bool value_when_not_explicitly_set,
				    const std::string& help,
				    const char* catagory,
				    bool option_is_simple = true)
    { return findBoolValue(key, value, value_when_not_explicitly_set, help, 
			   option_is_simple, std::string(catagory)); }
    
    bool assertNoOptions() const;

    void printUsage(std::ostream& stream, 
		    bool print_only_simple_options = false,
		    unsigned USAGE_OPTION_MAX = 30) const;

    const std::vector<OptionRecord>& 
    getOptionRecords() const { return fOptionRecords; }

    bool wasOptionRegistered(const std::string& key) const;
    bool wasOptionFound(const std::string& key) const;
    const OptionRecord& getOptionRecord(const std::string& key) const;

    void addCatagory(const std::string& catagory, const std::string& text = "")
    { fCatagories.push_back(CatagoryHelp(catagory, text)); }

  private:
    struct OptionInfo
    {
      std::string                               fVal;
      bool                                      fValAttached;
      std::vector<unsigned>                     fArgBlank;
      std::vector<unsigned>                     fValBlank;
      OptionInfo(): fVal(), fValAttached(false), fArgBlank(), fValBlank() {}
    };

    struct OptionHelp
    {
      std::string                               fOption;
      bool                                      fArgRequired;
      std::string                               fArgTypeName;
      std::string                               fText;
      bool                                      fIsSimple;
      std::string                               fCatagory;
      OptionHelp(const std::string option, bool arg_required, 
		 const std::string arg_type_name, const std::string& text,
		 bool option_is_simple, const std::string& catagory):
	fOption(option), fArgRequired(arg_required), 
	fArgTypeName(arg_type_name), fText(text), 
	fIsSimple(option_is_simple), fCatagory(catagory) { }
    };

    struct CatagoryHelp
    {
      std::string                               fCatagory;
      std::string                               fText;
      CatagoryHelp(const std::string& catagory, const std::string& text):
	fCatagory(catagory), fText(text) { }
    };

    VSOptions();
    VSOptions(const VSOptions&);
    const VSOptions& operator =(const VSOptions&);
    
    FindStatus doFind(const std::string& key, std::string* value,
		      const std::string& help, bool option_is_simple,
		      const std::string& type, bool attached_val_only,
		      const std::string& catagory);

    std::vector<OptionRecord>::const_iterator 
    returnOptionRecord(const std::string& key) const;

    std::map<std::string, OptionInfo> fOptions;
    std::string                       fArg0;
    std::vector<char*>                fArgs;
    int&                              fArgC;
    char**                            fArgV;

    bool                              fPrintDefault;
    std::list<OptionHelp>             fUsageInfo;
    std::list<CatagoryHelp>           fCatagories;

    std::vector<OptionRecord>         fOptionRecords;
  };

  template<typename T> inline bool 
  VSOptionValueExtractor<T>::extract(const std::string& buffer, T& value)
  {
    return(VSDatumConverter<T>::fromString(value,buffer.c_str()));
  }

  template<typename T> inline std::string
  VSOptionValueExtractor<T>::typeName()
  {
    return(VSDatumConverter<T>::typeName());
  }
  
  inline VSOptions::FindStatus VSOptions::
  find(const std::string& key, const std::string& help,
       bool option_is_simple, const std::string& catagory)
  {
    FindStatus status = doFind(key,0,help,option_is_simple,"",false,catagory);

    OptionRecord record;
    record.key = key;
    record.status = status;
    record.val_requested = false;
    record.val = "";
    fOptionRecords.push_back(record);

    return status;
  }
  
  template<typename T> inline VSOptions::FindStatus VSOptions::
  findWithValue(const std::string& key, T& value, const std::string& help,
		bool option_is_simple, bool attached_val_only,
		const std::string& catagory)
  {
    std::string val_string;
    std::string type = VSOptionValueExtractor<T>::typeName();
    if(fPrintDefault)
      {
	std::string default_value;
	VSDataConverter::toString(default_value,value,true);
	if(default_value.empty())default_value="\"\"";
	type+="=";
	type+=default_value;
      }
    FindStatus status = doFind(key, &val_string, help, option_is_simple, 
			       type, attached_val_only, catagory);

    if((status!=FS_NOT_FOUND)&&(status!=FS_FOUND_BUT_WITHOUT_VALUE))
      {
	if(VSOptionValueExtractor<T>::extract(val_string,value))
	  status = FS_FOUND;
	else
	  status = FS_FOUND_BUT_COULD_NOT_CONVERT_VALUE;
      }
    
    OptionRecord record;
    record.key           = key;
    record.status        = status;
    record.val_requested = true;
    VSDataConverter::toString(record.val,value,true);
    fOptionRecords.push_back(record);

    return status;
  }

  inline VSOptions::FindStatus VSOptions::
  findBoolValue(const std::string& key, bool& value, 
		bool value_when_not_explicitly_set, const std::string& help, 
		bool option_is_simple, const std::string& catagory)
  {
    FindStatus fs = 
      findWithValue(key, value, help, option_is_simple, true, catagory);
    if(fs == FS_FOUND_BUT_WITHOUT_VALUE)value = value_when_not_explicitly_set;
    return fs;
  }
  
#ifndef _OAWG
}
#endif

#endif // VSOPTIONS_HPP
