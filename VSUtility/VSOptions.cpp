//-*-mode:c++; mode:font-lock;-*-

// This file is shared between the Simulations and OAWG. To avoid a
// fork in development, please coordinate changes with Stephen Fegan.
// email: sfegan@astro.ucla.edu

// Old header for consistency with Simulations package

/*! \file VSOptions.cpp
  Class for processing of command line options

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

// New header for consistency with OAWG package

/**
 * \class VSOptions
 * \ingroup common
 * \brief This is a one-line description of this cpp file.
 *
 * Here is a tedious verbose multi-line description of all
 * the details of the code, more than you would
 * ever want to read. Generally, all the important documentation
 * goes in the .cpp files.
 *
 * Original Author: Stephen Fegan
 * $Author: sfegan $
 * $Date: 2008/03/14 01:41:01 $
 * $Revision: 1.11 $
 * $Tag$
 *
 **/

// These deinitions tell the makefile which library the cpp file
// should be included in
// VA_LIBRARY_TAG: libSP24common.a
// VA_LIBRARY_TAG: libSP24commonLite.a

#include<iostream>
#include<cctype>

#include"VSOptions.hpp"

#ifndef _OAWG
using namespace VERITAS;
#endif

VSOptions::VSOptions(int& argc, char**argv, bool print_default):
  fOptions(), fArg0(*argv), fArgs(), fArgC(argc), fArgV(argv),
  fPrintDefault(print_default), fUsageInfo(), fCatagories(), fOptionRecords()
{
  for(int i=0;i<=argc;i++)fArgs.push_back(argv[i]);
  
  for(unsigned iarg=0;iarg<fArgs.size();iarg++)
    {
      if(fArgs[iarg] == 0)continue;
      std::string arg(fArgs[iarg]);

      if((arg.length()==0)||(arg[0]!='-'))continue;

      if((arg.length()==1)||((arg.length()==2)&&(arg[1]=='-')))
	{
	  fArgs[iarg] = 0;
	  break;
	}
      std::string::size_type i_key_start = arg[1]=='-' ? 2 : 1;
      std::string::size_type i_key_end = arg.find('=');

      std::string i_key;
      if(i_key_end == std::string::npos)i_key = arg.substr(i_key_start);
      else i_key = arg.substr(i_key_start, i_key_end - i_key_start);

      if(i_key.empty())continue;
      
      fOptions[i_key].fArgBlank.push_back(iarg);
 
      if(i_key_end != std::string::npos)
	{
	  i_key_end++;
	  fOptions[i_key].fValAttached = true;
	  if(arg.length() > i_key_end)
	    fOptions[i_key].fVal = arg.substr(i_key_end);
	  else fOptions[i_key].fVal.erase();
	}
      else
	{
	  unsigned jarg = iarg+1;
	  if((jarg < fArgs.size())&&(fArgs[jarg] != 0)&&(fArgs[jarg][0]!='-'))
	    {
	      fOptions[i_key].fVal = std::string(fArgs[jarg]);
	      fOptions[i_key].fValBlank.push_back(jarg);
	    }
	}
    }
}

VSOptions::~VSOptions()
{
  // nothing to see here
}

void VSOptions::addOption(const std::string& opt, bool override)
{
  if((override)||(fOptions.find(opt) == fOptions.end()))
    {
      fOptions[opt].fVal         = "";
      fOptions[opt].fValAttached = false;
    }
}

void VSOptions::
addOptionWithValue(const std::string& opt, const std::string& val, 
		   bool override)
{
  if((override)||(fOptions.find(opt) == fOptions.end()))
    {
      fOptions[opt].fVal         = val;
      fOptions[opt].fValAttached = true;
    }
}

VSOptions::FindStatus 
VSOptions::doFind(const std::string& key, std::string* value,
		  const std::string& help, bool option_is_simple,
		  const std::string& type, bool attached_val_only,
		  const std::string& catagory)
{
  OptionHelp option(key, value!=0, type, help, option_is_simple, catagory);
  fUsageInfo.push_back(option);

  std::map<std::string, OptionInfo>::iterator f = fOptions.find(key);
  if(f == fOptions.end())return FS_NOT_FOUND;
  
  for(std::vector<unsigned>::iterator i=(*f).second.fArgBlank.begin();
      i!=(*f).second.fArgBlank.end(); i++)
    fArgs[*i] = 0;

  VSOptions::FindStatus status = FS_NOT_FOUND;

  if((value == 0)&&((*f).second.fValAttached))
    status = FS_FOUND_WITH_UNDESIRED_VALUE;
  else if(value == 0)
    status = FS_FOUND;
  else if((!(*f).second.fValAttached)
	  &&((attached_val_only) || ((*f).second.fValBlank.empty())))
    {
      vsassert((attached_val_only)||(*f).second.fVal.empty());
      status = FS_FOUND_BUT_WITHOUT_VALUE;
    }
  else
    {
      for(std::vector<unsigned>::iterator i=(*f).second.fValBlank.begin();
	  i!=(*f).second.fValBlank.end(); i++)
	fArgs[*i] = 0;
      *value = (*f).second.fVal;
      status = FS_FOUND;
    }
  
  if(status != FS_NOT_FOUND)
    {
      fArgC=0;
      for(std::vector<char*>::iterator i = fArgs.begin(); i!=fArgs.end(); i++)
	if(*i)*(fArgV+(fArgC++)) = *i;
      fArgV[fArgC]=0;
    }

  return status;
}

std::vector<VSOptions::OptionRecord>::const_iterator 
VSOptions::returnOptionRecord(const std::string& key) const
{
  std::vector<OptionRecord>::const_iterator ior = fOptionRecords.begin();
  while((ior != fOptionRecords.end())&&(ior->key != key))ior++;
  return ior;
}

bool VSOptions::wasOptionRegistered(const std::string& key) const
{
  std::vector<OptionRecord>::const_iterator ior = returnOptionRecord(key);
  return ior != fOptionRecords.end();
}

bool VSOptions::wasOptionFound(const std::string& key) const
{
  std::vector<OptionRecord>::const_iterator ior = returnOptionRecord(key);
  vsassert(ior != fOptionRecords.end());
  return ior->status != FS_NOT_FOUND;
}

const VSOptions::OptionRecord& 
VSOptions::getOptionRecord(const std::string& key) const
{
  std::vector<OptionRecord>::const_iterator ior = returnOptionRecord(key);
  vsassert(ior != fOptionRecords.end());
  return *ior;
}

bool VSOptions::assertNoOptions() const
{
  std::map<std::string, OptionInfo>::const_iterator opt = fOptions.begin();
  while(opt!=fOptions.end())
    {
      for(std::vector<unsigned>::const_iterator i=
	    (*opt).second.fArgBlank.begin();
	  i!=(*opt).second.fArgBlank.end(); i++)
	if(fArgs[*i] != 0)return false;
      opt++;
    }
  return true;
}

#define MAXIMIZE(x,y) (x) = ((x)>(y))?(y):(x)

#define USAGE_INDENT       1
#define USAGE_SPACING      2
//#define USAGE_OPTION_MAX   30
#define USAGE_WIDTH        80

static void format(std::list<std::string>& lines, const std::string& text,
		   unsigned width)
{
  std::vector<std::string> tokens;
  std::string token;
  for(std::string::const_iterator ichar = text.begin(); 
      ichar != text.end(); ichar++)
    {
      if(isspace(*ichar))
	{
	  if(!token.empty())
	    {
	      tokens.push_back(token);
	      token.erase();
	    }
	}
      else
	token+=*ichar;
    }

  if(!token.empty())tokens.push_back(token);

  std::string oneline;
  for(std::vector<std::string>::const_iterator itoken = tokens.begin();
      itoken != tokens.end(); itoken++)
    {
      if(oneline.length() + itoken->length() + 1 > width)
	{
	  lines.push_back(oneline);
	  oneline.clear();
	}
      else if(!oneline.empty())
	{
	  oneline += ' ';
	}
	      
      oneline += *itoken;
    }
  
  if(!oneline.empty())lines.push_back(oneline);
  return;
}

void VSOptions::printUsage(std::ostream& stream, 
			   bool print_only_simple_options,
			   unsigned USAGE_OPTION_MAX) const
{
  std::map<std::string, std::list<std::string> > cat_map;
  std::list<std::string> cat_list;
  cat_map[std::string()] = std::list<std::string>();
  cat_list.push_back(std::string());
  
  for(std::list<CatagoryHelp>::const_iterator icat = fCatagories.begin();
      icat != fCatagories.end(); icat++)
    {
      if(cat_map.find(icat->fCatagory) == cat_map.end())
	cat_list.push_back(icat->fCatagory);
      std::list<std::string> text;
      format(text,icat->fText,USAGE_WIDTH);
      if(!text.empty())
	{
	  if(!cat_map[icat->fCatagory].empty())
	    cat_map[icat->fCatagory].push_back(std::string());
	  for(std::list<std::string>::const_iterator itext=text.begin();
	      itext!=text.end(); itext++)
	    cat_map[icat->fCatagory].push_back(*itext);
	}
    }

  std::map<std::string, std::list<std::list<OptionHelp>::const_iterator> > 
    cat_opt;

  unsigned max_option_length = 0;
  for(std::list<OptionHelp>::const_iterator iopt = fUsageInfo.begin();
      iopt != fUsageInfo.end(); iopt++)
    {
      if((print_only_simple_options)&&(!iopt->fIsSimple))continue;
      std::string cat = iopt->fCatagory;
      if(cat_map.find(cat) == cat_map.end())cat=std::string();
      cat_opt[cat].push_back(iopt);
      unsigned len = USAGE_INDENT+1+iopt->fOption.length()+USAGE_SPACING;
      if(iopt->fArgRequired)len += 3+iopt->fArgTypeName.length();
      if(len > max_option_length)max_option_length=len;
    }

  unsigned option_print_width = max_option_length;
  MAXIMIZE(option_print_width,USAGE_OPTION_MAX-USAGE_SPACING);

  bool first_cat = true;
  for(std::list<std::string>::const_iterator icat = cat_list.begin();
      icat != cat_list.end(); icat++)
    if(cat_opt.find(*icat) != cat_opt.end())
      {
	if(!first_cat)stream << '\n';
	else first_cat=false;

	if(!cat_map[*icat].empty())
	  {
	    for(std::list<std::string>::const_iterator iline = 
		  cat_map[*icat].begin(); iline != cat_map[*icat].end(); 
		iline++)
	      stream << *iline << '\n';
	    stream << '\n';
	  }

	for(std::list<std::list<OptionHelp>::const_iterator>::const_iterator
	      iiopt = cat_opt[*icat].begin(); iiopt != cat_opt[*icat].end();
	    iiopt++)
	  {
	    std::list<OptionHelp>::const_iterator iopt = *iiopt;

	    if((print_only_simple_options)&&(!iopt->fIsSimple))continue;

	    std::string linestart;
	    linestart = std::string(USAGE_INDENT,' ') + '-' + iopt->fOption;
	    if(iopt->fArgRequired)
	      linestart += 
		std::string("={") + iopt->fArgTypeName + std::string("}");

	    std::list<std::string> lines;
	    format(lines, iopt->fText, 
		   USAGE_WIDTH-option_print_width-USAGE_SPACING);

	    if(lines.empty())
	      {
		stream << linestart << '\n';
		continue;
	      }

	    if(linestart.length()>option_print_width)
	      {
		stream << linestart << '\n';
		linestart = std::string(option_print_width+USAGE_SPACING,' ');
	      }
	    else
	      linestart += 
		std::string(option_print_width
			    - linestart.length()+USAGE_SPACING,' ');

	    for(std::list<std::string>::const_iterator iline = lines.begin();
		iline != lines.end(); iline++)
	      {
		stream << linestart;
		if(iline == lines.begin())
		  linestart = 
		    std::string(option_print_width+USAGE_SPACING,' ');
		stream << *iline << '\n';
	      }
	  }  
      }
}

#ifdef TEST_MAIN
int main(int argc, char** argv)
{
  VSOptions options(argc, argv);
  bool b = false;
  double d = 0;
  std::string s = "";
  std::cout << options.findWithValue("double",d) << std::endl;
  std::cout << options.findWithValue("string",s) << std::endl;
  std::cout << options.findWithValue("bool",b) << std::endl;
  std::cout << options.find("bloop") << std::endl;
  std::cout << std::endl;
  std::cout << b << std::endl;
  std::cout << d << std::endl;
  std::cout << s << std::endl;
  for(int i=0;i<argc;i++){ if(i!=0)std::cout << ' '; std::cout << argv[i]; }
  std::cout << std::endl;
}
#endif
