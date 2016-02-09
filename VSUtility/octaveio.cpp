/*! \file octaveio.cpp

  Manipulate files written in Octave H5 format

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/15/2005
*/

#include <iostream>
#include<iterator>
#include<algorithm>
#include <string>

#include <VSOptions.hpp>
#include <VSOctaveIO.hpp>

using namespace VERITAS;

void ls(VSOctaveH5ReaderStruct* s, const std::string& name, 
	bool long_list, bool cat_ls, bool recurse,
	std::ostream& stream = std::cout, const std::string& prefix = "");
void ls(VSOctaveH5ReaderCellArray* c, 
	bool long_list, bool cat_ls, bool recurse,
	std::ostream& stream = std::cout, const std::string& prefix = "");
void cat(VSOctaveH5ReaderStruct* s, const std::string& name,
	 std::ostream& stream = std::cout);

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] filename [variable]" 
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);

  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  unsigned modes = 0;
  bool mode_ls = false;
  bool mode_cat = false;
  bool recursive = false;
  bool long_list = false;
  bool cat_ls = false;

  if(options.find("ls","List variable fields [default].")
     != VSOptions::FS_NOT_FOUND)mode_ls=true, modes++;
  if(options.find("cat","Dump data from variable.")
     != VSOptions::FS_NOT_FOUND)mode_cat=true, modes++;
  
  if(options.find("r","Recursive list.")
     != VSOptions::FS_NOT_FOUND)recursive=true;
  if(options.find("l","Long listing.")
     != VSOptions::FS_NOT_FOUND)long_list=true;
  if(options.find("v","Print values in listing.")
     != VSOptions::FS_NOT_FOUND)cat_ls=true;
  
  if(modes>1)
    {
      std::cerr << "Only one of \"-ls\" or \"-cat\" can be specified."
		<< std::endl;
      usage(progname,options,std::cerr);
      exit(EXIT_FAILURE);
    }

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  if(argc < 1)
    {
      std::cerr << progname << ": need at least 1 argument, got " 
		<< argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::vector< std::string > filenames;
  std::string variable;

  while(argc)
    {
      if(argc == 1 && filenames.size() != 0) variable = *argv;
      else filenames.push_back(*argv);
      argv++,argc--;
    }

  VSOctaveH5Reader* reader = 0;
  try
    {
      const unsigned nfile = filenames.size();
      for(unsigned ifile = 0; ifile < nfile; ifile++)
	{
	  reader = new VSOctaveH5Reader(filenames[ifile]);
	  if(mode_cat)cat(reader, variable, std::cout);
	  else ls(reader, variable, long_list, cat_ls, recursive, std::cout); 
	  delete reader, reader = NULL;
	}
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << e.message() << std::endl;
    }

  try
    {
      delete reader;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << e.message() << std::endl;
    }

  return EXIT_SUCCESS;
}

template<typename T> 
void cat_scalar(VSOctaveH5ReaderStruct* s, const std::string& name, 
		 std::ostream& stream)
{
  T x;
  if(s->readScalar(name,x))
    stream << VSDataConverter::toString(x,true) << '\n';
  else
    stream << name << ": could not read scalar" << std::endl;
}

template<typename T> 
void cat_scalar(VSOctaveH5ReaderCellArray* c, unsigned irow, unsigned icol,
		 std::ostream& stream)
{
  T x;
  if(c->readScalar(irow,icol,x))
    stream << VSDataConverter::toString(x,true) << '\n';
  else
    stream << irow << ',' << icol << ": could not read scalar" << std::endl;
}

template<typename T> 
void cat_vector(VSOctaveH5ReaderStruct* s, const std::string& name, 
		std::ostream& stream)
{
  std::vector<T> v;
  if(s->readVector(name,v))
    {
      typename std::vector<T>::const_iterator i = v.begin();
      while(i != v.end())
	{
	  stream << VSDataConverter::toString(*i,true) << '\n';
	  i++;
	}
    }
  else
    stream << name << ": could not read vector" << std::endl;
}

template<typename T> 
void cat_matrix(VSOctaveH5ReaderStruct* s, const std::string& name, 
		std::ostream& stream)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  T *m = 0;
  if(s->readMatrix(name,nrow,ncol,m))
    {
      for(unsigned irow=0;irow<nrow;irow++)
	{
	  if(ncol)
	    stream << VSDataConverter::toString(m[irow*ncol],true);
	  for(unsigned icol=1;icol<ncol;icol++)
	    stream << ' ' << VSDataConverter::toString(m[irow*ncol+icol],true);
	  stream << '\n';
	}
      delete[] m;
    }
  else
    stream << name << ": could not read matrix" << std::endl;
}


void ls(VSOctaveH5ReaderStruct* s, const std::string& name, 
	bool long_list, bool cat_ls, bool recurse,
	std::ostream& stream, const std::string& prefix)
{
  if(name.empty())
    {
      std::vector<std::string> fieldnames = s->variables();
      for(std::vector<std::string>::iterator ifield = fieldnames.begin();
	  ifield != fieldnames.end(); ifield++)
	{
	  std::string var;
	  if(!prefix.empty())var = prefix+std::string(".");
	  var += *ifield;
	  if(long_list)
	    {
	      unsigned vrows;
	      unsigned vcols;
	      s->dimensions(*ifield,vrows,vcols);
	      stream
		<< std::left
		<< std::setw(6) << vrows << ' '
		<< std::setw(6) << vcols << ' '
		<< std::setw(9) << s->variableType(*ifield) << ' '
		<< std::setw(7) << s->elementType(*ifield) << ' ';
	    }
	  if(cat_ls)
	    {
	      unsigned w = (var.length()/20+1)*20;
	      stream << std::left << std::setw(w) << var << ' ';
	      if(s->isString(*ifield))
		{
		  std::string str;
		  s->readString(*ifield,str);
		  stream << str << '\n';
		}
	      else if(s->isScalar(*ifield))
		{
		  std::string eltype = s->elementType(*ifield);
		  if(eltype == "bool")
		    cat_scalar<bool>(s, *ifield, stream);
		  else if(eltype == "uint8")
		    cat_scalar<uint8_t>(s, *ifield, stream);
		  else if(eltype == "uint16")
		    cat_scalar<uint16_t>(s, *ifield, stream);
		  else if(eltype == "uint32")
		    cat_scalar<uint32_t>(s, *ifield, stream);
		  else if(eltype == "uint64")
		    cat_scalar<uint64_t>(s, *ifield, stream);
		  else if(eltype == "int8")
		    cat_scalar<int8_t>(s, *ifield, stream);
		  else if(eltype == "int16")
		    cat_scalar<int16_t>(s, *ifield, stream);
		  else if(eltype == "int32")
		    cat_scalar<int32_t>(s, *ifield, stream);
		  else if(eltype == "int64")
		    cat_scalar<int64_t>(s, *ifield, stream);
		  else 
		    cat_scalar<double>(s, *ifield, stream);
		}
#if 0
	      else if(s->isVector(*ifield))
		{
		  std::string eltype = s->elementType(*ifield);
		  if(eltype == "bool")
		    cat_vector<bool>(s, *ifield, stream);
		  else if(eltype == "uint8")
		    cat_vector<uint8_t>(s, *ifield, stream);
		  else if(eltype == "uint16")
		    cat_vector<uint16_t>(s, *ifield, stream);
		  else if(eltype == "uint32")
		    cat_vector<uint32_t>(s, *ifield, stream);
		  else if(eltype == "uint64")
		    cat_vector<uint64_t>(s, *ifield, stream);
		  else if(eltype == "int8")
		    cat_vector<int8_t>(s, *ifield, stream);
		  else if(eltype == "int16")
		    cat_vector<int16_t>(s, *ifield, stream);
		  else if(eltype == "int32")
		    cat_vector<int32_t>(s, *ifield, stream);
		  else if(eltype == "int64")
		    cat_vector<int64_t>(s, *ifield, stream);
		  else 
		    cat_vector<double>(s, *ifield, stream);
		}
#endif
	      else
		{
		  stream << '\n';
		}
	    }
	  else
	    {
	      stream << std::left << var << '\n';
	    }
	  if(recurse)
	    {
	      if(s->isStruct(*ifield))
		{
		  VSOctaveH5ReaderStruct* rs = s->readStruct(*ifield);
		  ls(rs,"",long_list,cat_ls,recurse,stream,var);
		}
	      else if(s->isCellArray(*ifield))
		{
		  VSOctaveH5ReaderCellArray* rc = s->readCellArray(*ifield);
		  ls(rc,long_list,cat_ls,recurse,stream,var);
		}
	    }
	}
    }
  else
    {
      if(!s->isValid(name))
	{
	  stream << name << ": not found" << std::endl;
	  return;
	}

      std::string var;
      if(!prefix.empty())var = prefix+std::string(".");
      var += name;

      if(long_list)
	{
	  unsigned vrows;
	  unsigned vcols;
	  s->dimensions(name,vrows,vcols);
	  stream
	    << std::left
	    << std::setw(6) << vrows << ' '
	    << std::setw(6) << vcols << ' '
	    << std::setw(9) << s->variableType(name) << ' '
	    << std::setw(7) << s->elementType(name) << ' ';
	}
      stream << var << '\n';
      if(s->isStruct(name))
	{
	  VSOctaveH5ReaderStruct* rs = s->readStruct(name);
	  ls(rs,"",long_list,cat_ls,recurse,stream,var);
	}
      else if(s->isCellArray(name))
	{
	  VSOctaveH5ReaderCellArray* rc = s->readCellArray(name);
	  ls(rc,long_list,cat_ls,recurse,stream,var);
	}
    }
}

void ls(VSOctaveH5ReaderCellArray* c, 
	bool long_list, bool cat_ls, bool recurse, 
	std::ostream& stream, const std::string& prefix)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  c->dimensions(nrow,ncol);
  for(unsigned irow=0;irow<nrow;irow++)
    for(unsigned icol=0;icol<ncol;icol++)
      {
	std::ostringstream varstream;
	if(nrow==1)varstream << '{' << icol << '}';
	else if(ncol==1)varstream << '{' << irow << '}';
	else varstream << '{' << irow << ',' << icol << '}';
	std::string var;
	if(!prefix.empty())var = prefix;
	var += varstream.str();
	if(long_list)
	  {
	    unsigned vrows;
	    unsigned vcols;
	    c->dimensions(irow,icol,vrows,vcols);
	    stream
	      << std::left
	      << std::setw(6) << vrows << ' '
	      << std::setw(6) << vcols << ' '
	      << std::setw(9) << c->variableType(irow,icol) << ' '
	      << std::setw(7) << c->elementType(irow,icol) << ' ';
	  }
	if(cat_ls)
	  {
	    stream << std::left << std::setw(25) << var.c_str() << ' ';
	    if(c->isString(irow,icol))
	      {
		std::string str;
		c->readString(irow,icol,str);
		stream << str << '\n';
	      }
	    else if(c->isScalar(irow,icol))
	      {
		std::string eltype = c->elementType(irow,icol);
		if(eltype == "bool")
		  cat_scalar<bool>(c, irow, icol, stream);
		else if(eltype == "uint8")
		  cat_scalar<uint8_t>(c, irow, icol, stream);
		else if(eltype == "uint16")
		  cat_scalar<uint16_t>(c, irow, icol, stream);
		else if(eltype == "uint32")
		  cat_scalar<uint32_t>(c, irow, icol, stream);
		else if(eltype == "uint64")
		  cat_scalar<uint64_t>(c, irow, icol, stream);
		else if(eltype == "int8")
		  cat_scalar<int8_t>(c, irow, icol, stream);
		else if(eltype == "int16")
		  cat_scalar<int16_t>(c, irow, icol, stream);
		else if(eltype == "int32")
		  cat_scalar<int32_t>(c, irow, icol, stream);
		else if(eltype == "int64")
		  cat_scalar<int64_t>(c, irow, icol, stream);
		else 
		  cat_scalar<double>(c, irow, icol, stream);
	      }
#if 0
	    else if(s->isVector(name))
	      {
		std::string eltype = s->elementType(name);
		if(eltype == "bool")
		  cat_vector<bool>(c, irow, icol, stream);
		else if(eltype == "uint8")
		  cat_vector<uint8_t>(c, irow, icol, stream);
		else if(eltype == "uint16")
		  cat_vector<uint16_t>(c, irow, icol, stream);
		else if(eltype == "uint32")
		  cat_vector<uint32_t>(c, irow, icol, stream);
		else if(eltype == "uint64")
		  cat_vector<uint64_t>(c, irow, icol, stream);
		else if(eltype == "int8")
		  cat_vector<int8_t>(c, irow, icol, stream);
		else if(eltype == "int16")
		  cat_vector<int16_t>(c, irow, icol, stream);
		else if(eltype == "int32")
		  cat_vector<int32_t>(c, irow, icol, stream);
		else if(eltype == "int64")
		  cat_vector<int64_t>(c, irow, icol, stream);
		else 
		  cat_vector<double>(c, irow, icol, stream);
	      }
#endif
	    else
	      {
		stream << '\n';
	      }
	  }
	else
	  {
	    stream << std::left << var << '\n';
	  }
	if(recurse)
	  {
	    if(c->isStruct(irow,icol))
	      {
		VSOctaveH5ReaderStruct* rs = c->readStruct(irow,icol);
		ls(rs,"",long_list,cat_ls,recurse,stream,var);
	      }
	    else if(c->isCellArray(irow,icol))
	      {
		VSOctaveH5ReaderCellArray* rc = c->readCellArray(irow,icol);
		ls(rc,long_list,cat_ls,recurse,stream,var);
	      }
	  }
      }
}

void cat(VSOctaveH5ReaderStruct* s, const std::string& name,
	 std::ostream& stream)
{
  if(name.empty())
    stream << "No variable name given" << std::endl;    
  else if(!s->isValid(name))
    stream << name << ": not found" << std::endl;
  else if(s->isStruct(name))
    stream << name << ": is structure" << std::endl;
  else if(s->isCellArray(name))
    stream << name << ": is cell array" << std::endl;
  else if(s->isEmpty(name))
    stream << name << ": empty" << std::endl;
  else if(s->isString(name))
    {
      std::string str;
      s->readString(name,str);
      stream << str << '\n';
    }
  else if(s->isScalar(name))
    {
      std::string eltype = s->elementType(name);
      if(eltype == "bool")cat_scalar<bool>(s, name, stream);
      else if(eltype == "uint8")cat_scalar<uint8_t>(s, name, stream);
      else if(eltype == "uint16")cat_scalar<uint16_t>(s, name, stream);
      else if(eltype == "uint32")cat_scalar<uint32_t>(s, name, stream);
      else if(eltype == "uint64")cat_scalar<uint64_t>(s, name, stream);
      else if(eltype == "int8")cat_scalar<int8_t>(s, name, stream);
      else if(eltype == "int16")cat_scalar<int16_t>(s, name, stream);
      else if(eltype == "int32")cat_scalar<int32_t>(s, name, stream);
      else if(eltype == "int64")cat_scalar<int64_t>(s, name, stream);
      else cat_scalar<double>(s, name, stream);
    }
  else if(s->isVector(name))
    {
      std::string eltype = s->elementType(name);
      if(eltype == "bool")cat_vector<bool>(s, name, stream);
      else if(eltype == "uint8")cat_vector<uint8_t>(s, name, stream);
      else if(eltype == "uint16")cat_vector<uint16_t>(s, name, stream);
      else if(eltype == "uint32")cat_vector<uint32_t>(s, name, stream);
      else if(eltype == "uint64")cat_vector<uint64_t>(s, name, stream);
      else if(eltype == "int8")cat_vector<int8_t>(s, name, stream);
      else if(eltype == "int16")cat_vector<int16_t>(s, name, stream);
      else if(eltype == "int32")cat_vector<int32_t>(s, name, stream);
      else if(eltype == "int64")cat_vector<int64_t>(s, name, stream);
      else cat_vector<double>(s, name, stream);
    }
  else if(s->isMatrix(name))
    {
      std::string eltype = s->elementType(name);
      if(eltype == "bool")cat_matrix<bool>(s, name, stream);
      else if(eltype == "uint8")cat_matrix<uint8_t>(s, name, stream);
      else if(eltype == "uint16")cat_matrix<uint16_t>(s, name, stream);
      else if(eltype == "uint32")cat_matrix<uint32_t>(s, name, stream);
      else if(eltype == "uint64")cat_matrix<uint64_t>(s, name, stream);
      else if(eltype == "int8")cat_matrix<int8_t>(s, name, stream);
      else if(eltype == "int16")cat_matrix<int16_t>(s, name, stream);
      else if(eltype == "int32")cat_matrix<int32_t>(s, name, stream);
      else if(eltype == "int64")cat_matrix<int64_t>(s, name, stream);
      else cat_matrix<double>(s, name, stream);
    }
}
