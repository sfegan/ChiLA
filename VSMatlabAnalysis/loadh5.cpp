//-*-mode:c++; mode:font-lock;-*-

/*! \file loadh5.cpp

  Load Octave HDF5 file into MATLAB

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       09/21/2005
*/

// $Id: loadh5.cpp,v 1.4 2009/12/16 00:38:17 matthew Exp $

// Compile with: mex -cxx -I${HOME}/ChiLA/VSUtility -L${HOME}/ChiLA/VSUtility -lVSUtility -lhdf5 -lz loadh5.cpp

// Changes:
// *) MATLAB is column ordered, so must flip cols,rows - at present a matrix
//    will be imported transposed from the C++ code - SJF 2006-09-22
// *) If invoked with no output parameters the variables will be loaded into
//    the caller's address space - SJF 2006-09-22

// On gamma5: g++ -march=nocona -O3 -g -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -DUSEALLOCA -fPIC -I/usr/include/mysql -I../VSCommon -I../VSUtility -I.  -I/veritas/matlab/extern/include -DMATLAB_MEX_FILE -fPIC -fno-omit-frame-pointer -ansi -D_GNU_SOURCE -pthread -O -DNDEBUG -fexceptions -shared -Wl,--version-script,/veritas/matlab/extern/lib/glnx86/mexFunction.map -Wl,-rpath-link,/veritas/matlab/extern/../bin/glnx86 -L/veritas/matlab/extern/../bin/glnx86 -o loadh5.mexglx loadh5.cpp /veritas/matlab/extern/src/mexversion.c ../VSUtility/VSOctaveH5Reader.cpp  ../VSUtility/VSOctaveH5Writer.cpp ../VSUtility/VSOctaveIO.cpp ../VSUtility/VSFileUtility.cpp ../VSUtility/VSDataConverter.cpp   -lhdf5 -lz -lmx -lmex -lmat -lm

#include <string>
#include <stdint.h>

#include <mex.h>

#include <VSFileUtility.hpp>
#include <VSOctaveH5Reader.hpp>

using namespace VERITAS;

template<typename T> class ClassID { public: static mxClassID get() { return 0; } };

#define DEFCLASSID(T,x) \
template<> class ClassID<T> { public: static mxClassID get() { return x; } }

DEFCLASSID(int8_t, mxINT8_CLASS);
DEFCLASSID(uint8_t, mxUINT8_CLASS);
DEFCLASSID(int16_t, mxINT16_CLASS);
DEFCLASSID(uint16_t, mxUINT16_CLASS);
DEFCLASSID(int32_t, mxINT32_CLASS);
DEFCLASSID(uint32_t, mxUINT32_CLASS);
DEFCLASSID(int64_t, mxINT64_CLASS);
DEFCLASSID(uint64_t, mxUINT64_CLASS);
DEFCLASSID(float, mxSINGLE_CLASS);
DEFCLASSID(double, mxDOUBLE_CLASS);

static const char* newStrCpy(const std::string& s)
{
  char* cs = static_cast<char*>(mxCalloc(s.length()+1, sizeof(*cs)));
  strcpy(cs, s.c_str());
  return cs;
}

// ============================================================================
// CELL LOADERS FOR PRIMATIVES
// ============================================================================

mxArray* 
loadLogicalScalar(VSOctaveH5ReaderCellArray* c, unsigned row, unsigned col)
{
  unsigned val = 0;
  c->readScalar(row, col, val);
  return mxCreateLogicalScalar(val);
}

mxArray* 
loadLogicalMatrix(VSOctaveH5ReaderCellArray* c, unsigned row, unsigned col)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  c->dimensions(row,col,nrow,ncol);
  mxArray* mv = mxCreateLogicalMatrix(/*nrow,ncol*/ncol,nrow);
  c->readMatrix(row,col,static_cast<mxLogical*>(mxGetLogicals(mv)));
  return mv;
}

template<typename T> mxArray* 
loadRealScalar(VSOctaveH5ReaderCellArray* c, unsigned row, unsigned col)
{
  mxArray* mv = mxCreateNumericMatrix(1,1,ClassID<T>::get(),mxREAL);
  c->readScalar(row,col,*static_cast<T*>(mxGetData(mv)));
  return mv;
}

template<typename T> mxArray* 
loadRealMatrix(VSOctaveH5ReaderCellArray* c, unsigned row, unsigned col)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  c->dimensions(row,col,nrow,ncol);
  mxArray* mv = 
    mxCreateNumericMatrix(/*nrow,ncol*/ncol,nrow,ClassID<T>::get(),mxREAL);
  c->readMatrix(row,col,static_cast<T*>(mxGetData(mv)));
  return mv;
}

// ============================================================================
// STRUCT LOADERS FOR PRIMATIVES
// ============================================================================

mxArray* loadLogicalScalar(VSOctaveH5ReaderStruct* s, const std::string& fn)
{
  unsigned val = 0;
  s->readScalar(fn, val);
  return mxCreateLogicalScalar(val);
}

mxArray* loadLogicalMatrix(VSOctaveH5ReaderStruct* s, const std::string& fn)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  s->dimensions(fn,nrow,ncol);
  mxArray* mv = mxCreateLogicalMatrix(/*nrow,ncol*/ncol,nrow);
  s->readMatrix(fn,static_cast<mxLogical*>(mxGetLogicals(mv)));
  return mv;
}

template<typename T> mxArray* loadRealScalar(VSOctaveH5ReaderStruct* s, 
					     const std::string& fn)
{
  mxArray* mv = mxCreateNumericMatrix(1,1,ClassID<T>::get(),mxREAL);
  s->readScalar(fn,*static_cast<T*>(mxGetData(mv)));
  return mv;
}

template<typename T> mxArray* loadRealMatrix(VSOctaveH5ReaderStruct* s, 
					     const std::string& fn)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  s->dimensions(fn,nrow,ncol);
  mxArray* mv = 
    mxCreateNumericMatrix(/*nrow,ncol*/ncol,nrow,ClassID<T>::get(),mxREAL);
  s->readMatrix(fn,static_cast<T*>(mxGetData(mv)));
  return mv;
}

// ============================================================================
// CELL LOADER
// ============================================================================

mxArray* loadStruct(VSOctaveH5ReaderStruct* s,
		    bool load_to_matlab_workspace = false);

mxArray* loadCell(VSOctaveH5ReaderCellArray* c)
{
  unsigned nrow = 0;
  unsigned ncol = 0;
  c->dimensions(nrow,ncol);

  // Allocate the cell array
  mxArray* cv = mxCreateCellMatrix(nrow, ncol);
  
  for(unsigned irow=0;irow<nrow;irow++)
    for(unsigned icol=0;icol<ncol;icol++)
      {
	const unsigned ifield = icol*nrow+irow;
	mxArray* ev = 0;

	if(c->isStruct(irow, icol))
	{
	  VSOctaveH5ReaderStruct* rs = c->readStruct(irow, icol);
	  ev=loadStruct(rs);
	  delete rs;
	}
      else if(c->isCellArray(irow, icol))
	{
	  VSOctaveH5ReaderCellArray* rc = c->readCellArray(irow, icol);
	  ev=loadCell(rc);
	  delete rc;
	}
      else if(c->isString(irow, icol))
	{
	  std::string str;
	  c->readString(irow, icol, str);
	  ev=mxCreateString(str.c_str());
	}
      else if(c->isScalar(irow, icol))
	{
	  std::string eltype = c->elementType(irow, icol);
	  if(eltype == "bool")ev=loadLogicalScalar(c, irow, icol);
	  else if(eltype == "uint8")ev=loadRealScalar<uint8_t>(c, irow, icol);
	  else if(eltype == "uint16")ev=loadRealScalar<uint16_t>(c, irow,icol);
	  else if(eltype == "uint32")ev=loadRealScalar<uint32_t>(c, irow,icol);
	  else if(eltype == "uint64")ev=loadRealScalar<uint64_t>(c, irow,icol);
	  else if(eltype == "int8")ev=loadRealScalar<int8_t>(c, irow, icol);
	  else if(eltype == "int16")ev=loadRealScalar<int16_t>(c, irow, icol);
	  else if(eltype == "int32")ev=loadRealScalar<int32_t>(c, irow, icol);
	  else if(eltype == "int64")ev=loadRealScalar<int64_t>(c, irow, icol);
	  else ev=loadRealScalar<double>(c, irow, icol);
	}
      else if(c->isMatrix(irow, icol))
	{
	  std::string eltype = c->elementType(irow, icol);
	  if(eltype == "bool")ev=loadLogicalMatrix(c, irow, icol);
	  else if(eltype == "uint8")ev=loadRealMatrix<uint8_t>(c, irow,icol);
	  else if(eltype == "uint16")ev=loadRealMatrix<uint16_t>(c, irow,icol);
	  else if(eltype == "uint32")ev=loadRealMatrix<uint32_t>(c, irow,icol);
	  else if(eltype == "uint64")ev=loadRealMatrix<uint64_t>(c, irow,icol);
	  else if(eltype == "int8")ev=loadRealMatrix<int8_t>(c, irow, icol);
	  else if(eltype == "int16")ev=loadRealMatrix<int16_t>(c, irow, icol);
	  else if(eltype == "int32")ev=loadRealMatrix<int32_t>(c, irow, icol);
	  else if(eltype == "int64")ev=loadRealMatrix<int64_t>(c, irow, icol);
	  else ev=loadRealMatrix<double>(c, irow, icol);
	}
      else
	{
	  ev=mxCreateDoubleScalar(double(ifield));
	}

	/* Set each field in output structure */
	mxSetFieldByNumber(cv, 0, ifield, ev);
      }

  return cv;
}

// ============================================================================
// STRUCT LOADER
// ============================================================================

mxArray* loadStruct(VSOctaveH5ReaderStruct* s, 
		    bool load_to_matlab_workspace)
{
  std::vector<std::string> fieldnames = s->variables();
  unsigned nfield = fieldnames.size();

  mxArray* sv = 0;

  if(!load_to_matlab_workspace)
    {
      // Make copy of the field names into C-type string for Matlab
      const char** fnames = 
	static_cast<const char**>(mxCalloc(nfield, sizeof(*fnames)));
      for(unsigned ifield=0;ifield<nfield;ifield++)
	fnames[ifield] = newStrCpy(fieldnames[ifield]);
      
      // Allocate the struct
      sv = mxCreateStructMatrix(1, 1, nfield, fnames);
      if(sv==0)mexErrMsgTxt("Could not allocate structure.");
      
      // Delete C-type strings and array
      for(unsigned ifield=0;ifield<nfield;ifield++)
	mxFree(static_cast<void*>(const_cast<char*>(fnames[ifield])));
      mxFree(static_cast<void*>(fnames));
    }

  // Loop over fields
  for(unsigned ifield=0;ifield<nfield;ifield++)
    {
      const std::string& fn(fieldnames[ifield]);
      mxArray* ev = 0;

      if(s->isStruct(fn))
	{
	  VSOctaveH5ReaderStruct* rs = s->readStruct(fn);
	  ev=loadStruct(rs);	
	  delete rs;
	}
      else if(s->isCellArray(fn))
	{
	  VSOctaveH5ReaderCellArray* rc = s->readCellArray(fn);
	  ev=loadCell(rc);
	  delete rc;
	}
      else if(s->isString(fn))
	{
	  std::string str;
	  s->readString(fn, str);
	  ev=mxCreateString(str.c_str());
	}
      else if(s->isScalar(fn))
	{
	  std::string eltype = s->elementType(fn);
	  if(eltype == "bool")ev=loadLogicalScalar(s, fn);
	  else if(eltype == "uint8")ev=loadRealScalar<uint8_t>(s, fn);
	  else if(eltype == "uint16")ev=loadRealScalar<uint16_t>(s, fn);
	  else if(eltype == "uint32")ev=loadRealScalar<uint32_t>(s, fn);
	  else if(eltype == "uint64")ev=loadRealScalar<uint64_t>(s, fn);
	  else if(eltype == "int8")ev=loadRealScalar<int8_t>(s, fn);
	  else if(eltype == "int16")ev=loadRealScalar<int16_t>(s, fn);
	  else if(eltype == "int32")ev=loadRealScalar<int32_t>(s, fn);
	  else if(eltype == "int64")ev=loadRealScalar<int64_t>(s, fn);
	  else ev=loadRealScalar<double>(s, fn);
	}
      else if(s->isMatrix(fn))
	{
	  std::string eltype = s->elementType(fn);
	  if(eltype == "bool")ev=loadLogicalMatrix(s, fn);
	  else if(eltype == "uint8")ev=loadRealMatrix<uint8_t>(s, fn);
	  else if(eltype == "uint16")ev=loadRealMatrix<uint16_t>(s, fn);
	  else if(eltype == "uint32")ev=loadRealMatrix<uint32_t>(s, fn);
	  else if(eltype == "uint64")ev=loadRealMatrix<uint64_t>(s, fn);
	  else if(eltype == "int8")ev=loadRealMatrix<int8_t>(s, fn);
	  else if(eltype == "int16")ev=loadRealMatrix<int16_t>(s, fn);
	  else if(eltype == "int32")ev=loadRealMatrix<int32_t>(s, fn);
	  else if(eltype == "int64")ev=loadRealMatrix<int64_t>(s, fn);
	  else ev=loadRealMatrix<double>(s, fn);
	}
      else
	{
	  ev=mxCreateDoubleScalar(double(ifield));
	}

      /* Set each field in output structure */
      if(!load_to_matlab_workspace)
	mxSetFieldByNumber(sv, 0, ifield, ev);
      else
	mexPutVariable("caller", fn.c_str(), ev);
    }
    
  return sv;
}

/* The gateway routine.  */
extern "C" {

  void mexFunction(int nlhs, mxArray *plhs[],
		   int nrhs, const mxArray *prhs[])
  {
    int status = 0; // general usage later

    if (nrhs < 1)
      mexErrMsgTxt("One input required.");
    else if (nlhs > 1)
      mexErrMsgTxt("Too many output arguments.");
    else if (mxIsChar(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a string.");
    if (mxGetM(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a row vector.");
    
    /* Get the length of input string and copy it locally */
    unsigned filenamelen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
    char* filenamebuf = 
      static_cast<char*>(mxCalloc(filenamelen, sizeof(char)));
    status = mxGetString(prhs[0], filenamebuf, filenamelen);
    if (status != 0) 
      mexErrMsgTxt("Not enough space for filename string.");
    
    std::string filename(filenamebuf);
    VSFileUtility::expandFilename(filename);

    std::string path;

    if(nrhs == 2)
      {
	unsigned pathlen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
	char* pathbuf = static_cast<char*>(mxCalloc(pathlen, sizeof(char)));
	status = mxGetString(prhs[1], pathbuf, pathlen);
	if (status != 0) 
	  mexErrMsgTxt("Not enough space for path string.");

	path = std::string(pathbuf);
      }

    /* Construct the reader and load the data */
    VSOctaveH5Reader* reader = 0;
    try
      {
	reader = new VSOctaveH5Reader(filename);
	if(!path.empty() && !reader->isStruct(path))
	  {
	    mexErrMsgTxt(std::string("Bad path name: " + path).c_str());
	  }
	else if(!path.empty() && reader->isStruct(path))
	  {
	    VSOctaveH5ReaderStruct* s = reader->readStruct(path);
	    if(nlhs == 1)plhs[0] = loadStruct(s);
	    else loadStruct(s,true);
	    delete s;
	  }
	else
	  {
	    if(nlhs == 1)plhs[0] = loadStruct(reader);
	    else loadStruct(reader,true);
	  }
	delete reader;
      }
    catch(const VSOctaveH5Exception& e)
      {
	delete reader;
	mexErrMsgTxt(e.message().c_str());
      }
    H5close();
  }

} // extern "C"
