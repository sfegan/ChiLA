//-*-mode:c++; mode:font-lock;-*-

/*! \file VSH5Library.hpp

  Functions to implement library of objects indexed by a tree of "paths" in
  the octave H5 format 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/21/2007

  $Id: VSH5Library.hpp,v 1.8 2008/02/27 07:15:10 sfegan Exp $

*/

#ifndef VSH5LIBRARY_HPP
#define VSH5LIBRARY_HPP

#include<string>
#include<vector>
#include<map>
#include<vsassert>

#include<VSOctaveIO.hpp>

namespace VERITAS
{

  // ==========================================================================
  //
  // Classes relating to paths
  //
  // ==========================================================================

  class VSH5LibraryDefaultPath;
  class VSH5LibraryNamedPath;
  class VSH5LibraryIndexedPath;
  class VSH5LibraryRangedPath;

  class VSH5LibraryPathSet;

  class VSH5LibraryPath
  {
  public:
    VSH5LibraryPath(const std::string& _name):
      fName(_name) { vsassert(!_name.empty()); }

    virtual ~VSH5LibraryPath();

    virtual bool allowsDefault() const = 0;
    virtual std::string type() const = 0;
    const std::string& name() const { return fName; }

    virtual const VSH5LibraryDefaultPath* castToDefaultPath() const;
    virtual const VSH5LibraryNamedPath* castToNamedPath() const;
    virtual const VSH5LibraryIndexedPath* castToIndexedPath() const;
    virtual const VSH5LibraryRangedPath* castToRangedPath() const;

    virtual bool isValidForDefinition() const;
    virtual bool isValidForReading(const VSH5LibraryPath* val) const;
    virtual bool isValidForWriting(const VSH5LibraryPath* val) const;

    virtual VSH5LibraryPath* copy() const = 0;

    static VSH5LibraryPath* readPath(VSOctaveH5ReaderStruct* s);
    void writePath(VSOctaveH5WriterStruct* s) const;

    class Writer
    {
    public:
      Writer() { /* nothing to see here */ }
      virtual ~Writer();
      virtual VSOctaveH5WriterStruct* write(const VSH5LibraryPathSet* path)=0;
    };
    
    virtual VSH5LibraryPath::Writer* 
    getWriter(VSOctaveH5WriterStruct* s,
	      const VSH5LibraryPathSet* remainder) const = 0;

    virtual VSOctaveH5ReaderStruct* 
    read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const = 0;

    virtual VSH5LibraryPathSet* list(VSOctaveH5ReaderStruct* s) const = 0;

  protected:
    virtual void writePathDetails(VSOctaveH5WriterStruct* s) const;

    std::string fName;
  };

  class VSH5LibraryDefaultPath: public VSH5LibraryPath
  {
  public:
    VSH5LibraryDefaultPath(const std::string& _name):
      VSH5LibraryPath(_name) { /* nothing to see here */ }

    virtual ~VSH5LibraryDefaultPath();

    virtual bool allowsDefault() const;
    virtual std::string type() const;

    virtual const VSH5LibraryDefaultPath* castToDefaultPath() const;

    virtual bool isValidForDefinition() const;
    virtual bool isValidForReading(const VSH5LibraryPath* val) const;
    virtual bool isValidForWriting(const VSH5LibraryPath* val) const;

    virtual VSH5LibraryPath* copy() const;

    static std::string typeName() { return "default"; }

    virtual VSH5LibraryPath::Writer* 
    getWriter(VSOctaveH5WriterStruct* s,
	      const VSH5LibraryPathSet* remainder) const;

    virtual VSOctaveH5ReaderStruct* 
    read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const;

    virtual VSH5LibraryPathSet* list(VSOctaveH5ReaderStruct* s) const;
  };

  class VSH5LibraryNamedPath: public VSH5LibraryPath
  {
  public:
    VSH5LibraryNamedPath(const std::string& _name,
			 const std::string& named_value = ""):
      VSH5LibraryPath(_name), fNamedValue(named_value)
    { /* nothing to see here */ }

    const std::string& namedValue() const { return fNamedValue; }

    virtual ~VSH5LibraryNamedPath();

    virtual bool allowsDefault() const;
    virtual std::string type() const;

    virtual const VSH5LibraryNamedPath* castToNamedPath() const;

    virtual bool isValidForDefinition() const;
    virtual bool isValidForReading(const VSH5LibraryPath* val) const;
    virtual bool isValidForWriting(const VSH5LibraryPath* val) const;

    virtual VSH5LibraryPath* copy() const;

    static std::string typeName() { return "named"; }
    static VSH5LibraryNamedPath* readPath(const std::string& _name,
					VSOctaveH5ReaderStruct* s);

    class Writer: public VSH5LibraryPath::Writer
    {
    public:
      Writer(VSOctaveH5WriterStruct* s, const VSH5LibraryNamedPath* path,
	     const VSH5LibraryPathSet* remainder);
      virtual ~Writer();
      virtual VSOctaveH5WriterStruct* write(const VSH5LibraryPathSet* pathset);
    private:
      const VSH5LibraryNamedPath*                    fPath;
      const VSH5LibraryPathSet*                      fRemainder;
      VSOctaveH5WriterStruct*                        fBaseStruct;
      VSOctaveH5WriterStruct*                        fDataStruct;
      std::map<std::string,VSH5LibraryPath::Writer*> fData;
      VSH5LibraryPath::Writer*                       fDefault;
    };
    
    virtual VSH5LibraryPath::Writer* 
    getWriter(VSOctaveH5WriterStruct* s,
	      const VSH5LibraryPathSet* remainder) const;

    virtual VSOctaveH5ReaderStruct* 
    read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const;

    virtual VSH5LibraryPathSet* list(VSOctaveH5ReaderStruct* s) const;

  private:
    std::string fNamedValue;
  };

  class VSH5LibraryIndexedPath: public VSH5LibraryPath
  {
  public:
    VSH5LibraryIndexedPath(const std::string& _name,
			  unsigned indexed_value):
      VSH5LibraryPath(_name), fIndexedValue(indexed_value)
    { /* nothing to see here */ }

    unsigned indexedValue() const { return fIndexedValue; }

    virtual ~VSH5LibraryIndexedPath();

    virtual bool allowsDefault() const;
    virtual std::string type() const;

    virtual const VSH5LibraryIndexedPath* castToIndexedPath() const;

    virtual bool isValidForDefinition() const;
    virtual bool isValidForReading(const VSH5LibraryPath* val) const;
    virtual bool isValidForWriting(const VSH5LibraryPath* val) const;

    virtual VSH5LibraryPath* copy() const;

    static std::string typeName() { return "indexed"; }

    static VSH5LibraryIndexedPath* readPath(const std::string& _name,
					  VSOctaveH5ReaderStruct* s);

    class Writer: public VSH5LibraryPath::Writer
    {
    public:
      Writer(VSOctaveH5WriterStruct* s, const VSH5LibraryIndexedPath* path,
	     const VSH5LibraryPathSet* remainder);
      virtual ~Writer();
      virtual VSOctaveH5WriterStruct* write(const VSH5LibraryPathSet* pathset);
    private:
      const VSH5LibraryIndexedPath*                  fPath;
      const VSH5LibraryPathSet*                      fRemainder;
      VSOctaveH5WriterStruct*                        fBaseStruct;
      VSOctaveH5WriterCellVector*                    fDataCellVector;
      std::vector<VSH5LibraryPath::Writer*>          fData;
      VSH5LibraryPath::Writer*                       fDefault;
    };
    
    virtual VSH5LibraryPath::Writer* 
    getWriter(VSOctaveH5WriterStruct* s, 
	      const VSH5LibraryPathSet* remainder) const;

    virtual VSOctaveH5ReaderStruct* 
    read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const;

    virtual VSH5LibraryPathSet* list(VSOctaveH5ReaderStruct* s) const;

  private:
    virtual void writePathDetails(VSOctaveH5WriterStruct* s) const;

    unsigned fIndexedValue;
  };

  class VSH5LibraryRangedPath: public VSH5LibraryPath
  {
  public:
    VSH5LibraryRangedPath(const std::string& _name,
			  double ranged_value_lo, double ranged_value_hi):
      VSH5LibraryPath(_name), 
      fRangedValueLo(ranged_value_lo), fRangedValueHi(ranged_value_hi)
    { vsassert(ranged_value_lo <= ranged_value_hi); }

    VSH5LibraryRangedPath(const std::string& _name,
			 double ranged_value = 0):
      VSH5LibraryPath(_name), 
      fRangedValueLo(ranged_value), fRangedValueHi(ranged_value)
    { /* nothing to see here */ }

    double rangedValueLo() const { return fRangedValueLo; }
    double rangedValueHi() const { return fRangedValueHi; }

    virtual ~VSH5LibraryRangedPath();

    virtual bool allowsDefault() const;
    virtual std::string type() const;

    virtual const VSH5LibraryRangedPath* castToRangedPath() const;

    virtual bool isValidForDefinition() const;
    virtual bool isValidForReading(const VSH5LibraryPath* val) const;
    virtual bool isValidForWriting(const VSH5LibraryPath* val) const;

    virtual VSH5LibraryPath* copy() const;

    static std::string typeName() { return "ranged"; }
    static VSH5LibraryRangedPath* readPath(const std::string& _name,
					 VSOctaveH5ReaderStruct* s);

    class Writer: public VSH5LibraryPath::Writer
    {
    public:
      Writer(VSOctaveH5WriterStruct* s, const VSH5LibraryRangedPath* path,
	     const VSH5LibraryPathSet* remainder);
      virtual ~Writer();
      virtual VSOctaveH5WriterStruct* write(const VSH5LibraryPathSet* pathset);
      typedef std::pair<double,double> key_type;
    private:
      const VSH5LibraryRangedPath*                   fPath;
      const VSH5LibraryPathSet*                      fRemainder;
      VSOctaveH5WriterStruct*                        fBaseStruct;
      VSOctaveH5WriterExpandableCellVector*          fDataECellVector;
      VSOctaveH5WriterCompositeVector<key_type>*     fKeyList;
      std::map<key_type,VSH5LibraryPath::Writer*>    fData;
      VSH5LibraryPath::Writer*                       fDefault;
    };
    
    virtual VSH5LibraryPath::Writer* 
    getWriter(VSOctaveH5WriterStruct* s,
	      const VSH5LibraryPathSet* remainder) const;

    virtual VSOctaveH5ReaderStruct* 
    read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const;

    virtual VSH5LibraryPathSet* list(VSOctaveH5ReaderStruct* s) const;

  public:
    double fRangedValueLo;
    double fRangedValueHi;
  };

  // ==========================================================================
  //
  // Classes relating to path set
  //
  // ==========================================================================

  class VSH5LibraryPathSet
  {
  public:
    VSH5LibraryPathSet(unsigned path_count): 
      fPaths(path_count) { vsassert(path_count); }

    VSH5LibraryPathSet(const VSH5LibraryPathSet& o);
    VSH5LibraryPathSet& operator=(const VSH5LibraryPathSet& o);

    ~VSH5LibraryPathSet();

    unsigned size() const { return fPaths.size(); }
    bool empty() const { return fPaths.empty(); }

    void resize(unsigned path_count);

    bool validate() const;
    bool isValidForDefinition() const;

    void set(unsigned ipath, const VSH5LibraryPath* path)
    { delete fPaths[ipath]; fPaths[ipath] = path; }

    typedef std::vector<const VSH5LibraryPath*>::const_iterator const_iterator;
    const_iterator begin() const { return fPaths.begin(); }
    const_iterator end() const { return fPaths.end(); }

    const VSH5LibraryPath*const& operator[](unsigned ipath) const 
    { return fPaths[ipath]; }

    VSH5LibraryPathSet* copy() const;
    VSH5LibraryPathSet* shift() const;

    void writePathSet(VSOctaveH5WriterStruct* s, const std::string& name) 
      const;

    static VSH5LibraryPathSet* readPathSet(VSOctaveH5ReaderStruct* s, 
					   const std::string& name);

  private:
    std::vector<const VSH5LibraryPath*> fPaths;
  };

  // ==========================================================================
  //
  // Classes relating to library
  //
  // ==========================================================================

  class VSH5LibraryWriter
  {
  public:
    VSH5LibraryWriter(VSOctaveH5WriterStruct* base,
		      const VSH5LibraryPathSet& paths);

    ~VSH5LibraryWriter();

    VSOctaveH5WriterStruct* write(const VSH5LibraryPathSet& pathset);

  private:
    VSOctaveH5WriterStruct*  fBase;
    VSH5LibraryPathSet*      fPaths;
    VSH5LibraryPath::Writer* fWriter;
  };

  class VSH5LibraryReader
  {
  public:
    VSH5LibraryReader(VSOctaveH5ReaderStruct* base);
    ~VSH5LibraryReader();

    VSOctaveH5ReaderStruct* read(const VSH5LibraryPathSet& pathset) const;

    VSH5LibraryPathSet* list() const; // List the root directory entries
    VSH5LibraryPathSet* list(const VSH5LibraryPathSet& pathset) const;
    
  private:
    VSOctaveH5ReaderStruct*  fBase;
    VSH5LibraryPathSet*      fPaths;
  };

}

#endif // defined VSH5LIBRARY_HPP
