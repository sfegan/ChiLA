//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveH5Writer.hpp

  Classes to write to structured GNU/Octave HDF5 file. Data can be
  read directly into octave.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/20/2006
*/

#include"VSOctaveIO.hpp"

#ifndef VSOCTAVEH5WRITER_HPP
#define VSOCTAVEH5WRITER_HPP

#include<string>
#include<vector>
#include<set>
#include<map>
#include<vsassert>

#include<hdf5.h>

namespace VERITAS
{
  
  // ==========================================================================
  // Octave H5 writer base class
  // ==========================================================================

  class VSOctaveH5WriterStruct;
  class VSOctaveH5WriterCellVector;
  class VSOctaveH5WriterExpandableCellVector;

  class VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterBase(VSOctaveH5WriterBase* parent):
      m_parent(parent), m_interested_parties(), m_children()
    { /* nothing to see here */ }
    virtual ~VSOctaveH5WriterBase();
    virtual void flush();
    void registerInterestedParty(VSOctaveH5WriterBase* x);

    static unsigned defaultChunk() { return 1024; }
    static unsigned defaultCache() { return defaultChunk(); }
    static bool defaultCompress() { return s_default_compress; }
    static bool setDefaultCompress(bool c = true) 
    { bool old=s_default_compress; s_default_compress=c; return old; }
    
  protected:
    template<typename T>
    void makeAttribute(hid_t id, const std::string& name, const T& value);
    void makeEmptyMatrix(hid_t gid, unsigned rows, unsigned cols);

    virtual void notifyOfChildDestruction(VSOctaveH5WriterBase* child);

    VSOctaveH5WriterBase*           m_parent;
    std::set<VSOctaveH5WriterBase*> m_interested_parties;
    std::set<VSOctaveH5WriterBase*> m_children;
  private:
    VSOctaveH5WriterBase(const VSOctaveH5WriterBase&);
    VSOctaveH5WriterBase& operator=(const VSOctaveH5WriterBase&);

    static bool s_default_compress;
  };

  // ==========================================================================
  // Octave H5 writer helpers
  // ==========================================================================

  template<typename T> class VSOctaveH5WriterLowLevelHelper
  {
  public:
    static inline herr_t
    write_scalar(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
		 const T& x);
    static inline herr_t
    write_vector(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
		 const std::vector<T>& v);
    static inline herr_t
    write_matrix(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
		 unsigned n, const T* m);
  };

  template<> class VSOctaveH5WriterLowLevelHelper<bool>
  {
  public:
    static inline herr_t
    write_scalar(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
		 const bool& x);
    static inline herr_t
    write_vector(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
		 const std::vector<bool>& v);
    static inline herr_t
    write_matrix(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
		 unsigned n, const bool* m);
  };

  template<typename T> class VSOctaveH5WriterHelper
  {
  public:
    static inline bool write(VSOctaveH5WriterStruct* s,
			     const std::string& name, const T& x);
  };
  
  template<> class VSOctaveH5WriterHelper<std::string>
  {
  public:
    static inline bool write(VSOctaveH5WriterStruct* s,
			     const std::string& name, const std::string& x);
  };
  
  template<typename T> class VSOctaveH5WriterHelper<std::vector<T> >
  {
  public:
    static inline bool write(VSOctaveH5WriterStruct* s,
			     const std::string& name, const std::vector<T>& x);
  };

  template<> class VSOctaveH5WriterHelper<std::vector<std::string> >
  {
  public:
    static inline bool write(VSOctaveH5WriterStruct* s,
			     const std::string& name, 
			     const std::vector<std::string>& x);
  };

  // ==========================================================================
  // Octave H5 writer expandable vector
  // ==========================================================================

  template<typename T> 
  class VSOctaveH5WriterVector: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterVector(hid_t gid, VSOctaveH5WriterBase* parent,
			   unsigned cachesize = defaultCache(), 
			   unsigned chunksize = defaultChunk(),
			   bool compress = defaultCompress());
    virtual ~VSOctaveH5WriterVector();
    virtual void flush();
    bool append(const T& x);
    unsigned rows() const { return m_rows; }
  private:
    VSOctaveH5WriterVector(const VSOctaveH5WriterVector&);
    VSOctaveH5WriterVector& operator=(const VSOctaveH5WriterVector&);

    bool write(const T* ptr, unsigned n);

    hid_t      m_gid;
    hid_t      m_tid;
    hid_t      m_did;
    unsigned   m_rows;
    unsigned   m_cachesize;
    unsigned   m_cacheoccupancy;
    T*         m_cache;
    unsigned   m_chunksize;
    bool       m_compress;
  };
  
  // ==========================================================================
  // Octave H5 writer expandable table (matrix)
  // ==========================================================================

  template<typename T> 
  class VSOctaveH5WriterTable: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterTable(hid_t gid, unsigned cols,
			  VSOctaveH5WriterBase* parent, 
			  unsigned cachesize = defaultCache(), 
			  unsigned chunksize = defaultChunk(),
			  bool compress = defaultCompress());
    virtual ~VSOctaveH5WriterTable();
    bool append(const std::vector<T>& row);
    unsigned rows() const { return m_rows; }
  private:
    VSOctaveH5WriterTable(const VSOctaveH5WriterTable&);
    VSOctaveH5WriterTable& operator=(const VSOctaveH5WriterTable&);

    hid_t m_gid;
    hid_t m_tid;
    hid_t m_did;
    unsigned m_cols;
    unsigned m_rows;
    unsigned m_cachesize;
    unsigned m_chunksize;
    bool     m_compress;
  };

  // ==========================================================================
  // Octave H5 writer expandable vector for composites
  // ==========================================================================

  class VSOH5WCVElementWriter
  {
  public:
    virtual ~VSOH5WCVElementWriter();
    virtual bool write(const void* data,
		       const VSOctaveH5CompositeDefinition::Member& field) = 0;
  };

  template<typename T> class VSOH5WCVEWTemplate: public VSOH5WCVElementWriter
  {
  public:
    VSOH5WCVEWTemplate(VSOctaveH5WriterVector<T>* w): m_writer(w) { }
    virtual ~VSOH5WCVEWTemplate() {  }
    virtual bool write(const void* data,
		       const VSOctaveH5CompositeDefinition::Member& field)
    {
      return 
	m_writer->append(*((const T*)((const char*)(data) + field.offset))); 
    }
  private:
    VSOH5WCVEWTemplate(const VSOH5WCVEWTemplate&);
    VSOH5WCVEWTemplate& operator=(const VSOH5WCVEWTemplate&);

    VSOctaveH5WriterVector<T>* m_writer;
  };

  class VSOH5WCVEWString: public VSOH5WCVElementWriter
  {
  public:
    VSOH5WCVEWString(VSOctaveH5WriterExpandableCellVector* ecv): m_ecv(ecv) { }
    virtual ~VSOH5WCVEWString();
    virtual bool write(const void* data,
		       const VSOctaveH5CompositeDefinition::Member& field);
  private:
    VSOH5WCVEWString(const VSOH5WCVEWString&);
    VSOH5WCVEWString& operator=(const VSOH5WCVEWString&);

    VSOctaveH5WriterExpandableCellVector* m_ecv;
  };

  template<typename T> 
  class VSOctaveH5WriterCompositeVector: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterCompositeVector(VSOctaveH5WriterBase* parent,
				    VSOctaveH5WriterStruct* s,
				    bool struct_is_my_child,
				    const std::string& prefix = "",
				    unsigned cachesize = defaultCache(), 
				    unsigned chunksize = defaultChunk(),
				    bool compress = defaultCompress());
    virtual ~VSOctaveH5WriterCompositeVector();

    virtual void flush();
    bool append(const T& x);
    unsigned rows() const { return m_rows; }

  protected:
    virtual void notifyOfChildDestruction(VSOctaveH5WriterBase* child);

  private:
    VSOctaveH5WriterCompositeVector(const VSOctaveH5WriterCompositeVector&);
    VSOctaveH5WriterCompositeVector& 
    operator=(const VSOctaveH5WriterCompositeVector&);

    struct MemWriter
    {
      MemWriter(const VSOctaveH5CompositeDefinition::Member& _field,
		VSOH5WCVElementWriter* _writer = 0,
		VSOctaveH5WriterBase* _actual_writer = 0):
	field(_field), writer(_writer), actual_writer(_actual_writer) { }
      MemWriter(const MemWriter& o):
	field(o.field), writer(o.writer), actual_writer(o.actual_writer) { }
      MemWriter& operator=(const MemWriter& o)
      {
	field         = o.field;
	writer        = o.writer;
	actual_writer = o.actual_writer;
	return *this;
      }

      VSOctaveH5CompositeDefinition::Member field;
      VSOH5WCVElementWriter*                writer;
      VSOctaveH5WriterBase*                 actual_writer;
    };

    std::vector<MemWriter>   m_members;
    unsigned                 m_rows;
  };

  // ==========================================================================
  // Octave H5 writer cell array - just an interface to struct
  // ==========================================================================

  class VSOctaveH5WriterCellArray: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterCellArray(hid_t gid, unsigned rows, unsigned cols,
			      VSOctaveH5WriterBase* parent);
    virtual ~VSOctaveH5WriterCellArray();
    virtual void flush();

    bool writeEmptyMatrix(unsigned row, unsigned col);

    bool writeString(unsigned row, unsigned col, const std::string& s);
    
    template<typename T>
    bool writeScalar(unsigned row, unsigned col, const T& x);
    
    template<typename T>
    bool writeVector(unsigned row, unsigned col, 
		     const std::vector<T>& v, 
		     VSOctaveH5VecOrientation orient = VO_ROW);
    
    template<typename T>
    bool writeMatrix(unsigned row, unsigned col, 
		     unsigned rows, unsigned cols, const T* m);
    
    template<typename T>
    bool write(unsigned row, unsigned col, const T& x);

    VSOctaveH5WriterStruct* writeStruct(unsigned row, unsigned col);
    
    VSOctaveH5WriterCellArray* 
    writeCellArray(unsigned row, unsigned col, unsigned rows, unsigned cols);

    VSOctaveH5WriterCellVector* 
    writeCellVector(unsigned row, unsigned col, unsigned nel,
		    VSOctaveH5VecOrientation orient = VO_ROW);

    VSOctaveH5WriterExpandableCellVector* 
    writeExpandableCellVector(unsigned row, unsigned col,
			      VSOctaveH5VecOrientation orient = VO_ROW);

    template<typename T> VSOctaveH5WriterVector<T>* 
    writeExpandableVector(unsigned row, unsigned col, 
			  unsigned cachesize = defaultCache(), 
			  unsigned chunksize = defaultChunk(),
			  bool compress = defaultCompress());
    
    template<typename T> VSOctaveH5WriterTable<T>* 
    writeExpandableTable(unsigned row, unsigned col, unsigned cols,
			 unsigned cachesize = defaultCache(), 
			 unsigned chunksize = defaultChunk(),
			 bool compress = defaultCompress());

    template<typename T> bool
    writeComposite(unsigned row, unsigned col, const T& x);

    template<typename T> bool 
    writeCompositeVector(unsigned row, unsigned col,
			 const std::vector<T>& x);

    template<typename T> VSOctaveH5WriterCompositeVector<T>*
    writeCompositeExpandableVector(unsigned row, unsigned col,
				   unsigned cachesize = defaultCache(), 
				   unsigned chunksize = defaultChunk(),
				   bool compress = defaultCompress());

  private:
    VSOctaveH5WriterCellArray(const VSOctaveH5WriterCellArray&);
    VSOctaveH5WriterCellArray& operator=(const VSOctaveH5WriterCellArray&);

    std::string name(unsigned row, unsigned col);

    VSOctaveH5WriterStruct* m_struct;
    unsigned m_rows;
    unsigned m_cols;
    unsigned m_name_precision;
    std::vector<bool> m_written;
  };

  // ==========================================================================
  // Octave H5 writer cell vector - just an interface to struct
  // ==========================================================================

  class VSOctaveH5WriterCellVector: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterCellVector(hid_t gid, unsigned nel,
			       VSOctaveH5WriterBase* parent);
    virtual ~VSOctaveH5WriterCellVector();
    virtual void flush();

    bool writeEmptyMatrix(unsigned iel);

    bool writeString(unsigned iel, const std::string& s);
    
    template<typename T> bool writeScalar(unsigned iel, const T& x);
    
    template<typename T>
    bool writeVector(unsigned iel, const std::vector<T>& v, 
		     VSOctaveH5VecOrientation orient = VO_ROW);
    
    template<typename T>
    bool writeMatrix(unsigned iel, unsigned rows, unsigned cols, const T* m);
    
    template<typename T>
    bool write(unsigned iel, const T& x);

    VSOctaveH5WriterStruct* writeStruct(unsigned iel);
    
    VSOctaveH5WriterCellArray* 
    writeCellArray(unsigned iel, unsigned rows, unsigned cols);

    VSOctaveH5WriterCellVector* 
    writeCellVector(unsigned iel, unsigned nel,
		    VSOctaveH5VecOrientation orient = VO_ROW);

    VSOctaveH5WriterExpandableCellVector* 
    writeExpandableCellVector(unsigned iel, 
			      VSOctaveH5VecOrientation orient = VO_ROW);

    template<typename T> VSOctaveH5WriterVector<T>* 
    writeExpandableVector(unsigned iel, 
			  unsigned cachesize = defaultCache(), 
			  unsigned chunksize = defaultChunk(),
			  bool compress = defaultCompress());
    
    template<typename T> VSOctaveH5WriterTable<T>* 
    writeExpandableTable(unsigned iel, unsigned cols,
			 unsigned cachesize = defaultCache(), 
			 unsigned chunksize = defaultChunk(),
			 bool compress = defaultCompress());

    template<typename T> bool
    writeComposite(unsigned iel, const T& x);

    template<typename T> bool 
    writeCompositeVector(unsigned iel, const std::vector<T>& x);

    template<typename T> VSOctaveH5WriterCompositeVector<T>*
    writeCompositeExpandableVector(unsigned iel,
				   unsigned cachesize = defaultCache(), 
				   unsigned chunksize = defaultChunk(),
				   bool compress = defaultCompress());

  private:
    VSOctaveH5WriterCellVector(const VSOctaveH5WriterCellVector&);
    VSOctaveH5WriterCellVector& operator=(const VSOctaveH5WriterCellVector&);

    std::string name(unsigned iel);

    VSOctaveH5WriterStruct* m_struct;
    unsigned m_nel;
    unsigned m_name_precision;
    std::vector<bool> m_written;
  };

  // ==========================================================================
  // Octave H5 writer expandable cell vector - just an interface to struct
  // ==========================================================================

  class VSOctaveH5WriterExpandableCellVector: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterExpandableCellVector(hid_t gid, 
					 VSOctaveH5WriterBase* parent,
					 VSOctaveH5VecOrientation orient);
    virtual ~VSOctaveH5WriterExpandableCellVector();
    virtual void flush();

    bool appendEmptyMatrix();

    bool appendString(const std::string& s);
    
    template<typename T> bool appendScalar(const T& x);
    
    template<typename T>
    bool appendVector(const std::vector<T>& v, 
		      VSOctaveH5VecOrientation orient = VO_ROW);
    
    template<typename T>
    bool appendMatrix(unsigned rows, unsigned cols, const T* m);
    
    template<typename T>
    bool append(const T& x);

    VSOctaveH5WriterStruct* appendStruct();
    
    VSOctaveH5WriterCellArray* appendCellArray(unsigned rows, unsigned cols);

    VSOctaveH5WriterCellVector* 
    appendCellVector(unsigned nel, VSOctaveH5VecOrientation orient = VO_ROW);

    VSOctaveH5WriterExpandableCellVector* 
    appendExpandableCellVector(VSOctaveH5VecOrientation orient = VO_ROW);

    template<typename T> VSOctaveH5WriterVector<T>* 
    appendExpandableVector(unsigned cachesize = defaultCache(), 
			   unsigned chunksize = defaultChunk(),
			   bool compress = defaultCompress());
    
    template<typename T> VSOctaveH5WriterTable<T>* 
    appendExpandableTable(unsigned cols,
			  unsigned cachesize = defaultCache(), 
			  unsigned chunksize = defaultChunk(),
			  bool compress = defaultCompress());

    template<typename T> bool appendComposite(const T& x);

    template<typename T> bool appendCompositeVector(const std::vector<T>& x);

    template<typename T> VSOctaveH5WriterCompositeVector<T>*
    appendCompositeExpandableVector(unsigned cachesize = defaultCache(), 
				    unsigned chunksize = defaultChunk(),
				    bool compress = defaultCompress());

  private:
    VSOctaveH5WriterExpandableCellVector
    (const VSOctaveH5WriterExpandableCellVector&);
    
    VSOctaveH5WriterExpandableCellVector& operator=
    (const VSOctaveH5WriterExpandableCellVector&);

    std::string name();

    unsigned                  m_nel;
    VSOctaveH5VecOrientation  m_orient;
    hid_t                     m_gid;
    hid_t                     m_value_gid;
    VSOctaveH5WriterStruct*   m_struct;
  };

  // ==========================================================================
  // Octave H5 writer struct - this is where most of the action is
  // ==========================================================================

  class VSOctaveH5WriterStruct: public VSOctaveH5WriterBase
  {
  public:
    VSOctaveH5WriterStruct(const std::string& filename, bool overwrite=false,
			   const std::string& comment="");
    VSOctaveH5WriterStruct(hid_t gid, VSOctaveH5WriterBase* parent);
    virtual ~VSOctaveH5WriterStruct();
    virtual void flush();

    bool writeEmptyMatrix(const std::string& name);

    bool writeString(const std::string& name, const std::string& s);
    
    template<typename T>
    bool writeScalar(const std::string& name, const T& x);
    
    template<typename T>
    bool writeVector(const std::string& name, 
		     const std::vector<T>& v, 
		     VSOctaveH5VecOrientation orient = VO_ROW);
    
    template<typename T>
    bool writeMatrix(const std::string& name, 
		     unsigned rows, unsigned cols, const T* m);
    
    template<typename T>
    bool write(const std::string& name, const T& x);

    VSOctaveH5WriterStruct* writeStruct(const std::string& name);
    
    VSOctaveH5WriterCellArray* 
    writeCellArray(const std::string& name, unsigned rows, unsigned cols);

    template<typename T> bool
    writeStructCellVector(const std::string& name, const std::vector<T>& x);
			   
    VSOctaveH5WriterCellVector* 
    writeCellVector(const std::string& name, unsigned nel,
		    VSOctaveH5VecOrientation orient = VO_ROW);

    VSOctaveH5WriterExpandableCellVector* 
    writeExpandableCellVector(const std::string& name,
			      VSOctaveH5VecOrientation orient = VO_ROW);

    template<typename T> VSOctaveH5WriterVector<T>* 
    writeExpandableVector(const std::string& name,
			  unsigned cachesize = defaultCache(), 
			  unsigned chunksize = defaultChunk(),
			  bool compress = defaultCompress());
    
    template<typename T> VSOctaveH5WriterTable<T>* 
    writeExpandableTable(const std::string& name, unsigned cols, 
			 unsigned cachesize = defaultCache(), 
			 unsigned chunksize = defaultChunk(),
			 bool compress = defaultCompress());

    template<typename T> bool
    writeCompositeHere(const T& x, const std::string& prefix = "");

    template<typename T> bool 
    writeCompositeVectorHere(const std::vector<T>& x,
			     const std::string& prefix = "");

    template<typename T> VSOctaveH5WriterCompositeVector<T>*
    writeCompositeExpandableVectorHere(const std::string& prefix = "",
				       unsigned cachesize = defaultCache(), 
				       unsigned chunksize = defaultChunk(),
				       bool compress = defaultCompress());

    template<typename T> bool
    writeComposite(const std::string& name, const T& x);

    template<typename T> bool
    writeCompositeVector(const std::string& name, const std::vector<T>& x);

    template<typename T> VSOctaveH5WriterCompositeVector<T>*
    writeCompositeExpandableVector(const std::string& name,
				   unsigned cachesize = defaultCache(), 
				   unsigned chunksize = defaultChunk(),
				   bool compress = defaultCompress());

  protected:
    virtual void notifyOfChildDestruction(VSOctaveH5WriterBase* child);

  private:
    VSOctaveH5WriterStruct(const VSOctaveH5WriterStruct&);
    VSOctaveH5WriterStruct& operator=(const VSOctaveH5WriterStruct&);

    hid_t makeVariable(const std::string& name, 
		       const std::string& type, const std::string& eltype="");
    void registerChild(VSOctaveH5WriterBase* child, hid_t gid);

    bool m_my_file;
    hid_t m_gid;
    std::set<std::string> m_vars;
    std::map<VSOctaveH5WriterBase*, hid_t> m_gid_map;
  };

  typedef VSOctaveH5WriterStruct VSOctaveH5Writer;

  // **************************************************************************
  // **************************************************************************

  // FUNCTION DEFINITIONS

  // **************************************************************************
  // **************************************************************************

  // ==========================================================================
  // BASE
  // ==========================================================================

  template<typename T> void VSOctaveH5WriterBase::
  makeAttribute(hid_t id, const std::string& name, const T& value)
  {
    hid_t tid = VSOctaveH5Type<T>::h5_type();
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate(id, name.c_str(), tid, sid, H5P_DEFAULT);
    H5Awrite(aid, tid, &value);
    H5Aclose(aid);
    H5Sclose(sid);
    H5Tclose(tid);
  }

  // ==========================================================================
  // HELPER
  // ==========================================================================

  template<typename T> inline herr_t VSOctaveH5WriterLowLevelHelper<T>::
  write_scalar(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	       const T& x)
  {
    return H5Dwrite(did, tid, core_sid, disk_sid, H5P_DEFAULT, &x);
  }

  template<typename T> inline herr_t VSOctaveH5WriterLowLevelHelper<T>::
  write_vector(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	       const std::vector<T>& v)
  {
    return H5Dwrite(did, tid, core_sid, disk_sid, H5P_DEFAULT, &v[0]);
  }

  template<typename T> inline herr_t VSOctaveH5WriterLowLevelHelper<T>::
  write_matrix(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	       unsigned n, const T* m)
  {
    return H5Dwrite(did, tid, core_sid, disk_sid, H5P_DEFAULT, m);
  }

  herr_t inline VSOctaveH5WriterLowLevelHelper<bool>::
  write_scalar(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	       const bool& x)
  {
    unsigned y = x?1:0;
    return H5Dwrite(did, tid, core_sid, disk_sid, H5P_DEFAULT, &y);
  }

  herr_t inline VSOctaveH5WriterLowLevelHelper<bool>::
  write_vector(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	       const std::vector<bool>& v)
  {
    unsigned nv = v.size();
    std::vector<unsigned> u(nv);
    for(unsigned iv=0;iv<nv;iv++)u[iv]=v[iv]?1:0;
    return H5Dwrite(did, tid, core_sid, disk_sid, H5P_DEFAULT, &u[0]);
  }

  herr_t inline VSOctaveH5WriterLowLevelHelper<bool>::
  write_matrix(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	       unsigned n, const bool* m)
  {
    std::vector<unsigned> v(n);
    for(unsigned im=0;im<n;im++)v[im]=m[im]?1:0;
    return H5Dwrite(did, tid, core_sid, disk_sid, H5P_DEFAULT, &v[0]);
  }

  template<typename T> bool inline VSOctaveH5WriterHelper<T>::
  write(VSOctaveH5WriterStruct* s, const std::string& name, const T& x)
  {
    return s->writeScalar(name,x);
  };
  
  bool inline VSOctaveH5WriterHelper<std::string>::
  write(VSOctaveH5WriterStruct* s,
	const std::string& name, const std::string& x)
  {
    return s->writeString(name,x);
  };
  
  template<typename T> bool inline VSOctaveH5WriterHelper<std::vector<T> >::
  write(VSOctaveH5WriterStruct* s, 
	const std::string& name, const std::vector<T>& x)
  {
    return s->writeVector(name,x);
  }

  bool inline VSOctaveH5WriterHelper<std::vector<std::string> >::
  write(VSOctaveH5WriterStruct* s, 
	const std::string& name, const std::vector<std::string>& x)
  {
    unsigned nstring = x.size();
    VSOctaveH5WriterCellVector* c = s->writeCellVector(name,nstring);
    if(c==0)return false;
    for(unsigned istring=0;istring<nstring;istring++)
      if(!c->writeString(istring,x[istring]))
	{
	  delete c;
	  return false;
	}
    delete c;
    return true;
  }

  // ==========================================================================
  // VECTOR
  // ==========================================================================
  
  template<typename T> VSOctaveH5WriterVector<T>::
  VSOctaveH5WriterVector(hid_t gid, VSOctaveH5WriterBase* parent,
			 unsigned cachesize, unsigned chunksize, bool compress)
    : VSOctaveH5WriterBase(parent),
      m_gid(gid), m_tid(), m_did(), m_rows(0),
      m_cachesize(cachesize), m_cacheoccupancy(), m_cache(),
      m_chunksize(chunksize), m_compress(compress)
  {
    if(cachesize)m_cache = new T[cachesize];
  }

  template<typename T> VSOctaveH5WriterVector<T>::
  ~VSOctaveH5WriterVector()
  {
    if(m_rows == 0)
      {
	makeEmptyMatrix(m_gid,1,0);
      }
    else
      {
	if(m_cacheoccupancy)write(m_cache,m_cacheoccupancy);
	H5Dclose(m_did);
	H5Tclose(m_tid);
      }
    if(m_cachesize)delete[] m_cache;
  }

  template<typename T> void VSOctaveH5WriterVector<T>::
  VSOctaveH5WriterVector::flush()
  {
    if(m_cacheoccupancy)
      {
	write(m_cache,m_cacheoccupancy);
	m_cacheoccupancy=0;
	H5Fflush(m_gid,H5F_SCOPE_LOCAL);
      }
  }

  template<typename T> bool VSOctaveH5WriterVector<T>::
  append(const T& x)
  {
    bool status = true;
    m_rows++;
    if(m_cachesize == 0)
      status=write(&x,1);
    else
      {
	m_cache[m_cacheoccupancy++]=x;
	if(m_cacheoccupancy==m_cachesize)
	  {
	    status=write(m_cache,m_cacheoccupancy);
	    m_cacheoccupancy=0;
	  }
      }
    return status;
  }

  template<typename T> bool VSOctaveH5WriterVector<T>::
  write(const T* ptr, unsigned n)
  {
    unsigned nwritten = m_rows-n;

    hsize_t count[2] = { n, 1 };
    h5_select_hyperslab_t start[2] = { nwritten, 0 };
    hsize_t dims[2] = { m_rows, 1 };

    hid_t disk_sid;
    if(nwritten == 0)
      {
	hsize_t maxdims[2] = { H5S_UNLIMITED, 1 };
	m_tid = VSOctaveH5Type<T>::h5_type();
	disk_sid = H5Screate_simple(2, dims, maxdims);
    
	hsize_t chunk[2] = { m_chunksize, 1 };
	hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(pid, H5D_CHUNKED);
	H5Pset_chunk(pid, 2, chunk);
	if(m_compress)
	  {
	    H5Pset_filter(pid, H5Z_FILTER_DEFLATE, H5Z_FLAG_OPTIONAL, 0, 0);
	    H5Pset_deflate(pid, 6);
	  }

	m_did = H5Dcreate(m_gid, OCTAVE_H5_EL_VALUE, m_tid, disk_sid, pid);
	if(m_did<0)
	  {
	    VSOctaveH5Exception 
	      e(std::string("Could not create dataset: value"));
	    throw(e);
	  }
    
	H5Pclose(pid);
      }
    else
      {
	if(H5Dextend(m_did,dims)<0)
	  {
	    VSOctaveH5Exception e(std::string("Could not extend dataset"));
	    throw(e);
	  }
	disk_sid = H5Dget_space(m_did);
      }

    H5Sselect_hyperslab(disk_sid, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t core_sid = H5Screate_simple(2, count, NULL);

    if(VSOctaveH5WriterLowLevelHelper<T>::
       write_matrix(m_did, m_tid, core_sid, disk_sid, n, ptr) < 0)
      {
	VSOctaveH5Exception 
	  e(std::string("Could not write to dataset: value"));
	throw(e);
      }

    H5Sclose(disk_sid);
    H5Sclose(core_sid);
    return true;
  }

  // ==========================================================================
  // TABLE
  // ==========================================================================

  template<typename T> VSOctaveH5WriterTable<T>::
  VSOctaveH5WriterTable(hid_t gid, unsigned cols,
			VSOctaveH5WriterBase* parent, 
			unsigned cachesize, unsigned chunksize, bool compress):
    VSOctaveH5WriterBase(parent), 
    m_gid(gid), m_tid(), m_did(), m_cols(cols), m_rows(0), 
    m_cachesize(cachesize), m_chunksize(chunksize), m_compress(compress)
  {
    // nothing to see here -- creation of HDF elements deferred to append()
  }

  template<typename T> VSOctaveH5WriterTable<T>::
  VSOctaveH5WriterTable::~VSOctaveH5WriterTable()
  {
    if(m_rows == 0)
      {
	makeEmptyMatrix(m_gid,m_cols,0); /* meaning of cols is wrong! */
      }
    else
      {
	H5Dclose(m_did);
	H5Tclose(m_tid);
      }
  }

  template<typename T> bool VSOctaveH5WriterTable<T>::
  append(const std::vector<T>& x)
  {
    if(x.size() != m_cols)return false;

    hsize_t count[2] = { 1, m_cols };
    hsize_t start[2] = { m_rows, 0 };

    m_rows++;
    hsize_t dims[2] = { m_rows, m_cols };

    hid_t disk_sid;
    if(m_rows == 1)
      {
	hsize_t maxdims[2] = { H5S_UNLIMITED, m_cols };
	m_tid = VSOctaveH5Type<T>::h5_type();
	disk_sid = H5Screate_simple(2, dims, maxdims);
    
	hsize_t chunk[2] = { m_chunksize, 1 };
	hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_layout(pid, H5D_CHUNKED);
	H5Pset_chunk(pid, 2, chunk);
	if(m_compress)
	  {
	    H5Pset_filter(pid, H5Z_FILTER_DEFLATE, H5Z_FLAG_OPTIONAL, 0, 0);
	    H5Pset_deflate(pid, 6);
	  }

	m_did = H5Dcreate(m_gid, OCTAVE_H5_EL_VALUE, m_tid, disk_sid, pid);
	if(m_did<0)
	  {
	    VSOctaveH5Exception 
	      e(std::string("Could not create dataset: value"));
	    throw(e);
	  }
    
	H5Pclose(pid);
      }
    else
      {
	if(H5Dextend(m_did,dims)<0)
	  {
	    VSOctaveH5Exception e(std::string("Could not extend dataset"));
	    throw(e);
	  }
	disk_sid = H5Dget_space(m_did);
      }

    H5Sselect_hyperslab(disk_sid, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t core_sid = H5Screate_simple(2, count, NULL);

    if(VSOctaveH5WriterLowLevelHelper<T>::
       write_matrix(m_did, m_tid, core_sid, disk_sid, m_cols, &x[0]) < 0)
      {
	VSOctaveH5Exception 
	  e(std::string("Could not write to dataset: value"));
	throw(e);
      }

    H5Sclose(disk_sid);
    H5Sclose(core_sid);

    return true;
  }

  // ==========================================================================
  // COMPOSITE VECTOR
  // ==========================================================================

#define __CASENEW(t)							\
  if(im->type == VSOctaveH5Type<t>::int_type())				\
    {									\
      VSOctaveH5WriterVector<t>* vec =					\
	s->writeExpandableVector<t>(im->name, cachesize,chunksize,compress); \
      vec->registerInterestedParty(this);				\
      m.writer = new VSOH5WCVEWTemplate<t>(vec);			\
      m.actual_writer = vec;						\
      m_members.push_back(m);						\
      continue;								\
    }

  template<typename T> VSOctaveH5WriterCompositeVector<T>::
  VSOctaveH5WriterCompositeVector(VSOctaveH5WriterBase* parent,
				  VSOctaveH5WriterStruct* s,
				  bool struct_is_my_child,
				  const std::string& prefix,
				  unsigned cachesize, unsigned chunksize,
				  bool compress):
    VSOctaveH5WriterBase(parent), m_members(), m_rows()
  {
    VSOctaveH5CompositeDefinition c(prefix);
    VSOctaveH5CompositeCompose<T>::compose(c);
    const std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
    for(std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator 
	  im = m.begin(); im!=m.end(); im++)
      {
	MemWriter m(*im,0);
	__CASENEW(bool);
	__CASENEW(int8_t);
	__CASENEW(int16_t);
	__CASENEW(int32_t);
	__CASENEW(int64_t);
	__CASENEW(uint8_t);
	__CASENEW(uint16_t);
	__CASENEW(uint32_t);
	__CASENEW(uint64_t);
	__CASENEW(float);
	__CASENEW(double);
	if(im->type == OIT_STRING)
	  {
	    VSOctaveH5WriterExpandableCellVector* cv =
	      s->writeExpandableCellVector(im->name);
	    cv->registerInterestedParty(this);
	    m.writer = new VSOH5WCVEWString(cv);
	    m.actual_writer = cv;
	    m_members.push_back(m);
	    continue;
	  }

	vsassert(0);
      }
    if(struct_is_my_child)m_children.insert(s);
  }
#undef __CASENEW

  template<typename T> VSOctaveH5WriterCompositeVector<T>::
  ~VSOctaveH5WriterCompositeVector()
  {
    for(typename std::vector<MemWriter>::const_iterator im = m_members.begin();
	im!=m_members.end(); im++)
      delete im->actual_writer;

    if(!m_children.empty())
      {
	delete *m_children.begin();
	m_children.erase(*m_children.begin());
      }
    vsassert(m_children.empty());
  }
  
  template<typename T> void VSOctaveH5WriterCompositeVector<T>::flush()
  {
    for(typename std::vector<MemWriter>::const_iterator im = m_members.begin();
	im!=m_members.end(); im++)
      im->actual_writer->flush();
  }

  template<typename T> bool VSOctaveH5WriterCompositeVector<T>::
  append(const T& x)
  {
    bool status = true;
    for(typename std::vector<MemWriter>::const_iterator im = m_members.begin();
	im!=m_members.end(); im++)
      status &= im->writer->write(&x, im->field);
    m_rows++;
    return status;
  }

  template<typename T> void VSOctaveH5WriterCompositeVector<T>::
  notifyOfChildDestruction(VSOctaveH5WriterBase* child)
  {
    for(typename std::vector<MemWriter>::iterator im = m_members.begin(); 
	im!=m_members.end(); im++)
      if(child == im->actual_writer)
	{
	  delete im->writer;
	  im->writer = 0;
	  im->actual_writer = 0;
	  return;
	}
    vsassert(0);
  }

  // ==========================================================================
  // CELL ARRAY
  // ==========================================================================

  template<typename T> bool VSOctaveH5WriterCellArray::
  writeScalar(unsigned row, unsigned col, const T& x)
  {
    return m_struct->writeScalar(name(row,col), x);
  }
    
  template<typename T> bool VSOctaveH5WriterCellArray::
  writeVector(unsigned row, unsigned col, 
	      const std::vector<T>& v, VSOctaveH5VecOrientation orient)
  {
    return m_struct->writeVector(name(row,col), v, orient);
  }
    
  template<typename T> bool VSOctaveH5WriterCellArray::
  writeMatrix(unsigned row, unsigned col, 
	      unsigned rows, unsigned cols, const T* m)
  {
    return m_struct->writeMatrix(name(row,col), rows, cols, m);
  }

  template<typename T> bool VSOctaveH5WriterCellArray::
  write(unsigned row, unsigned col, const T& x)
  {
    return m_struct->write(name(row,col), x);
  }
    
  template<typename T> VSOctaveH5WriterVector<T>* 
  VSOctaveH5WriterCellArray::
  writeExpandableVector(unsigned row, unsigned col, 
			unsigned cachesize, unsigned chunksize, bool compress)
  {
    return m_struct->
      writeExpandableVector<T>(name(row,col),cachesize,chunksize,compress);
  }
    
  template<typename T> VSOctaveH5WriterTable<T>* 
  VSOctaveH5WriterCellArray::
  writeExpandableTable(unsigned row, unsigned col, unsigned cols, 
		       unsigned cachesize, unsigned chunksize, bool compress)
  {
    return m_struct->
      writeExpandableTable<T>(name(row,col), cols,
			      cachesize, chunksize, compress);
  }

  template<typename T> bool VSOctaveH5WriterCellArray::
  writeComposite(unsigned row, unsigned col, const T& x)
  {
    return m_struct->writeComposite(name(row,col),x);
  }

  template<typename T> bool VSOctaveH5WriterCellArray::
  writeCompositeVector(unsigned row, unsigned col,
		       const std::vector<T>& x)
  {
    return m_struct->writeCompositeVector(name(row,col),x);
  }

  template<typename T> VSOctaveH5WriterCompositeVector<T>*
  VSOctaveH5WriterCellArray::
  writeCompositeExpandableVector(unsigned row, unsigned col,
				 unsigned cachesize, unsigned chunksize,
				 bool compress)
  {
    return m_struct->
      writeCompositeExpandableVector<T>(name(row,col), 
					cachesize, chunksize, compress);
  }

  // ==========================================================================
  // CELL VECTOR
  // ==========================================================================

  template<typename T> bool VSOctaveH5WriterCellVector::
  writeScalar(unsigned iel, const T& x)
  {
    return m_struct->writeScalar(name(iel), x);
  }
    
  template<typename T> bool VSOctaveH5WriterCellVector::
  writeVector(unsigned iel, 
	      const std::vector<T>& v, VSOctaveH5VecOrientation orient)
  {
    return m_struct->writeVector(name(iel), v, orient);
  }
    
  template<typename T> bool VSOctaveH5WriterCellVector::
  writeMatrix(unsigned iel, 
	      unsigned rows, unsigned cols, const T* m)
  {
    return m_struct->writeMatrix(name(iel), rows, cols, m);
  }

  template<typename T> bool VSOctaveH5WriterCellVector::
  write(unsigned iel, const T& x)
  {
    return m_struct->write(name(iel), x);
  }
    
  template<typename T> VSOctaveH5WriterVector<T>* 
  VSOctaveH5WriterCellVector::
  writeExpandableVector(unsigned iel, 
			unsigned cachesize, unsigned chunksize, bool compress)
  {
    return m_struct->
      writeExpandableVector<T>(name(iel),cachesize,chunksize,compress);
  }
    
  template<typename T> VSOctaveH5WriterTable<T>* 
  VSOctaveH5WriterCellVector::
  writeExpandableTable(unsigned iel, unsigned cols, 
		       unsigned cachesize, unsigned chunksize, bool compress)
  {
    return m_struct->
      writeExpandableTable<T>(name(iel), cols, cachesize, chunksize, compress);
  }

  template<typename T> bool VSOctaveH5WriterCellVector::
  writeComposite(unsigned iel, const T& x)
  {
    return m_struct->writeComposite(name(iel),x);
  }

  template<typename T> bool VSOctaveH5WriterCellVector::
  writeCompositeVector(unsigned iel, const std::vector<T>& x)
  {
    return m_struct->writeCompositeVector(name(iel),x);
  }

  template<typename T> VSOctaveH5WriterCompositeVector<T>*
  VSOctaveH5WriterCellVector::
  writeCompositeExpandableVector(unsigned iel,
				 unsigned cachesize, unsigned chunksize, 
				 bool compress)
  {
    return m_struct->
      writeCompositeExpandableVector<T>(name(iel), 
					cachesize, chunksize, compress);
  }

  // ==========================================================================
  // EXPANDABLE CELL VECTOR
  // ==========================================================================

  template<typename T> bool VSOctaveH5WriterExpandableCellVector::
  appendScalar(const T& x)
  {
    return m_struct->writeScalar(name(), x);
  }
    
  template<typename T> bool VSOctaveH5WriterExpandableCellVector::
  appendVector(const std::vector<T>& v, VSOctaveH5VecOrientation orient)
  {
    return m_struct->writeVector(name(), v, orient);
  }
    
  template<typename T> bool VSOctaveH5WriterExpandableCellVector::
  appendMatrix(unsigned rows, unsigned cols, const T* m)
  {
    return m_struct->writeMatrix(name(), rows, cols, m);
  }

  template<typename T> bool VSOctaveH5WriterExpandableCellVector::
  append(const T& x)
  {
    return m_struct->write(name(), x);
  }
    
  template<typename T> VSOctaveH5WriterVector<T>* 
  VSOctaveH5WriterExpandableCellVector::
  appendExpandableVector(unsigned cachesize, unsigned chunksize, bool compress)
  {
    return m_struct->
      writeExpandableVector<T>(name(),cachesize,chunksize, compress);
  }
    
  template<typename T> VSOctaveH5WriterTable<T>* 
  VSOctaveH5WriterExpandableCellVector::
  appendExpandableTable(unsigned cols, 
			unsigned cachesize, unsigned chunksize, bool compress)
  {
    return m_struct->
      writeExpandableTable<T>(name(), cols, cachesize, chunksize, compress);
  }

  template<typename T> bool VSOctaveH5WriterExpandableCellVector::
  appendComposite(const T& x)
  {
    return m_struct->writeComposite(name(),x);
  }

  template<typename T> bool VSOctaveH5WriterExpandableCellVector::
  appendCompositeVector(const std::vector<T>& x)
  {
    return m_struct->writeCompositeVector(name(),x);
  }

  template<typename T> VSOctaveH5WriterCompositeVector<T>*
  VSOctaveH5WriterExpandableCellVector::
  appendCompositeExpandableVector(unsigned cachesize, unsigned chunksize,
				  bool compress)
  {
    return m_struct->
      writeCompositeExpandableVector<T>(name(), 
					cachesize, chunksize, compress);
  }

  // ==========================================================================
  // STRUCT
  // ==========================================================================

  template<typename T> bool VSOctaveH5WriterStruct::
  writeScalar(const std::string& name, const T& x)
  {
    hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_SCALAR,
			     VSOctaveH5Type<T>::oct_type());
  
    hid_t tid = VSOctaveH5Type<T>::h5_type();
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t did = H5Dcreate(gid, OCTAVE_H5_EL_VALUE, tid, sid, H5P_DEFAULT);
    if(did<0)
      {
	VSOctaveH5Exception 
	  e(name+std::string(": could not create dataset: value"));
	throw(e);
      }
    
    if(VSOctaveH5WriterLowLevelHelper<T>::
       write_scalar(did, tid, H5S_ALL, H5S_ALL, x) < 0)
      {
	VSOctaveH5Exception 
	  e(name+std::string(": could not write to dataset: value"));
	throw(e);
      }

    H5Dclose(did);
    H5Sclose(sid);
    H5Tclose(tid);
    H5Gclose(gid);

    return true;
  }
  
  template<typename T> bool VSOctaveH5WriterStruct::
  VSOctaveH5WriterStruct::writeVector(const std::string& name, 
				      const std::vector<T>& v, 
				      VSOctaveH5VecOrientation orient)
  {
    hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_MATRIX,
			     VSOctaveH5Type<T>::oct_type());
  
    if(v.size() == 0)
      {
	if(orient == VO_COLUMN)makeEmptyMatrix(gid,0,1);
	else makeEmptyMatrix(gid,1,0);
	H5Gclose(gid);
	return true;
      }

    hid_t tid = VSOctaveH5Type<T>::h5_type();
    hsize_t dims[2] = { 1, 1 };
    if(orient == VO_COLUMN)dims[1] = v.size();
    else dims[0] = v.size();
    hid_t sid = H5Screate_simple(2, dims, NULL);
    hid_t did = H5Dcreate(gid, OCTAVE_H5_EL_VALUE, tid, sid, H5P_DEFAULT);
    if(did<0)
      {
	VSOctaveH5Exception 
	  e(name+std::string(": could not create dataset: value"));
	throw(e);
      }

    if(VSOctaveH5WriterLowLevelHelper<T>::
       write_vector(did, tid, H5S_ALL, H5S_ALL, v) < 0)
      {
	VSOctaveH5Exception 
	  e(name+std::string(": could not write to dataset: value"));
	throw(e);
      }

    H5Dclose(did);
    H5Sclose(sid);
    H5Tclose(tid);
    H5Gclose(gid);

    return true;
  }

  template<typename T> bool VSOctaveH5WriterStruct::
  writeMatrix(const std::string& name, 
	      unsigned rows, unsigned cols, const T* m)
  {
    hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_MATRIX,
			     VSOctaveH5Type<T>::oct_type());
  
    if((rows == 0)||(cols == 0))
      {
	// This is messed up - I would hate to totally go to the
	// FORTRAN standard - so I have to transpose the matrix, so that
	// a nr x nc matric in C++ becomes transposed in octave
	makeEmptyMatrix(gid,cols,rows); 
	H5Gclose(gid);
	return true;
      }

    hid_t tid = VSOctaveH5Type<T>::h5_type();
    hsize_t dims[2] = { rows, cols };
    hid_t sid = H5Screate_simple(2, dims, NULL);
    hid_t did = H5Dcreate(gid, OCTAVE_H5_EL_VALUE, tid, sid, H5P_DEFAULT);
    if(did<0)
      {
	VSOctaveH5Exception 
	  e(name+std::string(": could not create dataset: value"));
	throw(e);
      }
    
    if(VSOctaveH5WriterLowLevelHelper<T>::
       write_matrix(did, tid, H5S_ALL, H5S_ALL, rows*cols, m) < 0)
      {
	VSOctaveH5Exception 
	  e(name+std::string(": could not write to dataset: value"));
	throw(e);
      }

    H5Dclose(did);
    H5Sclose(sid);
    H5Tclose(tid);
    H5Gclose(gid);

    return true;
  }

  template<typename T> bool VSOctaveH5WriterStruct::
  write(const std::string& name, const T& x)
  {
    return VSOctaveH5WriterHelper<T>(this,name,x);
  }

  template<typename T> VSOctaveH5WriterVector<T>* VSOctaveH5WriterStruct::
  writeExpandableVector(const std::string& name, 
			unsigned cachesize, unsigned chunksize, bool compress)
  {
    hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_MATRIX, 
			     VSOctaveH5Type<T>::oct_type());
  
    VSOctaveH5WriterVector<T>* v = 
      new VSOctaveH5WriterVector<T>(gid,this,cachesize,chunksize,compress);
    registerChild(v,gid);
    
    return v;
  }

  template<typename T> VSOctaveH5WriterTable<T>*  VSOctaveH5WriterStruct::
  writeExpandableTable(const std::string& name, unsigned cols, 
		       unsigned cachesize, unsigned chunksize, bool compress)
  {
    hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_MATRIX, 
			     VSOctaveH5Type<T>::oct_type());
  
    VSOctaveH5WriterTable<T>* m = 
      new VSOctaveH5WriterTable<T>(gid,cols,this,cachesize,chunksize,compress);
    registerChild(m,gid);
    
    return m;
  }

#define __CASESCALAR(t,x,im)						\
  if(im->type == VSOctaveH5Type<t>::int_type())				\
    {									\
      status &= writeScalar(im->name,*((t*)((char*)(&x) + im->offset))); \
      continue;								\
    }  
  
  template<typename T> bool VSOctaveH5WriterStruct::
  writeCompositeHere(const T& x, const std::string& prefix)
  {
    bool status = true;
    VSOctaveH5CompositeDefinition c(prefix);
    VSOctaveH5CompositeCompose<T>::compose(c);
    const std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
    for(std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator 
	  im = m.begin(); im!=m.end(); im++)
      {
	__CASESCALAR(bool, x, im);
	__CASESCALAR(int8_t, x, im);
	__CASESCALAR(int16_t, x, im);
	__CASESCALAR(int32_t, x, im);
	__CASESCALAR(int64_t, x, im);
	__CASESCALAR(uint8_t, x, im);
	__CASESCALAR(uint16_t, x, im);
	__CASESCALAR(uint32_t, x, im);
	__CASESCALAR(uint64_t, x, im);
	__CASESCALAR(float, x, im);
	__CASESCALAR(double, x, im);
	if(im->type == OIT_STRING)
	  {
	    status &=
	      writeString(im->name,*((std::string*)((char*)(&x)+im->offset)));
	    continue;
	  }
	vsassert(0);
      }
    return status;
  }
#undef __CASESCALAR
  
#define __CASEVECTOR(t,x,im)						\
  if(im->type == VSOctaveH5Type<t>::int_type())				\
    {									\
      unsigned nel = x.size();						\
      std::vector<t> v(nel);						\
      for(unsigned iel=0;iel<nel;iel++)					\
	v[iel] = *((t*)((char*)(&x[iel]) + im->offset));		\
      status &= writeVector(im->name,v);				\
      continue;								\
    }
  
  template<typename T> bool VSOctaveH5WriterStruct::
  writeStructCellVector(const std::string& name, const std::vector<T>& x)
  {
    bool status = true;
    VSOctaveH5WriterCellVector* wc = writeCellVector(name, x.size());
    for(unsigned i = 0; i < x.size(); i++)
      {
	VSOctaveH5WriterStruct* ws = wc->writeStruct(i);
	x[i].save(ws);
	delete ws;
      }
    delete wc;
    return status;
  }

  template<typename T> bool VSOctaveH5WriterStruct::
  writeCompositeVectorHere(const std::vector<T>& x, const std::string& prefix)
  {
    bool status = true;
    VSOctaveH5CompositeDefinition c(prefix);
    VSOctaveH5CompositeCompose<T>::compose(c);
    const std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
    for(std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator 
	  im = m.begin(); im!=m.end(); im++)
      {
	__CASEVECTOR(bool, x, im);
	__CASEVECTOR(int8_t, x, im);
	__CASEVECTOR(int16_t, x, im);
	__CASEVECTOR(int32_t, x, im);
	__CASEVECTOR(int64_t, x, im);
	__CASEVECTOR(uint8_t, x, im);
	__CASEVECTOR(uint16_t, x, im);
	__CASEVECTOR(uint32_t, x, im);
	__CASEVECTOR(uint64_t, x, im);
	__CASEVECTOR(float, x, im);
	__CASEVECTOR(double, x, im);
	if(im->type == OIT_STRING)
	  {
	    unsigned nel = x.size();
	    VSOctaveH5WriterCellVector* cv = writeCellVector(im->name,nel);
	    if(cv)
	      {
		for(unsigned iel=0;iel<nel;iel++)
		  status &= cv->writeString(iel,
			       *((std::string*)((char*)(&x[iel])+im->offset)));
		delete cv;
	      }
	    continue;
	  }
	vsassert(0);
      }
    return status;
  }
#undef __CASEVECTOR

  template<typename T> 
  VSOctaveH5WriterCompositeVector<T>* VSOctaveH5WriterStruct::
  writeCompositeExpandableVectorHere(const std::string& prefix,
				     unsigned cachesize, unsigned chunksize, 
				     bool compress)
  {
    VSOctaveH5WriterCompositeVector<T>* c = 
      new VSOctaveH5WriterCompositeVector<T>(this, this, false, prefix,
					     cachesize,chunksize,compress);
    registerChild(c,-1);
    return c;
  }

  template<typename T> bool VSOctaveH5WriterStruct::
  writeComposite(const std::string& name, const T& x)
  {
    VSOctaveH5WriterStruct* s = writeStruct(name);
    bool status = s->writeCompositeHere(x);
    delete s;
    return status;
  }
  
  template<typename T> bool VSOctaveH5WriterStruct::
  writeCompositeVector(const std::string& name, const std::vector<T>& x)
  {
    VSOctaveH5WriterStruct* s = writeStruct(name);
    bool status = s->writeCompositeVectorHere(x);
    delete s;
    return status;
  }
    
  template<typename T> VSOctaveH5WriterCompositeVector<T>*
  VSOctaveH5WriterStruct::
  writeCompositeExpandableVector(const std::string& name,
				 unsigned cachesize, unsigned chunksize,
				 bool compress)
  {
    hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_STRUCT);
    hid_t value_gid = H5Gcreate(gid, OCTAVE_H5_EL_VALUE, 0);
    if(value_gid<0)
      {
	VSOctaveH5Exception 
	  e(std::string("Could not create group: " OCTAVE_H5_EL_VALUE));
	throw(e);
      }
  
    VSOctaveH5WriterStruct* s = new VSOctaveH5WriterStruct(value_gid,0);
    VSOctaveH5WriterCompositeVector<T>* c = 
      new VSOctaveH5WriterCompositeVector<T>(this, s, true, "",
					     cachesize,chunksize,compress);
    registerChild(c,gid);
    return c;
  }

} // namespace VERITAS

#endif // VSOCTAVEH5WRITER_HPP
