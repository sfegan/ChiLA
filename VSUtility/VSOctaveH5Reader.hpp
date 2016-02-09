//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveH5Reader.hpp

  Classes to read from structured GNU/Octave HDF5 file. Data can be
  written directly from octave.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/12/2006
*/

#include"VSOctaveIO.hpp"

#ifndef VSOCTAVEH5READER_HPP
#define VSOCTAVEH5READER_HPP

#include<vsassert>
#include<string>
#include<vector>
#include<set>
#include<map>
#include<limits>
#include<stdexcept>

#include<hdf5.h>

#include"VSOctaveH5Writer.hpp"

namespace VERITAS
{
  
  // ==========================================================================
  // Octave H5 reader base class
  // ==========================================================================

  class VSOctaveH5ReaderStruct;
  class VSOctaveH5ReaderCellVector;

  class VSOctaveH5ReaderBase
  {
  public:
    typedef VSOctaveH5CompositeDefinition::MemberSubset MemberSubset;

    VSOctaveH5ReaderBase(VSOctaveH5ReaderBase* parent):
      m_parent(parent), m_children() { /* nothing to see here */ }
    virtual ~VSOctaveH5ReaderBase();

  protected:

#define __CPC(x) x(o.x)
#define __CPA(x) x = o.x;

    class Variable
    {
    public:
      Variable(): 
	name(), parent(), gid(), attr(), type(), eltype(), reader() 
      { /* nothing to see here */ }

      Variable(const Variable& o): 
	__CPC(name),
	__CPC(parent),
	__CPC(gid),
	__CPC(attr),
	__CPC(type),
	__CPC(eltype),
	__CPC(reader)
      { /* nothing to see here */ }

      Variable& operator= (const Variable& o)
      {
	__CPA(name);
	__CPA(parent);
	__CPA(gid);
	__CPA(attr);
	__CPA(type);
	__CPA(eltype);
	__CPA(reader);
	return *this;
      }

      static herr_t iterateOverAtr(hid_t gid, const char*  name, void *object);
      void registerAttribute(const char* name);
      void readAttributes();
      template<typename T> bool getAttribute(const std::string& name,
					     T& x);

      std::string name;
      VSOctaveH5ReaderStruct* parent;
      hid_t gid;
      std::set<std::string> attr;
      std::string type;
      std::string eltype;
      VSOctaveH5ReaderBase* reader;
    };

#undef __CPC
#undef __CPA

    bool getEmpty(Variable* var, unsigned& rows, unsigned& cols);
    bool getCString(hid_t gid, const std::string& name, std::string& s);

    void registerChild(VSOctaveH5ReaderBase *reader);

    virtual void notifyOfChildDestruction(VSOctaveH5ReaderBase* child);

    VSOctaveH5ReaderBase* m_parent;
    std::set<VSOctaveH5ReaderBase*> m_children;
  private:
    VSOctaveH5ReaderBase(const VSOctaveH5ReaderBase&);
    VSOctaveH5ReaderBase& operator=(const VSOctaveH5ReaderBase&);
  };

  // ==========================================================================
  // Octave H5 reader helper -- needed to handle stupid boolean behaviour
  // ==========================================================================

  template<typename T> class VSOctaveH5ReaderActual
  {
  public:
    static inline herr_t 
    read(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	 unsigned n, T* ptr);
    static inline herr_t 
    read_scalar(hid_t did, hid_t tid, T& x);
    static inline herr_t 
    read_vector(hid_t did, hid_t tid, std::vector<T>& v);
    static inline herr_t 
    read_matrix(hid_t did, hid_t tid, unsigned n, T* m);
  };

  template<> class VSOctaveH5ReaderActual<bool>
  {
  public:
    static inline herr_t 
    read(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid, 
	 unsigned n, bool* ptr);
    static inline herr_t 
    read_scalar(hid_t did, hid_t tid, bool& x);
    static inline herr_t 
    read_vector(hid_t did, hid_t tid, std::vector<bool>& v);
    static inline herr_t 
    read_matrix(hid_t did, hid_t tid, unsigned n, bool* m);
  };

  // ==========================================================================
  // Octave H5 reader giving element-by-element access to a vector
  // ==========================================================================

  class VSOctaveH5ReaderMatrix;

  template<typename T> 
  class VSOctaveH5ReaderVector: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderVector(VSOctaveH5ReaderMatrix* matrix,
			   unsigned cachesize=1024);
    virtual ~VSOctaveH5ReaderVector();
    bool element(T& x, unsigned index);
    inline unsigned rows() const;
    T at(unsigned index) 
    { T x; if(!element(x,index))throw std::out_of_range(__PRETTY_FUNCTION__); 
      return x; }
    T operator[] (unsigned index) { T x; element(x,index); return x; }

  private:
    VSOctaveH5ReaderVector(const VSOctaveH5ReaderVector&);
    VSOctaveH5ReaderVector& operator=(const VSOctaveH5ReaderVector&);

    bool read(T* ptr, hsize_t n, hsize_t index);
    VSOctaveH5ReaderMatrix* m_matrix;
    hid_t                   m_tid;
    unsigned                m_cachesize;
    unsigned                m_cachestart;
    T*                      m_cache;
  };
  
  // ==========================================================================
  // Octave H5 reader giving row-by-row access to a table (matrix)
  // ==========================================================================

  template<typename T> 
  class VSOctaveH5ReaderTable: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderTable(hid_t gid, unsigned cols,
			  VSOctaveH5ReaderBase* parent, 
			  unsigned cachesize=1024);
    virtual ~VSOctaveH5ReaderTable();
    bool row(std::vector<T>& row, unsigned index);
    unsigned rows() const { return m_rows; }
    unsigned cols() const { return m_cols; }
    std::vector<T> operator[] (unsigned index) 
    { std::vector<T> x; row(x,index); return x; }
  private:
    VSOctaveH5ReaderTable(const VSOctaveH5ReaderTable&);
    VSOctaveH5ReaderTable& operator=(const VSOctaveH5ReaderTable&);

    bool read(T* ptr, unsigned n, unsigned index);

    hid_t m_gid;
    hid_t m_tid;
    hid_t m_did;
    unsigned m_cols;
    unsigned m_rows;
    unsigned m_cachesize;
  };

  // ==========================================================================
  // Octave H5 reader giving element-by-element access to a composite vector
  // ==========================================================================

  class VSOH5RCVElementReader
  {
  public:
    virtual ~VSOH5RCVElementReader();
    virtual bool read(void* data, unsigned index,
		      const VSOctaveH5CompositeDefinition::Member& field) = 0;
    virtual void releaseReader() = 0;
    virtual bool isNullReader() const;
  };

  class VSOH5RCVERString: public VSOH5RCVElementReader
  {
  public:
    VSOH5RCVERString(VSOctaveH5ReaderCellVector* cv):
      VSOH5RCVElementReader(), m_cv(cv) { }
    virtual ~VSOH5RCVERString();
    virtual bool read(void* data, unsigned index,
		      const VSOctaveH5CompositeDefinition::Member& field);
    virtual void releaseReader();
  private:
    VSOctaveH5ReaderCellVector* m_cv;
  };

  template<typename T> class VSOH5RCVERNull: public VSOH5RCVElementReader
  {
  public:
    VSOH5RCVERNull(): VSOH5RCVElementReader() { }
    virtual ~VSOH5RCVERNull() { }
    virtual bool read(void* data, unsigned index,
		      const VSOctaveH5CompositeDefinition::Member& field)
    { *((T*)((char*)(data)+field.offset)) = T(); return true; }
    virtual void releaseReader() { }
    virtual bool isNullReader() const { return true; }
  private:
    VSOH5RCVERNull(const VSOH5RCVERNull&);
    VSOH5RCVERNull& operator=(const VSOH5RCVERNull&);
  };

  template<typename T> class VSOH5RCVERTemplate: public VSOH5RCVElementReader
  {
  public:
    VSOH5RCVERTemplate(VSOctaveH5ReaderVector<T>* r): 
      VSOH5RCVElementReader(), m_reader(r) { }
    virtual ~VSOH5RCVERTemplate() { delete m_reader; }
    virtual bool read(void* data, unsigned index,
		      const VSOctaveH5CompositeDefinition::Member& field)
    {
      if((m_reader)
	 &&(m_reader->element(*((T*)((char*)(data)+field.offset)), index)))
	return true;
      *((T*)((char*)(data)+field.offset)) = T();
      return false;
    }
    virtual void releaseReader() { m_reader=0; }
  private:
    VSOH5RCVERTemplate(const VSOH5RCVERTemplate&);
    VSOH5RCVERTemplate& operator=(const VSOH5RCVERTemplate&);

    VSOctaveH5ReaderVector<T>* m_reader;
  };

  template<typename T> 
  class VSOctaveH5ReaderCompositeVector: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderCompositeVector(VSOctaveH5ReaderStruct* s,
				    const std::string& prefix,
				    const MemberSubset& member_subset = MemberSubset(),
				    unsigned cachesize=1024);
    virtual ~VSOctaveH5ReaderCompositeVector();
    bool element(T& x, unsigned index);
    inline unsigned rows() const { return m_rows; }
    T at(unsigned index) 
    { T x; if(!element(x,index))throw std::out_of_range(__PRETTY_FUNCTION__); 
      return x; }
    T operator[] (unsigned index) { T x; element(x,index); return x; }
    bool hasElement(const std::string& name) const;

  private:
    VSOctaveH5ReaderCompositeVector(const VSOctaveH5ReaderCompositeVector&);
    VSOctaveH5ReaderCompositeVector& 
    operator=(const VSOctaveH5ReaderCompositeVector&);

    struct MemReader
    {
      MemReader(const VSOctaveH5CompositeDefinition::Member& _field,
		VSOH5RCVElementReader* _reader = 0):
	field(_field), reader(_reader) { }
      MemReader(const MemReader& o):
	field(o.field), reader(o.reader) { }
      MemReader& operator=(const MemReader& o)
      { field  = o.field; reader = o.reader; return *this; }

      VSOctaveH5CompositeDefinition::Member field;
      VSOH5RCVElementReader*                reader;
    };

  std::vector<MemReader>   m_members;
    unsigned                 m_rows;
  };  

  // ==========================================================================
  // Octave H5 reader for scalar, vector, matrix and string
  // ==========================================================================

  class VSOctaveH5ReaderMatrix: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderMatrix(Variable* var);
    virtual ~VSOctaveH5ReaderMatrix();

    bool dimensions(unsigned& rows, unsigned& cols) 
    { rows=m_rows; cols=m_cols; return true; }

    bool isEmpty() { return (m_rows==0)||(m_cols==0); }
    bool isScalar() { return (m_rows==1)&&(m_cols==1); }
    bool isVector() { return (m_rows==1)||(m_cols==1); }

    bool readString(std::string& s);
    
    template<typename T> bool readScalar(T& x);
    
    template<typename T> bool readVector(std::vector<T>& v);
    
    template<typename T>
    inline bool readMatrix(unsigned& rows, unsigned& cols, T*& m);

    template<typename T> bool readMatrix(T* m);

    template<typename T> VSOctaveH5ReaderVector<T>* 
    readVectorElements(unsigned cachesize=1024);
    
  private:
    VSOctaveH5ReaderMatrix(const VSOctaveH5ReaderMatrix&);
    VSOctaveH5ReaderMatrix& operator=(const VSOctaveH5ReaderMatrix&);

    template<typename T> friend class VSOctaveH5ReaderVector;
    template<typename T> friend class VSOctaveH5ReaderTable;

    hid_t m_did;
    unsigned m_cols;
    unsigned m_rows;
  };

  // ==========================================================================
  // Octave H5 writer cell array - just an interface to struct
  // ==========================================================================

  class VSOctaveH5ReaderCellArray: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderCellArray(Variable* var);
    virtual ~VSOctaveH5ReaderCellArray();

    bool isEmpty() { return (m_rows==0)||(m_cols==0); }

    bool dimensions(unsigned& rows, unsigned& cols) 
    { rows=m_rows; cols=m_cols; return true; }

    std::string variableType(unsigned row, unsigned col);
    std::string elementType(unsigned row, unsigned col);

    bool isValid(unsigned row, unsigned col);
    bool isEmpty(unsigned row, unsigned col);
    bool isScalar(unsigned row, unsigned col);
    bool isVector(unsigned row, unsigned col);
    bool isMatrix(unsigned row, unsigned col);
    bool isString(unsigned row, unsigned col);
    bool isCellArray(unsigned row, unsigned col);
    bool isCellVector(unsigned row, unsigned col);
    bool isStruct(unsigned row, unsigned col);

    bool dimensions(unsigned row, unsigned col,
		    unsigned& rows, unsigned& cols);

    bool readString(unsigned row, unsigned col, std::string& s);
    
    template<typename T>
    bool readScalar(unsigned row, unsigned col, T& x);
    
    template<typename T>
    bool readVector(unsigned row, unsigned col, std::vector<T>& v);
    
    template<typename T>
    bool readMatrix(unsigned row, unsigned col, 
		    unsigned& rows, unsigned& cols, T*& m);

    template<typename T>
    bool readMatrix(unsigned row, unsigned col, T* m);
    
    VSOctaveH5ReaderStruct* readStruct(unsigned row, unsigned col);
    
    VSOctaveH5ReaderCellArray* readCellArray(unsigned row, unsigned col);

    VSOctaveH5ReaderCellVector* readCellVector(unsigned row, unsigned col);
    
    template<typename T> VSOctaveH5ReaderVector<T>* 
    readVectorElements(unsigned row, unsigned col, unsigned cachesize=1024);
    
    //template<typename T> VSOctaveH5ReaderTable<T>* 
    //readTableRows(unsigned row, unsigned col, unsigned cachesize=1024);

    template<typename T> bool
    readComposite(unsigned row, unsigned col, T& x, 
		  const MemberSubset& member_subset = MemberSubset());

    template<typename T> bool 
    readCompositeVector(unsigned row, unsigned col, std::vector<T>& x,
			const MemberSubset& member_subset = MemberSubset());

    template<typename T> VSOctaveH5ReaderCompositeVector<T>*
    readCompositeVectorElements(unsigned row, unsigned col,
			const MemberSubset& member_subset = MemberSubset(), 
				unsigned cachesize=1024);

    bool copyTo(unsigned srow, unsigned scol,
		VSOctaveH5WriterStruct* s, std::string dname);
    bool copyTo(unsigned srow, unsigned scol,
		VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol);
    bool copyTo(unsigned srow, unsigned scol,
		VSOctaveH5WriterCellVector* c, unsigned diel);
    bool copyTo(unsigned srow, unsigned scol,
		VSOctaveH5WriterExpandableCellVector* c);

  private:
    VSOctaveH5ReaderCellArray(const VSOctaveH5ReaderCellArray&);
    VSOctaveH5ReaderCellArray& operator=(const VSOctaveH5ReaderCellArray&);

    friend class VSOctaveH5ReaderStruct; // so that it can get m_struct

    std::string name(unsigned row, unsigned col);

    VSOctaveH5ReaderStruct* m_struct;
    unsigned m_rows;
    unsigned m_cols;
    unsigned m_name_precision;
  };

  // ==========================================================================
  // Octave H5 writer cell vector - just an interface to struct
  // ==========================================================================

  class VSOctaveH5ReaderCellVector: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderCellVector(VSOctaveH5ReaderCellArray* cell);
    virtual ~VSOctaveH5ReaderCellVector();

    bool isEmpty() { return m_nel==0; }

    unsigned dimensions() { return m_nel; }
    //bool dimensions(unsigned& nel) { nel=m_nel; return true; }

    std::string variableType(unsigned iel);
    std::string elementType(unsigned iel);

    bool isValid(unsigned iel);
    bool isEmpty(unsigned iel);
    bool isScalar(unsigned iel);
    bool isVector(unsigned iel);
    bool isMatrix(unsigned iel);
    bool isString(unsigned iel);
    bool isCellArray(unsigned iel);
    bool isCellVector(unsigned iel);
    bool isStruct(unsigned iel);

    bool dimensions(unsigned iel, unsigned& rows, unsigned& cols);

    bool readString(unsigned iel, std::string& s);
    
    template<typename T> bool readScalar(unsigned iel, T& x);
    
    template<typename T> bool readVector(unsigned iel, std::vector<T>& v);
    
    template<typename T> bool 
    readMatrix(unsigned iel, unsigned& rows, unsigned& cols, T*& m);

    template<typename T> bool readMatrix(unsigned iel, T* m);
    
    VSOctaveH5ReaderStruct* readStruct(unsigned iel);
    
    VSOctaveH5ReaderCellArray* readCellArray(unsigned iel);

    VSOctaveH5ReaderCellVector* readCellVector(unsigned iel);
    
    template<typename T> VSOctaveH5ReaderVector<T>* 
    readVectorElements(unsigned iel, unsigned cachesize=1024);
    
    //template<typename T> VSOctaveH5ReaderTable<T>* 
    //readTableRows(unsigned iel, unsigned cachesize=1024);

    template<typename T> bool readComposite(unsigned iel, T& x, 
			const MemberSubset& member_subset = MemberSubset());

    template<typename T> bool 
    readCompositeVector(unsigned iel, std::vector<T>& x,
			const MemberSubset& member_subset = MemberSubset());

    template<typename T> VSOctaveH5ReaderCompositeVector<T>*
    readCompositeVectorElements(unsigned iel,
			const MemberSubset& member_subset = MemberSubset(), 
				unsigned cachesize=1024);

  private:
    VSOctaveH5ReaderCellVector(const VSOctaveH5ReaderCellVector&);
    VSOctaveH5ReaderCellVector& operator=(const VSOctaveH5ReaderCellVector&);

    unsigned row(unsigned iel) { return (m_column)?iel:0; }
    unsigned col(unsigned iel) { return (m_column)?0:iel; }

    VSOctaveH5ReaderCellArray* m_cell;
    unsigned m_nel;
    bool m_column;
  };

  // ==========================================================================
  // Octave H5 reader struct - this is where most of the action is
  // ==========================================================================

  class VSOctaveH5ReaderStruct: public VSOctaveH5ReaderBase
  {
  public:
    VSOctaveH5ReaderStruct(const std::string& filename);
    VSOctaveH5ReaderStruct(Variable* var);
    virtual ~VSOctaveH5ReaderStruct();

    std::vector<std::string> variables() const;
    std::vector<std::string> fieldNames() const { return variables(); }

    bool isEmpty() { return m_vars.empty(); }
    bool dimensions(unsigned& rows, unsigned& cols)
    { rows=m_vars.size(); cols=1; return true; }

    bool hasVariable(const std::string& name);
    std::string variableType(const std::string& name);
    std::string elementType(const std::string& name);

    bool isValid(const std::string& name);
    bool isEmpty(const std::string& name);
    bool isScalar(const std::string& name);
    bool isVector(const std::string& name);
    bool isMatrix(const std::string& name);
    bool isString(const std::string& name);
    bool isCellArray(const std::string& name);
    bool isCellVector(const std::string& name);
    bool isStruct(const std::string& name);

    bool dimensions(const std::string& name, unsigned& rows, unsigned& cols);

    bool readString(const std::string& name, std::string& s);
    
    template<typename T>
    bool readScalar(const std::string& name, T& x);
    
    template<typename T>
    bool readVector(const std::string& name, std::vector<T>& v);
    
    template<typename T>
    bool readMatrix(const std::string& name, 
		    unsigned& rows, unsigned& cols, T*& m);

    template<typename T>
    bool readMatrix(const std::string& name, T* m);
    
    VSOctaveH5ReaderStruct* readStruct(const std::string& name);
    
    VSOctaveH5ReaderCellArray* readCellArray(const std::string& name);
    
    template<typename T> 
    bool readStructCellVector(const std::string& name, std::vector<T>& x);

    VSOctaveH5ReaderCellVector* readCellVector(const std::string& name);

    template<typename T> VSOctaveH5ReaderVector<T>* 
    readVectorElements(const std::string& name, unsigned cachesize=1024);
    
    //template<typename T> VSOctaveH5ReaderTable<T>* 
    //readTableRows(const std::string& name, unsigned cachesize=1024);

    template<typename T> bool
    readCompositeHere(T& x, const std::string& prefix = std::string(), 
		      const MemberSubset& member_subset = MemberSubset());

    template<typename T> bool 
    readCompositeVectorHere(std::vector<T>& x,
			    const std::string& prefix = std::string(),
		      const MemberSubset& member_subset = MemberSubset());

    template<typename T> VSOctaveH5ReaderCompositeVector<T>*
    readCompositeVectorElementsHere(const std::string& prefix = std::string(),
		      const MemberSubset& member_subset = MemberSubset(),
    				    unsigned cachesize=1024);

    template<typename T> bool
    readComposite(const std::string& name, T& x, 
		      const MemberSubset& member_subset = MemberSubset());

    template<typename T> bool 
    readCompositeVector(const std::string& name, std::vector<T>& x,
		      const MemberSubset& member_subset = MemberSubset());

    template<typename T> VSOctaveH5ReaderCompositeVector<T>*
    readCompositeVectorElements(const std::string& name,
		      const MemberSubset& member_subset = MemberSubset(),
				unsigned cachesize=1024);

    static bool isHDF5(const std::string& filename);

    bool copyTo(const std::string& sname,
		VSOctaveH5WriterStruct* s, std::string dname = "");
    bool copyTo(const std::string& sname,
		VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol);
    bool copyTo(const std::string& sname,
		VSOctaveH5WriterCellVector* c, unsigned diel);
    bool copyTo(const std::string& sname,
		VSOctaveH5WriterExpandableCellVector* c);

    bool copyAllTo(VSOctaveH5WriterStruct* s);

  protected:
    virtual void notifyOfChildDestruction(VSOctaveH5ReaderBase* child);

  private:
    VSOctaveH5ReaderStruct(const VSOctaveH5ReaderStruct&);
    VSOctaveH5ReaderStruct& operator=(const VSOctaveH5ReaderStruct&);

    void readGroup();
    static herr_t iterateOverGrp(hid_t gid, const char*  name, void *object);
    void registerVariable(const char* name);
    void registerAssistant(VSOctaveH5ReaderBase *reader);

    Variable* resolveVariable(const std::string& name);
    static VSOctaveH5ReaderBase* loadVariable(Variable* var);

    template<typename T> bool copyRealScalar(const std::string& sname,
		   VSOctaveH5WriterStruct* s, const std::string dname);
    template<typename T> bool copyRealMatrix(const std::string& sname,
		   VSOctaveH5WriterStruct* s, const std::string dname);
    template<typename T> bool copyRealScalar(const std::string& sname,
		VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol);
    template<typename T> bool copyRealMatrix(const std::string& sname,
		VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol);
    template<typename T> bool copyRealScalar(const std::string& sname,
		VSOctaveH5WriterCellVector* c, unsigned diel);
    template<typename T> bool copyRealMatrix(const std::string& sname,
		VSOctaveH5WriterCellVector* c, unsigned diel);
    template<typename T> bool copyRealScalar(const std::string& sname,
		VSOctaveH5WriterExpandableCellVector* c);
    template<typename T> bool copyRealMatrix(const std::string& sname,
		VSOctaveH5WriterExpandableCellVector* c);

    bool                               m_my_file;
    hid_t                              m_gid;
    std::map<std::string, Variable*>   m_vars;
    hsize_t                            m_index;
    std::set<VSOctaveH5ReaderBase*>    m_assistants;
  };

  typedef VSOctaveH5ReaderStruct VSOctaveH5Reader;

  // **************************************************************************
  // **************************************************************************
  //
  // FUNCTION DEFINITIONS
  //
  // **************************************************************************
  // **************************************************************************

  // ==========================================================================
  //
  // BASE / VARIABLE
  //
  // ==========================================================================

  template<typename T> bool VSOctaveH5ReaderBase::Variable::
  getAttribute(const std::string& name, T& x)
  {
    if(attr.find(name) == attr.end())return false;
    hid_t tid = VSOctaveH5Type<T>::h5_type();
    hid_t aid = H5Aopen_name(gid, name.c_str());
    if((aid<0)||(H5Aread(aid, tid, &x)<0))
      {
	H5Tclose(tid);
	return false;
      }
    H5Aclose(aid);
    H5Tclose(tid);
    return true;
  }

  // ==========================================================================
  //
  // READER HELPER
  //
  // ==========================================================================

  template<typename T> inline herr_t VSOctaveH5ReaderActual<T>::
  read(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid,
       unsigned n, T* ptr)
  {
    return H5Dread(did, tid, core_sid, disk_sid, H5P_DEFAULT, ptr);
  }

  template<typename T> inline herr_t VSOctaveH5ReaderActual<T>::
  read_scalar(hid_t did, hid_t tid, T& x)
  {
    return H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &x);
  }

  template<typename T> inline herr_t VSOctaveH5ReaderActual<T>::
  read_vector(hid_t did, hid_t tid, std::vector<T>& v)
  {
    return H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*v.begin()));
  }
   
  template<typename T> inline herr_t VSOctaveH5ReaderActual<T>::
  read_matrix(hid_t did, hid_t tid, unsigned n, T* m)
  {
    return H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, m);
  }
  
  inline herr_t VSOctaveH5ReaderActual<bool>::
  read(hid_t did, hid_t tid, hid_t core_sid, hid_t disk_sid,
       unsigned n, bool* ptr)
  {
    std::vector<unsigned> v(n);
    herr_t status = 
      H5Dread(did, tid, core_sid, disk_sid, H5P_DEFAULT, &(*v.begin()));
    for(unsigned i=0;i<n;i++)ptr[i]=v[i];
    return status;
  }

  inline herr_t VSOctaveH5ReaderActual<bool>::
  read_scalar(hid_t did, hid_t tid, bool& x)
  {
    unsigned y;
    herr_t status = H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &y);
    x=y;
    return status;
  }

  inline herr_t VSOctaveH5ReaderActual<bool>::
  read_vector(hid_t did, hid_t tid, std::vector<bool>& v)
  {
    unsigned n = v.size();
    std::vector<unsigned> u(n);
    herr_t status = 
      H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*u.begin()));
    for(unsigned i=0;i<n;i++)v[i]=u[i];
    return status;
  }
   
  inline herr_t VSOctaveH5ReaderActual<bool>::
  read_matrix(hid_t did, hid_t tid, unsigned n, bool* m)
  {
    std::vector<unsigned> v(n);
    herr_t status = 
      H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*v.begin()));
    for(unsigned i=0;i<n;i++)m[i]=v[i];
    return status;
  }

  // ==========================================================================
  //
  // VECTOR
  //
  // ==========================================================================

  template<typename T> VSOctaveH5ReaderVector<T>::
  VSOctaveH5ReaderVector(VSOctaveH5ReaderMatrix* matrix,
			 unsigned cachesize):
    VSOctaveH5ReaderBase(matrix), 
    m_matrix(matrix), m_tid(VSOctaveH5Type<T>::h5_type()),
    m_cachesize(cachesize), m_cachestart(matrix->m_rows),
    m_cache(new T[cachesize])
  {
    // nothing to see here
  }
  
  template<typename T> VSOctaveH5ReaderVector<T>::~VSOctaveH5ReaderVector()
  {
    H5Tclose(m_tid);
    delete[] m_cache;
  }
  
  template<typename T> unsigned VSOctaveH5ReaderVector<T>::rows() const
  { 
    if(m_matrix->m_rows==1)return m_matrix->m_cols;
    else return m_matrix->m_rows; 
  }
  
  template<typename T> bool VSOctaveH5ReaderVector<T>::
  element(T& x, unsigned index)
  {
    unsigned nrows = rows();
    if(index>=nrows)return false;
    if(m_cachesize<=1)return read(&x,1,index);
    if(index<m_cachestart)
      {
	// assume we are iterating in reverse
	if(index<m_cachesize)m_cachestart=0;
	else m_cachestart=index-m_cachesize+1;
	unsigned n = m_cachesize;
	if(m_cachestart+n>nrows)n=nrows-m_cachestart;
	if(!read(m_cache,n,m_cachestart))return false;
      }
    else if(index>=m_cachestart+m_cachesize)
      {
	// assume we are iterating forwards
	m_cachestart=index;
	unsigned n = m_cachesize;
	if(m_cachestart+n>nrows)n=nrows-m_cachestart;
	if(!read(m_cache,n,m_cachestart))return false;	
      }
    x = m_cache[index-m_cachestart];
    return true;
  }

  template<typename T> bool VSOctaveH5ReaderVector<T>::
  read(T* ptr, hsize_t n, hsize_t index)
  {
    hid_t core_sid = H5Screate_simple(1,&n,NULL);
    hid_t disk_sid = H5Dget_space(m_matrix->m_did);
    if(H5Sget_simple_extent_ndims(disk_sid) == 1)
      {
	h5_select_hyperslab_t start = index;
	if(H5Sselect_hyperslab(disk_sid,H5S_SELECT_SET,
			       &start,NULL,&n,NULL) < 0)
	  {
	    VSOctaveH5Exception e(std::string("Could not select hyperslab"));
	    throw(e);
	  }
      }
    else
      {
	h5_select_hyperslab_t start[2] = { index, 0 };
	hsize_t count[2] = { n, 1 };
	if(m_matrix->m_rows == 1)
	  {
	    start[0] = 0, start[1] = index;
	    count[0] = 1, count[1] = n;
	  }
	
	if(H5Sselect_hyperslab(disk_sid,H5S_SELECT_SET,
			       start,NULL,count,NULL) < 0)
	  {
	    VSOctaveH5Exception e(std::string("Could not select hyperslab"));
	    throw(e);
	  }
      }
    if(VSOctaveH5ReaderActual<T>::
       read(m_matrix->m_did,m_tid,core_sid,disk_sid,n,ptr) < 0)
      {
	VSOctaveH5Exception 
	  e(std::string("Could not read dataset for hyperslab"));
	throw(e);
      }
    H5Sclose(disk_sid);
    H5Sclose(core_sid);
    return true;
  }

  // ==========================================================================
  //
  // COMPOSITE VECTOR
  //
  // ==========================================================================


#define __CASENEW(t)							\
  if(im->type == VSOctaveH5Type<t>::int_type())				\
    {									\
      if(s->isVector(im->name))						\
	{								\
	  VSOctaveH5ReaderVector<t>* vec =				\
	    s->readVectorElements<t>(im->name, cachesize);		\
	  m.reader = new VSOH5RCVERTemplate<t>(vec);			\
	  if(vec)							\
	    {								\
	      unsigned rows = vec->rows();				\
	      if((m_members.empty())||(rows>m_rows))m_rows=rows;	\
	    }								\
	}								\
      else								\
	{								\
	  m.reader = new VSOH5RCVERNull<t>;				\
	}								\
      m_members.push_back(m);						\
      continue;								\
    }

  template<typename T> VSOctaveH5ReaderCompositeVector<T>::
  VSOctaveH5ReaderCompositeVector(VSOctaveH5ReaderStruct* s,
				  const std::string& prefix,
				  const MemberSubset& member_subset,
				  unsigned cachesize):
    VSOctaveH5ReaderBase(s), m_members(), m_rows()
  {
    VSOctaveH5CompositeDefinition c(prefix, member_subset);
    VSOctaveH5CompositeCompose<T>::compose(c);
    const std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
    for(std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator 
	  im = m.begin(); im!=m.end(); im++)
      {
	MemReader m(*im);
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
	    if(s->isCellVector(im->name))
	      {
		VSOctaveH5ReaderCellVector* cv = s->readCellVector(im->name);
		m.reader = new VSOH5RCVERString(cv);
		if(cv)
		  {
		    unsigned rows = cv->dimensions();
		    if((m_members.empty())||(rows>m_rows))m_rows=rows;
		  }
	      }
	    else
	      {
		m.reader = new VSOH5RCVERNull<std::string>;
	      }
	    m_members.push_back(m);
	    continue;
	  }
	
	vsassert(0);
      }
  }
#undef __CASENEW

  template<typename T> VSOctaveH5ReaderCompositeVector<T>::
  ~VSOctaveH5ReaderCompositeVector()
  {
    for(typename std::vector<MemReader>::const_iterator im = 
	  m_members.begin(); im!=m_members.end(); im++)
      delete im->reader;
  }
  
  template<typename T> bool VSOctaveH5ReaderCompositeVector<T>::
  element(T& x, unsigned index)
  {
    if(index>=m_rows)return false;
    bool status = true;
    for(typename std::vector<MemReader>::const_iterator im = 
	  m_members.begin(); im!=m_members.end(); im++)
      status &= im->reader->read(&x, index, im->field);
    return status;
  }

  template<typename T> bool VSOctaveH5ReaderCompositeVector<T>::
  hasElement(const std::string& name) const
  { 
    for(typename std::vector<MemReader>::const_iterator im = 
	  m_members.begin(); im!=m_members.end(); im++)
      if(im->field.name == name)return !im->reader->isNullReader();
    return false;
  }

  // ==========================================================================
  //
  // MATRIX
  //
  // ==========================================================================

  template<typename T>
  bool VSOctaveH5ReaderMatrix::readScalar(T& x)
  {
    if(!isScalar())return false;
    hid_t core_tid = VSOctaveH5Type<T>::h5_type();
    if(VSOctaveH5ReaderActual<T>::read_scalar(m_did,core_tid,x) < 0)
      {
	VSOctaveH5Exception e(std::string("Could not read from dataset"));
	throw(e);
      }
    H5Tclose(core_tid);
    return true;
  }
    
  template<typename T>
  bool VSOctaveH5ReaderMatrix::readVector(std::vector<T>& v)
  {
    if(!isVector())return false;
    if(isEmpty()){ v.clear(); return true; }
    if(m_rows!=1)v.resize(m_rows);
    else v.resize(m_cols);
    hid_t core_tid = VSOctaveH5Type<T>::h5_type();
    if(VSOctaveH5ReaderActual<T>::read_vector(m_did,core_tid,v) < 0)
      {
	VSOctaveH5Exception e(std::string("Could not read from dataset"));
	throw(e);
      }
    H5Tclose(core_tid);
    return true;
  }
    
  template<typename T>
  inline bool VSOctaveH5ReaderMatrix::
  readMatrix(unsigned& rows, unsigned& cols, T*& m)
  {
    rows = m_rows;
    cols = m_cols;
    if(isEmpty()) { m = 0; return true; }
    m = new T[rows*cols];
    return readMatrix(m);
  }

  template<typename T>
  bool VSOctaveH5ReaderMatrix::readMatrix(T* m)
  {
    if(isEmpty()) { return true; }
    hid_t core_tid = VSOctaveH5Type<T>::h5_type();
    if(VSOctaveH5ReaderActual<T>::
       read_matrix(m_did,core_tid,m_rows*m_cols,m) < 0)
      {
	VSOctaveH5Exception e(std::string("Could not read from dataset"));
	throw(e);
      }
    H5Tclose(core_tid);
    return true;
  }

  template<typename T> VSOctaveH5ReaderVector<T>* 
  VSOctaveH5ReaderMatrix::readVectorElements(unsigned cachesize)
  {
    if(!isVector())return 0;
    VSOctaveH5ReaderVector<T>* vec = 
      new VSOctaveH5ReaderVector<T>(this,cachesize);
    m_children.insert(vec);
    return vec;
  }

  // ==========================================================================
  //
  // CELL ARRAY
  //
  // ==========================================================================

  template<typename T> bool VSOctaveH5ReaderCellArray::
  readScalar(unsigned row, unsigned col, T& x)
  {
    return m_struct->readScalar(name(row,col),x);
  }
  
  template<typename T> bool VSOctaveH5ReaderCellArray::
  readVector(unsigned row, unsigned col, std::vector<T>& v)
  {
    return m_struct->readVector(name(row,col),v);
  }
  
  template<typename T> bool VSOctaveH5ReaderCellArray::
  readMatrix(unsigned row, unsigned col, unsigned& rows, unsigned& cols, T*& m)
  {
    return m_struct->readMatrix(name(row,col),rows,cols,m);
  }
  
  template<typename T> bool VSOctaveH5ReaderCellArray::
  readMatrix(unsigned row, unsigned col, T* m)
  {
    return m_struct->readMatrix(name(row,col),m);
  }

  template<typename T> VSOctaveH5ReaderVector<T>* VSOctaveH5ReaderCellArray::
  readVectorElements(unsigned row, unsigned col, unsigned cachesize)    
  {
    return m_struct->readVectorElements<T>(name(row,col),cachesize);
  }

  template<typename T> bool VSOctaveH5ReaderCellArray::
  readComposite(unsigned row, unsigned col, T& x, 
		const MemberSubset& member_subset)
  {
    return m_struct->readComposite(name(row,col),x,member_subset);
  }

  template<typename T> bool VSOctaveH5ReaderCellArray::
  readCompositeVector(unsigned row, unsigned col, std::vector<T>& x,
		      const MemberSubset& member_subset)
  {
    return m_struct->readCompositeVector(name(row,col),x,member_subset);
  }

  template<typename T> VSOctaveH5ReaderCompositeVector<T>*
  VSOctaveH5ReaderCellArray::
  readCompositeVectorElements(unsigned row, unsigned col,
			      const MemberSubset& member_subset, 
			      unsigned cachesize)
  {
    return m_struct->readCompositeVectorElements<T>(name(row,col),
						    member_subset,cachesize);
  }

  // ==========================================================================
  //
  // CELL VECTOR
  //
  // ==========================================================================

  template<typename T> bool VSOctaveH5ReaderCellVector::
  readScalar(unsigned iel, T& x)
  {
    return m_cell->readScalar(row(iel),col(iel),x);
  }
  
  template<typename T> bool VSOctaveH5ReaderCellVector::
  readVector(unsigned iel, std::vector<T>& v)
  {
    return m_cell->readVector(row(iel),col(iel),v);
  }
  
  template<typename T> bool VSOctaveH5ReaderCellVector::
  readMatrix(unsigned iel, unsigned& rows, unsigned& cols, T*& m)
  {
    return m_cell->readMatrix(row(iel),col(iel),rows,cols,m);
  }
  
  template<typename T> bool VSOctaveH5ReaderCellVector::
  readMatrix(unsigned iel, T* m)
  {
    return m_cell->readMatrix(row(iel),col(iel),m);
  }

  template<typename T> VSOctaveH5ReaderVector<T>* VSOctaveH5ReaderCellVector::
  readVectorElements(unsigned iel, unsigned cachesize)    
  {
    return m_cell->readVectorElements<T>(row(iel),col(iel),cachesize);
  }

  template<typename T> bool VSOctaveH5ReaderCellVector::
  readComposite(unsigned iel, T& x,
		const MemberSubset& member_subset)
  {
    return m_cell->readComposite(row(iel),col(iel),x,member_subset);
  }

  template<typename T> bool VSOctaveH5ReaderCellVector::
  readCompositeVector(unsigned iel, std::vector<T>& x,
		      const MemberSubset& member_subset)
  {
    return m_cell->readCompositeVector(row(iel),col(iel),x,member_subset);
  }

  template<typename T> VSOctaveH5ReaderCompositeVector<T>*
  VSOctaveH5ReaderCellVector::
  readCompositeVectorElements(unsigned iel, const MemberSubset& member_subset, 
			      unsigned cachesize)
  {
    return m_cell->readCompositeVectorElements<T>(row(iel),col(iel),
						  member_subset,cachesize);
  }

  // ==========================================================================
  //
  // STRUCT
  //
  // ==========================================================================
  
  template<typename T> bool VSOctaveH5ReaderStruct::
  readScalar(const std::string& name, T& x)
  {
    Variable* var = resolveVariable(name);
    if((var == 0)
       ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	  (var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
    VSOctaveH5ReaderBase* reader = loadVariable(var);
    VSOctaveH5ReaderMatrix* mtx = 
      dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
    vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
    return mtx->readScalar(x);
  }
  
  template<typename T> bool VSOctaveH5ReaderStruct::
  readVector(const std::string& name, std::vector<T>& v)
  {
    Variable* var = resolveVariable(name);
    if((var == 0)
       ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	  (var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
    VSOctaveH5ReaderBase* reader = loadVariable(var);
    VSOctaveH5ReaderMatrix* mtx = 
      dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
    vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
    return mtx->readVector(v);
  }
  
  template<typename T> bool VSOctaveH5ReaderStruct::
  readMatrix(const std::string& name, unsigned& rows, unsigned& cols, T*& m)
  {
    Variable* var = resolveVariable(name);
    if((var == 0)
       ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	  (var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
    VSOctaveH5ReaderBase* reader = loadVariable(var);
    VSOctaveH5ReaderMatrix* mtx = 
      dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
    vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
    return mtx->readMatrix(rows,cols,m);
  }
  
  template<typename T> bool VSOctaveH5ReaderStruct::
  readMatrix(const std::string& name, T* m)
  {
    Variable* var = resolveVariable(name);
    if((var == 0)
       ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	  (var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
    VSOctaveH5ReaderBase* reader = loadVariable(var);
    VSOctaveH5ReaderMatrix* mtx = 
      dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
    vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
    return mtx->readMatrix(m);
  }

  template<typename T> VSOctaveH5ReaderVector<T>* VSOctaveH5ReaderStruct::
  readVectorElements(const std::string& name, unsigned cachesize)
  {
    Variable* var = resolveVariable(name);
    if((var == 0)
       ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	  (var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
    VSOctaveH5ReaderBase* reader = loadVariable(var);
    VSOctaveH5ReaderMatrix* mtx = 
      dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
    vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
    return mtx->readVectorElements<T>(cachesize);
  }

#define __CASESCALAR(t,x,im)						\
  if(im->type == VSOctaveH5Type<t>::int_type())				\
    {									\
      status &= readScalar(im->name,*((t*)((char*)(&x) + im->offset))); \
      continue;								\
    }  

  template<typename T> bool VSOctaveH5ReaderStruct::
  readCompositeHere(T& x, const std::string& prefix,
		    const MemberSubset& member_subset)
  {
    bool status = true;
    VSOctaveH5CompositeDefinition c(prefix,member_subset);
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
	      readString(im->name,*((std::string*)((char*)(&x) + im->offset)));
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
      std::vector<t> v;							\
      status &= readVector(im->name,v);					\
      unsigned nel = v.size();						\
      if(nel>x.size())x.resize(nel);					\
      for(unsigned iel=0;iel<nel;iel++)					\
	*((t*)((char*)(&x[iel]) + im->offset)) = v[iel];		\
      continue;								\
    }

  template<typename T> bool VSOctaveH5ReaderStruct::
  readStructCellVector(const std::string& name, std::vector<T>& x)
  {
    bool status = true;
    VSOctaveH5ReaderCellVector* wc = readCellVector(name);
    const unsigned n = wc->dimensions();
    x.resize(n);
    for(unsigned i = 0; i < n; i++)
      {
	VSOctaveH5ReaderStruct* ws = wc->readStruct(i);
	x[i].load(ws);
	delete ws;
      }
    delete wc;
    return status;
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  readCompositeVectorHere(std::vector<T>& x, const std::string& prefix,
			  const MemberSubset& member_subset)
  {
    bool status = true;
    VSOctaveH5CompositeDefinition c(prefix,member_subset);
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
	    VSOctaveH5ReaderCellVector* cv = readCellVector(im->name);	    
	    if(cv)
	      {
		unsigned nel = cv->dimensions();
		if(nel>x.size())x.resize(nel);
		for(unsigned iel=0;iel<nel;iel++)
		  status &= cv->readString(iel,
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

  template<typename T> VSOctaveH5ReaderCompositeVector<T>* 
  VSOctaveH5ReaderStruct::
  readCompositeVectorElementsHere(const std::string& prefix,
				  const MemberSubset& member_subset,
				  unsigned cachesize)
  {
    VSOctaveH5ReaderCompositeVector<T>* v =
      new VSOctaveH5ReaderCompositeVector<T>(this, prefix, 
					     member_subset, cachesize);
    registerAssistant(v);
    return v;
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  readComposite(const std::string& name, T& x,
		const MemberSubset& member_subset)
  {
    VSOctaveH5ReaderStruct* s = readStruct(name);
    if(s==0)return false;
    return s->readCompositeHere(x, std::string(), member_subset);
  }
  
  template<typename T> bool VSOctaveH5ReaderStruct::
  readCompositeVector(const std::string& name, std::vector<T>& x,
		      const MemberSubset& member_subset)
  {
    VSOctaveH5ReaderStruct* s = readStruct(name);
    if(s==0)return false;
    return s->readCompositeVectorHere(x,std::string(),member_subset);
  }
  
  template<typename T> VSOctaveH5ReaderCompositeVector<T>* 
  VSOctaveH5ReaderStruct::
  readCompositeVectorElements(const std::string& name,
			      const MemberSubset& member_subset,
			      unsigned cachesize)
  {
    VSOctaveH5ReaderStruct* s = readStruct(name);
    if(s==0)return false;
    return s->readCompositeVectorElementsHere<T>(std::string(),member_subset,
						 cachesize);
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealScalar(const std::string& sname,
		 VSOctaveH5WriterStruct* s, const std::string dname)
  {
    T x;
    if(!readScalar(sname, x))return false;
    return s->writeScalar(dname, x);
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealScalar(const std::string& sname,
		 VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol)
  {
    T x;
    if(!readScalar(sname, x))return false;
    return c->writeScalar(dirow, dicol, x);
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealScalar(const std::string& sname,
		 VSOctaveH5WriterCellVector* c, unsigned diel)
  {
    T x;
    if(!readScalar(sname, x))return false;
    return c->writeScalar(diel, x);
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealScalar(const std::string& sname,
		 VSOctaveH5WriterExpandableCellVector* c)
  {
    T x;
    if(!readScalar(sname, x))return false;
    return c->appendScalar(x);
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealMatrix(const std::string& sname,
		 VSOctaveH5WriterStruct* s, const std::string dname)
  {
    unsigned nrow = 0;
    unsigned ncol = 0;
    T* m;
    if(!readMatrix(sname, nrow, ncol, m))return false;
    bool _ret = s->writeMatrix(dname, nrow, ncol, m);
    delete[] m;
    return _ret;
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealMatrix(const std::string& sname,
		 VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol)
  {
    unsigned nrow = 0;
    unsigned ncol = 0;
    T* m;
    if(!readMatrix(sname, nrow, ncol, m))return false;
    bool _ret = c->writeMatrix(dirow, dicol, nrow, ncol, m);
    delete[] m;
    return _ret;
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealMatrix(const std::string& sname,
		 VSOctaveH5WriterCellVector* c, unsigned diel)
  {
    unsigned nrow = 0;
    unsigned ncol = 0;
    T* m;
    if(!readMatrix(sname, nrow, ncol, m))return false;
    bool _ret = c->writeMatrix(diel, nrow, ncol, m);
    delete[] m;
    return _ret;
  }

  template<typename T> bool VSOctaveH5ReaderStruct::
  copyRealMatrix(const std::string& sname,
		 VSOctaveH5WriterExpandableCellVector* c)
  {
    unsigned nrow = 0;
    unsigned ncol = 0;
    T* m;
    if(!readMatrix(sname, nrow, ncol, m))return false;
    bool _ret = c->appendMatrix(nrow, ncol, m);
    delete[] m;
    return _ret;
  }

} // namespace VERITAS

#endif // VSOCTAVEH5READER_HPP
