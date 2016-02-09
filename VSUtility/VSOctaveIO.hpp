//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveIO.hpp

  Class to write to GNU/Octave file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/20/2006
*/

#ifndef VSOCTAVEIO_HPP
#define VSOCTAVEIO_HPP

#include<string>
#include<vector>
#include<set>
#include<cstddef>
#include<stdint.h>

#include<hdf5.h>

#include<VSDataConverter.hpp>

namespace VERITAS
{
  
  // ==========================================================================
  // Exception
  // ==========================================================================

  class VSOctaveH5Exception
  {
  public:
    VSOctaveH5Exception(const std::string& _message): m_message(_message) { }
    virtual ~VSOctaveH5Exception();
    std::string message() const { return m_message; }
  private:
    std::string m_message;
  };

  // ==========================================================================
  // Type : Generate HDF5 type for variable
  // ==========================================================================
  
  enum VSOctaveH5InternalType { OIT_BOOL, 
				OIT_INT8,
				OIT_INT16,
				OIT_INT32,
				OIT_INT64,
				OIT_UINT8,
				OIT_UINT16,
				OIT_UINT32,
				OIT_UINT64,
				OIT_FLOAT,
				OIT_DOUBLE,
				OIT_LONG_DOUBLE,
				OIT_STRING };

  template<typename T> class VSOctaveH5Type
  {
  public:
    static VSOctaveH5InternalType int_type();
    static std::string oct_type();
    static hid_t h5_type();
  };

#define DEFINE_H5OCTTYPE(T,i,s,x)			   \
  template<> class VSOctaveH5Type<T>			   \
  {							   \
  public:						   \
    static VSOctaveH5InternalType int_type() { return i; } \
    static std::string oct_type() { return s; }		   \
    static hid_t h5_type(){ return x; }			   \
  }
  
  DEFINE_H5OCTTYPE(bool, OIT_BOOL, "bool", H5Tcopy(H5T_NATIVE_UINT));

  DEFINE_H5OCTTYPE(int8_t, OIT_INT8, "int8", 
		   H5Tget_native_type(H5T_STD_I8LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(int16_t, OIT_INT16, "int16", 
		   H5Tget_native_type(H5T_STD_I16LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(int32_t, OIT_INT32, "int32", 
		   H5Tget_native_type(H5T_STD_I32LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(int64_t, OIT_INT64, "int64", 
		   H5Tget_native_type(H5T_STD_I64LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(uint8_t, OIT_UINT8, "uint8", 
		   H5Tget_native_type(H5T_STD_U8LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(uint16_t, OIT_UINT16, "uint16", 
		   H5Tget_native_type(H5T_STD_U16LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(uint32_t, OIT_UINT32, "uint32", 
		   H5Tget_native_type(H5T_STD_U32LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(uint64_t, OIT_UINT64, "uint64", 
		   H5Tget_native_type(H5T_STD_U64LE,H5T_DIR_ASCEND));

  DEFINE_H5OCTTYPE(float, OIT_FLOAT, "", H5Tcopy(H5T_NATIVE_FLOAT));

  DEFINE_H5OCTTYPE(double, OIT_DOUBLE, "", H5Tcopy(H5T_NATIVE_DOUBLE));

#if __WORDSIZE == 64 && defined __GLIBC_HAVE_LONG_LONG
  DEFINE_H5OCTTYPE(long long, OIT_INT64, "int64", 
		   H5Tget_native_type(H5T_STD_I64LE,H5T_DIR_ASCEND));
  DEFINE_H5OCTTYPE(unsigned long long, OIT_UINT64, "uint64", 
		   H5Tget_native_type(H5T_STD_U64LE,H5T_DIR_ASCEND));
#endif

#if __WORDSIZE == 32
  DEFINE_H5OCTTYPE(long, OIT_INT32, "int32", 
		   H5Tget_native_type(H5T_STD_I32LE,H5T_DIR_ASCEND));
  DEFINE_H5OCTTYPE(unsigned long, OIT_UINT32, "uint32", 
		   H5Tget_native_type(H5T_STD_U32LE,H5T_DIR_ASCEND));
#endif  

#if 0 // probably not supported by octave
  DEFINE_H5OCTTYPE(long double, OIT_LONG_DOUBLE, 
		   "", H5Tcopy(H5T_NATIVE_LDOUBLE));
#endif

#define DEFINE_H5ENUM(T)						\
  DEFINE_H5OCTTYPE(T, OIT_INT32, "int32",				\
		   H5Tget_native_type(H5T_STD_I32LE,H5T_DIR_ASCEND));

  enum VSOctaveH5VecOrientation { VO_ROW, VO_COLUMN };

#define OCTAVE_H5_EL_TYPE              "type"
#define OCTAVE_H5_EL_VALUE             "value"
#define OCTAVE_H5_EL_DIMS              "dims"

#define OCTAVE_H5_TYPE_CELL            "cell"
#define OCTAVE_H5_TYPE_STRUCT          "struct"
#define OCTAVE_H5_TYPE_STRING          "sq_string"
#define OCTAVE_H5_TYPE_MATRIX          "matrix"
#define OCTAVE_H5_TYPE_SCALAR          "scalar"

#define OCTAVE_H5_ATTR_NEW_FORMAT      "OCTAVE_NEW_FORMAT"
#define OCTAVE_H5_ATTR_EMPTY_MATRIX    "OCTAVE_EMPTY_MATRIX"

  class VSOctaveH5Common
  {
  public:
    static unsigned cellPrecision(unsigned nel);
  };

  // ==========================================================================
  // Octave H5 composite
  // ==========================================================================

  template<typename T> 
  class VSOctaveH5CompositeMemberType
  {
  public:
    static VSOctaveH5InternalType int_type() 
    { 
      return VSOctaveH5Type<T>::int_type();
    }
  };

  template<> class VSOctaveH5Type<std::string>
  {
  public:
    static VSOctaveH5InternalType int_type() { return OIT_STRING; }
  };

  class VSOctaveH5CompositeDefinition
  {
  public:
    class Member
    {
    public:
      Member(const std::string& _name, VSOctaveH5InternalType _type,
	     size_t _offset): name(_name), type(_type), offset(_offset) { }
      std::string            name;
      VSOctaveH5InternalType type;
      size_t                 offset;
    };
    
    typedef std::set<std::string> MemberSubset;

    VSOctaveH5CompositeDefinition(const std::string& prefix = "",
			   const MemberSubset member_subset = MemberSubset()): 
      m_prefix(), m_members(), m_member_subset(member_subset)
    {
      if(!prefix.empty())m_prefix = prefix + std::string("_");
    }

    ~VSOctaveH5CompositeDefinition() { }

    template<typename T> 
    void addMember(const std::string& name,size_t offset)
    {
      if(!isInSubset(name))return;
      VSOctaveH5InternalType t = VSOctaveH5CompositeMemberType<T>::int_type();
      m_members.push_back(Member(m_prefix+name, t, offset));
    }

    template<typename T> 
    void addMember(const T* __attribute((unused)) unused_arg,
		   const std::string& name, size_t offset)
    { 
      addMember<T>(name,offset);
    }

    template<typename T> 
    void addSimpleComposite(const std::string& name,size_t offset);

    template<typename T> 
    void addSimpleComposite(const T* __attribute((unused)) unused_arg,
			    const std::string& name, size_t offset)
    { 
      addSimpleComposite<T>(name,offset);
    }

    const std::vector<Member>& members() const { return m_members; }

  private:
    std::string           m_prefix;
    std::vector<Member>   m_members;
    MemberSubset          m_member_subset;

    bool isInSubset(const std::string& s)
    { return m_member_subset.empty() 
	|| m_member_subset.find(s)!=m_member_subset.end(); }
  };

#define H5_OFFSETOF(T,M) ((size_t)(&((T*)1000)->M)-1000)
#define H5_POINTERTO(T,M) (&((T*)1000)->M)

#define H5_ADDMEMBER(c,TYPE,MEMBER)					\
  (c.addMember(H5_POINTERTO(TYPE,MEMBER), #MEMBER, H5_OFFSETOF(TYPE,MEMBER)))
    
#define H5_ADDSIMPLECOMPOSITE(c,TYPE,MEMBER)			\
  (c.addSimpleComposite(H5_POINTERTO(TYPE,MEMBER), #MEMBER,	\
			H5_OFFSETOF(TYPE,MEMBER)))
  
#define H5_ADDNAMEDMEMBER(c,TYPE,MEMBER,NAME)				\
  (c.addMember(H5_POINTERTO(TYPE,MEMBER), NAME, H5_OFFSETOF(TYPE,MEMBER)))
  
#define H5_ADDNAMEDSIMPLECOMPOSITE(c,TYPE,MEMBER,NAME)		\
  (c.addSimpleComposite(H5_POINTERTO(TYPE,MEMBER), NAME,	\
			H5_OFFSETOF(TYPE,MEMBER)))

#if (H5_VERS_MAJOR>1)||(H5_VERS_MINOR>6)||(H5_VERS_RELEASE>2)
  typedef hsize_t h5_select_hyperslab_t;
#else
  typedef hssize_t h5_select_hyperslab_t;
#endif

  template<typename T> class VSOctaveH5CompositeCompose
  {
  public:
    static void compose(VSOctaveH5CompositeDefinition& c)
    {
      return T::_compose(c);
    }
    
    static VSOctaveH5CompositeDefinition compose()
    {
      VSOctaveH5CompositeDefinition c; compose(c); return c;
    }
  };

  template<typename T1, typename T2> 
  class VSOctaveH5CompositeCompose<std::pair<T1, T2> >
  {
  public:
    typedef std::pair<T1,T2> pair_type;

    static void compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDNAMEDMEMBER(c,pair_type,first,"first");
      H5_ADDNAMEDMEMBER(c,pair_type,second,"second");
      return;
    }
    
    static VSOctaveH5CompositeDefinition compose()
    {
      VSOctaveH5CompositeDefinition c; compose(c); return c;
    }
  };

  template<typename T1, typename T2, typename T3> 
  class VSOctaveH5CompositeCompose<triple<T1, T2, T3> >
  {
  public:
    typedef triple<T1,T2,T3> triple_type;

    static void compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDNAMEDMEMBER(c,triple_type,first,"first");
      H5_ADDNAMEDMEMBER(c,triple_type,second,"second");
      H5_ADDNAMEDMEMBER(c,triple_type,third,"third");
      return;
    }
    
    static VSOctaveH5CompositeDefinition compose()
    {
      VSOctaveH5CompositeDefinition c; compose(c); return c;
    }
  };

  template<typename T> void VSOctaveH5CompositeDefinition::
  addSimpleComposite(const std::string& name,size_t offset)
  {
    if(!isInSubset(name))return;
    VSOctaveH5CompositeDefinition c(name);
    VSOctaveH5CompositeCompose<T>::compose(c);
    const std::vector<VSOctaveH5CompositeDefinition::Member> m = c.members();
    for(std::vector<VSOctaveH5CompositeDefinition::Member>::const_iterator 
	  im = m.begin(); im!=m.end(); im++)
      m_members.push_back(Member(im->name, im->type, offset+im->offset));
  }
  
} // namespace VERITAS

#include "VSOctaveH5Writer.hpp"
#include "VSOctaveH5Reader.hpp"

#endif // VSOCTAVEIO_HPP
