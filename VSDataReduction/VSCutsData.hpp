//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutsData.hpp

  Data structures for array and scope cuts

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/08/2007

  $Id: VSCutsData.hpp,v 3.5 2008/11/24 02:02:02 matthew Exp $

*/

#ifndef VSCUTSDATA_HPP
#define VSCUTSDATA_HPP

#include<vector>

#include<VSOctaveIO.hpp>

namespace VERITAS
{

  struct VSScopeCutsDatum
  {
    VSScopeCutsDatum():
      passed_scope_cut()
    { /* nothing to see here */ }

    bool                           passed_scope_cut;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSScopeCutsDatum,passed_scope_cut);
    }
  };

  struct VSArrayCutsDatum
  {
    VSArrayCutsDatum(const std::vector<unsigned>& nchan = 
		     std::vector<unsigned>()):
      passed_cuts(),
      passed_array_cut(), npassed_scope_cut(), scope(nchan.size())
    { 
      for(unsigned iscope=0;iscope<nchan.size();iscope++)
	if(nchan[iscope])scope[iscope]=new VSScopeCutsDatum;
    }

    bool                           passed_cuts;
    bool                           passed_array_cut;
    unsigned                       npassed_scope_cut;
    std::vector<VSScopeCutsDatum*> scope;

    void add(VSArrayCutsDatum &datum)
    {
      passed_cuts &= datum.passed_cuts;
      passed_array_cut &= datum.passed_array_cut;
      const unsigned nscope = scope.size();
      vsassert(datum.scope.size() == nscope);
      npassed_scope_cut = 0;
      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  scope[iscope]->passed_scope_cut &= 
	    datum.scope[iscope]->passed_scope_cut;
	  if(scope[iscope]->passed_scope_cut)
	    npassed_scope_cut++;
	}
    }

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSArrayCutsDatum,passed_cuts);
      H5_ADDMEMBER(c,VSArrayCutsDatum,passed_array_cut);
      H5_ADDMEMBER(c,VSArrayCutsDatum,npassed_scope_cut);
    }
  };

  class VSArrayCutsWriter
  {
  public:
    VSArrayCutsWriter(VSOctaveH5WriterStruct* s,
		      const std::vector<unsigned>& nchan);
    virtual ~VSArrayCutsWriter();

    bool append(const VSArrayCutsDatum& x);

  private:
    VSArrayCutsWriter(const VSArrayCutsWriter&);
    VSArrayCutsWriter& operator=(const VSArrayCutsWriter&);

    typedef VSOctaveH5WriterCompositeVector<VSArrayCutsDatum> ArrayCutsWriter;
    typedef VSOctaveH5WriterCompositeVector<VSScopeCutsDatum> ScopeCutsWriter;

    ArrayCutsWriter*                              m_array_writer;
    std::vector<ScopeCutsWriter*>                 m_scope_writer;
  };

  class VSArrayCutsReader
  {
  public:
    VSArrayCutsReader(VSOctaveH5ReaderStruct* s);
    VSArrayCutsReader();
    virtual ~VSArrayCutsReader();

    bool element(VSArrayCutsDatum& x, unsigned index);

    inline unsigned rows() const { return m_array_reader->rows(); }

    VSArrayCutsDatum at(unsigned index) 
    { 
      VSArrayCutsDatum x; 
      if(!element(x,index))throw std::out_of_range(__PRETTY_FUNCTION__); 
      return x; 
    }

    VSArrayCutsDatum operator[] (unsigned index) 
    {
      VSArrayCutsDatum x; 
      element(x,index); 
      return x;
    }

    static bool loadAllArrayCuts(VSOctaveH5ReaderStruct* s, 
				 std::vector<VSArrayCutsDatum>& x);

  private:
    VSArrayCutsReader(const VSArrayCutsReader&);
    VSArrayCutsReader& operator=(const VSArrayCutsReader&);

    typedef VSOctaveH5ReaderCompositeVector<VSArrayCutsDatum> ArrayCutsReader;
    typedef VSOctaveH5ReaderCompositeVector<VSScopeCutsDatum> ScopeCutsReader;

    ArrayCutsReader*                              m_array_reader;
    std::vector<ScopeCutsReader*>                 m_scope_reader;
  };

}


#endif // defined VSCUTSDATA_HPP
